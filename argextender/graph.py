"""Next-round targets from the assembly graph (strategy B-graph).

1. Anchors. `Bandage reduce GRAPH OUT --scope aroundblast --query SEEDS --distance 0 --ifilter 90
   --alfilter 100` BLASTs the seeds against the segments of SPAdes' `strain_graph.gfa`; the segments with a
   hit (identity >= 90%, alignment >= 100 bp) are the ARG anchors.
2. Components. The GFA becomes an undirected graph (segments = nodes, links = edges, orientation ignored).
   Every connected component holding an anchor is kept whole; the others are dropped.
3. Targets. Each kept segment is written unchanged as its own record (`>seg_<name>_arg` for anchors,
   `>seg_<name>_component` for the rest), never joined. Order: anchors by name, then the other segments by
   link distance from the nearest anchor, then by name.

If the graph is missing or unusable, or has no anchor, the round falls back to the selected scaffolds as
targets and the reason is recorded.
"""

from __future__ import annotations

import os
import re
import resource
import shutil
import subprocess
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

import networkx as nx

from argextender import ArgExtenderError
from argextender.spades import Assembly

BANDAGE_FILTERS = ("--ifilter", "90", "--alfilter", "100")
NO_HITS = "To draw the graph around BLAST hits, you must first conduct a BLAST search."
OVERLAP_RE = re.compile(r"^(\d+)M$")


class GraphError(ArgExtenderError):
    pass


@dataclass(frozen=True)
class Link:
    a: str
    b: str


@dataclass
class Graph:
    segments: dict[str, str]  # name -> sequence
    links: list[Link]


def parse_gfa(path: Path) -> Graph:
    """Segments and links of a GFA 1 file (SPAdes); GraphError on malformed records."""
    segments: dict[str, str] = {}
    links: list[Link] = []
    for n, line in enumerate(path.read_text().splitlines(), start=1):
        if not line.strip() or line.startswith("#"):
            continue
        f = line.split("\t")
        if f[0] == "S":
            if len(f) < 3 or not f[1] or f[1] in segments:
                raise GraphError(f"{path}:{n}: bad or duplicate segment record")
            segments[f[1]] = f[2]
        elif f[0] == "L":
            if len(f) < 6 or f[2] not in "+-" or f[4] not in "+-" or OVERLAP_RE.match(f[5]) is None:
                raise GraphError(f"{path}:{n}: bad link record")
            links.append(Link(f[1], f[3]))
        elif f[0] not in ("H", "P"):
            raise GraphError(f"{path}:{n}: unsupported GFA record type '{f[0]}'")
    missing = sorted({x for link in links for x in (link.a, link.b) if x not in segments})
    if missing:
        raise GraphError(f"{path}: links refer to missing segments {missing[:5]}")
    return Graph(segments, links)


def natural(name: str) -> tuple[int, str]:
    """Segment names in numeric order (SPAdes names are numbers), other names after them."""
    return (int(name), "") if name.isdigit() else (1 << 30, name)


def select_components(graph: Graph, anchors: Iterable[str]) -> list[tuple[str, str]]:
    """(segment, role) of every segment in a connected component holding an anchor, in target order."""
    anchor_set = set(anchors)
    unknown = sorted(anchor_set - set(graph.segments))
    if unknown:
        raise ArgExtenderError(f"anchors not in the graph: {', '.join(unknown[:5])}")
    g = nx.Graph()
    g.add_nodes_from(graph.segments)
    g.add_edges_from((link.a, link.b) for link in graph.links)
    chosen: set[str] = set().union(*(c for c in nx.connected_components(g) if c & anchor_set))
    dist = nx.multi_source_dijkstra_path_length(g, anchor_set)
    rest = sorted(chosen - anchor_set, key=lambda s: (int(dist[s]), natural(s)))
    return [*((s, "arg") for s in sorted(anchor_set, key=natural)), *((s, "component") for s in rest)]


def _no_core() -> None:
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))


@dataclass
class Bandage:
    exe: str
    log: Path | None = None

    def anchors(self, graph: Path, seeds: Path, work: Path) -> list[str]:
        """Names of the segments with a filtered BLAST hit of a seed ([] if Bandage finds no hit at all)."""
        scratch = work / "bandage_run"
        shutil.rmtree(scratch, ignore_errors=True)
        (scratch / "runtime").mkdir(parents=True, mode=0o700)
        (scratch / "tmp").mkdir()
        out = work / "bandage_hits.gfa"
        out.unlink(missing_ok=True)
        cmd = [
            self.exe,
            "reduce",
            str(graph),
            str(out),
            "--scope",
            "aroundblast",
            "--query",
            str(seeds),
            "--distance",
            "0",
            *BANDAGE_FILTERS,
        ]
        env = {
            **os.environ,
            "QT_QPA_PLATFORM": "offscreen",  # a Qt application, run headless
            "XDG_RUNTIME_DIR": str(scratch / "runtime"),
            "TMPDIR": str(scratch / "tmp"),
        }
        # Bandage segfaults when no BLAST hit remains; no core dump for that
        proc = subprocess.run(
            cmd, cwd=scratch, capture_output=True, text=True, check=False, env=env, preexec_fn=_no_core
        )
        if self.log:
            with self.log.open("a") as log:
                log.write("Running cmd:\t" + " ".join(cmd) + "\n" + proc.stdout + proc.stderr)
        shutil.rmtree(scratch, ignore_errors=True)
        if proc.returncode != 0:
            if NO_HITS in proc.stdout + proc.stderr and not out.exists():
                return []
            raise ArgExtenderError(f"Bandage failed (exit {proc.returncode}): {' '.join(cmd)}")
        if not out.is_file():
            raise ArgExtenderError(f"Bandage wrote no graph: {' '.join(cmd)}")
        return [line.split("\t")[1] for line in out.read_text().splitlines() if line.startswith("S\t")]


def graph_targets(
    assembly: Assembly | None, selected: Path, seeds: Path, bandage: Bandage, out: Path
) -> tuple[str, str | None]:
    """Write the next targets from the strain graph to `out`; returns (source, fallback reason)."""
    if assembly is None or assembly.mode is None:
        return _fallback(selected, out, "no assembly in this round")
    graph_file = assembly.strain_graph
    if not graph_file.is_file() or graph_file.stat().st_size == 0:
        return _fallback(selected, out, f"strain_graph.gfa missing (SPAdes mode {assembly.mode})")
    try:
        graph = parse_gfa(graph_file)
    except GraphError as exc:
        return _fallback(selected, out, f"unusable graph: {exc}")
    names = bandage.anchors(graph_file, seeds, out.parent)
    unknown = sorted(set(names) - set(graph.segments))
    if unknown:
        raise ArgExtenderError(f"Bandage named segments not in the graph: {', '.join(unknown[:5])}")
    if not names:
        return _fallback(selected, out, "no recruitment anchor in the graph")
    chosen = select_components(graph, names)
    out.write_text("".join(f">seg_{s}_{role}\n{graph.segments[s]}\n" for s, role in chosen))
    return "strain-graph", None


def _fallback(selected: Path, out: Path, reason: str) -> tuple[str, str]:
    shutil.copyfile(selected, out)
    return "scaffolds (fallback)", reason
