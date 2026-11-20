"""ARGextender: extend the genetic context of antimicrobial resistance genes (ARGs) by iterative read
recruitment and targeted assembly.

One round: KMA recruits the reads that map to the current targets (k-mer mapping only, `-sasm`), both mates
of every recruited pair join the pool, and on assembly rounds SPAdes assembles the pool. The scaffolds with
a hit to a seed are selected, the ARGs are placed on them (a separate KMA `-sam` call) to measure the flanks,
and the next targets are built. The pool is reset after every assembly.

Strategies (one loop; A is B with an assembly after every round):
  A-scaffold  (default) next targets = the selected scaffolds; assembly every round
  B-graph     next targets = Bandage anchors + their connected components in strain_graph.gfa; assembly
              only every N rounds (--assemble-every) and after the last round; between assemblies the
              pool's reads are added to the targets

Stopping, checked after each assembly: --rounds (hard limit), --min-growth (the flanks around the contigs' best
ARGs grew by less than X % per round since the previous assembly), --max-flank (at least --max-flank-share % of the
contigs, default 50 = the median, have both flanks of their best ARG >= L bp), and empty targets.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import math
import os
import platform
import resource
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from argextender import ArgExtenderError, __version__
from argextender.flanks import Hit, flank_total, parse_sam, share_reached, write_hits
from argextender.graph import Bandage, graph_targets
from argextender.kma import (
    Kma,
    detected_references,
    fasta_stats,
    frag_query_names,
    match_asm,
    read_fasta,
    read_mapstat,
    select_seeds,
    write_fasta,
)
from argextender.reads import Reads, parse_internal_id
from argextender.spades import Assembly, Spades

STRATEGIES = ("A-scaffold", "B-graph")
DEFAULT_ASSEMBLE_EVERY = 3  # B-graph
VERSION_ARGS = {"kma": ["-v"], "spades.py": ["--version"], "Bandage": ["--version"], "blastn": ["-version"]}
OUTPUT_SUFFIXES = (".fasta", ".gfa", ".frag.gz", ".frag_raw.gz")
ROUND_COLUMNS = (
    "round",
    "assembled",
    "target_records",
    "target_bp",
    "reads_recruited",
    "pool_reads",
    "spades_mode",
    "selected_contigs",
    "selected_bp",
    "arg_contigs",
    "flank_bp",
    "max_flank_reached_pct",
    "flank_growth_pct",
    "next_target_source",
    "fallback_reason",
    "stopped_by",
    "seconds",
)


@dataclass(frozen=True)
class Config:
    r1: Path
    r2: Path
    singletons: Path | None
    mapstat: Path
    panres: Path
    out: Path
    strategy: str = "A-scaffold"
    threads: int = 1
    rounds: int = 30
    min_growth: float = 2.0  # percent per round; 0 disables
    max_flank: int = 0  # bp; 0 disables
    max_flank_share: float = 50.0  # percent of the contigs; 50 = the median
    assemble_every: int = 1
    keep_work: bool = False

    @property
    def graph_targets(self) -> bool:
        return self.strategy == "B-graph"

    @property
    def work(self) -> Path:
        return Path(f"{self.out}.work")


# --- stopping ------------------------------------------------------------------------------------------


def flank_growth(flank: int, previous: int | None, rounds_between: int) -> float | None:
    """Growth per round of the flank total since the previous assembly: (F / F_prev)^(1/n) - 1.

    F = the summed up- and downstream flanks of every selected contig's best ARG. n = rounds since the previous
    assembly (1 under A), so B's cycles give a per-round value.
    """
    if previous is None:
        return None
    if previous == 0:
        return math.inf if flank > 0 else 0.0
    return (flank / previous) ** (1 / rounds_between) - 1


def stop_reasons(cfg: Config, k: int, growth: float | None, reached: float, targets_empty: bool) -> list[str]:
    """Why the run stops after assembly round k (empty: it continues); reached = share of contigs at L."""
    reasons = []
    if k >= cfg.rounds:
        reasons.append("max_rounds")
    if cfg.min_growth > 0 and growth is not None and growth * 100 < cfg.min_growth:
        reasons.append("no_growth")
    if cfg.max_flank > 0 and reached * 100 >= cfg.max_flank_share:
        reasons.append("max_flank")
    if targets_empty:
        reasons.append("empty_target")
    return reasons


# --- tools ---------------------------------------------------------------------------------------------


def find_tool(name: str) -> str:
    path = shutil.which(name)
    if path is None:
        raise ArgExtenderError(f"Required tool '{name}' was not found on PATH")
    return path


def tool_version(name: str, path: str) -> str:
    env = {**os.environ, "QT_QPA_PLATFORM": "offscreen"} if name == "Bandage" else None
    proc = subprocess.run([path, *VERSION_ARGS[name]], capture_output=True, text=True, check=False, env=env)
    text = (proc.stdout or proc.stderr).strip()
    return text.splitlines()[0] if text else ""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1 << 20):
            digest.update(block)
    return digest.hexdigest()


# --- the run -------------------------------------------------------------------------------------------


class Run:
    def __init__(self, cfg: Config) -> None:
        self.cfg = cfg
        self.work = cfg.work.resolve()
        names = ["kma", "spades.py"] + (["Bandage", "blastn"] if cfg.graph_targets else [])
        self.tools = {n: {"path": find_tool(n)} for n in names}
        for n, t in self.tools.items():
            t["version"] = tool_version(n, t["path"])
        if cfg.graph_targets:
            find_tool("makeblastdb")  # Bandage calls it
        log = self.work / "tools.log"
        self.kma = Kma(self.tools["kma"]["path"], cfg.threads, log)
        self.spades = Spades(self.tools["spades.py"]["path"], cfg.threads, log)
        self.bandage = Bandage(self.tools["Bandage"]["path"], log) if cfg.graph_targets else None
        self.rows: list[dict[str, Any]] = []
        self.seed_lengths: dict[str, int] = {}

    def execute(self, argv: list[str]) -> dict[str, Any]:
        cfg = self.cfg
        t0 = time.perf_counter()
        started = dt.datetime.now(dt.UTC).isoformat()
        for path in (cfg.r1, cfg.r2, cfg.mapstat, cfg.panres, *([cfg.singletons] if cfg.singletons else [])):
            if not path.is_file():
                raise ArgExtenderError(f"Input file not found: {path}")
        seeds = select_seeds(detected_references(cfg.mapstat), cfg.panres)
        if self.work.exists():
            raise ArgExtenderError(f"Work directory {self.work} already exists; remove it or choose another --out")
        self.work.mkdir(parents=True)
        cfg.out.parent.mkdir(parents=True, exist_ok=True)
        seeds_fasta = self.work / "seeds.fasta"
        write_fasta(seeds, seeds_fasta)
        self.seed_lengths = {h.split()[0]: len(s) for h, s in seeds}
        reads, internal = Reads.prepare(cfg.r1, cfg.r2, cfg.singletons, self.work / "reads")
        stopped_by, out_round = "no_detected_references", None
        if seeds:
            (self.work / "db").mkdir()
            self.kma.index(seeds_fasta, self.work / "db" / "selection")
            stopped_by, out_round = self.loop(reads, internal, seeds_fasta)
        outputs = self.write_outputs(out_round)
        self.write_rounds()
        summary = self.summary(argv, started, t0, reads, seeds, stopped_by, out_round, outputs)
        Path(f"{cfg.out}.summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        if not cfg.keep_work:
            shutil.rmtree(self.work)
        return summary

    def loop(self, reads: Reads, internal: dict[str, Path], seeds_fasta: Path) -> tuple[str, int | None]:
        """All rounds; returns (stopped_by, the round the outputs come from)."""
        cfg = self.cfg
        targets = seeds_fasta  # the targets built from the last assembly (the seeds before the first)
        pool_pairs: set[int] = set()
        pool_singles: set[int] = set()
        previous_flank: int | None = None  # the flank total of the previous assembly
        last_assembly = 0
        out_round = None
        for k in range(1, cfg.rounds + 1):
            t = time.perf_counter()
            rd = self.work / f"round_{k:03d}"
            rd.mkdir()
            # templates: the targets, plus the pool's reads between assemblies (B-graph; empty under A)
            with (rd / "target.fasta").open("w") as handle:
                handle.write(targets.read_text())
                reads.write_templates(pool_pairs, pool_singles, handle)
            self.kma.index_targets(rd / "target.fasta", rd / "target_db")
            self.kma.recruit([internal[m] for m in ("1", "2", "s") if m in internal], rd / "target_db", rd / "recruit")
            frag = rd / "recruit.frag.gz"
            direct = {parse_internal_id(n) for n in frag_query_names(frag)} if frag.is_file() else set()
            pool_pairs |= {i for i, m in direct if m in ("1", "2")}  # both mates of every recruited pair
            pool_singles |= {i for i, m in direct if m == "s"}
            n_t, bp_t = fasta_stats(rd / "target.fasta")
            row: dict[str, Any] = {c: "" for c in ROUND_COLUMNS}
            row.update(
                round=k,
                assembled=0,
                target_records=n_t,
                target_bp=bp_t,
                reads_recruited=len(direct),
                pool_reads=2 * len(pool_pairs) + len(pool_singles),
            )
            stop: list[str] = []
            if k % cfg.assemble_every == 0 or k == cfg.rounds:
                assembly, hits, source, reason = self.assemble(rd, reads, pool_pairs, pool_singles, seeds_fasta)
                flank = flank_total(hits)
                growth = flank_growth(flank, previous_flank, k - last_assembly)
                targets = rd / "next_targets.fasta"
                reached = share_reached(hits, cfg.max_flank)
                stop = stop_reasons(cfg, k, growth, reached, fasta_stats(targets)[0] == 0)
                n_sel, bp_sel = fasta_stats(rd / "selected.fasta")
                row.update(
                    assembled=int(assembly is not None and assembly.mode is not None),
                    spades_mode=(assembly.mode or "") if assembly else "",
                    selected_contigs=n_sel,
                    selected_bp=bp_sel,
                    arg_contigs=sum(h.best_for_contig for h in hits),
                    flank_bp=flank,
                    max_flank_reached_pct=round(reached * 100, 1) if cfg.max_flank else "",
                    flank_growth_pct="" if growth is None else round(growth * 100, 3),
                    next_target_source=source,
                    fallback_reason=reason or "",
                )
                if assembly is not None and assembly.graph.is_file() and (rd / "select.kma.frag.gz").is_file():
                    out_round = k  # an assembly graph and a selection run: the outputs can come from here
                previous_flank, last_assembly = flank, k
                pool_pairs, pool_singles = set(), set()
            row.update(stopped_by="+".join(stop), seconds=round(time.perf_counter() - t, 1))
            self.rows.append(row)
            if stop:
                return row["stopped_by"], out_round
        return "max_rounds", out_round  # not reached: the last round is always an assembly round

    def assemble(
        self, rd: Path, reads: Reads, pairs: set[int], singles: set[int], seeds_fasta: Path
    ) -> tuple[Assembly | None, list[Hit], str, str | None]:
        """Assemble the pool, select the scaffolds with a seed hit, place the ARGs, build the next targets."""
        selected, next_targets = rd / "selected.fasta", rd / "next_targets.fasta"
        assembly = None
        selected.write_text("")
        if pairs or singles:
            p1, p2, ps = rd / "pool_1.fq", rd / "pool_2.fq", rd / "pool_s.fq"
            reads.write_pool(pairs, singles, p1, p2, ps)
            assembly = self.spades.assemble(p1, p2, ps, rd / "spades")
            scaffolds = assembly.scaffolds
            if scaffolds.is_file() and scaffolds.stat().st_size > 0:
                self.kma.select(scaffolds, self.work / "db" / "selection", rd / "select.kma")
                write_fasta(match_asm(rd / "select.kma.frag.gz"), selected)
        hits = self.place(selected, rd / "place.kma")
        write_hits(hits, rd / "arg_hits.tsv")
        if self.bandage is None:
            shutil.copyfile(selected, next_targets)
            return assembly, hits, "scaffolds", None
        source, reason = graph_targets(assembly, selected, seeds_fasta, self.bandage, next_targets)
        return assembly, hits, source, reason

    def place(self, contigs: Path, out: Path) -> list[Hit]:
        """The seeds' positions and flanks on the contigs (KMA -sam)."""
        if contigs.stat().st_size == 0:
            return []
        sam = self.kma.place(contigs, self.work / "db" / "selection", out)
        lengths = {h.split()[0]: len(s) for h, s in read_fasta(contigs)}
        return parse_sam(sam, lengths, self.seed_lengths)

    # --- outputs -----------------------------------------------------------------------------------------

    def write_outputs(self, out_round: int | None) -> dict[str, Any]:
        """The outputs, all from one round: the last with an assembly graph and a selection run."""
        cfg = self.cfg
        dest = {s: Path(f"{cfg.out}{s}") for s in (*OUTPUT_SUFFIXES, ".arg_hits.tsv")}
        if out_round is None:
            for path in dest.values():
                path.write_bytes(b"")
            write_hits([], dest[".arg_hits.tsv"])
            return {s: str(p) for s, p in dest.items()}
        rd = self.work / f"round_{out_round:03d}"
        shutil.copyfile(rd / "selected.fasta", dest[".fasta"])
        shutil.copyfile(rd / "spades" / "assembly_graph_with_scaffolds.gfa", dest[".gfa"])
        shutil.copyfile(rd / "arg_hits.tsv", dest[".arg_hits.tsv"])
        frag, raw = rd / "select.kma.frag.gz", rd / "select.kma.frag_raw.gz"
        if dest[".fasta"].stat().st_size > 0:  # the final contigs against the complete PanRes
            ann = self.work / "annotation"
            ann.mkdir()
            self.kma.index(cfg.panres.resolve(), ann / "db")
            self.kma.select(dest[".fasta"], ann / "db", ann / "final")
            frag, raw = ann / "final.frag.gz", ann / "final.frag_raw.gz"
        shutil.copyfile(frag, dest[".frag.gz"])
        if raw.is_file():
            shutil.copyfile(raw, dest[".frag_raw.gz"])
        else:
            dest[".frag_raw.gz"].write_bytes(b"")
        return {s: str(p) for s, p in dest.items()}

    def write_rounds(self) -> None:
        with Path(f"{self.cfg.out}.rounds.tsv").open("w") as handle:
            handle.write("\t".join(ROUND_COLUMNS) + "\n")
            handle.writelines("\t".join(str(row[c]) for c in ROUND_COLUMNS) + "\n" for row in self.rows)

    def summary(
        self,
        argv: list[str],
        started: str,
        t0: float,
        reads: Reads,
        seeds: list[tuple[str, str]],
        stopped_by: str,
        out_round: int | None,
        outputs: dict[str, Any],
    ) -> dict[str, Any]:
        cfg = self.cfg
        own = resource.getrusage(resource.RUSAGE_SELF)
        tools = resource.getrusage(resource.RUSAGE_CHILDREN)
        files = {"r1": cfg.r1, "r2": cfg.r2, "mapstat": cfg.mapstat, "panres": cfg.panres}
        if cfg.singletons:
            files["singletons"] = cfg.singletons
        return {
            "argextender_version": __version__,
            "command": argv,
            "started_utc": started,
            "parameters": {
                "strategy": cfg.strategy,
                "rounds": cfg.rounds,
                "min_growth_pct": cfg.min_growth,
                "max_flank_bp": cfg.max_flank,
                "max_flank_share_pct": cfg.max_flank_share,
                "assemble_every": cfg.assemble_every,
                "threads": cfg.threads,
            },
            "tools": self.tools,
            "inputs": {
                "files": {k: {"path": str(p.resolve()), "bytes": p.stat().st_size} for k, p in files.items()},
                "mapstat_sha256": sha256_file(cfg.mapstat),
                "pairs": reads.pairs,
                "singletons": reads.singletons,
                "seeds": len(seeds),
                "seed_source": read_mapstat(cfg.mapstat)[0],
            },
            "rounds": len(self.rows),
            "assembly_rounds": [r["round"] for r in self.rows if r["spades_mode"]],
            "stopped_by": stopped_by,
            "outputs_from_round": out_round,
            "fallbacks": [
                {"round": r["round"], "reason": r["fallback_reason"]} for r in self.rows if r["fallback_reason"]
            ],
            "outputs": outputs,
            "resources": {
                "wall_seconds": round(time.perf_counter() - t0, 1),
                "cpu_seconds": round(own.ru_utime + own.ru_stime + tools.ru_utime + tools.ru_stime, 1),
                "cpu_seconds_tools": round(tools.ru_utime + tools.ru_stime, 1),
                "peak_rss_kb": max(own.ru_maxrss, tools.ru_maxrss),
            },
            "host": {"node": platform.node(), "python": platform.python_version()},
        }


# --- command line --------------------------------------------------------------------------------------


def parse_args(argv: list[str]) -> Config:
    p = argparse.ArgumentParser(
        prog="argextender", description=__doc__.split("\n\n")[0], formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p.add_argument("--r1", type=Path, required=True, help="forward reads (FASTQ, optionally gzipped)")
    p.add_argument("--r2", type=Path, required=True, help="reverse reads")
    p.add_argument("--singletons", type=Path, help="unpaired reads")
    p.add_argument("--mapstat", type=Path, required=True, help="ARGprofiler's unfiltered KMA <sample>.mapstat")
    p.add_argument("--panres", type=Path, required=True, help="PanRes gene FASTA (the seeds' sequences)")
    p.add_argument("--strategy", choices=STRATEGIES, default="A-scaffold")
    p.add_argument("--out", type=Path, required=True, help="output prefix")
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--rounds", type=int, default=30, help="maximum number of recruitment rounds")
    p.add_argument(
        "--min-growth",
        type=float,
        default=2.0,
        help="stop when the flanks around the contigs' best ARGs grow by less than this percentage per round since "
        "the previous assembly (0: off)",
    )
    p.add_argument(
        "--max-flank",
        type=int,
        default=0,
        help="stop when --max-flank-share %% of the contigs have both flanks of their best ARG >= this (bp; 0: off)",
    )
    p.add_argument(
        "--max-flank-share",
        type=float,
        default=50.0,
        help="percentage of the contigs that must reach --max-flank (50: the median; 100: all)",
    )
    p.add_argument(
        "--assemble-every", type=int, help=f"B-graph: assemble every N rounds (default {DEFAULT_ASSEMBLE_EVERY})"
    )
    p.add_argument("--keep-work", action="store_true", help="keep PREFIX.work/ (the round folders)")
    p.add_argument("--version", action="version", version=f"argextender {__version__}")
    a = p.parse_args(argv)
    if a.assemble_every is not None and a.strategy != "B-graph":
        p.error("--assemble-every applies to B-graph only")
    every = (a.assemble_every or DEFAULT_ASSEMBLE_EVERY) if a.strategy == "B-graph" else 1
    if a.threads < 1 or a.rounds < 1 or every < 1:
        p.error("--threads, --rounds and --assemble-every must be at least 1")
    if a.min_growth < 0 or a.max_flank < 0:
        p.error("--min-growth and --max-flank must be >= 0")
    if not 0 < a.max_flank_share <= 100:
        p.error("--max-flank-share must be in (0, 100]")
    return Config(
        a.r1,
        a.r2,
        a.singletons,
        a.mapstat,
        a.panres,
        a.out,
        a.strategy,
        a.threads,
        a.rounds,
        a.min_growth,
        a.max_flank,
        a.max_flank_share,
        every,
        a.keep_work,
    )


def main(argv: list[str] | None = None) -> None:
    argv = sys.argv[1:] if argv is None else argv
    cfg = parse_args(argv)
    try:
        summary = Run(cfg).execute(["argextender", *argv])
    except ArgExtenderError as exc:
        sys.exit(f"argextender: error: {exc}")
    print(f"argextender: {summary['rounds']} rounds, stopped by {summary['stopped_by']}; outputs {cfg.out}.*")


if __name__ == "__main__":
    main()
