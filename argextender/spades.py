"""SPAdes 3.15.5: `spades.py --only-assembler -k 21,33,55,87 -o OUT -t T` plus the input-dependent mode.

Paired reads -> `--meta -1 R1 -2 R2 [-s S]`; only singletons -> `--sc -s S`; no reads -> not run.
After the run only the files in KEPT are kept, to save space.
"""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from argextender import ArgExtenderError

FLAGS = ("--only-assembler", "-k", "21,33,55,87")
KEPT = ("scaffolds.fasta", "contigs.fasta", "assembly_graph_with_scaffolds.gfa", "strain_graph.gfa", "spades.log")


def nonempty(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0


@dataclass(frozen=True)
class Assembly:
    directory: Path
    mode: str | None  # "--meta" (paired), "--sc" (singletons only) or None (no reads: not run)

    @property
    def scaffolds(self) -> Path:
        return self.directory / "scaffolds.fasta"

    @property
    def graph(self) -> Path:
        """The final assembly graph (the run's .gfa output)."""
        return self.directory / "assembly_graph_with_scaffolds.gfa"

    @property
    def strain_graph(self) -> Path:
        """The graph B-graph takes their next targets from."""
        return self.directory / "strain_graph.gfa"


@dataclass
class Spades:
    exe: str
    threads: int = 1
    log: Path | None = None

    def command(self, r1: Path, r2: Path, singletons: Path, out: Path) -> tuple[list[str] | None, str | None]:
        cmd = [self.exe, *FLAGS, "-o", str(out), "-t", str(self.threads)]
        has_s, has_pair = nonempty(singletons), nonempty(r1)
        if has_pair and has_s:
            return [*cmd, "--meta", "-1", str(r1), "-2", str(r2), "-s", str(singletons)], "--meta"
        if has_pair:
            return [*cmd, "--meta", "-1", str(r1), "-2", str(r2)], "--meta"
        if has_s:
            return [*cmd, "--sc", "-s", str(singletons)], "--sc"
        return None, None

    def assemble(self, r1: Path, r2: Path, singletons: Path, out: Path) -> Assembly:
        if out.exists():
            shutil.rmtree(out)
        cmd, mode = self.command(r1, r2, singletons, out)
        if cmd is None:
            return Assembly(out, None)
        with self.log.open("a") if self.log else open("/dev/null", "w") as log:
            log.write("Running cmd:\t" + " ".join(cmd) + "\n")
            log.flush()
            proc = subprocess.run(cmd, cwd=out.parent, stdout=log, stderr=subprocess.STDOUT, check=False)
        if proc.returncode != 0:
            raise ArgExtenderError(f"SPAdes failed (exit {proc.returncode}): {' '.join(cmd)}")
        for child in list(out.iterdir()) if out.is_dir() else []:
            if child.name not in KEPT:
                shutil.rmtree(child) if child.is_dir() else child.unlink()
        return Assembly(out, mode)
