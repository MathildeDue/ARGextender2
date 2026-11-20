"""Input reads: validation, internal read IDs for KMA, and writing the pool for SPAdes.

Reads are tracked by their position in the input files. Each record gets an internal ID `x<index>_<mate>`
(mate `1`, `2` or `s` for singletons) in a renamed copy of the inputs that KMA maps. KMA's output names then
identify the exact record and its mate, whatever the original naming style (`/1 /2`, Illumina comments,
SRA). Everything written for SPAdes uses the original records, unchanged.

Only paths and counts are kept in memory. Records are read from the input files whenever they are needed,
by record position, so mates stay paired; a pass stops after the highest requested record.
"""

from __future__ import annotations

import gzip
import re
from collections.abc import Collection, Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import IO

from argextender import ArgExtenderError

INTERNAL_RE = re.compile(r"^x(\d+)_([12s])$")


@dataclass(frozen=True)
class FastqRecord:
    header: str  # the header line without the leading "@" and the newline
    seq: str
    plus: str  # the separator line without the newline
    qual: str

    @property
    def name(self) -> str:
        """First whitespace-separated token of the header."""
        return self.header.split(maxsplit=1)[0] if self.header.strip() else ""

    def text(self) -> str:
        return f"@{self.header}\n{self.seq}\n{self.plus}\n{self.qual}\n"


def open_text(path: Path) -> IO[str]:
    with path.open("rb") as probe:
        gz = probe.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if gz else path.open()


def read_fastq(path: Path) -> Iterator[FastqRecord]:
    """Records of a FASTQ file (plain or gzip); ArgExtenderError on malformed records."""
    with open_text(path) as handle:
        n = 0
        while True:
            header = handle.readline()
            if not header:
                return
            if header.strip() == "":
                continue
            seq, plus, qual = handle.readline(), handle.readline(), handle.readline()
            n += 1
            if not header.startswith("@") or not plus.startswith("+") or not qual:
                raise ArgExtenderError(f"{path}: malformed FASTQ record {n}")
            seq, qual = seq.rstrip("\n"), qual.rstrip("\n")
            if len(seq) != len(qual):
                raise ArgExtenderError(f"{path}: record {n} has {len(seq)} bases but {len(qual)} quality values")
            yield FastqRecord(header[1:].rstrip("\n"), seq, plus.rstrip("\n"), qual)


def normalise_name(name: str) -> str:
    """Mate-independent read name: the first token without a trailing /1 or /2."""
    return re.sub(r"/[12]$", "", name)


def internal_id(index: int, mate: str) -> str:
    return f"x{index}_{mate}"


def parse_internal_id(name: str) -> tuple[int, str]:
    m = INTERNAL_RE.match(name)
    if m is None:
        raise ArgExtenderError(f"unexpected read identifier in KMA output: '{name}'")
    return int(m.group(1)), m.group(2)


class Reads:
    """The input reads: R1/R2 pairs and optional singletons."""

    def __init__(self, r1: Path, r2: Path, s: Path | None, pairs: int, singletons: int) -> None:
        self.r1, self.r2 = r1, r2
        self.s = s if s is not None and s.stat().st_size > 0 else None
        self.pairs, self.singletons = pairs, singletons

    @classmethod
    def prepare(cls, r1: Path, r2: Path, s: Path | None, copies: Path) -> tuple[Reads, dict[str, Path]]:
        """Validate R1/R2 and write the internal copies for KMA, in one pass.

        R1 and R2 must have the same number of records, and the i-th records must carry the same normalised
        name. An empty or missing singletons file is allowed. Returns the reads and the internal copies by
        mate ("1", "2" and, if there are singletons, "s").
        """
        copies.mkdir(parents=True, exist_ok=True)
        out = {"1": copies / "reads_internal_1.fastq", "2": copies / "reads_internal_2.fastq"}
        n = 0
        with out["1"].open("w") as h1, out["2"].open("w") as h2:
            it1, it2 = read_fastq(r1), read_fastq(r2)
            while True:
                x, y = next(it1, None), next(it2, None)
                if x is None or y is None:
                    if x is not None or y is not None:
                        a = n + (x is not None) + sum(1 for _ in it1)
                        b = n + (y is not None) + sum(1 for _ in it2)
                        raise ArgExtenderError(f"R1 has {a} records but R2 has {b} ({r1}, {r2})")
                    break
                if normalise_name(x.name) != normalise_name(y.name):
                    raise ArgExtenderError(f"R1 and R2 are out of sync at record {n + 1}: '{x.name}' vs '{y.name}'")
                h1.write(f"@{internal_id(n, '1')}\n{x.seq}\n+\n{x.qual}\n")
                h2.write(f"@{internal_id(n, '2')}\n{y.seq}\n+\n{y.qual}\n")
                n += 1
        singles = 0
        if s is not None and s.stat().st_size > 0:
            out["s"] = copies / "reads_internal_s.fastq"
            with out["s"].open("w") as hs:
                for r in read_fastq(s):
                    hs.write(f"@{internal_id(singles, 's')}\n{r.seq}\n+\n{r.qual}\n")
                    singles += 1
            if singles == 0:  # no copy for a singletons file without records
                out.pop("s").unlink()
        return cls(r1, r2, s, n, singles), out

    def _scan(self, wanted_pairs: set[int], wanted_singles: set[int]) -> Iterator[tuple[str, int, FastqRecord]]:
        """Records of the wanted pairs (mates "1" and "2") and singletons, in record order.

        The files were validated by `prepare`, so records that are not wanted are only stepped over.
        """
        if wanted_pairs:
            last = max(wanted_pairs)
            with open_text(self.r1) as h1, open_text(self.r2) as h2:
                for i, (x, y) in enumerate(zip(_raw(h1), _raw(h2), strict=False)):
                    if i in wanted_pairs:
                        yield "1", i, _record(x)
                        yield "2", i, _record(y)
                    if i >= last:
                        break
        if wanted_singles and self.s is not None:
            last = max(wanted_singles)
            with open_text(self.s) as hs:
                for i, r in enumerate(_raw(hs)):
                    if i in wanted_singles:
                        yield "s", i, _record(r)
                    if i >= last:
                        break

    def write_pool(
        self, pairs: Collection[int], singles: Collection[int], out_1: Path, out_2: Path, out_s: Path
    ) -> None:
        """The original records of the given pairs (both mates) and singletons, in record order."""
        with out_1.open("w") as h1, out_2.open("w") as h2, out_s.open("w") as hs:
            for m, _, r in self._scan(set(pairs), set(singles)):
                (h1 if m == "1" else h2 if m == "2" else hs).write(r.text())

    def write_templates(self, pairs: Collection[int], singles: Collection[int], handle: IO[str]) -> None:
        """Reads as KMA templates (strategy B): `>x<i>_1`, `>x<i>_2` per pair in index order, then singletons."""
        for m, i, r in self._scan(set(pairs), set(singles)):
            handle.write(f">{internal_id(i, m)}\n{r.seq}\n")


def _raw(handle: IO[str]) -> Iterator[tuple[str, str, str, str]]:
    """The four raw lines of each record, skipping blank lines between records as read_fastq does."""
    readline = handle.readline
    while True:
        header = readline()
        if not header:
            return
        if header.strip() == "":
            continue
        yield header, readline(), readline(), readline()


def _record(lines: tuple[str, str, str, str]) -> FastqRecord:
    header, seq, plus, qual = lines
    return FastqRecord(header[1:].rstrip("\n"), seq.rstrip("\n"), plus.rstrip("\n"), qual.rstrip("\n"))
