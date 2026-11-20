"""ARG placement on the selected contigs and flank lengths, from KMA's SAM output.

After each assembly the selected contigs are mapped against the seeds (`Kma.place`: the selection flags plus
`-sam`, a separate call). In each mapped SAM record the reference is the gene and the CIGAR is in gene
orientation, so on either strand a leading soft clip is the upstream flank and a trailing one the
downstream flank. A hit is kept if identity (`=` / (`=` + `X` + `I` + `D`)) >= 0.90 and it covers at least
min(100, gene length) gene bases (`=` + `X` + `D`).

Per contig, the best ARG is the hit with the highest identity x seed coverage, then coverage, then seed
name. The flank total of a round is the sum, over contigs, of upstream + downstream of the contig's best ARG.

Classes: `complete` (seed coverage >= 0.90); `partial-truncated` (the contig ends inside the gene);
`partial-internal` (the rest).
"""

from __future__ import annotations

import re
from dataclasses import dataclass, replace
from pathlib import Path

MIN_IDENTITY = 0.90
MIN_GENE_BP = 100
COMPLETE_COVERAGE = 0.90
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
COLUMNS = (
    "contig",
    "contig_length",
    "seed",
    "seed_length",
    "identity",
    "seed_coverage",
    "seed_start",
    "seed_end",
    "contig_start",
    "contig_end",
    "strand",
    "upstream_flank_bp",
    "downstream_flank_bp",
    "class",
    "best_for_contig",
)


@dataclass(frozen=True)
class Hit:
    contig: str
    contig_length: int
    seed: str
    seed_length: int
    identity: float
    seed_coverage: float
    seed_start: int  # 1-based, inclusive, gene coordinates
    seed_end: int
    contig_start: int  # 1-based, inclusive, forward contig coordinates
    contig_end: int
    strand: str
    upstream_flank_bp: int  # in gene orientation
    downstream_flank_bp: int
    hit_class: str
    best_for_contig: bool = False

    @property
    def flank_bp(self) -> int:
        return self.upstream_flank_bp + self.downstream_flank_bp


def classify(coverage: float, seed_start: int, seed_end: int, seed_length: int, up: int, down: int) -> str:
    if coverage >= COMPLETE_COVERAGE:
        return "complete"
    truncated = (seed_start > 1 and up == 0) or (seed_end < seed_length and down == 0)
    return "partial-truncated" if truncated else "partial-internal"


def parse_sam(sam: Path, contig_lengths: dict[str, int], seed_lengths: dict[str, int]) -> list[Hit]:
    """Kept hits of a KMA SAM, with each contig's best hit marked; one hit per distinct placement."""
    hits: dict[tuple[str, str, int, int], Hit] = {}
    for line in sam.read_text().splitlines():
        if not line or line.startswith("@"):
            continue
        f = line.split("\t")
        flag = int(f[1])
        if flag & 4 or f[5] == "*":
            continue
        contig, seed = f[0].split()[0], f[2].split()[0]
        ops = [(int(n), op) for n, op in CIGAR_RE.findall(f[5])]
        up = ops[0][0] if ops[0][1] in "SH" else 0
        down = ops[-1][0] if ops[-1][1] in "SH" else 0
        matches = sum(n for n, op in ops if op in "=M")
        aligned = sum(n for n, op in ops if op in "=XMID")
        gene_bp = sum(n for n, op in ops if op in "=XMD")
        identity = matches / aligned if aligned else 0.0
        seed_length = seed_lengths.get(seed, gene_bp)
        if identity < MIN_IDENTITY or gene_bp < min(MIN_GENE_BP, seed_length):
            continue
        n = contig_lengths[contig]
        reverse = bool(flag & 16)
        start = (down if reverse else up) + 1  # forward contig coordinates
        end = n - (up if reverse else down)
        seed_start = int(f[3])
        seed_end = seed_start + gene_bp - 1
        coverage = gene_bp / seed_length
        hit = Hit(
            contig,
            n,
            seed,
            seed_length,
            round(identity, 4),
            round(coverage, 4),
            seed_start,
            seed_end,
            start,
            end,
            "-" if reverse else "+",
            up,
            down,
            classify(coverage, seed_start, seed_end, seed_length, up, down),
        )
        hits.setdefault((contig, seed, start, end), hit)
    best: dict[str, tuple[str, str, int, int]] = {}  # contig -> key of its best hit
    for key, h in hits.items():
        if h.contig not in best or rank(h) > rank(hits[best[h.contig]]):
            best[h.contig] = key
    chosen = set(best.values())
    return [replace(h, best_for_contig=True) if key in chosen else h for key, h in hits.items()]


def rank(h: Hit) -> tuple[float, float, str]:
    """A contig's best ARG: highest identity x seed coverage, then coverage, then seed name."""
    return (h.identity * h.seed_coverage, h.seed_coverage, h.seed)


def best_hits(hits: list[Hit]) -> list[Hit]:
    return [h for h in hits if h.best_for_contig]


def flank_total(hits: list[Hit]) -> int:
    """Sum over contigs of upstream + downstream flank of the contig's best ARG."""
    return sum(h.flank_bp for h in best_hits(hits))


def share_reached(hits: list[Hit], length: int) -> float:
    """Share of the contigs whose best ARG has both flanks >= length (0 without hits)."""
    sides = [min(h.upstream_flank_bp, h.downstream_flank_bp) for h in best_hits(hits)]
    return sum(s >= length for s in sides) / len(sides) if sides else 0.0


def write_hits(hits: list[Hit], path: Path) -> None:
    def cell(h: Hit, column: str) -> str:
        v = getattr(h, "hit_class" if column == "class" else column)
        return str(int(v)) if isinstance(v, bool) else str(v)

    with path.open("w") as handle:
        handle.write("\t".join(COLUMNS) + "\n")
        for h in hits:
            handle.write("\t".join(cell(h, c) for c in COLUMNS) + "\n")
