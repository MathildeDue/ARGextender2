"""KMA: seeds from ARGprofiler's mapstat, indexing, read recruitment and contig selection.

Seeds. ARGprofiler maps the reads against the complete PanRes with KMA; every PanRes gene that receives
reads has a row in the unfiltered `<sample>.mapstat`. All of them are seeds: their PanRes sequences are the
first round's targets and the selection database. The filtered `<sample>.mapstat.filtered` is refused.

Commands (flags fixed; the benchmarked settings):
  selection / annotation index  kma index -i FASTA -o DB -m 14
  round target index            kma index -i FASTA -o DB -k 25
  recruitment                   kma -i READS... -o OUT -tmp -t_db DB -mem_mode -1t1 -mrc 0.25 -ID 50 -na -nc -ef
                                    -sasm -t T
  selection                     kma -i CONTIGS -o OUT -tmp -t_db DB -ID 50 -na -nc -a -proxi -0.7 -t T
  ARG placement                 the selection command plus -sam (a separate call: -sam changes what the
                                selection reports, so it is never added to the selection call itself)
"""

from __future__ import annotations

import gzip
import subprocess
from dataclasses import dataclass
from pathlib import Path

from argextender import ArgExtenderError

INDEX_FLAGS = ("-m", "14")
TARGET_INDEX_FLAGS = ("-k", "25")
RECRUIT_FLAGS = ("-mem_mode", "-1t1", "-mrc", "0.25", "-ID", "50", "-na", "-nc", "-ef", "-sasm")
SELECT_FLAGS = ("-ID", "50", "-na", "-nc", "-a", "-proxi", "-0.7")

# columns that ARGprofiler's mapstatFilters.R appends to its filtered output
FILTER_COLUMNS = (
    "meanDepthCovered",
    "spuriosCovRatio",
    "FragReadRatio",
    "meanMapScore",
    "readConsensRefIdentity",
    "readRefIdentity",
    "depthCV",
    "propCovered",
)


@dataclass
class Kma:
    """Runs KMA; every command is appended to `log` with the tool's output."""

    exe: str
    threads: int = 1
    log: Path | None = None

    def _run(self, args: list[str], cwd: Path, stdout: Path | None = None) -> None:
        cmd = [self.exe, *args]
        with self.log.open("a") if self.log else open("/dev/null", "w") as log:
            log.write("Running cmd:\t" + " ".join(cmd) + "\n")
            log.flush()
            if stdout is None:
                proc = subprocess.run(cmd, cwd=cwd, stdout=log, stderr=subprocess.STDOUT, check=False)
            else:
                with stdout.open("w") as out:
                    proc = subprocess.run(cmd, cwd=cwd, stdout=out, stderr=log, check=False)
        if proc.returncode != 0:
            raise ArgExtenderError(f"KMA failed (exit {proc.returncode}): {' '.join(cmd)}")

    def index(self, fasta: Path, db: Path, flags: tuple[str, ...] = INDEX_FLAGS) -> None:
        self._run(["index", "-i", str(fasta), "-o", str(db), *flags], db.parent)

    def index_targets(self, fasta: Path, db: Path) -> None:
        self.index(fasta, db, TARGET_INDEX_FLAGS)

    def recruit(self, reads: list[Path], db: Path, out: Path) -> None:
        args = ["-i", *map(str, reads), "-o", str(out), "-tmp", "-t_db", str(db), *RECRUIT_FLAGS]
        self._run([*args, "-t", str(self.threads)], out.parent)

    def select(self, contigs: Path, db: Path, out: Path) -> None:
        """Writes OUT.frag.gz, OUT.frag_raw.gz and OUT.res."""
        args = ["-i", str(contigs), "-o", str(out), "-tmp", "-t_db", str(db), *SELECT_FLAGS]
        self._run([*args, "-t", str(self.threads)], out.parent)

    def place(self, contigs: Path, db: Path, out: Path) -> Path:
        """The selection command plus -sam; the SAM (KMA's stdout) goes to OUT.sam."""
        sam = Path(f"{out}.sam")
        args = ["-i", str(contigs), "-o", str(out), "-tmp", "-t_db", str(db), *SELECT_FLAGS]
        self._run([*args, "-t", str(self.threads), "-sam"], out.parent, stdout=sam)
        return sam


def frag_query_names(frag: Path) -> set[str]:
    """First token of column 7 (the query name) of every `.frag.gz` line."""
    names = set()
    with gzip.open(frag, "rt") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 7 and fields[6]:
                names.add(fields[6].split(" ")[0])
    return names


def match_asm(frag: Path) -> list[tuple[str, str]]:
    """Selected contigs from a selection `.frag.gz`: (name, sequence) from columns 7 and 1, the first
    occurrence of each name, in file order."""
    seen: dict[str, str] = {}
    with gzip.open(frag, "rt") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 7 and fields[6] not in seen:
                seen[fields[6]] = fields[0]
    return list(seen.items())


def write_fasta(records: list[tuple[str, str]], path: Path) -> None:
    path.write_text("".join(f">{name}\n{seq}\n" for name, seq in records))


def read_fasta(path: Path) -> list[tuple[str, str]]:
    """(header line without '>', sequence) of every record, sequences joined, in file order."""
    records: list[tuple[str, list[str]]] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            records.append((line[1:].strip(), []))
        elif records and line.strip():
            records[-1][1].append(line.strip())
    return [(h, "".join(s)) for h, s in records]


def fasta_stats(path: Path) -> tuple[int, int]:
    """Number of records and total bases of a FASTA file; (0, 0) if it is missing."""
    if not path.is_file():
        return 0, 0
    n = bp = 0
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                n += 1
            else:
                bp += len(line.strip())
    return n, bp


def read_mapstat(mapstat: Path) -> tuple[dict[str, str], list[dict[str, str]]]:
    """KMA's `##` header fields and the data rows (keyed by the column header) of a mapstat."""
    meta: dict[str, str] = {}
    columns: list[str] = []
    rows: list[dict[str, str]] = []
    for line in mapstat.read_text().splitlines():
        if line.startswith("##"):
            key, _, value = line[2:].strip().partition("\t")
            meta[key.strip()] = value.strip()
        elif line.startswith("#"):
            columns = [c.strip() for c in line[1:].split("\t")]
        elif line.strip():
            rows.append(dict(zip(columns or ["refSequence"], line.split("\t"), strict=False)))
    filtered = [c for c in FILTER_COLUMNS if c in columns]
    if filtered:
        raise ArgExtenderError(
            f"{mapstat} is a filtered mapstat (it has the columns {', '.join(filtered)} added by mapstatFilters.R); "
            "give the unfiltered KMA <sample>.mapstat"
        )
    return meta, rows


def detected_references(mapstat: Path) -> list[str]:
    """Reference names of the data rows of an unfiltered mapstat, in file order, without duplicates."""
    names: list[str] = []
    for row in read_mapstat(mapstat)[1]:
        name = row.get("refSequence", "").strip()
        if name and name not in names:
            names.append(name)
    return names


def select_seeds(names: list[str], panres: Path) -> list[tuple[str, str]]:
    """The PanRes records of the detected references, in detection order.

    A name matches a record whose full header equals it (KMA template names are the full header line),
    or failing that whose first header token does.
    """
    records = read_fasta(panres)
    by_header = {h: (h, s) for h, s in records}
    by_token: dict[str, tuple[str, str]] = {}
    for h, s in records:
        by_token.setdefault(h.split()[0], (h, s))
    missing = [n for n in names if n not in by_header and n not in by_token]
    if missing:
        shown = ", ".join(missing[:10]) + (" ..." if len(missing) > 10 else "")
        raise ArgExtenderError(f"{len(missing)} detected reference(s) not found in {panres}: {shown}")
    return [by_header.get(n) or by_token[n] for n in names]
