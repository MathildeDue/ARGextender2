"""End-to-end runs on the synthetic test data with the real tools (KMA, SPAdes, Bandage + BLAST+)."""

import csv
import json
import shutil
from pathlib import Path

import pytest

from argextender.argextender import main

DATA = Path(__file__).parent / "data"
TOOLS = ("kma", "spades.py", "Bandage", "blastn", "makeblastdb")
pytestmark = [
    pytest.mark.tools,
    pytest.mark.skipif(any(shutil.which(t) is None for t in TOOLS), reason="KMA, SPAdes or Bandage not on PATH"),
]


def run(out: Path, strategy: str, *extra: str, mapstat: Path = DATA / "test.mapstat") -> dict:
    main(
        [
            "--r1",
            str(DATA / "reads_1.fastq"),
            "--r2",
            str(DATA / "reads_2.fastq"),
            "--singletons",
            str(DATA / "reads_singleton.fastq"),
            "--mapstat",
            str(mapstat),
            "--panres",
            str(DATA / "panres_subset20.fasta"),
            "--strategy",
            strategy,
            "--out",
            str(out),
            "--threads",
            "2",
            *extra,
        ]
    )
    return json.loads(Path(f"{out}.summary.json").read_text())


def rounds(out: Path) -> list[dict[str, str]]:
    with Path(f"{out}.rounds.tsv").open() as h:
        return list(csv.DictReader(h, delimiter="\t"))


@pytest.mark.parametrize("strategy", ["A-scaffold", "B-graph"])
def test_fixed_rounds(tmp_path: Path, strategy: str) -> None:
    out = tmp_path / "fx"
    summary = run(out, strategy, "--rounds", "6", "--min-growth", "0")
    assert summary["stopped_by"] == "max_rounds" and summary["rounds"] == 6
    expected = [1, 2, 3, 4, 5, 6] if strategy == "A-scaffold" else [3, 6]
    assert summary["assembly_rounds"] == expected
    for suffix in (".fasta", ".gfa", ".frag.gz", ".arg_hits.tsv"):
        assert Path(f"{out}{suffix}").stat().st_size > 0
    assert not Path(f"{out}.work").exists()  # removed without --keep-work
    with Path(f"{out}.arg_hits.tsv").open() as h:
        best = [r for r in csv.DictReader(h, delimiter="\t") if r["best_for_contig"] == "1"]
    assert [(r["seed"], r["class"]) for r in best] == [("pan_27", "complete")]
    flanks = [int(r["flank_bp"]) for r in rounds(out) if r["flank_bp"]]
    assert flanks == sorted(flanks) and flanks[-1] > flanks[0]  # measured after every assembly; the flanks grow
    assert flanks[-1] == sum(int(r["upstream_flank_bp"]) + int(r["downstream_flank_bp"]) for r in best)
    if strategy == "B-graph":
        assert {r["next_target_source"] for r in rounds(out) if r["assembled"] == "1"} == {"strain-graph"}


def test_no_growth_stops_early(tmp_path: Path) -> None:
    out = tmp_path / "fx"
    summary = run(out, "A-scaffold", "--rounds", "6", "--min-growth", "1000")
    assert summary["stopped_by"] == "no_growth" and summary["rounds"] == 2  # round 1 has no earlier pool
    growth = [r["flank_growth_pct"] for r in rounds(out)]
    assert growth[0] == "" and float(growth[1]) < 1000


def test_max_flank_stops(tmp_path: Path) -> None:
    out = tmp_path / "fx"
    summary = run(out, "A-scaffold", "--rounds", "10", "--min-growth", "0", "--max-flank", "400")
    assert summary["stopped_by"] == "max_flank" and summary["rounds"] < 10
    assert all(r["max_flank_reached_pct"] for r in rounds(out))
    with Path(f"{out}.arg_hits.tsv").open() as h:
        (best,) = [r for r in csv.DictReader(h, delimiter="\t") if r["best_for_contig"] == "1"]
    assert min(int(best["upstream_flank_bp"]), int(best["downstream_flank_bp"])) >= 400


def test_no_detected_references(tmp_path: Path) -> None:
    empty = tmp_path / "empty.mapstat"
    empty.write_text("## method\tKMA\n# refSequence\treadCount\n")
    out = tmp_path / "fx"
    summary = run(out, "B-graph", mapstat=empty)
    assert summary["stopped_by"] == "no_detected_references" and summary["rounds"] == 0
    assert Path(f"{out}.fasta").stat().st_size == 0
