import gzip
from pathlib import Path

import pytest

from argextender import ArgExtenderError
from argextender.graph import graph_targets, parse_gfa, select_components
from argextender.kma import detected_references, match_asm, read_mapstat, select_seeds
from argextender.reads import Reads, parse_internal_id


def fastq(path: Path, names: list[str]) -> Path:
    path.write_text("".join(f"@{n} extra\n{'ACGT' * (i + 1)}\n+\n{'I' * 4 * (i + 1)}\n" for i, n in enumerate(names)))
    return path


def test_reads_prepare_pool_and_templates(tmp_path: Path) -> None:
    r1 = fastq(tmp_path / "r1.fq", ["a/1", "b/1", "c/1"])
    r2 = fastq(tmp_path / "r2.fq", ["a/2", "b/2", "c/2"])
    s = fastq(tmp_path / "s.fq", ["x"])
    reads, internal = Reads.prepare(r1, r2, s, tmp_path / "copies")
    assert (reads.pairs, reads.singletons) == (3, 1)
    assert internal["1"].read_text().splitlines()[0] == "@x0_1"
    reads.write_pool({2, 0}, {0}, *(tmp_path / f"p{m}.fq" for m in "12s"))
    assert (tmp_path / "p1.fq").read_text().splitlines()[::4] == ["@a/1 extra", "@c/1 extra"]
    assert (tmp_path / "ps.fq").read_text().startswith("@x extra\n")
    with (tmp_path / "t.fa").open("w") as h:
        reads.write_templates({1}, set(), h)
    assert (tmp_path / "t.fa").read_text() == ">x1_1\nACGTACGT\n>x1_2\nACGTACGT\n"
    assert parse_internal_id("x12_s") == (12, "s")


def test_reads_out_of_sync(tmp_path: Path) -> None:
    r1 = fastq(tmp_path / "r1.fq", ["a/1", "b/1"])
    with pytest.raises(ArgExtenderError, match="out of sync"):
        Reads.prepare(r1, fastq(tmp_path / "r2.fq", ["a/2", "z/2"]), None, tmp_path / "c")
    with pytest.raises(ArgExtenderError, match="records but R2"):
        Reads.prepare(r1, fastq(tmp_path / "r3.fq", ["a/2"]), None, tmp_path / "d")


def test_seeds_from_mapstat(tmp_path: Path) -> None:
    m = tmp_path / "s.mapstat"
    m.write_text("## method\tKMA\n# refSequence\treadCount\ng2\t5\ng1 desc\t3\ng2\t1\n")
    assert detected_references(m) == ["g2", "g1 desc"]
    fa = tmp_path / "p.fa"
    fa.write_text(">g1 desc\nAAA\n>g2\nCC\nGG\n>g3\nT\n")
    assert select_seeds(detected_references(m), fa) == [("g2", "CCGG"), ("g1 desc", "AAA")]
    with pytest.raises(ArgExtenderError, match="not found"):
        select_seeds(["g9"], fa)
    m.write_text("# refSequence\treadCount\tpropCovered\ng1\t3\t0.9\n")
    with pytest.raises(ArgExtenderError, match="filtered mapstat"):
        read_mapstat(m)


def test_match_asm_first_occurrence(tmp_path: Path) -> None:
    frag = tmp_path / "s.frag.gz"
    with gzip.open(frag, "wt") as h:
        h.write("AAA\t1\t9\t0\t3\tg1\tc2\nAAA\t1\t9\t0\t3\tg2\tc2\nCCC\t1\t9\t0\t3\tg1\tc1\n")
    assert match_asm(frag) == [("c2", "AAA"), ("c1", "CCC")]


GFA = "S\t1\tAAAA\nS\t2\tCCCC\nS\t3\tGGGG\nS\t10\tTTTT\nS\t4\tACAC\nL\t1\t+\t2\t+\t2M\nL\t2\t+\t10\t-\t2M\nL\t3\t+\t4\t+\t2M\n"


def test_components_order(tmp_path: Path) -> None:
    (tmp_path / "g.gfa").write_text(GFA)
    graph = parse_gfa(tmp_path / "g.gfa")
    # anchor 2: its component {1, 2, 10}; anchors first, then by link distance, then by name
    assert select_components(graph, ["2"]) == [("2", "arg"), ("1", "component"), ("10", "component")]
    assert [s for s, _ in select_components(graph, ["10", "4"])] == ["4", "10", "2", "3", "1"]


def test_graph_targets_fallback_without_assembly(tmp_path: Path) -> None:
    sel = tmp_path / "selected.fasta"
    sel.write_text(">c1\nAAAA\n")
    source, reason = graph_targets(None, sel, tmp_path / "seeds.fa", None, tmp_path / "next.fa")  # type: ignore[arg-type]
    assert (source, reason) == ("scaffolds (fallback)", "no assembly in this round")
    assert (tmp_path / "next.fa").read_text() == ">c1\nAAAA\n"
