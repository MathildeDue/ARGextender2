from pathlib import Path

from argextender.flanks import flank_total, parse_sam, share_reached, write_hits


def sam(tmp_path: Path, *records: str) -> Path:
    path = tmp_path / "x.sam"
    path.write_text("@HD\tVN:1.0\n" + "".join(r + "\n" for r in records))
    return path


def rec(contig: str, flag: int, seed: str, pos: int, cigar: str) -> str:
    return f"{contig}\t{flag}\t{seed}\t{pos}\t255\t{cigar}\t*\t0\t0\tACGT\t*\tAS:i:1"


def test_forward_hit_flanks_and_coordinates(tmp_path: Path) -> None:
    # contig of 1000 bp: 200 bp, then gene bases 1-500 of a 500 bp seed, then 300 bp
    hits = parse_sam(sam(tmp_path, rec("c1", 0, "g", 1, "200S500=300S")), {"c1": 1000}, {"g": 500})
    (h,) = hits
    assert (h.upstream_flank_bp, h.downstream_flank_bp) == (200, 300)
    assert (h.contig_start, h.contig_end, h.strand) == (201, 700, "+")
    assert (h.seed_start, h.seed_end, h.seed_coverage, h.hit_class) == (1, 500, 1.0, "complete")
    assert h.best_for_contig


def test_reverse_hit_flanks_follow_the_gene(tmp_path: Path) -> None:
    # reverse strand: the CIGAR is in gene orientation, so the leading clip is still upstream
    (h,) = parse_sam(sam(tmp_path, rec("c1", 16, "g", 1, "100S500=400S")), {"c1": 1000}, {"g": 500})
    assert (h.upstream_flank_bp, h.downstream_flank_bp) == (100, 400)
    assert (h.contig_start, h.contig_end, h.strand) == (401, 900, "-")


def test_filters_identity_length_and_unmapped(tmp_path: Path) -> None:
    records = (
        rec("c1", 4, "*", 0, "*"),  # unmapped
        rec("c2", 0, "g", 1, "50=50X"),  # identity 0.5
        rec("c3", 0, "g", 1, "80="),  # 80 gene bp < 100
        rec("c4", 0, "short", 1, "60="),  # a 60 bp gene fully covered: kept
    )
    hits = parse_sam(sam(tmp_path, *records), {f"c{i}": 1000 for i in range(1, 5)}, {"g": 500, "short": 60})
    assert [h.contig for h in hits] == ["c4"]


def test_identity_counts_indels_and_classes(tmp_path: Path) -> None:
    # gene bases 101-400 of 500: 300 gene bp (= X D), aligned 290 = + 5 X + 5 I + 5 D = 305
    (h,) = parse_sam(sam(tmp_path, rec("c1", 0, "g", 101, "290=5X5I5D")), {"c1": 300}, {"g": 500})
    assert h.identity == round(290 / 305, 4)
    assert (h.seed_start, h.seed_end) == (101, 400)
    assert h.hit_class == "partial-truncated"  # no flank on either side: the contig ends inside the gene
    (h2,) = parse_sam(sam(tmp_path, rec("c1", 0, "g", 101, "50S300=50S")), {"c1": 400}, {"g": 500})
    assert h2.hit_class == "partial-internal"


def test_best_per_contig_and_totals(tmp_path: Path) -> None:
    records = (
        rec("c1", 0, "g1", 1, "100S500=400S"),  # complete
        rec("c1", 0, "g2", 1, "700S200=100S"),  # partial: coverage 0.4
        rec("c2", 0, "g1", 1, "10S500="),
    )
    hits = parse_sam(sam(tmp_path, *records), {"c1": 1000, "c2": 510}, {"g1": 500, "g2": 500})
    best = {h.contig: h.seed for h in hits if h.best_for_contig}
    assert best == {"c1": "g1", "c2": "g1"}
    assert flank_total(hits) == (100 + 400) + (10 + 0)
    assert share_reached(hits, 100) == 0.5  # c1: min(100, 400); c2: min(10, 0)
    assert share_reached(hits, 0) == 1.0 and share_reached(hits, 101) == 0.0 and share_reached([], 1) == 0.0


def test_duplicate_records_kept_once_and_table(tmp_path: Path) -> None:
    r = rec("c1", 0, "g", 1, "200S500=300S")
    hits = parse_sam(sam(tmp_path, r, r), {"c1": 1000}, {"g": 500})
    assert len(hits) == 1
    write_hits(hits, tmp_path / "t.tsv")
    header, row = (tmp_path / "t.tsv").read_text().splitlines()
    assert header.split("\t")[-2:] == ["class", "best_for_contig"]
    assert row.split("\t")[-2:] == ["complete", "1"]
