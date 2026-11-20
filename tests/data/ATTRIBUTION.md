# Test data

Synthetic data for correctness tests only, never for scientific claims. Copied from `tests/fixtures/data/` of
ARGprofiler2 (the ARGextender 2.0 development repository), where `make_fixture.py` generated it.

- `reads_1.fastq`, `reads_2.fastq`, `reads_singleton.fastq`: error-free 2 x 150 bp reads (insert 270 +/- 20 bp,
  ~30x) from a random 30 kb replicon carrying one ARG, `pan_27` (blaCTX-M-15), at 14,001-14,876 (+), as in
  `truth.tsv`; 200 singletons.
- `panres_subset20.fasta`: 20 unmodified PanRes gene sequences (`pan_27` and 19 decoys) from
  `panres2_filtered_gene_sequences.fasta` (PanRes v2.0 without E. coli homologs) in `prerequisites/manual_panres/`
  of https://github.com/MathildeDue/ARGprofiler_human_microbiome (Apache License 2.0). PanRes is published on
  Zenodo (concept DOI 10.5281/zenodo.8055115) under CC BY 4.0.
- `test.mapstat`: a synthetic unfiltered KMA mapstat listing `pan_27` and one decoy (`pan_1`) as detected; the
  counts carry no meaning. It stands in for ARGprofiler's detection step.
