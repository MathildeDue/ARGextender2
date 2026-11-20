# ARGextender2

ARGextender extends the genetic context of antimicrobial resistance genes (ARGs) found by
[ARGprofiler](https://github.com/MathildeDue/ARGprofiler_human_microbiome). Starting from the ARGs detected in a
sample, it repeatedly recruits the reads that map to the current targets (KMA), assembles them (SPAdes), keeps
the scaffolds that carry an ARG and uses them, or the assembly graph around them, as the next round's targets.
With each round the contigs around the ARGs get longer.

## Usage

```
argextender --r1 R1.fq.gz --r2 R2.fq.gz [--singletons S.fq.gz] \
            --mapstat SAMPLE.mapstat --panres panres_genes.fa \
            --out PREFIX [--strategy {A-scaffold,B-graph}] [--threads 8]
```

`--mapstat` is ARGprofiler's unfiltered KMA `<sample>.mapstat` against PanRes: every gene listed there is a
seed. `python -m argextender` works as well as the `argextender` command.

### Strategies

| Strategy | Next round's targets | Assembly |
|---|---|---|
| `A-scaffold` (default) | the selected scaffolds | after every round |
| `B-graph` | the ARG segments of SPAdes' `strain_graph.gfa` (Bandage) and their whole connected components; between assemblies also the reads recruited so far | every `--assemble-every` rounds (default 3) and after the last round |

Both recruit by k-mer mapping only (`kma -sasm`, k = 25).

### Stopping

Checked after each assembly; the run stops when any rule applies.

| Option | Default | Stops when |
|---|---|---|
| `--rounds N` | 30 | N rounds are done (hard limit) |
| `--min-growth X` | 2 | the flanks around the contigs' best ARGs grew by less than X % per round since the previous assembly (`0` turns it off) |
| `--max-flank L` | 0 (off) | at least `--max-flank-share` % of the contigs have both flanks of their best ARG ≥ L bp |
| `--max-flank-share P` | 50 (the median) | the share of contigs `--max-flank` needs (100: all contigs) |

The run also stops when no targets remain. For a fixed number of rounds use `--rounds N --min-growth 0`.

The flanks are measured after every assembly: a separate KMA call (`-sam`) places the seeds on the selected
scaffolds, and the flank total is the summed up- and downstream flank of every scaffold's best ARG. The original
pipeline (`targetAsm.pl`) stopped when the reads recruited in a round did not grow; the read pool turned out to
be a noisy measure (it can shrink while the flanks still grow), so the flanks are used instead, with a
percentage threshold. Under `B-graph` the growth between two assemblies is converted to a per-round value,
(F / F_prev)^(1/N) − 1. The defaults are provisional.

### Outputs

| File | Content |
|---|---|
| `PREFIX.fasta` | the selected scaffolds of the last assembly |
| `PREFIX.gfa` | that assembly's graph (`assembly_graph_with_scaffolds.gfa`) |
| `PREFIX.frag.gz`, `PREFIX.frag_raw.gz` | the scaffolds mapped against the complete PanRes (KMA) |
| `PREFIX.arg_hits.tsv` | where each ARG sits on each final scaffold, its flanks and class; `best_for_contig` marks each scaffold's best ARG |
| `PREFIX.rounds.tsv` | one row per round: targets, recruited reads, assembly, flank total, flank growth, share at L (with `--max-flank`), stop reason |
| `PREFIX.summary.json` | command, parameters, tool versions, inputs, stop reason, run time |

The round folders go in `PREFIX.work/`, which is removed at the end unless `--keep-work` is given.

## Requirements

KMA 1.6.18, SPAdes 3.15.5, Bandage 0.9.0 with BLAST+ 2.17.0 (`B-graph` only), Python ≥ 3.11 with
networkx, all on `PATH` (`environment.yml`). These are the versions ARGextender was benchmarked with.

## Tests

```
pip install -e .[dev]
pytest
```

The end-to-end tests (`tests/test_pipeline.py`) run the real tools on a small synthetic data set
(`tests/data/`) and are skipped if a tool is missing.
