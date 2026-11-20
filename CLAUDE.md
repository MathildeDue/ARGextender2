# ARGextender2: notes for development sessions

- **What this is.** The minimal final version of ARGextender. The development, benchmarks and decisions live in
  `/hpc/projects/2026/unseenbio/ARGprofiler2`. Specification: `docs/planning/14_final_version_spec.md` there.
  That repository stays as the record of the benchmark results and is not changed from here.
- **Scope.** A-scaffold (default) and B-graph (option; decided 2026-11-20, A-graph dropped), mapping-only
  recruitment (`kma -sasm`, k = 25) and the stopping rules (`--rounds`, `--min-growth` on the flank growth
  (KMA `-sam` placement after every assembly), `--max-flank`). No evaluation code, no legacy parity, no checkpoints/resume. New features only
  if they improve performance, and only after the user agrees.
- **Pinned tools.** KMA 1.6.18, SPAdes 3.15.5, Bandage 0.9.0 + BLAST+ 2.17.0, networkx 3.4.2. Do not change
  versions or KMA/SPAdes/Bandage flags without discussion: the benchmark results depend on them.
- **Environment on the HPC.** Use the existing envs (never install into shared conda envs):
  `E=/hpc/projects/2026/unseenbio/ARGprofiler2-envs; export PATH=$E/argextender-dev/bin:$PATH:$E/argextender-bandage/bin`
  The Bandage env goes at the *end* of PATH (its own python must not shadow the dev env's). The dev env has the
  old `argextender` installed editable: run from this repository's root without `python -I`, so this package is
  imported.
- **Tests.** `python -m pytest -q --basetemp=/hpc/projects/2026/unseenbio/ARGprofiler2/runs/dev/argx2_pytest -p no:cacheprovider`
  (compute-node `/tmp` is noexec). Run directly inside the user's interactive allocation; never nested `srun`.
- **Compute.** All real runs go through Slurm. Never submit batch jobs without the user's explicit approval;
  `sbatch --test-only` is fine. Run outputs go under `ARGprofiler2/runs/`, never into this repository.
- **Git.** Commit only when the user asks; never push unless asked. No data or run outputs in the repository
  (tests use `tests/data/`, synthetic).
