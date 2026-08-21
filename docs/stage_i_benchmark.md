# Stage I solver benchmark

`stage_i_benchmark.py` compares the complete preferred adaptive multi-SEM
`SparseFSpec::spectra(paramsNew)` call in separate Release builds.  It uses the
ex3/PREM-style type-4 all-mode attenuating input, fixes the OpenMP and BLAS
thread settings, performs exactly one untimed warm-up per backend/fmax, and
alternates which backend is launched first at each fmax.

The final grid is `fmax = 5, 10, 20, 35, 55 mHz`, with three measured
repetitions per backend and grid point.  The driver writes machine-readable TSV
timings and JSON reproduction metadata.  `stage_i_plot.py` writes a single-axis
two-line median wall-time figure and a JSON table containing timing spread and
`Eigen median / LAPACK median` speedups.

Before the full grid, a representative `fmax=5 mHz` pilot used the actual ex3
template `lmax=750`, one warm-up, and one measured repetition per backend.  It
confirmed finite `6 x 16385` outputs with consistent norms, complete-`spectra()`
timing, fair backend order/thread controls, and practical runtime: Eigen took
`10.388261539 s` and LAPACK took `10.102468484 s`.  Its exact raw data and run
metadata are retained as `results/stage_i/stage_i_pilot_timings.tsv` and
`results/stage_i/stage_i_pilot_metadata.json`.

The PR1 LAPACK path intentionally measures the current substitution boundary:

```text
Eigen sparse frequency matrix -> LAPACK band packing -> zgbtrf/zgbtrs
```

Thus the reported LAPACK time includes sparse-to-band packing overhead; it does
not measure the later direct-band-assembly optimization.

Example final run (after both Release builds exist):

```bash
python3 scripts/stage_i_benchmark.py \
  --eigen-build /tmp/dspecm1d-stage-i-eigen \
  --lapack-build /tmp/dspecm1d-stage-i-lapack \
  --output-dir results/stage_i \
  --threads 2 --repetitions 3

python3 scripts/stage_i_plot.py \
  results/stage_i/stage_i_timings.tsv \
  --output results/stage_i/stage_i_wall_time.png \
  --summary results/stage_i/stage_i_summary.json
```

The retained Stage I run is in:

- `results/stage_i/stage_i_timings.tsv`: 30 raw measured calls;
- `results/stage_i/stage_i_metadata.json`: machine, build, thread, and
  segmented-execution provenance;
- `results/stage_i/stage_i_summary.json`: medians, min/max spread, and speedup;
- `results/stage_i/stage_i_wall_time.png`: the primary one-axis timing figure.
