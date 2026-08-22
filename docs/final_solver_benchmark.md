# Final production solver benchmark

This campaign compares the two user-facing adaptive multi-SEM backends from
commit `d08fc5a063816d9e4a8f4851a2c40b9e65019b63` using one unprofiled GCC 15.2
Release executable.  The executable was built with
`DSPECM1D_ENABLE_LAPACK_BAND_SOLVER=ON`; each run selected either
`SolverBackend::EigenSparseLU` or `SolverBackend::LapackBandLU` at runtime.

Both selections used the same ex3/PREM all-mode attenuating calculation,
`lmax=750`, source, receivers, output grid, adaptive chunking, and numerical
settings.  Every point performed one untimed warm-up per backend before timing
the complete preferred `spectra()` call.  Backend launch order alternated
across points.  `OMP_DYNAMIC=FALSE` and the OpenBLAS, MKL, and BLIS thread
counts were fixed to one.  The machine was an AMD EPYC 9334 server with two
32-core sockets.

The LAPACK implementation measured here is the accepted direct-band path:
the sparse SEM base operators are converted to band storage once per
degree/chunk, then each frequency forms the active band matrix directly before
`zgbtrf`/`zgbtrs`.  It is not the earlier per-frequency sparse-pack experiment.

## Primary 32-thread result

`Speedup` is `t_Eigen / t_LAPACK`; values below one mean Eigen is faster.  The
wall-time reduction is `100 * (1 - t_LAPACK / t_Eigen)`, so a negative value
means LAPACK increased wall time.

| fmax (mHz) | Eigen min/median/max (s) | LAPACK min/median/max (s) | Speedup | LAPACK reduction |
|---:|---:|---:|---:|---:|
| 5  | 0.656 / 0.657 / 0.661 | 2.916 / 2.919 / 2.922 | 0.225 | -344.3% |
| 10 | 1.487 / 1.489 / 1.489 | 5.770 / 5.850 / 5.854 | 0.255 | -292.9% |
| 20 | 4.059 / 4.059 / 4.060 | 12.159 / 12.168 / 12.169 | 0.334 | -199.7% |
| 35 | 11.171 / 11.175 / 11.184 | 22.488 / 22.497 / 22.639 | 0.497 | -101.3% |
| 55 | 30.323 / 30.329 / 30.342 | 40.155 / 40.214 / 40.303 | 0.754 | -32.6% |

At the production thread count, Eigen SparseLU is faster at every measured
frequency.  The LAPACK disadvantage decreases as `fmax` rises, but this grid
alone does not establish the cause of that trend.

## Single-thread algorithmic check

These are one measured call per backend after warm-up, so they are sanity
measurements rather than distribution estimates.

| fmax (mHz) | Eigen (s) | LAPACK (s) | Speedup | LAPACK reduction |
|---:|---:|---:|---:|---:|
| 5  | 20.339 | 13.343 | 1.524 | 34.4% |
| 10 | 46.561 | 30.911 | 1.506 | 33.6% |
| 20 | 127.650 | 84.648 | 1.508 | 33.7% |

The direct-band backend has a consistent serial algorithmic advantage of
about 1.51x, or roughly 34% lower complete wall time, over these cheap points.

## Thread scaling at 35 mHz

The 10, 20, and 50-thread rows use two measured repetitions.  The 32-thread
row reuses the three-repetition primary measurement.

| Threads | Eigen median (s) | LAPACK median (s) | Speedup | LAPACK reduction |
|---:|---:|---:|---:|---:|
| 10 | 35.613 | 31.550 | 1.129 | 11.4% |
| 20 | 17.855 | 25.022 | 0.714 | -40.1% |
| 32 | 11.175 | 22.497 | 0.497 | -101.3% |
| 50 | 7.252 | 21.317 | 0.340 | -194.0% |

Between 10 and 50 threads, the measured complete-call wall-time improvement
is 4.91x for Eigen and 1.48x for LAPACK.  This is an observed scaling ratio,
not a parallel-efficiency claim.  It explains the practical crossover: LAPACK
is faster at one and ten threads, while Eigen is faster at 20 threads and
above on this machine.  The campaign does not isolate which OpenMP scheduling,
memory-traffic, or per-system factorization effects produce the different
scaling; that would require a separate profiling experiment.

## Numerical validation

Every result was finite and had shape `6 x 16385`.  Full output matrices, not
only norms, were compared after each paired point.

| Series | fmax (mHz) | Threads | Maximum absolute difference |
|---|---:|---:|---:|
| primary | 5  | 32 | 1.326e-16 |
| primary | 10 | 32 | 3.917e-16 |
| primary | 20 | 32 | 5.862e-16 |
| primary | 35 | 32 | 8.175e-16 |
| primary | 55 | 32 | 7.655e-16 |
| single-thread | 5  | 1 | 1.326e-16 |
| single-thread | 10 | 1 | 3.917e-16 |
| single-thread | 20 | 1 | 5.862e-16 |
| scaling | 35 | 10 | 8.175e-16 |
| scaling | 35 | 20 | 8.175e-16 |
| scaling | 35 | 50 | 8.175e-16 |

The largest discrepancy was `8.175450883859586e-16`, far below the enforced
`1e-9` tolerance.

## Reproduction and artifacts

The driver appends each completed point to its TSV files, flushes and syncs
them, and updates metadata atomically.  Re-running the same command and tag
skips completed points and discards only incomplete timing rows lacking a
paired comparison checkpoint.

- `final_production_timings.tsv`: 48 measured complete calls.
- `final_production_comparisons.tsv`: 11 full-output comparisons.
- `final_production_metadata.json`: machine, build, request, hashes, and
  completion state.
- `final_production_summary.json`: validated derived results.
- `final_production_primary.csv`, `final_production_single_thread.csv`,
  `final_production_scaling.csv`, and `final_production_numerical.csv`: compact
  machine-readable tables.
- `final_production_wall_time.png`: primary 32-thread wall-time figure.
- `final_production_thread_scaling.png`: 35 mHz thread-scaling figure.

The raw timing SHA-256 is
`cff2223db59efdff5d431cf6c9ffcd4b0b1cad1da28df13364c1794cdcc76803`.
The raw comparison SHA-256 is
`833bc4f401c1377c6ba031d24d3534bca2f2e9e9e591b8e6011f66178159b859`.
