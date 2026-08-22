# Experimental post-`zgbtrf` backsolve

This stage tests a thread-local Eigen/C++ replay of LAPACK 3.11.0
`ZGBTRS(TRANS='N')`.  The factorization and `IPIV` are still produced by
`zgbtrf`; the solve replays `ZSWAP`, `ZGERU`, and non-unit upper-band
back-substitution with Eigen row/block expressions.  It supports arbitrary
right-hand-side counts and the trailing active-column offset used by
`LapackBandSolver`.  No BLAS or LAPACK routine is called after factorization
on the experimental path.  The old `LAPACKE_zgbtrs` call remains locally in
`LapackBandSolver.h` as a reference implementation.

The translation was checked against
`Reference-LAPACK/lapack v3.11.0/SRC/zgbtrs.f`, whose SHA-256 is
`96e35c7c9ceb47aa7cc51bd21f16365a5f713b44ceb635802fa5eb1ac58ae550`.
This reference source is recorded separately from the actual runtime: the
probe loads Debian `liblapacke.so.3.11.0` plus the pthread OpenBLAS-provided
`liblapack.so.3`, `libopenblas.so.0`, and `libblas.so.3`.  Their resolved paths
and hashes are retained in the experiment metadata.

The focused probe is
`benchmarks/lapack_custom_backsolve_probe.cpp`.  It validates random and
deliberately pivoted complex band systems, single and multiple RHS, and full
and trailing systems against `LAPACKE_zgbtrs` and original-matrix residuals.
It then times identical pre-factorized systems at `n=17`, `256`, and `1554`
with one warm-up and one measured run at 1, 10, 32, and 50 OpenMP threads.
The analysis script validates the raw TSV and writes the JSON summary and
metadata.

All validation cases passed; two of six exercised row pivots.  The largest
normalized solve difference was below `7e-16`, and the largest normalized
residual was below `2e-16`.  At 32 threads the custom solve was about 102x
faster for `n=17`, 4x faster for `n=256`, and within a few percent of
`LAPACKE_zgbtrs` for `n=1554`.  The fmax=35 mHz / OMP=32 production pilot
used the universal experimental path and measured a 7.86 s LAPACK wall time
versus 22.50 s for the retained LAPACK result, with full-output discrepancy
`8.18e-16` against the Eigen reference.

The production evidence is retained from two separate executions, each with
one untimed warm-up and one measured run per backend.  The profiling execution
uses `stage_profile_solver` and is the primary source for wall time and solve
worker-time accounting (`custom_f35_t32_pilot_profiles.jsonl`).  The separate
`final_solver_benchmark` execution writes the complete final output comparison
(`custom_f35_t32_final_check.txt`); its wall time is retained only as a
cross-check because that harness runs both backends in one process and is not
the primary timing source.  Both exact commands, executable paths, and hashes
are recorded in the two metadata files.

The raw output, production profile, numerical check, summary, and metadata
are retained in `results/lapack_custom_backsolve/`.  Timings are exploratory
and use one warm-up plus one measured run; no size threshold or source-depth
policy is introduced.

## Final single-thread production gate

A final `fmax=20 mHz`, `OMP_NUM_THREADS=1` comparison used one warm-up and one
measured complete spectrum per backend.  The old `zgbtrs` executable was built
from accepted revision `8f7a5d3` with the same GCC 15.2 Release configuration
and the same parameter-file hash as the custom build.

| backend | wall (s) | factorization worker-s | solve worker-s |
|---|---:|---:|---:|
| LAPACK `zgbtrs` | 87.7334 | 44.3517 | 27.2282 |
| custom Eigen/C++ backsolve | 89.3311 | 44.3283 | 29.0430 |
| Eigen SparseLU control | 130.1250 | 75.9676 | 21.2108 |

Making the custom backsolve universal therefore costs 1.82% complete wall time
and 6.67% solve time in this serial workload; factorization changes by -0.05%.
The LAPACK workload counts are identical (983,155 systems and solves; 2,948,155
RHS), and the absolute difference between the retained output norms is
`3.9031e-18` (`3.8914e-18` after normalization).  This gate did not change the
implementation or add a size threshold.
