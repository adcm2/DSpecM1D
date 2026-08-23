# E1b Eigen NaturalOrdering integration

This stage changes only the adaptive multi-SEM Eigen SparseLU ordering. The
single-SEM `detail::SingleSemSolver` remains
`Eigen::COLAMDOrdering<int>`. The adaptive path has one affected alias,
`detail::EigenMultiSemSolver`, used by the radial `eigenSolver` and the
toroidal/spheroidal `eigenSolver1` worker-private instances. Both now use
`Eigen::NaturalOrdering<int>`.

No cadence, sparse matrix construction, compression, OpenMP ownership, LAPACK
path, custom backsolve, backend API, or physics code was changed. In
particular, `isSymmetric(true)` is not used by the production integration.

## Validation

The existing `PreferredSolverApiTests` suite passed all six tests. The focused
real-SEM validation was rerun with the rebuilt production source and covered
11 systems: radial, toroidal, spheroidal, an explicit fluid-solid spheroidal
case, and a high-resolution spheroidal case. All factorization and solve
status values were zero, every output was finite, and the maximum residual was
`1.375788140840027e-16`. The maximum relative solution discrepancy of
NaturalOrdering versus the COLAMD reference was
`2.848933726852794e-12`.

Mean focused timings for `isSymmetric(false)` were:

| ordering | analyze (ms) | factorize (ms) | solve (ms) |
|---|---:|---:|---:|
| COLAMD | 0.111623 | 0.726548 | 0.119949 |
| NaturalOrdering | 0.038110 | 0.544921 | 0.095636 |

## Required production pilots

Each pilot used one untimed warm-up and one measured complete `spectra()` call,
with fixed `OMP_DYNAMIC=FALSE`, `OPENBLAS_NUM_THREADS=1`,
`MKL_NUM_THREADS=1`, and `BLIS_NUM_THREADS=1`. The retained E1 COLAMD outputs
are the exact comparable controls; new NaturalOrdering outputs are in
`results/eigen_sparselu_e1b/`.

| pilot | COLAMD (s) | NaturalOrdering (s) | speedup |
|---|---:|---:|---:|
| 5 mHz, 1 thread | 20.264759730 | 19.320262735 | 1.04889x |
| 30 mHz, 32 threads | 8.362349297 | 7.937082845 | 1.05358x |

The complete-result norms differ by `7.1536e-18` at 5 mHz and
`1.1601e-17` at 30 mHz. The focused improvement therefore survives production
integration at both amended pilot points. These timings are diagnostic and
were not repeated.

## Reproduction

The NaturalOrdering binary was built from this worktree with GCC 15.2 and the
existing Release/LAPACK-enabled configuration:

```bash
cmake --build /tmp/dspecm1d-backend-on-gcc15-final \
  --target dspecm1d_solver_api_tests eigen_sparselu_e1 eigen_sparselu_e1_pilot -j4
```

The focused command was run with the environment above and PREM model
`/tmp/dspecm1d-backend-on-gcc15-final/data/models/prem.200.no.txt`. The pilot
commands used the matching `e1_pilot_f5.txt` and `e1_pilot_f30.txt` parameter
files. Executable and parameter hashes are retained in
`results/eigen_sparselu_e1b/eigen_sparselu_e1b_metadata.json`.

Stage E2, fixed-pattern/direct-value sparse `A(omega)` formation, is deferred.
