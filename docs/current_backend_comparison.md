# Current Eigen/LAPACK production-backend comparison

This uncommitted benchmark-only stage compares the two production backends at
committed revision `3fcb02df45169877b08cc24c8b11c331e6c96906`. No solver,
physics, backend API, threshold, ordering, matrix-formation, OpenMP, truncation,
or cadence code is changed.

## Method

The existing `final_solver_benchmark` executable selects `EigenSparseLU` or
`LapackBandLU` at runtime. Both backends therefore use the same executable,
PREM model, ex3 all-mode attenuating configuration, source/receivers, lmax 750,
frequency grid, damping, and output shape. The current Eigen path is
NaturalOrdering, default pivot threshold 1, fixed-pattern sparse A(omega), and
SparseLU. The current LAPACK path is direct band A(omega), `zgbtrf`, and the
local C++/Eigen backsolve; its profile records zero band packs.

Primary timings use GCC 15.2 Release (`-O3 -DNDEBUG`) and measure the complete
preferred adaptive multi-SEM `spectra()` call. Every backend/point has one
untimed warm-up and exactly one measured call. Backend order alternates:
Eigen-first at 5 mHz, LAPACK-first at 30 mHz, and Eigen-first at 100 mHz.
`OMP_DYNAMIC=FALSE`; OpenBLAS, MKL, and BLIS threads are fixed at one.

The accepted source has a pre-existing non-profiling guard around three
`patternRebuilt` declarations but not all uses. The unprofiled binary therefore
uses the retained declaration-only `/tmp` overlay in
`nonprofiling_compile_workaround.diff`. It changes no executed solver logic and
is identical for both runtime-selected backends. The separately instrumented
binary compiles the committed source directly. Profile category values are
summed OpenMP worker seconds and are diagnostic; unprofiled wall time is the
primary result.

An independent input audit found that the initial 5 and 30 mHz records used
`prem.200.no.txt` rather than the canonical `prem.200` used at 100 mHz. Those
invalid records were discarded and the two affected points were rerun. Every
retained timing and profile artifact now uses the same canonical `prem.200`;
the already-valid 100 mHz artifacts were left unchanged.

## Primary result

Speedup is `Eigen wall / LAPACK wall`; values greater than one mean LAPACK is
faster. “Less wall” is `(Eigen-LAPACK)/Eigen`.

| fmax | threads | Eigen wall | LAPACK wall | LAPACK speedup | LAPACK less wall |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 5 mHz | 1 | `16.836303676 s` | `13.696712328 s` | `1.22922x` | `18.65%` |
| 30 mHz | 32 | `6.917206552 s` | `5.610110935 s` | `1.23299x` | `18.90%` |
| 100 mHz | 32 | `107.134988079 s` | `85.849940865 s` | `1.24793x` | `19.87%` |

LAPACK wins all three points. The speedup stays within `1.229–1.248x`, so the
relative advantage does not change materially as fmax and maximum active
dimension grow. These are current-vs-current measurements; historical
intermediate timings are not used in these ratios.

## Kernel attribution

The following are summed worker seconds from separate profiling runs and do
not add to wall time.

| fmax/backend | A(omega) | analyze | factorization | solve | projection | unclassified |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 30 Eigen | `20.795` | `0.356` | `150.796` | `38.665` | `0.997` | `6.400` |
| 30 LAPACK | `12.761` | `0` | `93.196` | `60.336` | `0.981` | `13.010` |
| 100 Eigen | `296.681` | `2.357` | `2450.382` | `612.260` | `4.043` | `74.322` |
| 100 LAPACK | `206.139` | `0` | `1512.683` | `964.721` | `4.189` | `150.043` |

Eigen factorization costs `1.618x` the LAPACK factorization worker time at
30 mHz and `1.620x` at 100 mHz. Direct-band A(omega) is also cheaper than the
fixed sparse update. LAPACK's custom backsolve is slower than Eigen's solve
(`1.56–1.58x` worker time) and partly gives back the factorization gain.
Projection is effectively equal. The dominant remaining advantage is therefore
factorization, with matrix formation secondary—not solve.

The complete Eigen backend is about 19–20% slower, but the numerical
factorization kernel still has a substantial, stable gap at both parallel
points. A future experiment with an alternative sparse LU or structure-aware
factorizer remains justified; this stage does not begin one.

## Workload

| fmax | frequency systems | factorizations | solves | RHS | active n | kl/ku | output shape |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 5 | 244,663 | 244,663 | 244,663 | 733,663 | 17–354 | 4–15 | 6×16385 |
| 30 | 1,475,483 | 1,475,483 | 1,475,483 | 4,424,483 | 17–1362 | 4–15 | 6×16385 |
| 100 | 4,917,276 | 4,917,276 | 4,917,276 | 14,745,276 | 17–4302 | 4–15 | 6×16385 |

The matched profile counts validate identical scientific workloads per backend.
The 100 mHz point increases the largest active system from 1362 to 4302 while
retaining the same narrow bandwidth range.

## Numerical comparison

Every output was finite and had identical shape. The full matrix was compared
elementwise after the final measured call:

| fmax | Eigen norm | LAPACK norm | max absolute difference | difference/output norm |
| ---: | ---: | ---: | ---: | ---: |
| 5 | `1.3681986045269e-6` | `1.3681986045289e-6` | `2.6527e-17` | `1.9389e-11` |
| 30 | `1.6860737723921e-4` | `1.6860737723927e-4` | `4.8830e-16` | `2.8961e-12` |
| 100 | `2.1307195596634e-3` | `2.1307195596634e-3` | `1.3577e-15` | `6.3722e-13` |

The largest absolute discrepancy is `1.3577375468803944e-15` at 100 mHz,
well below the accepted `1e-9` tolerance.

## Reproduction

Exact build/run commands, binary/source/input hashes, thread settings, backend
order, and raw-output hashes are retained in
`results/current_backend_comparison/current_backend_metadata.json`.

The retained `current_f5.txt`, `current_f30.txt`, and `current_f100.txt` files
are the exact inputs passed to the executables.
For the recorded `/tmp/current-backend-comparison` layout, copy them to its
`data/params` directory and copy `data/models` beside that directory. The
temporary unprofiled include overlay is made by copying `DSpecM1D/src` to
`/tmp/current-backend-overlay/DSpecM1D/src` and applying the retained
declaration-only diff to that copied `FullSpecMultiSem.h`; the repository file
is not patched. The metadata then gives the exact compile, link, and invocation
commands.

Rebuild the derived tables and summary with:

```bash
python3 scripts/current_backend_comparison_analyze.py \
  results/current_backend_comparison
```
