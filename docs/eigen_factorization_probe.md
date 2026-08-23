# Eigen SparseLU numerical-factorization probe

This is a benchmark-only investigation.  It makes no production solver or
Eigen-header change and does not run a production spectrum.  The executable
was compiled from `benchmarks/eigen_factorization_probe.cpp` with GCC 15.2.0
and the same Release build flags as the accepted E2 build.

## Eigen source audit

The build uses Eigen 3.4.0 from
`/home/adcm2/Documents/c++/DSpecM3D/build-docs/_deps/eigen3-src`, source
revision `3147391d946bb4b6c68edd901f2add6ac1f31f8c`.  In
`Eigen/src/SparseLU/SparseLU.h`, `SparseLU` initializes
`m_diagpivotthresh` to `1.0` and `setPivotThreshold()` controls the diagonal
pivot acceptance threshold.  `pivotL()` selects a diagonal pivot when it is
large enough relative to the column maximum; otherwise it searches the
column for the largest acceptable entry and applies a row interchange.

The installed source defaults are:

| setting | default |
|---|---:|
| `panel_size` | 16 |
| `relax` | 1 |
| `maxsuper` | 128 |
| `rowblk` | 16 |
| `colblk` | 8 |
| `fillfactor` | 20 |

`rowblk` and `colblk` are not used by this `factorize()` implementation's
panel path (the temporary allocation passes `m` for the row block); they were
therefore not swept.  The probe uses a tiny derived type to set the protected
`panel_size`, `relax`, and `maxsuper` fields without modifying Eigen.

`analyzePattern()` computes the column ordering, permuted structural matrix,
elimination tree, and postorder.  Every `factorize()` then copies/permutates
the current values, allocates work arrays, relaxes supernodes, and repeats
panel DFS, panel updates, column DFS, column updates, pivot selection,
pruning, supernode construction, and final factor cleanup.  Pivoting and the
numeric updates necessarily depend on matrix values.  The panel/column DFS
and supernode decisions are also performed in `factorize()` in this Eigen
version, so they are not reusable through the public API even when the sparse
pattern is unchanged.

## Probe design

The seven real PREM SEM systems are radial small (`n=104`), toroidal small
(`n=37`), spheroidal small (`n=158`), spheroidal medium (`n=257`), spheroidal
large (`n=860`), fluid-solid spheroidal (`n=860`), plus spheroidal max
(`n=1610`, near the production maximum).  All use attenuation-enabled
`H-w^2P+chi Ha` matrices from the current E2-style sparse formation pathway.
The benchmark directly reuses the production
`SPARSESPEC::detail::FixedPatternFrequencyMatrix` helper: its union pattern
and component values are prepared outside the timed region and its `update()`
method fills only the compressed value buffer before each factorization.

The fixed frequencies are `0.18, 0.41, 0.67, 1.05`.  A cheap NaturalOrdering
scan over `0.10` through `1.20` selects the frequency with the largest
backward residual and the frequency with the largest row-pivot displacement
for each system; duplicates are removed.  The selected frequencies are
reported on stderr in `eigen_factorization_probe_scan.tsv` and retained with
the raw data.  This is a difficult-frequency proxy, not proof that a true
physical resonance was found: multiplier growth was not directly accessible
or measured.  Broader resonance/model validation is required before any
pivot-threshold change.

Each configuration performs one untimed warm-up, then one measured batch of
four factorization calls and four solve calls per system/frequency.  Results
contain per-factor times, factor/solve status, `nnz(L)`, `nnz(U)`, nonidentity
row-pivot counts, maximum pivot displacement, discrepancy against the
NaturalOrdering threshold-1 baseline, and normalized backward residual
`||Ax-b||inf/(||A||inf||x||inf+||b||inf)`.

Pivot thresholds tested were `1, 0.5, 0.1, 0.01, 0`.  Performance screening
used the default, panel sizes `1, 4, 8, 32`, `relax=4`, and
`maxsuper=64,256`, one factor at a time, followed by only two small
interactions (`panel=4,relax=4` and `panel=4,maxsuper=64`).  No setting was
integrated into production.

## Result and decision

All 600 measured rows were finite with zero factorization and solve status.
The maximum solution discrepancy from the threshold-1 baseline was
`4.288207415733035e-11`; the maximum normalized backward residual was
`1.6197161038352034e-16`.  The largest system had `n=1610`.

Relative to threshold 1, aggregate factorization speedups were approximately
`1.02x` at threshold 0.5, `1.07x` at 0.1, `1.21x` at 0.01, and `1.42x` at
0.0.  Threshold zero also reduced average fill and eliminated row pivots in
this selected set.  The performance-parameter screen was modest: panel 4
was about `1.14x`, panel 1/8 about `1.07x`, and the relax/maxsuper variants
were within about 3% of baseline.  The duplicate `pivot_1` and
`natural_default` configurations were deliberately measured independently;
their roughly 2.8% timing difference is an empirical run/order noise floor for
this compact batch.  Parameter results near that floor should therefore be
interpreted cautiously.

Decision outcome: **A — kernel candidate, no adoption in this stage.**
Threshold zero is a clear kernel candidate, but it is not adopted here: a
lower pivot threshold must be validated on a broader resonance and model
set, including production-level numerical comparisons, before changing the
production default.  The measured 1.42x factorization result would imply a
theoretical complete-worker speedup of about `1.25x` if the accepted E2
factorization share (`67.0996%`) transferred unchanged.  This is an upper
bound for this kernel experiment, not a production claim.

If threshold-zero validation fails on broader systems, the small parameter
sweeps do not justify further Eigen tuning.  The next investigation should
then be a structure-aware or alternative sparse LU implementation (KLU,
SuperLU, or UMFPACK) before considering modified Eigen internals.

## Reproduction

```bash
cd /space/adcm2/CleanBuilds/DSpecM1D
/opt/gcc-15.2.0/bin/g++ \
  -DDSPECM1D_ENABLE_PROFILING -DDSPECM1D_ENABLE_LAPACK_BAND_SOLVER \
  -DHAVE_LAPACK_CONFIG_H \
  -DLAPACK_COMPLEX_CPP \
  -I/tmp/dspecm1d-backend-on-gcc15-final/benchmarks \
  -I/space/adcm2/CleanBuilds/DSpecM1D \
  -I/home/adcm2/Documents/c++/DSpecM3D/build-docs/_deps/eigen3-src \
  -I/home/adcm2/Documents/c++/DSpecM3D/build-docs/_deps/fftwpp-src \
  -I/home/adcm2/Documents/c++/DSpecM3D/build-docs/_deps/numericconcepts-src/include \
  -I/home/adcm2/Documents/c++/DSpecM3D/build-docs/_deps/gaussquad-src \
  -I/home/adcm2/Documents/c++/DSpecM3D/build-docs/_deps/gshtrans-src \
  -I/home/adcm2/Documents/c++/DSpecM3D/build-docs/_deps/interpolation-src \
  -O3 -DNDEBUG -fopenmp -std=gnu++23 \
  benchmarks/eigen_factorization_probe.cpp -o /tmp/eigen_factorization_probe \
  /usr/lib/x86_64-linux-gnu/libfftw3.so \
  /usr/lib/x86_64-linux-gnu/libfftw3f.so \
  /usr/lib/x86_64-linux-gnu/libfftw3l.so \
  /opt/gcc-15.2.0/lib64/libgomp.so /lib/x86_64-linux-gnu/libpthread.a \
  /usr/lib/x86_64-linux-gnu/liblapacke.so -lm -ldl \
  /usr/lib/x86_64-linux-gnu/libopenblas.so
OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OPENBLAS_NUM_THREADS=1 \
MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 \
/tmp/eigen_factorization_probe /tmp/dspecm1d-backend-on-gcc15-final/data/models/prem.200.no.txt \
  > results/eigen_factorization_probe/eigen_factorization_probe_raw.tsv \
  2> results/eigen_factorization_probe/eigen_factorization_probe_scan.tsv
python3 scripts/eigen_factorization_probe_analyze.py \
  results/eigen_factorization_probe/eigen_factorization_probe_raw.tsv \
  results/eigen_factorization_probe/eigen_factorization_probe_summary.json
```
