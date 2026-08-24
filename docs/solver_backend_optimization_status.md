# Solver backend optimization stopping point

This document records the production architecture and the evidence-based
stopping decision at commit `90bffae47333fc93dac1a0a1dacc58aafa9a87fa` on
`dev/lapack-banded-pr1`. It describes the current source, not a proposed
backend or a historical experiment.

## Retained production architecture

The preferred adaptive multi-SEM API accepts
`SPARSESPEC::SolverBackend::{EigenSparseLU,LapackBandLU}` through
`InputParametersNew::setSolverBackend()`. `EigenSparseLU` is the default.
`DSPECM1D_ENABLE_LAPACK_BAND_SOLVER` is OFF by default; enabling it discovers
LAPACK/LAPACKE and adds `LapackBandLU` without removing Eigen. Without that
capability, selecting LAPACK through `InputParametersNew` throws
`std::invalid_argument`; the explicit legacy adaptive overload rejects the
same request with `std::runtime_error`.

Runtime backend selection applies to the adaptive multi-SEM path. The
single/reused-SEM overload remains Eigen-only (COLAMD ordering) and throws
`std::invalid_argument` if `LapackBandLU` is selected. Its caller-provided
`nskip` behavior is unchanged.

The stopping-point audit also restored ordinary non-profiling compilation by
keeping profiling-only `patternRebuilt` assignments inside their profiling
guards; this has no effect on solver behavior.

### EigenSparseLU

The adaptive path uses `Eigen::SparseLU` with `NaturalOrdering<int>`. No
production call changes the pivot threshold, so Eigen's robust default of 1.0
is retained. For each worker, SEM chunk, and active truncation, the E2 path
builds one compressed union pattern for H/P/(Ha), stores aligned coefficient
arrays, and forms each `A(omega)` by updating `valuePtr()` directly. It does
not insert, merge, or recompress sparse entries in the frequency hot path.

Radial systems call `compute()` at every frequency. Toroidal and spheroidal
systems use the existing reverse frequency traversal and the internally
derived `nskip = max(1, (i2-i1)/20)`: cadence points call `compute()`
(`analyzePattern` plus numerical factorization), while intermediate points
call `factorize()` against the established pattern. The existing
`allIndicesTor`/`allIndicesSph` results determine the trailing `ridx`; a
worker-private fixed-pattern object is rebuilt when that active truncation
changes. Each OpenMP worker owns its own SparseLU and fixed-pattern state, and
solves the existing arbitrary-column Eigen RHS before the unchanged receiver
projection.

### LapackBandLU

For each SEM/chunk and degree, the LAPACK path converts the sparse base H/P/Ha
operators once into exact column-major general-band storage. Every frequency
then forms `A(omega) = H - omega^2 P` (plus `chi Ha` with attenuation) directly
in a preallocated solver-owned band workspace. There is no per-frequency
sparse matrix construction or sparse-to-band packing.

Radial uses the complete band. Toroidal and spheroidal use the same adaptive
`ridx`/truncation sequence as Eigen and fill only the active trailing columns.
`factorize(ridx)` passes a column offset into the same workspace, so pivoted
`LAPACKE_zgbtrf` overwrites the active trailing matrix without a band-to-band
copy. Pivot and coefficient workspaces are reused. Each OpenMP worker owns its
solver state. The post-factorization solve is the local C++/Eigen replay of
the LAPACK `ZGBTRS` no-transpose operations and handles an arbitrary positive
number of RHS columns. `LAPACKE_zgbtrs` remains only as an unused private
reference implementation for validation; it is not the production solve hot
path.

## Current backend comparison

The accepted current-vs-current benchmark used the same GCC 15.2 Release
executable, PREM/ex3 all-mode attenuating physics, one warm-up and one measured
complete `spectra()` call per backend and point, fixed external BLAS threads,
and runtime backend selection.

| fmax / OpenMP threads | EigenSparseLU | LapackBandLU | Eigen/LAPACK |
| --- | ---: | ---: | ---: |
| 5 mHz / 1 | 16.8363 s | 13.6967 s | 1.2292x |
| 30 mHz / 32 | 6.9172 s | 5.6101 s | 1.2330x |
| 100 mHz / 32 | 107.1350 s | 85.8499 s | 1.2479x |

The largest retained full-output discrepancy is
`1.3577375468803944e-15`, below the accepted `1e-9` tolerance. The full method,
raw data, metadata, workload counts, and profiles are retained in
[`current_backend_comparison.md`](current_backend_comparison.md) and
`results/current_backend_comparison/`.

## Alternatives investigated and rejected

- Further `analyzePattern` caching was not justified: it was about 0.16% of
  summed worker time after E2. Evidence is in
  [`eigen_post_e2_profile.md`](eigen_post_e2_profile.md) and
  `results/eigen_post_e2_profile/`.
- Reduced Eigen pivot thresholds were faster, but were deliberately not
  adopted. The campaign retained threshold 1 because solver-error control and
  robustness near physical resonances outweigh the modest speedup. The 0.01
  robustness and error characterization are in
  [`eigen_pivot001_validation.md`](eigen_pivot001_validation.md),
  [`eigen_pivot_error_characterization.md`](eigen_pivot_error_characterization.md),
  and their correspondingly named result directories.
- Sequential SuperLU `SamePattern` was about 7.4% slower than Eigen overall;
  see [`superlu_factorization_probe.md`](superlu_factorization_probe.md) and
  `results/superlu_factorization_probe/`.
- Robust-natural UMFPACK was 22.8% slower than Eigen, and even its weaker
  default pivot policy did not win; see
  [`umfpack_factorization_probe.md`](umfpack_factorization_probe.md) and
  `results/umfpack_factorization_probe/`.
- Exact SEM static condensation was structurally valid and robust, but its
  complete intrinsic factor-and-solve kernel was effectively equal to LAPACK
  (`0.996x` LAPACK/condensed aggregate), not the required material win. See
  [`sem_structure_factorization_probe.md`](sem_structure_factorization_probe.md)
  and `results/sem_structure_factorization_probe/`.
- The generic sparse-LU library search is therefore stopped.

## Stopping decision

The direct frequency-domain backend optimization campaign is intentionally
stopped here. The retained design is a reliable runtime choice between the
robust/default sparse `EigenSparseLU` path and the faster `LapackBandLU` path
for adaptive multi-SEM calculations. Future major acceleration should use a
genuinely different algorithmic route rather than further tuning of these
direct LU pathways.
