# Level-1 Eigen SparseLU experiment

This benchmark is a diagnostic artifact only. It does not change the
production Eigen pathway. `benchmarks/eigen_sparselu_e1.cpp` builds real SEM
frequency matrices from PREM, uses one untimed warm-up and one measured run for
each of four public `Eigen::SparseLU` configurations, and writes separate
analysis, factorization, and solve timings.

## Reproduction

The source was built with GCC 15.2 and the LAPACK-enabled Release build. From
the repository root:

```bash
cmake --build /tmp/dspecm1d-backend-on-gcc15-final \
  --target eigen_sparselu_e1 -j4
OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OPENBLAS_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 \
  /tmp/dspecm1d-backend-on-gcc15-final/bin/eigen_sparselu_e1 \
  /tmp/dspecm1d-backend-on-gcc15-final/data/models/prem.200.no.txt \
  > results/eigen_sparselu_e1/eigen_sparselu_e1_raw.tsv
python3 scripts/eigen_sparselu_e1_analyze.py \
  --raw results/eigen_sparselu_e1/eigen_sparselu_e1_raw.tsv \
  --output-dir results/eigen_sparselu_e1 \
  --model /tmp/dspecm1d-backend-on-gcc15-final/data/models/prem.200.no.txt \
  --source benchmarks/eigen_sparselu_e1.cpp
```

The systems use PREM's fluid-solid mesh, three genuine SEM resolutions
(`maxstep` 0.15, 0.08 and 0.02), radial/toroidal/spheroidal operators, and a
separate large fluid-solid spheroidal row plus a high-resolution spheroidal
system (`maxstep` 0.01) near the observed production maximum of 1554 active
unknowns. The RHS has one, two, or four
columns for radial, toroidal, or spheroidal systems respectively. Every
configuration sees the same matrix and RHS for a given row.

The four configurations are:

1. `COLAMD`, `isSymmetric(false)` (the current default);
2. `COLAMD`, `isSymmetric(true)`;
3. `NaturalOrdering`, `isSymmetric(false)`;
4. `NaturalOrdering`, `isSymmetric(true)`.

`isSymmetric(true)` is only the SparseLU sparsity-pattern hint. The numerical
matrix remains complex symmetric, and SparseLU is used in all cases; no
Hermitian or symmetric numerical solver is selected.

## Current production audit

The audit was made against commit `d2d435d`.

Single-SEM (`FullSpecSingleSem.h`):

- radial constructs a new sparse `A(w)`, compresses it, and calls `solver.compute`
  at every frequency; no direct `analyzePattern` or factorize-only call is used;
  its active dimension and sparsity pattern are fixed for the reused SEM, but
  `compute` still repeats symbolic analysis at every frequency;
- toroidal and spheroidal construct/compress a new sparse `A(w)` at every
  frequency and call `factorizeOrCompute`; the first position and every `nskip`
  position call `compute`, while intermediate positions call `factorize`;
- `nskip` comes from `SpectraRunContext`, clamped to at least one;
- toroidal/spheroidal `ridx` and active dimension can change by frequency, so
  the active principal matrix and often its structure change;
- the radial loop does not use the cadence helper and analyzes repeatedly;
- toroidal/spheroidal recomputed `ridx` values can remain identical for adjacent
  frequencies, so a cadence `compute` can reanalyze an unchanged active
  pattern; between cadence points `factorize` correctly reuses the prior
  symbolic state;
- each OpenMP worker owns a private SparseLU object (`private(solver)` or
  `private(solver1)`), so symbolic state is not shared between workers. In the
  single-SEM path this is per worker and mode loop; in the adaptive path the
  radial worker pool owns one solver and the toroidal/spheroidal degree/chunk
  loops own another worker-private solver. There is no cross-worker or
  cross-degree symbolic reuse.

Adaptive multi-SEM (`FullSpecMultiSem.h`):

- radial constructs/compresses `A(w)` and calls `eigenSolver.compute` for every
  frequency;
  its reused SEM gives a fixed matrix dimension/pattern, but the current call
  nevertheless repeats `analyzePattern` for every frequency;
- toroidal and spheroidal call `factorizeOrCompute`, using an internally
  derived `nskip = max(1, (i2-i1)/20)`; `compute` therefore includes
  `analyzePattern + factorize` at cadence points and `factorize` alone is used
  between them;
- `ridx` and the active dimension change with frequency in toroidal and
  spheroidal paths, while each chunk reuses its SEM and base operators. A
  recomputed `ridx` may also be identical to the previous one, so cadence
  `compute` can reanalyze an unchanged active pattern; intermediate
  `factorize` calls reuse the worker's symbolic state;
- radial has one worker-private solver and toroidal/spheroidal have another
  worker-private solver (`private(eigenSolver)` / `private(eigenSolver1)`), so
  every worker independently repeats symbolic analysis when `compute` occurs;
- no production path explicitly caches a symbolic analysis across changed
  active patterns. In particular, sparse `A(w)` is rebuilt from expressions and
  compressed at each frequency before the solver call.

Within a toroidal/spheroidal adaptive loop, the worker-private solver is reused
through the reverse frequency traversal for a given degree and chunk. It is not
shared across OpenMP workers, degrees, or chunks, and a new `compute` cadence
restarts symbolic analysis even when the active pattern happens not to change.

Thus, for accounting, `compute(A)` is `analyzePattern(A)` plus
`factorize(A)`. The current code does not call `analyzePattern` directly in the
spectra loops.

The focused kernel measured aggregate COLAMD `analyzePattern` time of
approximately 1.23 ms versus 7.78 ms for factorization (15.8%), while
NaturalOrdering measured 0.423 ms versus 6.07 ms (7.0%). Symbolic analysis is
therefore not dominant, but it is not negligible when the radial path repeats
it for every frequency. Caching one symbolic analysis per unchanged
worker/pattern could remove roughly this measured 7--16% component from such a
workload. For toroidal/spheroidal paths, where analysis occurs only at cadence
points and `ridx` can change, the likely benefit is smaller and must be weighed
against pattern-change detection; no caching is implemented in E1.

## Results and decision

The raw TSV records `nnz(A)`, `kl`, `ku`, `nnz(L)`, `nnz(U)`, and the dense
LAPACK general-band storage envelope
`(2*kl + ku + 1)*n`. It also records sparse-factor/envelope and sparse-factor/A
ratios, separate timings, all three stage statuses, residual, and discrepancy
from the default COLAMD solution. The high-resolution spheroidal row contains
1610 unknowns, close to the observed production maximum of 1554. The JSON
summary is the machine-readable result. The envelope is a storage upper bound, not a claim that the current
Eigen pathway uses LAPACK storage.

Natural ordering substantially reduces symbolic-analysis time on these systems,
has lower mean factor fill, and lowers factorization and solve time in the
focused experiment. It produces a small (below the validation threshold)
solution difference relative to COLAMD. The symmetry hint does not provide a
consistent additional win. NaturalOrdering was therefore checked in exactly
the amended two-point production pilot: one warm-up and one measured run at
5 mHz/one thread and 30 mHz/32 threads. The pilot is diagnostic only; no
production source was changed.

Pilot provenance is retained in `eigen_sparselu_e1_metadata.json`: exact
executable paths and hashes, parameter paths and hashes, all thread-control
variables, commands, the common pilot-source hash, and the one-line
`natural_ordering_override.diff`. The default pilot is the committed
`FullSpecMultiSem.h` build. The NaturalOrdering pilot is a second compile-time
binary from the same source/configuration, using only a temporary copy whose
`EigenMultiSemSolver` alias changes from `COLAMDOrdering<int>` to
`NaturalOrdering<int>`; that override is not a committed production change.

Across the focused rows, the sparse factors occupy a smaller storage count than
the corresponding general-band envelope, but the envelope comparison alone does
not establish solve work: sparse factors can skip structural zeros and use
supernodes. The lower Eigen solve time is therefore consistent with its sparse
factor structure, while the focused experiment does not prove that structure is
the sole cause.

This experiment also leaves Level 2 deferred. Today `A(w)` is formed by sparse
Eigen expressions (`H - w*w*P`, plus attenuation), followed by compression.
The later Level-2 opportunity is a fixed compressed exact sparsity pattern with
direct numerical-value updates, avoiding repeated sparse-expression merge,
allocation, and compression. It is not implemented here.
