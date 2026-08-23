# UMFPACK factorization probe

This benchmark-only experiment compares current Eigen
`SparseLU<NaturalOrdering<int>>` with pivot threshold 1.0 against UMFPACK's
repeated complex-double numerical factorization. It does not change production
solver behavior, CMake dependencies, backend selection, or spectra code.

## Installation audit and provenance

No system UMFPACK package was registered in `pkg-config` or the loader cache,
and no header or library was found below `/usr` or `/opt`. A usable user-local
installation was subsequently found at `/home/adcm2/space/SuiteSparse` and is
the sole UMFPACK implementation used by this experiment:

- SuiteSparse_config 7.12.2 and UMFPACK 6.3.7;
- installed headers in `/home/adcm2/space/SuiteSparse/include/suitesparse`;
- installed libraries in `/home/adcm2/space/SuiteSparse/lib`;
- complex-double `umfpack_zi_*` and 64-bit-index `umfpack_zl_*` symbols;
- Release build with GCC 15.2.0, CHOLMOD enabled, and LP64 system OpenBLAS;
- source checkout commit `16178b1eca3a16ed44959451c3354062d66451fb`,
  tagged `v7.12.2-beta1`, whose commit message is `SuiteSparse 7.12.2`.

The source checkout contains build products and two unrelated modified SPQR GPU
headers. UMFPACK, AMD, and SuiteSparse_config sources are unmodified. Version,
defaults, and API claims below are confirmed from the installed 6.3.7 header
and its matching local UMFPACK source. Exact library hashes are retained in the
metadata. An official 7.12.3 benchmark-local build under `/tmp` was tested
during availability diagnosis but was not linked into the probe and contributes
no timings.

## Exact interface and policy

The probe uses the int32 complex-double interface:
`umfpack_zi_symbolic`, `umfpack_zi_numeric`, `umfpack_zi_solve`,
`umfpack_zi_get_lunz`, and the matching free routines. Eigen's compressed
column-major `StorageIndex` is `int`, so `outerIndexPtr()` and
`innerIndexPtr()` are passed directly with no index conversion. UMFPACK's
documented packed-complex form (`Az == nullptr`) is interleaved real/imaginary
doubles, which matches standard `std::complex<double>` array access. Thus the
values are also zero-copy. Conversion time is exactly zero; no conversion or
packing is hidden in numerical-factorization timing.

One Symbolic object is computed per system, ordering, and pivot policy, then
reused for both frequency matrices with the same checked CSC pattern. Each
measured `umfpack_zi_numeric` creates fresh Numeric factors; freeing the prior
Numeric object is outside the timed call. This matches UMFPACK's documented
Symbolic-reuse interface and the production fixed-pattern question.

The two orderings are exactly:

- `UMFPACK_ORDERING_NONE` (natural/no fill-reducing preordering);
- the installed default `UMFPACK_ORDERING_AMD`.

The default `UMFPACK_STRATEGY_AUTO` selected the symmetric strategy (value 3)
for every tested pattern. The installed header states that this strategy does
not refine the column preordering during numeric factorization. Natural order
therefore remains natural apart from UMFPACK's structural singleton/front
machinery; AMD is its normal default fill-reducing order.

The primary robust policy retains default row-sum scaling and AUTO strategy,
but sets both pivot controls to 1.0. In UMFPACK 6.3.7:

- the general threshold accepts a pivot only when its magnitude is at least
  that threshold times the largest entry in the column;
- under the symmetric strategy, the diagonal is first accepted against the
  symmetric threshold, otherwise the general threshold governs an
  off-diagonal pivot.

Consequently `PIVOT_TOLERANCE=1.0` and `SYM_PIVOT_TOLERANCE=1.0` require a
column-maximum pivot whether it is diagonal or off diagonal. This is the
conservative primary comparison, not a speed-tuned policy.

The cheap secondary diagnostic retains the installed defaults: general pivot
threshold 0.1 and symmetric diagonal threshold 0.001. It is numerically
validated separately and is not treated as a robust production recommendation.
The raw solve sets `UMFPACK_IRSTEP=0`; the separately timed refined solve uses
the installed default maximum of two refinement steps. Three right-hand sides
are solved individually because `umfpack_zi_solve` is a vector interface.

`Info` records the actual strategy, ordering, row scaling, off-diagonal pivot
count, flops, Numeric size, peak memory, reciprocal-condition estimate, and
refinement steps. `umfpack_zi_get_lunz` supplies actual L and U nonzero counts.
UMFPACK memory values are retained in its documented internal `Unit` measure,
not mislabeled as bytes.

## Physical matrices and timing

The accepted seven-system real-PREM set covers radial, toroidal,
spheroidal small/medium/large, fluid-solid spheroidal, and the previous
maximum-size system. Dimensions range from 37 to 1,626. Every matrix has
attenuation enabled and is formed by the accepted E2
`FixedPatternFrequencyMatrix` path:

```text
compressed H/P/Ha union pattern
-> direct valuePtr update of A(w)
```

Each system uses its accepted physical near-resonance production-grid point
and the adjacent grid point. Eigen, UMFPACK, and current pivoted LAPACK receive
the identical CSC values and deterministic three-column RHS. The full solution
and a three-location transfer sample are compared.

The retained run fixes all external BLAS libraries to one thread and uses one
untimed warm-up followed by one measured batch of 16 operations per matrix.
Symbolic/setup, numerical factorization, Eigen/UMFPACK solves, and LAPACK
reference validation are separate. An initial diagnostic run exposed first-use
library initialization in the first Symbolic timing; it was discarded, an
explicit Symbolic warm-up was added, and the complete short probe was rerun
once. Only the coherent rerun is retained.

## Results

The aggregate of 14 per-matrix mean times is:

| Pivot policy | Ordering | Eigen factor | UMFPACK Numeric | Eigen/UMFPACK | UMFPACK slower |
|---|---|---:|---:|---:|---:|
| robust 1.0/1.0 | natural | 0.00832263 s | 0.01021990 s | 0.81436x | 22.80% |
| robust 1.0/1.0 | default AMD | 0.00831913 s | 0.01221849 s | 0.68086x | 46.87% |
| default 0.1/0.001 | natural | 0.00830638 s | 0.00916507 s | 0.90631x | 10.34% |
| default 0.1/0.001 | default AMD | 0.00832175 s | 0.00971598 s | 0.85650x | 16.75% |

Natural ordering is the better UMFPACK result, but robust UMFPACK remains
22.8% slower than Eigen. Default, substantially weaker diagonal pivoting is
also slower. On the n=1,626 system, robust-natural Eigen is about 5.5% faster
at both frequencies. UMFPACK raw solves are roughly comparable to Eigen, while
default refinement took one step for each RHS and made the aggregate solve
about ten times more expensive without changing an already tiny residual in a
scientifically useful way.

One-time Symbolic cost summed over the seven structural patterns is
`0.00147924` s for robust natural order and `0.00134838` s for robust default
AMD. It is reported separately and never charged to repeated Numeric calls.

All 56 configuration rows report successful Symbolic/Numeric/solve statuses
and finite outputs. Maximum normalized backward residual is
`8.03517e-17` for Eigen, `8.62657e-17` for raw UMFPACK,
`8.30711e-17` for refined UMFPACK, and `1.14205e-16` for LAPACK. Maximum
UMFPACK/Eigen full-solution discrepancy is `1.48931e-8`, and maximum retained
transfer-sample discrepancy is `1.57866e-8`, both at the deliberately sensitive
small-spheroidal resonance. All solvers have tiny backward residuals there;
the difference is consistent with sensitivity and distinct robust pivots.

Decision: **C — robust UMFPACK is slower.** It does not close the current
Eigen/LAPACK factorization gap, and even its relaxed defaults do not beat Eigen.
Per the performance gate, no outer-OpenMP concurrency probe and no production
spectrum were run. The external sparse-LU library search should stop here.
Retain the current Eigen and LAPACK backends; any further factorization work
would need a structure-aware bespoke approach with a clearly justified scope.

## Reproduction

The exact compiler and link command is retained in metadata. Run the resulting
binary with:

```bash
LD_LIBRARY_PATH=/home/adcm2/space/SuiteSparse/lib:/opt/gcc-15.2.0/lib64 \
OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OPENBLAS_NUM_THREADS=1 \
MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 \
/tmp/umfpack_factorization_probe data/models/prem.200 \
  > results/umfpack_factorization_probe/umfpack_factorization_probe_raw.tsv \
  2> results/umfpack_factorization_probe/umfpack_factorization_probe.log

python3 scripts/umfpack_factorization_probe_analyze.py \
  results/umfpack_factorization_probe/umfpack_factorization_probe_raw.tsv \
  results/umfpack_factorization_probe/umfpack_factorization_probe_summary.json \
  results/umfpack_factorization_probe/umfpack_factorization_probe_summary.tsv
```
