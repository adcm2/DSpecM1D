# Sequential SuperLU factorization probe

This benchmark-only experiment compares the current Eigen
`SparseLU<NaturalOrdering<int>>` with pivot threshold 1.0 against sequential
SuperLU 7.0.1. It does not alter either production backend, CMake dependency
selection, or solver API.

## Library audit and interface semantics

No sequential SuperLU header, library, package-config file, or Debian package
was available in the searched system, `/opt`, user, or temporary locations.
The experiment therefore built the official stable `v7.0.1` tag, commit
`83a3c82b9e44c966ca2bc26991ffdb2bf63e1e2b`, under `/tmp`. The Release build
uses GCC 15.2.0, enables the complex-double (`z*`) routines, and links the
system OpenBLAS. Nothing was installed system-wide or added to production.

The exact 7.0.1 headers and `zgssvx.c`, `zgstrf.c`, and `sp_preorder.c` sources
establish these modes:

- `DOFACT` computes the factorization from scratch, including the column
  permutation and elimination tree.
- `SamePattern` reuses `perm_c` and the column elimination tree. `perm_r`
  remains an output, so `zgstrf` performs fresh threshold row pivoting.
- `SamePattern_SameRowPerm` additionally attempts to reuse the prior row
  permutation and factor storage. It was deliberately not tested because
  values change with frequency and fresh numerical pivoting is required.

The probe calls complex-double `zgstrf` and `zgstrs`, with `ColPerm=NATURAL`,
`DiagPivotThresh=1.0`, no equilibration, and no iterative refinement. SuperLU's
public `NATURAL` mode begins with identity column order; `sp_preorder` may then
apply the elimination-tree postorder required by SuperLU. No generic COLAMD or
minimum-degree heuristic is used.

## Matrix set and timing

Seven attenuation-enabled real PREM SEM systems cover radial, toroidal,
spheroidal small/medium/large, fluid-solid, and the prior maximum-size case.
Their dimensions range from 37 to 1,626. Each uses the accepted E2
`FixedPatternFrequencyMatrix`: one compressed union pattern followed by direct
`valuePtr()` updates. Two production-grid frequencies are used per system: the
nearest grid point to its catalogue mode and the adjacent point.

Eigen and SuperLU receive numerically identical CSC matrices and the same
three-column deterministic RHS. SuperLU's required `std::complex<double>` to
`doublecomplex` value copy is timed separately and excluded from numerical
factorization. Index conversion and `sp_preorder` are one-time structural
setup. The retained run uses one untimed warm-up and one measured 16-operation
batch per point, with external BLAS threads fixed to one.

## Results

The aggregate of the 14 per-matrix mean factorization times is:

| Kernel | summed mean time |
|---|---:|
| Eigen factorize | 0.00838122 s |
| SuperLU `SamePattern` | 0.00900316 s |

Thus Eigen/SuperLU is `0.93092x`: SuperLU is about 7.4% slower overall. The
per-system factor ratio ranges from `0.81980x` to `1.09222x`; the largest
system is also slower in SuperLU (`0.82975x`). SuperLU's first `DOFACT` is
slower still, as expected. This does not materially close the observed
Eigen-to-LAPACK production factorization gap.

All 14 Eigen, SuperLU, and pivoted LAPACK reference factorizations and solves
succeed and remain finite. Maximum normalized backward residual is `8.03517e-17`
for Eigen, `7.54285e-17` for SuperLU, and `1.14205e-16` for LAPACK. The maximum
full-solution discrepancy is `3.16709e-8` on the accepted near-resonance small
spheroidal peak. At that point Eigen/LAPACK differs by `3.76891e-9` and
SuperLU/LAPACK by `3.48290e-8`, while all three backward residuals are below
`9e-21`. The matrix is solution-sensitive, but each solver is backward stable;
the disagreement is consistent with different robust pivot choices rather
than factorization or solve failure. Since SuperLU also fails the performance
gate, this is not a production-adoption claim.

Decision: **B**. Sequential SuperLU `SamePattern` is not materially faster, so
the added dependency is not justified. Per the decision gate, no 32-thread
concurrency probe and no production spectra were run. A similarly small
UMFPACK probe is the next reasonable sparse-factorizer experiment, but is not
started here.

## Reproduction

Build SuperLU outside the repository:

```bash
cmake -S /tmp/superlu-7.0.1-src -B /tmp/superlu-7.0.1-build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=/tmp/superlu-7.0.1-install \
  -DCMAKE_C_COMPILER=/opt/gcc-15.2.0/bin/gcc \
  -DCMAKE_Fortran_COMPILER=/opt/gcc-15.2.0/bin/gfortran \
  -Denable_tests=OFF -Denable_examples=OFF -DBUILD_SHARED_LIBS=ON
cmake --build /tmp/superlu-7.0.1-build -j 8
cmake --install /tmp/superlu-7.0.1-build
```

Compile `scripts/superlu_factorization_probe.cpp` with the exact command in
`superlu_factorization_probe_metadata.json`, then run:

```bash
OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OPENBLAS_NUM_THREADS=1 \
MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 \
/tmp/superlu_factorization_probe data/models/prem.200 \
  > results/superlu_factorization_probe/superlu_factorization_probe_raw.tsv \
  2> results/superlu_factorization_probe/superlu_factorization_probe.log
python3 scripts/superlu_factorization_probe_analyze.py \
  results/superlu_factorization_probe/superlu_factorization_probe_raw.tsv \
  results/superlu_factorization_probe/superlu_factorization_probe_summary.json \
  results/superlu_factorization_probe/superlu_factorization_probe_summary.tsv
```
