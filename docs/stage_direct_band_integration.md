# Direct-band adaptive multi-SEM integration

The optional LAPACK production path now converts the existing sparse SEM
operators `H`, `P`, and `Ha` once per degree and frequency chunk into immutable
matrices with one common LAPACK general-band shape. Each frequency constructs
the active coefficient columns of the solver-owned band workspace directly,
then factorizes the same storage in place. The top `kl` LAPACK fill rows are
not cleared or included in the coefficient expression.

This stage deliberately retains sparse SEM assembly as the trusted source of
the base operators. It removes per-frequency sparse `A(w)` construction,
compression, bandwidth scanning, and sparse-to-band packing, but it does not
yet assemble `H`, `P`, or `Ha` directly in band form.

The Eigen SparseLU pathway remains available in every build. The
`DSPECM1D_ENABLE_LAPACK_BAND_SOLVER` option adds LAPACK capability without
replacing Eigen, and the preferred adaptive API defaults to Eigen:
`params.setSolverBackend(SPARSESPEC::SolverBackend::LapackBandLU)` selects the
direct-band path when that capability is enabled. The preferred reused/single-
SEM API remains Eigen-only for now and reports this limitation if the LAPACK
backend is requested.

## Performance pilot

The retained pilot compares the pre-integration sparse-pack LAPACK executable
with the direct-band executable at 20 mHz on the same machine. Both use the
same Release build, GCC 15.2, parameter file, two OpenMP threads, and serialized
BLAS settings. Each executable performs one untimed warm-up followed by three
measured complete `spectra()` calls.

| LAPACK pathway | Median wall time (s) | Dynamic matrix thread-time (s) | Band packing thread-time (s) | Band packs |
|---|---:|---:|---:|---:|
| Sparse pack | 64.423 | 24.964 | 21.166 | 983155 |
| Direct band | 44.380 | 6.138 | 0.000 | 0 |

The direct-band pilot is 1.452 times faster, a 31.1% wall-time reduction.
Combining the attenuated expression into one assignment reduced its candidate
dynamic-matrix median by a further 10.5%. Factorization and solve medians remain
stable, while the maximum measured output-norm difference is
`2.168404344971009e-19`. The raw selected timings and full reproducibility
metadata are retained under `results/direct_band_integration`.
