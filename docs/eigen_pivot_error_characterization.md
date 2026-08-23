# Eigen pivot-threshold numerical-error characterization

This is a benchmark-only characterization at base `0952f8474d43880d91484b93742c98ca01e8f4d5`.
It does not alter the production solver, run production spectra, or authorize a
production threshold change.

## Scope

`benchmarks/eigen_pivot_error_characterization.cpp` reuses the accepted
`eigen_pivot001_validation.cpp` physical stress machinery. It compares Eigen
`SparseLU<NaturalOrdering<int>>` thresholds **1.0, 0.1, and 0.01**, plus the
current pivoted LAPACK general-band factorization with custom C++ backsolve.
The six n=0 PREM catalogue systems are radial 0S0, toroidal 0T2/0T12,
spheroidal 0S2/0S12, and fluid-solid 0S40. Each uses the nearest physical
response-peak production-grid point and adjacent points, attenuation on/off,
and production epsilon, epsilon/10, and epsilon/100. This produces 396 rows:
22 frequency points × 6 epsilon/attenuation combinations × 3 thresholds.

Radial 0S0 and spheroidal 0S2 also retain a deliberately small ±2 grid stencil.
They were selected because radial 0S0 has the largest receiver/output peak and
spheroidal 0S2 contains the prior worst threshold-0.01 transfer discrepancy.
The response class is recorded as resonance, flank, antiresonance, or
elsewhere. The accepted ±4 threshold-1 scan used to select each physical peak
is retained separately.

Each row records factor/solve status and finiteness, normalized backward
residual, full solution discrepancy, receiver/output transfer infinity norm and
absolute/ordinary-relative errors, pivot fraction/displacement, maximum L/U/A,
U/A, and minimum diagonal U. It also records the actual complex receiver/output
entry with largest threshold-1 output magnitude (receiver row/RHS column and
candidate, threshold-1, and LAPACK real/imaginary values). The documented local
response scale is the maximum threshold-1 transfer norm over that case's ±2
stencil; local normalized error is absolute transfer error divided by that
scale, so a nearly-zero response point is not judged against an artificially
small denominator. All three complex backend values are evaluated at the same
threshold-1-selected receiver/RHS index; they are not independently selected
per backend.

## Results

All 396 Eigen and LAPACK factorizations/solves succeeded and were finite. The
largest errors and growth values were:

| diagnostic | threshold 1 | threshold 0.1 | threshold 0.01 |
| --- | ---: | ---: | ---: |
| max normalized residual | 1.5278e-16 | 1.5270e-16 | 1.5270e-16 |
| max relative solution error vs 1 | 0 | 1.3395e-7 | 1.0429e-7 |
| max relative transfer error vs 1 | 0 | 7.3614e-8 | 4.8085e-8 |
| max local-scale normalized transfer error | 0 | 7.3614e-8 | 4.8085e-8 |
| max off-diagonal |L| | 1.0000 | 9.9885 | 99.8652 |
| min diagonal |U| (across rows) | 378.035 | 104.454 | 104.454 |
| max |U|/|A| | 0.997738 | 0.998769 | 0.998949 |
| max pivot fraction | 0.98936 | 0.71059 | 0.51529 |

Threshold 0.01 is numerically indistinguishable from threshold 1 at the
receiver/output level for this focused set (about 0.05 ppm worst local-scaled
error; full-solution error about 0.10 ppm). This is an empirical statement
about these cases, not a general guarantee. Threshold 0.1 has lower L growth
than 0.01, so it is safer as a factor-growth rule, but it is **not** safer in
the observed output error: its maximum output and solution discrepancies are
larger here.

The prior approximately `4.81e-8` discrepancy is resolved exactly. It is
`spheroidal_0S2`, degree 2, frequency `0.30517578125` mHz (selected
resonance peak), attenuation on, epsilon/100, threshold 0.01, receiver row 0,
RHS column 2. The candidate complex transfer is
`-2.0233923965729995e-14 + 4.171064613164199e-15 i`; threshold 1 is
`-2.023392486560636e-14 + 4.1710650020532174e-15 i`; LAPACK is
`-2.0233924210756806e-14 + 4.171064720808044e-15 i`. The transfer norm is
`4.179125858345093e-14`, absolute matrix-inf error is
`2.0095227875749966e-21`, ordinary/local relative error is
`4.8084761287929855e-8`, and max L is `98.85627689126515` for that row.
The local ±2 scale is `4.179126055221622e-14`; the complex component error is
`9.803122684569987e-22`. It is therefore a resonance-peak rounding-level
response difference, not an unexplained LAPACK failure or a hidden large
physical-output error.

The ±2 stencils retain peak shape and location. For production epsilon with
attenuation on, radial 0S0 has threshold-1 transfer norms
`[5.5821e-13, 8.9854e-13, 1.1527e-12, 1.1209e-12, 9.8752e-13]` at offsets
`[-2,-1,0,+1,+2]`; spheroidal 0S2 has
`[2.4512e-15, 3.0088e-15, 3.4535e-15, 3.1297e-15, 2.4569e-15]`.
All three thresholds preserve the offset-0 peak location in these stencils.

## Conservative conclusion

The focused physical evidence supports treating 0.01 as output-indistinguishable
from 1 for these modes, while showing materially larger L multipliers (up to
about 100). Threshold 0.1 reduces that multiplier to about 10 but does not
improve observed output error. Neither result merits a production change from
this benchmark alone. Keep the production threshold unchanged; broader
conditioning, spectra, and timing evidence would be required before any
adoption decision.

Raw TSV, machine-readable summary/metadata, peak scan, and selected-mode list
are in `results/eigen_pivot_error_characterization/`. Re-run the analyzer with:

The exact validation build command was:

```bash
/opt/gcc-15.2.0/bin/g++ -DDSPECM1D_ENABLE_PROFILING -DDSPECM1D_ENABLE_LAPACK_BAND_SOLVER -DHAVE_LAPACK_CONFIG_H -DLAPACK_COMPLEX_CPP -I/tmp/dspecm1d-backend-on-gcc15-final/benchmarks -I/space/adcm2/CleanBuilds/DSpecM1D -I/home/adcm2/Documents/c++/DSpecM3D/build-docs/_deps/eigen3-src -I/home/adcm2/Documents/c++/DSpecM3D/build-docs/_deps/fftwpp-src -I/home/adcm2/Documents/c++/DSpecM3D/build-docs/_deps/numericconcepts-src/include -I/home/adcm2/Documents/c++/DSpecM3D/build-docs/_deps/gaussquad-src -I/home/adcm2/Documents/c++/DSpecM3D/build-docs/_deps/gshtrans-src -I/home/adcm2/Documents/c++/DSpecM3D/build-docs/_deps/interpolation-src -O2 -DNDEBUG -fopenmp -std=gnu++23 benchmarks/eigen_pivot_error_characterization.cpp -o /tmp/eigen_pivot_error_characterization /usr/lib/x86_64-linux-gnu/libfftw3.so /usr/lib/x86_64-linux-gnu/libfftw3f.so /usr/lib/x86_64-linux-gnu/libfftw3l.so /opt/gcc-15.2.0/lib64/libgomp.so /lib/x86_64-linux-gnu/libpthread.a /usr/lib/x86_64-linux-gnu/liblapacke.so -lm -ldl /usr/lib/x86_64-linux-gnu/libopenblas.so
```

The benchmark writes numerical rows to stdout and its physical peak scan to
stderr:

```bash
OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OPENBLAS_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 \
  /tmp/eigen_pivot_error_characterization data/models/prem.200.no.txt data/params/ex3.txt \
  > results/eigen_pivot_error_characterization/eigen_pivot_error_characterization_raw.tsv \
  2> results/eigen_pivot_error_characterization/eigen_pivot_error_characterization_peak_scan.tsv
```

Then run the analyzer:

```bash
python3 scripts/eigen_pivot_error_characterization_analyze.py \
  results/eigen_pivot_error_characterization/eigen_pivot_error_characterization_raw.tsv \
  results/eigen_pivot_error_characterization/eigen_pivot_error_characterization_summary.json
```
