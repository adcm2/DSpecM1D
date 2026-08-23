# Eigen SparseLU pivot-threshold 0.01 validation

This uncommitted validation tests Eigen 3.4.0 `SparseLU` with
`NaturalOrdering<int>` and `setPivotThreshold(0.01)`. It does not change the
production solver or its current threshold. It compares exactly threshold 1,
threshold 0.01, and the current pivoted LAPACK general-band factorization with
custom C++ backsolve.

## Physical stress set

The test deliberately reuses the accepted threshold-zero infrastructure and
does not perform another threshold sweep. Its catalogue-anchored PREM modes
are radial 0S0, toroidal 0T2/0T12, spheroidal 0S2/0S12, and fluid-solid
spheroidal 0S40. Dimensions span 81 through 1610. Every system uses the
nearest response-peak production-grid point and its two adjacent grid points,
attenuation both enabled and disabled, and production epsilon
`0.3229401443736645`, epsilon/10, and epsilon/100. The retained nine-point
threshold-one scans merely anchor each three-point set to its physical peak.

The RHS and receiver projections use the same production SEM source and all
available receiver/output quantities as the accepted threshold-zero test.
The raw 216 paired threshold rows record status, finiteness, normalized
backward residual, solution and transfer discrepancies, pivot count/fraction
and maximum displacement, U diagonal range, maximum off-diagonal L, maximum U,
maximum A, and `max|U|/max|A|`.

## Robustness result

All Eigen and LAPACK factorizations and solves succeeded, and all factors,
solutions, and outputs were finite. The important maxima are:

| diagnostic | threshold 0.01 | threshold 1 |
| --- | ---: | ---: |
| normalized backward residual | `1.5270140987379321e-16` | `1.5277938153621558e-16` |
| min `|diag(U)|` | `104.45362498194919` | `378.04835860160887` |
| max `|diag(U)|` | `5.6073327811816474e17` | `5.59749442525675e17` |
| max `|U|` | `8.889069633472128e17` | `8.878299165632896e17` |
| max `|U|/max|A|` | `0.9989493036290786` | `0.997738234636047` |
| max off-diagonal `|L|` | `99.86259004993806` | `0.999999998013822` |
| max nonidentity pivots | `331` | `1508` |
| max pivot fraction | `0.5152941176470588` | `0.9893617021276596` |
| max pivot displacement | `384` | `280` |

The threshold-0.01 L distribution has median `86.6839`, p95 `99.8582`, and
maximum `99.8626`; the corresponding threshold-one values are `0.999213`,
`1.000000`, and `1.000000`. This is substantial but controlled multiplier
growth at the scale implied by the 0.01 pivot rule: it remains below 100 over
the complete deliberately difficult set. It is qualitatively different from
the rejected threshold-zero maximum `5013.55`, which had no positive-threshold
bound. More importantly, the paired `max|U|/max|A|` ratio changes by at most
`1.001214x`; no isolated case has anomalous U growth. The worst L case is
spheroidal 0S2 at `0.274658203125` mHz, production epsilon, attenuation off.
The worst U/A case is the n=1610 fluid-solid 0S40 system at
`4.730224609375` mHz, epsilon/100, attenuation off.

The per-system maxima confirm that this conclusion is not hiding an isolated
resonance (`L001`/`L1` are the maximum L multipliers):

| system (n) | L001 | L1 | max U/A at 0.01 |
| --- | ---: | ---: | ---: |
| radial 0S0 (282) | `90.921` | `1.000` | `0.997738` |
| toroidal 0T2 (81) | `1.206` | `0.999` | `0.629071` |
| toroidal 0T12 (81) | `1.206` | `1.000` | `0.629074` |
| spheroidal 0S2 (425) | `99.863` | `1.000` | `0.997738` |
| spheroidal 0S12 (425) | `82.714` | `0.999` | `0.998429` |
| fluid-solid 0S40 (1610) | `96.304` | `1.000` | `0.998949` |

Maximum relative solution and receiver-output discrepancies against threshold
one are `1.0429308544890226e-7` and `4.8084761287929855e-8`; against LAPACK
they are `8.803423864866869e-8` and `2.532400820662189e-8`. These are about
0.1 ppm or less and scientifically negligible for this comparison. Taken with
the bounded L distribution and essentially unchanged U/A growth, threshold
0.01 **passes** this focused robustness gate.

## Minimal production pilots

Only the authorized production points were run, with GCC 15.2 Release flags,
identical physics and driver, `OMP_DYNAMIC=FALSE`, and all external BLAS thread
counts fixed to one. Each executable performs one untimed warm-up and one
measured run. Order was reversed at the parallel point.

| fmax | OpenMP threads | threshold 1 | threshold 0.01 | speedup |
| ---: | ---: | ---: | ---: | ---: |
| 5 mHz | 1 | `14.923160606 s` | `13.941919085 s` | `1.07038x` |
| 30 mHz | 32 | `6.538730759 s` | `5.961899304 s` | `1.09675x` |

Both configurations were checked against the same LAPACK result within each
point. By the triangle inequality, the candidate/default maximum-output
difference is bounded by `4.34e-17` at 5 mHz and `9.36e-16` at 30 mHz.

The accepted source currently has a pre-existing non-profiling compile guard
around `patternRebuilt` declarations but not all uses. To make the requested
non-profiling Release pilots without changing the repository, both temporary
build overlays contain the same declaration-only workaround. The baseline
overlay is retained as `nonprofiling_compile_workaround.diff`; the candidate
overlay is retained as `threshold001_pilot_overlay.diff` and differs further
only by a small solver subclass whose constructor calls
`setPivotThreshold(0.01)`. This workaround is identical between candidates
and therefore does not confound the threshold timing comparison.

## Reproduction

The physical validator is built with the same include and link set as the
accepted threshold-zero validator, changing only the source and output names:

```bash
/opt/gcc-15.2.0/bin/g++ -DDSPECM1D_ENABLE_PROFILING \
  -DDSPECM1D_ENABLE_LAPACK_BAND_SOLVER -DHAVE_LAPACK_CONFIG_H \
  -DLAPACK_COMPLEX_CPP -I/tmp/dspecm1d-backend-on-gcc15-final/benchmarks \
  -I/space/adcm2/CleanBuilds/DSpecM1D \
  -I/home/adcm2/Documents/c++/DSpecM3D/build-docs/_deps/eigen3-src \
  -I/home/adcm2/Documents/c++/DSpecM3D/build-docs/_deps/fftwpp-src \
  -I/home/adcm2/Documents/c++/DSpecM3D/build-docs/_deps/numericconcepts-src/include \
  -I/home/adcm2/Documents/c++/DSpecM3D/build-docs/_deps/gaussquad-src \
  -I/home/adcm2/Documents/c++/DSpecM3D/build-docs/_deps/gshtrans-src \
  -I/home/adcm2/Documents/c++/DSpecM3D/build-docs/_deps/interpolation-src \
  -O2 -DNDEBUG -fopenmp -std=gnu++23 benchmarks/eigen_pivot001_validation.cpp \
  -o /tmp/eigen_pivot001_validation \
  /usr/lib/x86_64-linux-gnu/libfftw3.so /usr/lib/x86_64-linux-gnu/libfftw3f.so \
  /usr/lib/x86_64-linux-gnu/libfftw3l.so /opt/gcc-15.2.0/lib64/libgomp.so \
  /lib/x86_64-linux-gnu/libpthread.a /usr/lib/x86_64-linux-gnu/liblapacke.so \
  -lm -ldl /usr/lib/x86_64-linux-gnu/libopenblas.so
OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OPENBLAS_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 \
  /tmp/eigen_pivot001_validation data/models/prem.200.no.txt data/params/ex3.txt \
  > results/eigen_pivot001_validation/eigen_pivot001_validation_raw.tsv \
  2> results/eigen_pivot001_validation/eigen_pivot001_validation_physical_peak_scan.tsv
python3 scripts/eigen_pivot001_validation_analyze.py \
  results/eigen_pivot001_validation/eigen_pivot001_validation_raw.tsv \
  results/eigen_pivot001_validation/eigen_pivot001_validation_summary.json PASS
```

The exact commands, binary/header hashes, parameter hashes, raw pilot outputs,
thread settings, ordering, and result summary are retained in
`eigen_pivot001_validation_metadata.json` and
`eigen_pivot001_validation_pilots.tsv`.
