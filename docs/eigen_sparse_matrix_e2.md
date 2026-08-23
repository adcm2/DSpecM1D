# E2: fixed-pattern Eigen sparse frequency matrices

This experiment targets matrix formation only in the adaptive multi-SEM Eigen
path.  It does not change ordering, solver cadence, `nskip`, truncation,
OpenMP decomposition, or the LAPACK backend.

## Source audit

For each SEM and angular degree, `FullSpecMultiSem.h` obtains sparse H, P and
Ha matrices from the SEM and casts them to complex Eigen sparse matrices.
Each base matrix is compressed once.  The previous Eigen frequency path then
formed a trailing sparse block (for toroidal and spheroidal modes), evaluated
`H - w*w*P`, optionally added `chi*Ha`, and compressed the result at every
frequency.  Those expressions allocate sparse temporaries, merge the
independent nonzero streams, and may reallocate while assembling the result.
Radial matrices are full-size and used for every frequency; toroidal and
spheroidal active dimensions follow `allIndicesTor`/`allIndicesSph` and the
existing internally derived cadence.

The E2 path constructs a compressed union pattern once per worker, active
truncation, and SEM chunk.  Its stored pattern is

```text
attenuation on:  pattern(H) union pattern(P) union pattern(Ha)
attenuation off: pattern(H) union pattern(P)
```

The coefficient arrays are copied during that setup only.  Frequency updates
write `valuePtr()` directly; no insertion, sparse merging, structural search,
or compression occurs in the hot path.  A changing `ridx` rebuilds the small
per-worker object for that active truncation.  The existing `compute()` /
`factorize()` cadence and solver ownership are unchanged.

## Kernel result

The benchmark uses GCC 15.2 Release code, one thread, one warm-up and one
measured batch of 64 frequencies.  It uses real PREM radial, toroidal,
spheroidal, and fluid-solid spheroidal SEM systems.  The raw rows are in
`results/eigen_sparse_matrix_e2/eigen_sparse_matrix_e2_raw.tsv`.

All rows retained identical compressed outer and inner index arrays.  The
direct path was 2.22--2.69 times faster without attenuation and 4.20--5.27
times faster with attenuation.  The off/on union counts are recorded
separately; for these systems they coincide.
The coefficient comparison is exact at the stored double values in the
benchmark, every explicit `analyzePattern` call completed, all subsequent
`factorize` and `solve` status values are zero for both pathways, and all
residuals are below `1.1e-16`. Eigen does not define `info()` as an analysis
status immediately after `analyzePattern()`, so success is confirmed by its
return followed by successful factorization.

## Production pilots

After the kernel result, the production path was tested at the two approved
points only: 5 mHz with one thread and 30 mHz with 32 threads.  Each backend
used one untimed warm-up and one measured `spectra()` call, with
`OMP_DYNAMIC=FALSE`, `OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, and
`BLIS_NUM_THREADS=1`.  The current NaturalOrdering path was retained as the
control.  The 5 mHz result was 14.999127575 s for E2 versus 17.257272809 s for
the control; the 30 mHz result was 6.537373321 s versus 7.518734013 s.  The
reported spectrum norms agree to the displayed precision.  The end-to-end
result is therefore approximately 13--15% faster at these points while
formation alone is 2--6 times faster; further production optimization is
deferred.

The implementation intentionally leaves future pattern caching beyond the
current per-worker/per-ridx workspace, and any Level-2 direct operator
assembly, out of scope.
