# Focused LAPACK parallel-scaling diagnostic

This diagnostic was run from revision `fdbc2a79ab05df46f3a8cead5bf5f221c9661755` without changing production solver code.  It uses the exact final benchmark executable and the existing profiling-capable build.  Each production profile point has one untimed warm-up followed by one measured `spectra()` call.

The exact executables were `/tmp/dspecm1d-final-benchmark-gcc15/bin/final_solver_benchmark` and `/tmp/dspecm1d-backend-on-gcc15-final/bin/stage_profile_solver`.  The two focused probes were built as `lapack_band_concurrency_probe` and `lapack_solve_allocation_probe` in `/tmp/dspecm1d-backend-on-gcc15-final/bin`.  Their absolute paths and SHA-256 hashes, along with every raw and normalized artifact path/hash, are retained in `lapack_scaling_diagnostic_metadata.json`.

## BLAS implementation and effective threads

Both executables load these SONAMEs, which resolve to:

```text
/lib/x86_64-linux-gnu/liblapacke.so.3 -> /usr/lib/x86_64-linux-gnu/liblapacke.so.3.11.0
/lib/x86_64-linux-gnu/libopenblas.so.0 -> /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.21.so
/lib/x86_64-linux-gnu/libblas.so.3 -> /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
/lib/x86_64-linux-gnu/liblapack.so.3 -> /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3
```

The loaded OpenBLAS is `0.3.21`, configured as:

```text
OpenBLAS 0.3.21 NO_LAPACKE DYNAMIC_ARCH NO_AFFINITY Cooperlake MAX_THREADS=64
```

`openblas_get_parallel()` returned `1` (pthread-threaded), `openblas_get_num_threads()` returned `1` under the benchmark controls, and `OPENBLAS_NUM_THREADS=1` was honored independently of `OMP_NUM_THREADS`.  A representative exact-executable run with `OMP_NUM_THREADS=32`, `OPENBLAS_NUM_THREADS=1`, `OMP_DYNAMIC=FALSE`, and the other BLAS thread variables set to one showed process `nlwp` changing from startup `1` to a steady `32`.  There were no hidden OpenBLAS worker threads.

## Category profile: fmax=35 mHz, lmax=750

The workload counts and numerical output were identical at both thread counts: 1,720,146 frequency systems, 1,720,146 factorizations, 1,720,146 solves, and 5,158,146 right-hand sides.  The observed dimensions were 17..1554 with `kl,ku` 4..15.

| OMP threads | wall (s) | total worker-s | base operators | dynamic A(w) | factorization | solve | projection | unclassified |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 | 31.5423 | 315.2025 | 3.1970 | 17.1423 | 125.7584 | 146.2970 | 1.3307 | 14.8266 |
| 32 | 22.7508 | 727.2588 | 3.4606 | 17.3710 | 126.5895 | 553.3806 | 1.5294 | 18.3056 |

Factorization worker time is flat (t32/t10 = 1.0066), while solve worker time grows 3.7826x.  Base-operator preparation, direct dynamic band construction, projection, and unclassified worker time also remain close.  The spheroidal solve contribution alone rises from 93.8931 to 346.1979 worker-s.  Thus the saturation is in the solve path, dominated by concurrent `zgbtrs` calls; factorization and the surrounding DSpecM1D work do not explain it.

## Focused band-call probes

`lapack_band_concurrency_probe` uses 128 independent `n=1554`, `kl=ku=15`, `nrhs=3` systems, matching the largest observed spheroidal dimensions and bandwidths.  It reports one warm-up and one measured run for each thread count.  The large-system calls scale well: `zgbtrs` is 0.0568 s at one thread, 0.0064 s at 10, 0.0022 s at 32, and 0.0020 s at 50.  `zgbtrf` shows the same healthy trend.

`lapack_solve_allocation_probe` uses 128 independent minimum-size systems (`n=17`, `kl=ku=4`, `nrhs=3`) and exactly 512 repeated solves per system.  This reproduces the small-call granularity problem:

| threads | allocate result (s) | reused result (s) |
|---:|---:|---:|
| 1 | 0.156624 | 0.155909 |
| 10 | 0.693127 | 0.695259 |
| 32 | 0.692078 | 0.694679 |
| 50 | 0.654029 | 0.649961 |

Allocation reuse changes the measured time by only a few percent and does not restore scaling.  The evidence points to the enormous number of small independent `zgbtrs` calls and their call-level/low-arithmetic-intensity overhead, not to Eigen `DenseMatrix result = rhs` allocation as the primary cause.

## Diagnosis and next step

OpenBLAS threading and oversubscription are not the cause.  The direct-band path exposes roughly 1.7 million small LAPACK solves to an outer OpenMP loop.  Large band calls scale, but minimum-size calls become slower as more workers contend for tiny calls, which matches the production solve-category collapse.  The next optimization should reduce call granularity—batch compatible right-hand sides/systems or otherwise fuse/reuse solve work—while preserving the one-thread LAPACK advantage.  Affinity/NUMA tuning is not warranted before that experiment because the focused call probe already reproduces the behavior without NUMA changes.

Raw and normalized machine-readable artifacts are in `results/lapack_scaling_diagnostic/`; the parser is `scripts/lapack_scaling_diagnostic.py`.  No full benchmark grid was rerun.
