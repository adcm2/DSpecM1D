# LAPACKE high-level versus work-interface experiment

This focused experiment starts from revision
`711eb5d36c15475093d4b18793f3995d89d2265b` and makes no production solver
change.  It tests whether the high-level LAPACKE wrapper, rather than the
underlying `zgbtrs`, explains the previously observed concurrent tiny-solve
scaling loss.

## Matching LAPACKE source

The installed `/usr/lib/x86_64-linux-gnu/liblapacke.so.3.11.0` comes from
Debian `lapack 3.11.0-2`.  The exact Debian source archive has SHA-256
`4b9ba79bfd4921ca820e83979db76ab3363155709444a787979e81c22285ffa9` and
matches the Reference LAPACK `v3.11.0` release.  Source-file hashes are retained
in `lapacke_interface_metadata.json`.

Inspection confirms that `LAPACKE_zgbtrs` validates the layout, calls
`LAPACKE_get_nancheck()`, scans the factored band with
`LAPACKE_zgb_nancheck()`, scans the RHS with `LAPACKE_zge_nancheck()`, and then
calls `LAPACKE_zgbtrs_work`.  For column-major storage the work interface calls
the underlying Fortran `zgbtrs` directly and only adjusts a negative `info`.
The factorization wrappers have the analogous structure: `LAPACKE_zgbtrf`
NaN-checks the band before `LAPACKE_zgbtrf_work`, whose column-major branch
calls Fortran `zgbtrf` directly.

## Method

`lapacke_interface_probe` compares the high-level and work interfaces using
128 independent systems at each of these observed/representative sizes:

| n | kl=ku | RHS | repeated solves/system |
|---:|---:|---:|---:|
| 17 | 4 | 3 | 512 |
| 256 | 8 | 3 | 32 |
| 1554 | 15 | 3 | 1 |

Each point uses one untimed warm-up and one measured run at 1, 10, 32, and 50
outer OpenMP threads.  OpenBLAS remains fixed at one thread.  Both solve
interfaces receive identical pre-factorized buffers and identical RHS data.
The secondary factorization comparison starts from identical band matrices.
All resulting solutions, factors, and pivots agree exactly (maximum absolute
difference zero).

## Results

`work speedup` is high-level time divided by work-interface time:

| n | threads | high-level solve (s) | work solve (s) | work speedup |
|---:|---:|---:|---:|---:|
| 17 | 1 | 0.155726 | 0.140854 | 1.106x |
| 17 | 10 | 0.694210 | 0.715662 | 0.970x |
| 17 | 32 | 0.693885 | 0.683851 | 1.015x |
| 17 | 50 | 0.671358 | 0.659630 | 1.018x |
| 256 | 1 | 0.174588 | 0.146884 | 1.189x |
| 256 | 10 | 0.049684 | 0.053072 | 0.936x |
| 256 | 32 | 0.043132 | 0.042601 | 1.012x |
| 256 | 50 | 0.041756 | 0.041981 | 0.995x |
| 1554 | 1 | 0.057096 | 0.065763 | 0.868x |
| 1554 | 10 | 0.007932 | 0.006912 | 1.147x |
| 1554 | 32 | 0.002794 | 0.002879 | 0.970x |
| 1554 | 50 | 0.010092 | 0.002729 | 3.698x |

The 50-thread large-call point shows sensitivity to the single exploratory
sample, but it does not affect the decision gate: large calls already scaled
well, while the production problem is millions of tiny calls.  At `n=17`,
the work interface changes concurrent throughput by only -3.0% to +1.8% and
does not remove the roughly fourfold slowdown relative to the one-thread call
stream.  Factorization work-interface differences are secondary and do not
alter the earlier conclusion that factorization scales well; the complete raw
table is retained in `results/lapacke_interface_experiment/`.

## Decision

`LAPACKE_zgbtrs_work` does not materially improve the concurrent tiny-system
solve path, so the conditional 35 mHz production pilot was not run and no
production call was changed.  The NaN checks have measurable serial overhead
at some sizes, but they are not the cause of the outer-concurrency collapse.

The next experiment, only after human approval, is the previously proposed
custom small-band-LU backsolve.  A size threshold should be measured rather
than assumed.  The optional production size histogram was not added because
the retained profiler records only the 17..1554 extrema; obtaining the bins
would require new per-system production instrumentation, outside this focused
interface experiment.
