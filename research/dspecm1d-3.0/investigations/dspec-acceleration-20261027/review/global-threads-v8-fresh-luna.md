# Fresh bounded review: combined global v8 OpenMP resource cleanup

**Result: PASS for static source/API and launch-plan review; runtime cleanup remains untested.** No SSH, build, or science run was performed.

The v7→v8 source diff is confined to `global_kernel_consumer.hpp` and the implementation status log. The new call occurs after shared grid and worker contexts are constructed and before `runOwnedTiles` starts the outer worker threads. It runs only when `gshThreads > 1`; a nonzero `omp_pause_resource_all(omp_pause_hard)` return throws before worker execution. The sidecar records the cleanup policy, attempted flag, and success flag; these are set true only on the path after a zero return. The existing runtime-thread-limit check still rejects an undersized OpenMP limit before setup.

The pinned local GCC 13 `omp.h` declares `omp_pause_resource_all(omp_pause_resource_t)` and `omp_pause_hard`; the parent independently confirmed the pinned remote GCC 13.2 header has the same declarations. The call is guarded by `_OPENMP`, and the v8 build configuration retains the OpenMP build setup. The one-thread/default path does not call the pause API and keeps its prior behavior. No numerical loops or transform calls changed in the v8 source delta, so arithmetic order and the exact comparison gates are unchanged. The prior resident-pool risk is addressed structurally: the main-thread grid-construction pool is hard-paused before worker threads can create their own OpenMP teams. This is a source-level control, not measured runtime proof.

The v8 build configuration uses a fresh `kernel-acceleration-20261027-global-threads-02` root, one profiled variant, and the staged combined-global-v8 source. Its build plan contains one job with an 8-thread build cap. Run-v3 contains 19 jobs; after substituting the new binary path, it is structurally identical to run-v2. All producer environments retain `OMP_NUM_THREADS=32`, `OMP_THREAD_LIMIT=32`, `OMP_DYNAMIC=FALSE`, `OMP_MAX_ACTIVE_LEVELS=1`, and the 40-thread plan cap; all nine comparisons retain `--exact`. The previous v2 review now carries the requested HOLD note: the active product of 32 alone did not bound resident OpenMP pools, and v2 was not cleared to launch.

## SHA-256

- `candidates/combined-global-v8/global/global_kernel_consumer.hpp`: `f544304569f0052dd9df73ff79045ee404be862813e45c9c9681a742e603ed3d`
- `candidates/combined-global-v8/implementation_status.md`: `17b5da3a7f6e7d9e34e2e6b4bc9e8105ff991fab87af2e14a5ee5f0085ba8abb`
- `configs/global-threads-v8-build.json`: `07b84c21cdc930146674b3b7fe04b7c2b1b336ac8747b2d73d825a4fb3f65ece`
- `configs/global-threads-v8-build-plan.json`: `18eb4f0d156275b3d2dc8006ccd4b2ad485331cbda56ac2208455e3be6bf9a91`
- `configs/global-threads-run-v3-plan.json`: `37b5736bf9ab82b82f5decaf4831766908f64e0a1d28fe5f8f56f692685f6ed5`
- local GCC 13 `omp.h`: `4c4add3bee649c66078df2031941826830afd495591dfc4fcc8735a6d0a0b2df`
- `candidates/combined-global-v7/global/global_kernel_consumer.hpp`: `5555a444b9ff0ec76555ad4fdec51c1dbb0d25ad06bfd8d776bb65c1c182dfd3`
