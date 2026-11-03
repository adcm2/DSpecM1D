# Fresh Luna review: combined global v7 thread controls

**Verdict: PASS for the reviewed v7 source and the retargeted run plan, with a scope note on dynamic-team metadata.** This was a bounded static review only; no source edits, SSH, build, or scientific run were performed.

## Findings

- The candidate keeps `--global-gsh-threads` defaulted to 1 and selects `GSHTrans::Execution::Sequential()` for that setting. Values above 1 fail closed when the pinned GSHTrans build lacks OpenMP support. The consumer checks the requested inner team against `omp_get_thread_limit()`, and the worker × inner-team product is bounded by 48 both before and after memory-based worker reduction.
- `GlobalGridPlan` passes the explicit team size to the GSHTrans execution policy. The reviewed transform path parallelizes independent theta/output rows, with per-thread workspaces; outer workers own disjoint side ranges. The inspected path introduces no cross-thread reduction or shared output-row writes. The runtime measurement-to-estimate policy remains intact.
- Scratch preflight conservatively charges each inner thread for generated Wigner TLS (zero in stored mode) and the retained FFT input/output buffers for the applicable batch shapes. Checked arithmetic multiplies that per-inner-thread amount by the requested team, adds it to per-worker workspace, and includes per-worker space in working-memory worker selection. It separately reports the per-worker and total TLS bounds.
- The build plan points at the v7 source, retains the reviewed OpenMP-enabled compiler command, and has a build-control cap of 8. The run plan points at the corresponding v7 build output, sets one outer worker and 32 GSH threads, exports `OMP_THREAD_LIMIT=32`, `OMP_DYNAMIC=FALSE`, and `OMP_MAX_ACTIVE_LEVELS=1`, and binds the run to CPUs 16–47. Its campaign cap is 40. With the documented app-main and supervisor allowance, the reviewed launch is about 34 threads, below that cap; the separate native cap of 4 gives the recorded combined ceiling of 44, below the phase-wide limit of 50.
- The recorded launch constraint requires the old cap-12 global-control queue to be terminal before the cap-40 inner-32 run. `continuation-current.json` records native cap 4 as completed. This review treats those as launch preconditions; it did not inspect external process state.

## Dynamic-team scope note

`omp_set_dynamic(0)` in the consumer's calling thread is not by itself a portable guarantee for later `std::thread` worker contexts: the OpenMP dynamic adjustment ICV is task-scoped. For the exact reviewed run plan, the process is started with `OMP_DYNAMIC=FALSE`, which initializes the runtime accordingly before any worker threads are created. Thus the plan's 32-thread execution is configured consistently. The driver's `global_consumer_openmp_dynamic_teams_disabled` metadata is an assertion derived from the requested GSH thread count, not an observed runtime query; outside this plan, a caller using `OMP_DYNAMIC=TRUE` could make that field inaccurate. For general-purpose metadata, set dynamic mode in each worker or record `omp_get_dynamic()` from the worker context.

## Hash binding

| Reviewed file | SHA-256 |
|---|---|
| `candidates/combined-global-v7/global/global_kernel_consumer.hpp` | `5555a444b9ff0ec76555ad4fdec51c1dbb0d25ad06bfd8d776bb65c1c182dfd3` |
| `candidates/combined-global-v7/global/global_gsh_tile.hpp` | `5d5f0b4965dd427f1d3f956e0f828f9dba80006a6bc1a60a038de8e24a3450b5` |
| `candidates/combined-global-v7/physical_driver.cpp` | `432ac75daa0e72f1bfed742240053dbca282065e64c84a984abc4c245c082388` |
| `candidates/combined-global-v7/implementation_status.md` | `ec9e555f11468905a55c70b3370dc05cc0424a1e381054d48f404904ecfcbefa` |
| `configs/global-threads-build.json` | `27335c9a87aa210829bfd6807587ccc34ed3226d99d2f2d0ee8c8316d8d70066` |
| `configs/global-threads-build-plan.json` | `858aee5be817e5abb9cbbf7b7d85178dc81a71298a48a03dfec08e2514e893e1` |
| `configs/global-threads-run-plan.json` | `6d372af5ef1199b87a11d6a3060d74d8e45054b4b924e08b66cac81bceddc3f8` |
| `continuation-current.json` | `e95ce263eca9160fec658e99dcd0e73084b3b618e2cf90b5150b5ea6edc512bb` |

