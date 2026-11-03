# Fresh static review: global GSH threads run plan v2

**Result: PASS for the bounded launch-plan delta.** No SSH, build, or science run was performed; actual process peaks remain unmeasured.

After removing the single added `OMP_NUM_THREADS=32` assignment from each producer command, the v2 plan is structurally identical to v1. It retains all 19 job IDs, inputs, outputs, limits, commands otherwise, and comparison arguments. All nine comparison jobs still request `--exact`; the checker keeps its existing `2e-13` non-exact gate as an additional reported metric. Each of the ten producer jobs requests one outer worker and 32 GSH threads (product 32), within the plan's aggregate `thread_cap` of 40. Producer environments retain `OMP_THREAD_LIMIT=32`, `OMP_DYNAMIC=FALSE`, `OMP_MAX_ACTIVE_LEVELS=1`, and the same CPU affinity. The only added per-job environment assignment is `OMP_NUM_THREADS=32`.

The environment choice matches the pinned GSHTrans API. `SphericalGrid` constructs stored Wigner values in its immutable implementation; `Wigner::ComputeAll` uses an OpenMP `parallel for` with no explicit `num_threads`, so it uses the runtime default established by `OMP_NUM_THREADS`. V1's supervisor environment sets that default to 1, which would serialize stored-grid setup even though later inverse transforms request 32 threads explicitly. V2's producer-level assignment overrides that supervisor default for the complete driver process. The pinned loop partitions independent `(upper-index, theta)` entries and writes their separate table blocks; this changes setup scheduling, not any reduction or the arithmetic order within an entry. The consumer still calls `Eigen::setNbThreads(1)`, and `stage_reviewed_run.py` sets OpenBLAS, MKL, BLIS, VECLIB, and NumExpr thread variables to 1 for the supervisor; these are inherited by producer jobs. No source or numerical kernel code changed in this delta.

The 40-thread cap is the configured supervisor bound, not an observed maximum. The plan requests 32 compute threads for each producer and leaves headroom under that bound; measured peak counts require the future run receipt. The plan remains marked diagnostic, and this static review makes no claim about speedup or numerical outcome.

## SHA-256

- `configs/global-threads-run-plan.json` (v1): `6d372af5ef1199b87a11d6a3060d74d8e45054b4b924e08b66cac81bceddc3f8`
- `configs/global-threads-run-v2-plan.json` (v2): `219ab197730e2b5aca17ea02e0b927dceb1fae39bba19c543a5388df7c0c25ea`
- `candidates/combined-global-v6/global/global_kernel_consumer.hpp`: `1a1600999680edf5f504a37471c80a66eca60db1727af9c4a7512b63ab3c8632`
- pinned `GSHTrans/SphericalGrid.hpp`: `5697c22bd34464449f419dc36cc0c4ca40f9c9c24b93a7693e611fe220e4cca0`
- pinned `GSHTrans/Wigner.hpp`: `1fd2197a9f37ee463b49ecf0c615300cf98ca53612362ec1c3787ad8397b62c4`
- `tools/stage_reviewed_run.py`: `55cb407ff609df44190600582f4d3f8b9c2ca9dc6ddd65429e6e6dc832c32848`

## Follow-up HOLD: resident OpenMP pools

A later independent Sol review identified a resource risk not bounded by the static `1 × 32` active-team product above: the main-thread Stored-Wigner setup can leave an OpenMP pool of 31 worker threads resident while the outer worker starts its own 31-thread pool. Peak resident threads can therefore approach 64, exceeding both the plan's 40-thread cap and the 50-thread host limit. The v2 plan is **not cleared to launch** until the new v8 implementation explicitly calls and checks `omp_pause_resource_all(omp_pause_hard)` between main-thread grid/context setup and worker execution, and that behavior is independently validated. The prior PASS applies only to active arithmetic teams and did not establish a resident-thread bound.
