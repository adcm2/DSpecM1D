# Fresh review: combined-cache-v2 memory and shared-grid policy

**Decision: PASS for this source-only review scope.** The v1 memory omissions are addressed in v2. The shared-grid/thread-local execution design is unchanged and matches the locally pinned dependency implementation. This PASS does not clear the separately known compile issues, authorize a build/science run, establish numerical agreement, or claim a speedup.

Reviewed SHA-256:

| File | SHA-256 |
|---|---|
| `candidates/combined-cache-v2/global/global_gsh_tile.hpp` | `c19f438af8ba56dadb84617a4f03302e52075d36fa6939b88146f7c392fe9eac` |
| `candidates/combined-cache-v2/global/global_kernel_consumer.hpp` | `3a84e1cfdf02e48ff48f8ffd08146e122894335117126ff3081781e8a1845b90` |
| `candidates/combined-cache-v2/global/pole_response_spool.hpp` | `5e8868439771cb502cb4e3005039f6986e9d0d6ad23a58b168a7a61c2c9152dc` |
| `candidates/combined-cache-v2/physical_driver.cpp` | `8b763d85bf1552809bdf44b6ed1a5e16af6c218c1f52ad07f30e6fe5c237ae74` |
| `candidates/combined-cache-v2/implementation_status.md` | `1e2b67cc029d1742080c5d6bef2cc6bd275dd00ab7ab4a8f2db8ee1ace251c6e` |
| `candidates/combined-cache-v2/engine_headers/Windowed3DFields.h` | `a66f74690d5655335598948d2b804c83d19eca6606481a26214389e0dff9feb0` |
| Local GSHTrans `GSHTrans/SphericalGrid.hpp` | `5697c22bd34464449f419dc36cc0c4ca40f9c9c24b93a7693e611fe220e4cca0` |
| Local FFTWpp `FFTWpp/src/Core.h` | `4e63807a85724f33de472c1bc5f9ceb14b39f3533edb8be667ea288a2ab6ee59` |
| Local FFTWpp `FFTWpp/src/Plan.h` | `524462848f733b5659793eba3c74ea4aa183032506922a2621012df65010ff27` |

Dependency source revisions are GSHTrans `8dd116df8df4222fe9b10a40c7213d34f2a274fb` and FFTWpp `cb152c0e85f3822493c2d9e85a56307ac7b4f1c0` from the local `DSpecM1D/build/3d-kernels-20261025/phase-c/_deps` checkout. No SSH was used.

## Review results

- **Memory accounting remediation: PASS.** V2 uses checked multiply/add helpers for the added bounds. It accounts for padded `sizeof(NativeSide)` and `sizeof(SpoolAudit)`, vector capacity slack, each `PoleResponseSpoolReader` object, the primary reader once as shared memory, and one reader per worker. It includes the in-flight `FrequencySideRecords` object and upper bounds for both decoded records' labels and six-complex strain arrays using the correct `5L-1` label count for local orders at degrees 0 through L. These resident and in-flight costs feed the shared and per-worker budgets before worker count selection. The 256 MiB maximum offset-validation allocation remains included and charged once because reader construction is serial. The configured 2 GiB-per-worker transform reserve remains separate and included. The manifest continues to call this an estimate and explicitly excludes unknown runtime/library/thread-stack peaks.
- **Shared grid and FFT plan concurrency: PASS.** The v2 `global_gsh_tile.hpp` is byte-identical to the reviewed v1 file. The locally inspected GSHTrans source documents immutable shared grid state and thread-local workspaces, plans, and scratch; its transform entry point is const. FFTWpp planner calls use a process-wide planner lock. Plans are still initialized lazily on first use in each worker, so construction is serialized safely by the dependency and charged to worker transform timing rather than the serial setup interval. Each worker owns its mutable plan, phase/norm/Wigner cache, reader, output stream, and counters; calls request sequential GSH execution.
- **Wigner table, phase cache, and high-L batching: PASS.** Stored-mode table bytes still use the five-order complete table formula and are charged once. Generated-mode scratch is charged per worker, bounded as a per-angle working set rather than all theta rows. The two-frame phase/norm cache and corresponding frame Wigners are charged per worker. Four same-spin channels are used through L=32; high-L uses two only if budget allows and otherwise falls to one. The chosen batch scratch participates in worker-count sizing.
- **Guard and execution attribution: PASS.** The spool still declares and validates all 53 contiguous records, including the out-of-band guard. Each worker reads/validates each side/bin record and counts an exact-zero `adjointMixing` as skipped; only nonzero weights increment `transformsExecuted` and call the global transform. For the stated producer case this represents 53 retained records, 52 transform executions, and one zero-mixing skip per side. Both counts are emitted in the sidecar and propagated to result metadata. This is auditable data; the checker must bind the per-bin zero to the declared guard before treating it as expected.
- **Inherited compile remediation: source evidence only.** V2's `Windowed3DFields.h` and spool header are byte-identical to combined-v4, carrying the visibility and stream-accessor fixes noted by the parent. The consumer has local checked angular-size helpers and no longer calls the private engine size helpers. This review did not compile or otherwise verify those fixes.

No code was changed. No build, remote access, or science run was performed.
