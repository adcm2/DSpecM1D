# Fresh review: combined-cache-v1 shared-grid and memory policy

**Decision: HOLD on resource-preflight acceptance; shared-grid/thread-safety design passes static inspection.** Read-only source review only; no build, science run, or SSH. The known combined-v3 compile problems called out by the parent are outside this review and must still be resolved in combined-v4.

Reviewed SHA-256:

| File | SHA-256 |
|---|---|
| `candidates/combined-cache-v1/global/global_gsh_tile.hpp` | `c19f438af8ba56dadb84617a4f03302e52075d36fa6939b88146f7c392fe9eac` |
| `candidates/combined-cache-v1/global/global_kernel_consumer.hpp` | `023d823a86c8c90c455fce68256c041e9eb90ddeee4308c0fdf0e000aa8d2f37` |
| `candidates/combined-cache-v1/global/pole_response_spool.hpp` | `39a8e1bd6de5d52257fdfcb365f68fb09c9609439696dbba3e9c81b7d3605ea4` |
| `candidates/combined-cache-v1/implementation_status.md` | `818c204b5ce5059dc36c477ef92d40174d3f6789a564c758a22cb9712c3c431e` |
| Local GSHTrans `GSHTrans/SphericalGrid.hpp` | `5697c22bd34464449f419dc36cc0c4ca40f9c9c24b93a7693e611fe220e4cca0` |
| Local FFTWpp `FFTWpp/src/Core.h` | `4e63807a85724f33de472c1bc5f9ceb14b39f3533edb8be667ea288a2ab6ee59` |
| Local FFTWpp `FFTWpp/src/Plan.h` | `524462848f733b5659793eba3c74ea4aa183032506922a2621012df65010ff27` |

The local dependency checkout used for the API inspection is under `DSpecM1D/build/3d-kernels-20261025/phase-c/_deps`; its GSHTrans revision is `8dd116df8df4222fe9b10a40c7213d34f2a274fb` and FFTWpp revision is `cb152c0e85f3822493c2d9e85a56307ac7b4f1c0`.

## Findings

- **Shared grid and per-worker transforms: PASS.** The inspected `SphericalGrid` API documents `shared_ptr<const Impl>` and immutable grid data. `InverseTransformation` is const; Wigner-generation scratch, coefficient scratch, transform buffers, and workspace/plan caches are `thread_local`. Each outer worker owns its own `GlobalGridPlan`, reader, and output stream; those workers call the shared grid with `Execution::Sequential()`, so no library-internal transform parallel layer is requested. The two-frame Wigner/phase/norm cache lives on the worker-owned plan.
- **FFTW plan creation: PASS for safety; lazy, not prebuilt.** GSHTrans creates each worker's transform plan on that worker's first inverse call and retains it in that worker's thread-local cache. FFTWpp's planner entry points take a shared `PlannerLock`, serializing planner operations across workers. Thus plans need not be built on the main thread for correctness. First-use plan cost is correctly included in worker transform time and excluded from serial grid/context setup; the manifest calls this out. Do not describe this as all plans being constructed in the serial setup phase.
- **Stored Wigner table and high-L batching: PASS as bounded accounting logic.** The formula `(2L+1)*(5*(L+1)^2-10)*sizeof(double)` matches the complete stored table for five upper-index values `n=-2..2` and all orders, and it is charged once as shared memory. Generated-mode scratch is charged per worker as a per-angle estimate rather than for every theta row; its helper includes a conservative all-`n` Wigner and coefficient-table extent, which bounds the GSH generated-mode per-thread buffers. Low-L uses four same-spin channels; high-L uses two only if the post-shared, post-reserve budget allows it, otherwise one. Worker count is reduced by the resulting per-worker bound and capped by requested workers, selected sides, and 48.
- **Memory preflight: HOLD.** `readerPersistentBytes` counts one side/audit payload per worker using packed field sums, but each `PoleResponseSpoolReader` stores vectors of `NativeSide` and `SpoolAudit` objects whose `sizeof` includes alignment padding. Also, the initial `metadata` reader remains alive while all worker readers are created and the transforms run, but the estimate counts only the worker readers. Finally, each worker's `read()` constructs a `FrequencySideRecords` value with two decoded pole records: each owns label and six-complex-value arrays proportional to up to `5*(L+1)` entries. These live during GSH evaluation and are not included in `knownPerWorkerBytes`. The configured 2 GiB transform reserve does not replace accounting for known allocations. Price the resident readers using `sizeof(NativeSide)` and `sizeof(SpoolAudit)` and account for the initial reader plus workers; add a checked bound for the two decoded records per active worker. Keep the 256 MiB validation bound and explicit per-worker transform reserve, and continue labeling the result as an estimate rather than an RSS guarantee.

## Scope boundary

This review covers only the shared angular grid, plan/cache ownership, phase cache, high-L channel batch selection, and working-memory estimate. It does not clear the separate known compile errors, verify numerical equivalence, or establish performance. No code was changed.
