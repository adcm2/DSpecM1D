# Short final Sol review: global run gates

2026-10-27. **HOLD for one checker/supervisor integration issue.** Numerical/layout gate review and local data tests otherwise pass. No SSH, compilation, science execution, or implementation edits.

## Required correction

`check_global_pair.py::_unchanged()` currently rejects a changed whole `state.json` hash. The unchanged supervisor rewrites the running checker job's watchdog every two seconds. A comparison spanning such an update can fail although the completed reference/candidate job records and all verified artifacts are unchanged.

An independent bounded pure-data counterexample verified a completed job, changed only another running job's watchdog elapsed time, and verified the completed job again. Its job object and verified artifact map remained identical, but `_unchanged()` raised `completed-job output or state hash changed during comparison`. This is an unrelated-state false failure, not numerical evidence.

Compare completed job objects and their verified artifact maps before/after, as the strict pair checker already does; keep whole-state hashes informational. Add the narrowly targeted unrelated-watchdog regression and obtain fresh review before launch. Preserve rejection of changes to the completed jobs or artifacts. No other correction is required by this review.

## Checks passed within this boundary

- All 16 supplied data-only tests passed independently. An additional end-to-end synthetic receipt harness used saved actual native-side metadata, angle-major direct arrays over all 179 positive sides, and side-major global arrays over the plan's eight selected sides. Direct→global and w1→w4 comparisons passed; refreshed-receipt kernel corruption failed both numeric and exact gates. These are synthetic gate tests, not solver results.
- The native identity map correctly selects and transposes the direct reference rather than requiring 179 and 8 side arrays to have equal shapes. The direct comparison retains the fixed per-channel `2e-13` scaled gate. Parallel comparison requires byte-exact kernel/trace/spectrum files and exact energy/qdot scalars, side order, and actual 1-versus-4 workers for eight sides. Worker elapsed/stage sums close against per-worker durations and remain separate from wall time.
- The fixed physics reports 52 output bins, while capture retains 53 including the zero-q guard. Frequency audits require contiguous identity, finite complex weights/damping, and residual ≤`1e-12`; sparse-factorization attribution requires `53*(2L+1)`. Physical metadata and independently recomputed energy identity remain gated. The existing mathematical/source review is unchanged.
- Independently asserted all 12 plan jobs: two original-v6 direct references, four global candidates and six checks at L=3/12. Binary at command index 5 matches declared input; candidates use combined build03, cache off/spool on, symbolic env 0, selected sides `0,28,59,100,149,154,172,173`, requested workers 1/4, CPU set 4–15, 32 GiB working cap and 64 GiB job cap. Physical arguments and forward-only-off method match references. Each direct gate is immediate; the parallel gate follows both, before advancing degree. Any nonzero gate stops later jobs under the existing supervisor.
- Reconstructed `stage_reviewed_run.py` selection confirms all staged checker/plan input paths exist in its selected package beneath `global-run-v1`. The helper uses exclusive staging/archive writes, verifies manifests and reviewed tools, and checks each declared physical-driver binary against successful unchanged-input build-command receipts before detached launch. It was not executed. Actual remote/build/resource state remains parent-owned launch evidence.

Reviewed hashes: checker `719f3581c24a1e9ebf4cd609ca6dbba1c67ec9fae09c5a8f3e1d160fe3dad74a`; tests `9c3a7c9719d273ddfae93dae18870426649cd2ef3b61c2017fad2b7f4d041fff`; plan `1e242190dbdf35a05f6c210f2afd8bd949db76af03bcbb0f67e576289908b361`; stage helper `55cb407ff609df44190600582f4d3f8b9c2ca9dc6ddd65429e6e6dc832c32848`.

No numerical correctness, measured speedup, scientific completion, or figure approval is established. Resolve only the identified HOLD and re-review the resulting delta.
