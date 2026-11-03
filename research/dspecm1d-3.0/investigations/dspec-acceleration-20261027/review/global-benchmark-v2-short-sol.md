# Short Sol follow-up: global benchmark API repair

2026-10-27. **PASS for the targeted API repair and the existing diagnostic benchmark plan.** This clears only the call-signature HOLD in `global-benchmark-short-sol.md`, following the fresh v2 Luna checker review. No science execution, SSH, compilation, or implementation edits.

Independent comparison against the preserved, hash-identical v1 checker confirms the only checker change is replacement of the invalid four-argument per-job call with `after = spool._unchanged(roots, jobs, receipts)`. This matches the actual imported helper, rereads both completed receipts, compares full completed job records and verified artifact maps, and returns final receipts for report hashes. Numerical comparisons and all other gates are unchanged.

Independently ran all ten current helper tests and all sixteen imported global-pair tests: **26 passed**. The new test executes `assess(..., exact=True)` end-to-end using complete synthetic job files, qualified completed receipts and matching generated/stored outputs; it now reaches PASS with exact files and zero kernel difference. Imported tests retain unrelated-watchdog acceptance and completed-job/artifact mutation rejection. Preserved v1 checker/tests retain their original reviewed hashes.

The 23-job benchmark plan hash is unchanged. Prior assertions for 14 producers, nine exact gates, matching degree/native-side/storage/worker invocations, final04 binaries, CPU set 4–15, cap 12 and bounded memory combinations carry forward. Requested workers do not exceed selected sides in this plan, so the separate generic worker/side compatibility HOLD remains inapplicable to these exact invocations, assuming memory remains nonbinding. Remote machine/build/resource checks remain parent-owned launch evidence.

Scope remains internal generated/stored/parallel diagnostics only. Exact internal parity does not clear the unresolved original-direct/global `2e-13` cross-route FAIL, establish scientific acceptance, justify a full179-side volume or authorize figures. No tolerance is relaxed and no measured numerical success or speedup is claimed here.

Reviewed SHA-256:

- Helper: `9890e033bfab6f9bd7a4c3bfdf42cab005fdddb97fa248cbd107d0d82c1a21af`
- Tests: `9308a8dbeaf252fb0d7cbbf8cc18c51a2bf7e75e57b0e175b4bbd743fb88562b`
- Shared global checker: `340efdf0c60b80eea23fa69d642595cf7b989112d9310ed792b21ebe4b36d51a`
- Plan: `9cc5c7a67f50139e169e2be0e39872ec54b1faf4f1fb6bef989d081af2e202a1`
