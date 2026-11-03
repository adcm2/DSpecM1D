# Short Sol follow-up: global state-recheck HOLD resolved

2026-10-27. **PASS for the targeted correction to the previously reviewed global run gate.** This clears only the unrelated-state false-failure HOLD in `global-run-short-sol.md`, following the fresh Luna v4 review. No implementation changes, SSH, build, or science execution.

`_unchanged()` now compares the full completed job record and its verified artifact map before/after comparison. Whole-state SHA-256 values remain informational. An unrelated running checker watchdog update therefore cannot invalidate unchanged completed evidence; completed job mutations and changed receipted artifacts remain rejected.

Independently reran the exact targeted regression test successfully. Re-exercised the original counterexample using a separate running-job watchdog record: changing only that watchdog now passes while state digests differ. Adding a field to the completed job record (while preserving its qualified status) still fails the new full-record comparison. This directly verifies both sides of the corrected boundary. The fresh Luna review reports all 16 data-only tests passing; no broader new test scope was added here.

`global-run-v2-plan.json` remains byte-identical at its prior reviewed hash. All prior numerical/layout/worker, native-subset, 52-output/53-spool, exact parallel parity, failure-stop and staging findings carry forward unchanged from `global-run-short-sol.md`. This is the existing build03/combined-v4 plan gate; it does not silently extend the checker to new final-candidate storage/timing schemas or to different run plans. Final-candidate build-only review is separately recorded in `final-build-short-sol.md`.

Reviewed SHA-256:

- `check_global_pair.py`: `340efdf0c60b80eea23fa69d642595cf7b989112d9310ed792b21ebe4b36d51a`
- `test_check_global_pair.py`: `1e425551cc878159a1142432ae008162f122b522ddcf5ff5b35f914916f09846`
- `global-run-v2-plan.json`: `1e242190dbdf35a05f6c210f2afd8bd949db76af03bcbb0f67e576289908b361`

This resolves the bounded checker integration HOLD, not compilation, runtime numerical correctness, speedup, scientific completion, or figure approval. Successful receipted build and current launch/resource checks remain parent-owned execution gates.
