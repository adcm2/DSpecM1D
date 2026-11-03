# Targeted fresh review: global cache state recheck

**Verdict: PASS for the declared benchmark plan and the targeted state-recheck fix.** The bounded worker/side combinations and the 23-job plan remain as reviewed in `global-benchmark-fresh-luna.md`; the separate generic worker/side HOLD remains applicable outside this plan.

Reviewed `tools/check_global_cache.py` SHA-256 `9890e033bfab6f9bd7a4c3bfdf42cab005fdddb97fa248cbd107d0d82c1a21af` and `tools/test_check_global_cache.py` SHA-256 `9308a8dbeaf252fb0d7cbbf8cc18c51a2bf7e75e57b0e175b4bbd743fb88562b`. The preserved v1 checker and tests remain byte-identical to their reviewed hashes (`303d7c1e…c04d63c9` and `1df4c165…1cc9dd3b`).

The only checker source change is at `assess()` completion: it now obtains the final receipts through `spool._unchanged(roots, jobs, receipts)`. That helper compares the completed-job records and verified artifact maps, allowing unrelated state-file heartbeat changes while still rejecting completed-job or receipted-artifact mutations. It returns the final receipt tuples used for the informational state and artifact hashes in the report. This corrects the prior duplicated per-job call without changing comparison tolerances or scientific gates.

Verification: `python3 -m unittest -v test_check_global_cache.py` passed all 10 tests, including a completed-receipt end-to-end `assess(..., exact=True)` fixture. `python3 -m unittest -v test_check_global_pair.py` passed all 16 tests, including unrelated-watchdog acceptance and completed-job/artifact mutation rejection in the shared `_unchanged` helper. The bounded 23-job plan remains compatible: every producer requests no more workers than selected sides, and every exact pair matches degree, selected sides, worker count, and storage mode.

No implementation files were changed during this review. No build, SSH access, or science run was performed.
