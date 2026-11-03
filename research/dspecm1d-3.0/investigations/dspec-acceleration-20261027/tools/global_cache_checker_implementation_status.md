# Global cache/storage pair checker status

## Change log

- Added the separate read-only `check_global_cache.py` for completed full-grid generated/stored and worker-count comparisons. It reuses completed-job receipt, physical-result, global-manifest, audit, timing-rollup, and numeric helpers without changing the reviewed spool-pair checker.
- Bound receipt command flags and explicit CLI worker/storage expectations to the result and manifest, including the selected-side order, full-grid dimensions, output bytes, sampled/full angular counts, Wigner allocation, batch size, estimated memory, and executed/skipped frequency-side counters. Requires one to eight selected sides, 180 total physical sides, 179 positive-radius rows, SEM degree 6, 52 output bins, and 53 spool audit records.
- Added bounded-memory per-channel kernel comparison and trace/spectrum/energy/qdot comparisons using the existing `2e-13` threshold. `--exact` additionally requires exact kernel/trace/spectrum hashes and exact energy/qdot scalars. Completed job objects and verified artifact maps are rechecked after comparisons; the broader state hash remains informational through the shared receipt helper.
- Added producer-shaped synthetic tests for supported degrees, generated/stored modes at matching worker counts, generated comparisons across 1/2/4 workers, worker/storage provenance mismatches, producer-counter mismatches, native identity/metadata, frequency and q audits, energy identity, and channel-level numeric regressions.
- Preserved the pre-integration-fix checker and test bytes in `preserved-global-cache-v1/`. Corrected the final stability recheck to call the frozen helper's three-argument API, and added a completed-job end-to-end `assess()` fixture covering receipt hashes, result/profile/manifest, comparisons, exact mode, and the final stability recheck.

## Validation

- `python3 -m unittest -v test_check_global_cache.py`: 10 tests passed, including end-to-end `assess()`.
- `python3 -m py_compile check_global_cache.py test_check_global_cache.py`: passed.
- No completed-job artifacts were available for a real receipt-to-output comparison. No build, SSH, or scientific execution was performed.
