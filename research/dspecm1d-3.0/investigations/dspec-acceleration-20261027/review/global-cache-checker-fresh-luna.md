# Fresh Luna review: global cache checker

**Verdict: HOLD for generic worker/side combinations.**

Reviewed `tools/check_global_cache.py` at SHA-256 `303d7c1ef67ef3758305595db681bb758448d600f83b99ad8a51f771c04d63c9` and `tools/test_check_global_cache.py` at SHA-256 `1df4c165d8ff13de00b858dded0e39a6c5f6bb9c119fe47bd0b6248a1cc9dd3b`. Ran `python3 -m unittest -v test_check_global_cache.py`: all 9 tests passed.

The checker correctly restricts the degree matrix to L=3/12/128/384/1536; binds the requested storage mode and worker request to the receipted command, result, and manifest; validates complete-grid structure, native-side tables, all frequency records (including guard), profile/factorization counts, finite kernel values, energy/qdot consistency, worker timing rollups, and result/manifest transform counters; compares trace, spectrum, kernel, energy, and qdot under the shared 2e-13 gate; and makes `--exact` an additional byte/scalar equality gate. It reuses the existing global manifest and immutable completed-job validators.

**Blocking compatibility finding:** `validate_invocation` at lines 68–72 requires manifest `workers` and timing `actual_workers` to equal the explicit requested count. The producer in `candidates/combined-v4/global/global_kernel_consumer.hpp` line 153 (and cache-v2 equivalent around line 200) chooses `min(requested workers, selected sides, 48, workers-by-memory)`. Thus, for a valid request of 8 workers with 1–7 selected sides, actual workers are capped by side count (or memory) and the checker rejects the completed run. This conflicts with the generic 1–8-side matrix. Validate requested workers against command/result as now, and validate actual workers against the producer cap (for the current memory-nonbinding plan, `min(requested, selected_count)`); if memory-capped cases are intended later, represent and validate the cap explicitly. Add a regression using requested=8 and selected sides=2 with actual=2.

No implementation files were changed. This was a static/data-only review; no SSH, build, or science run was performed.
