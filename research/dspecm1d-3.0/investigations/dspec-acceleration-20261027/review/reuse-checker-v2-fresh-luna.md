# Reuse checker: fresh Luna review, v2

**Result: PASS for the standalone data checker.** This reviews `tools/check_reuse.py` and its call into the existing strict pair comparator only. It does not approve candidate source, build integration, or a scientific run. No SSH, build, or science run was performed.

The two fallback fields now must each be explicit JSON `false`; missing, `null`, and true values fail. The fixed case requires exactly 52 frequency bins, so the formula resolves to 53 capture misses/factorizations and 52 cache hits per active family, `2L+1`. The exact cache entry, miss, preflight-entry, and hit counters must match; retained bytes must be positive and within both preflight estimate and budget; receiver-basis reconstruction is bounded by `1e-12`; factorization counts must close against symbolic analysis/reuse, with observed reuse required when symbolic mode is on.

Parity still goes through unchanged `check_pair.assess()`: the `2e-13` comparison limit and `PHYSICAL` metadata checks remain strict. Candidate `forward_only_baseline` must be false, and parity requires the reference to match. The candidate result is re-read only after parity passes and checked against the already verified artifact SHA-256 before attribution. This supports the documented method of running both reference and candidate with `--forward-only off`; resulting timings must remain labeled as this route when compared with prior forward-only-on measurements.

The parent reports pure-data checks passed for valid reuse, missing/null fallback flags, altered bin count, and missing hits. I did not rerun them. The remaining global/module compile integration review is separate and pending; this data-only PASS grants no public-source or build approval.
