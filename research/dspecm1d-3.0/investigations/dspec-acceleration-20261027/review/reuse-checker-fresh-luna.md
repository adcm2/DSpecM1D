# Reuse checker and experiment method: fresh Luna review

**Result: HOLD for final checker acceptance; experiment method is sound.** Data-only static review; no build, SSH, or science run.

## Findings

- **Pass — parity and route controls.** `check_reuse.py` calls the unchanged `check_pair.assess()` with its fixed `2e-13` per-column/channel gate and unchanged `PHYSICAL` field comparisons. It adds no parity exceptions or relaxed thresholds. The candidate must report `forward_only_baseline == false`; because that field remains in `check_pair.PHYSICAL`, the reference must report false as well. Running frozen original-v6 and candidate with `--forward-only off` makes the captured-unit baseline methodology explicit and keeps this field matched. This is a valid new parity fixture, but its timings must be labeled as the forward-only-off route; compare performance only to a reference run with that same option, and keep prior forward-only-on timings as a separate route.
- **Pass — expected reuse attribution.** With `frequency_bins == 52` and degree `L`, `families = 2L+1` yields 53 capture misses/factorizations and 52 kernel reuse hits per family. The checker derives `(frequency_bins+1)*(2L+1)` misses and `frequency_bins*(2L+1)` hits, requires cache entry/miss/preflight counts to agree, requires factorization count to equal capture misses, and closes symbolic analyses plus reuses against all factor calls. Requested symbolic-on also requires at least one actual reuse.
- **Pass — evidence binding and numeric gates.** `assess()` keeps the strict comparator and before/after receipt verification. The candidate `result.json` is re-read once and its bytes are hashed against the already-verified manifest before reuse metadata is parsed. The checker requires positive retained bytes within both estimate and budget, and basis reconstruction no greater than `1e-12`. No cache-use or fallback condition is added to or substituted for the parity gate.
- **Hold — fallback booleans are not fail-closed.** `if result.get('response_reuse_capacity_fallback') or result.get('response_reuse_preflight_fallback')` accepts missing keys and `null` as “no fallback.” Require each field to be explicitly `False` (and ideally validate JSON boolean types) so absent/unknown evidence cannot satisfy the no-fallback gate. Other primary route fields already use exact `is True`/`is False` checks.

## Method gates

1. Run both frozen original-v6 reference and candidate with `--forward-only off`; keep all `check_pair.PHYSICAL` fields and its `2e-13` tolerance unchanged.
2. For the intended 52-bin case, assert `frequency_bins == 52` in the run plan or checker; then the dynamic formulas resolve to exactly `53*(2L+1)` captures/factorizations and `52*(2L+1)` cache hits. Record the exact tested `L` in the run plan.
3. Require explicit false fallback metadata, actual cache use, no spool, cache bytes within the corrected estimate and budget, basis residual <=`1e-12`, and complete symbolic count closure before accepting attribution.
4. Report timings against the matching forward-only-off reference as a separate route from previous forward-only-on measurements.
