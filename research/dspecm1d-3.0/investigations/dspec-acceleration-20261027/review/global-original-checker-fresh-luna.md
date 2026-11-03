# Fresh targeted review: direct-original to global-spool checker

**Result: PASS for the checker wrapper and synthetic receipt fixture.** The fixture test passed (`test_full_direct_179_to_global_8_through_reviewed_assess`). No build, SSH, or science run was performed.

`check_global_original.py` delegates the comparison to the reviewed `check_global_pair.assess` with `probe64=False` and `global_reference=False`, so this is a complete-grid direct-original comparison, not a probe or global-to-global shortcut. The reviewed checker requires the original 179 positive-radius sides, validates the candidate's selected native-side identities as a subset, and compares the selected kernel data using the existing `check_pair.LIMIT` of `2e-13`. The synthetic receipt fixture maps 8 selected sides from the 179-side direct output and exercises receipt/output hashes and the full compare path; it passes.

The sole compatibility shim applies the existing `_final_timing_adapter` to an in-memory manifest copy. It requires `serial_grid_and_context_setup_wall_s` and checks that every original timing field is preserved and the compatibility value exactly equals the source. The reviewed timing validator is restored in `finally`, including if assessment raises. Provenance in the resulting report includes the reviewed checker helper SHA. The wrapper does not replace comparison logic or alter the tolerance. This establishes the adapter path and fixture behavior only; the current Estimate L3/L12 artifacts still require assessment against actual receipted outputs before their numerical status can be stated.

## SHA-256

- `tools/check_global_original.py`: `8e5b3792cce309b5e6570a3e7d946f3a0599e933b79c7ec449bb04804e58d56f`
- `tools/test_check_global_original.py`: `37ecbeb90ac66c21018f8d3c24f705eee22899949e54faba542f0fa80d89d2df`
- `tools/check_global_pair.py` (reviewed helper): `340efdf0c60b80eea23fa69d642595cf7b989112d9310ed792b21ebe4b36d51a`
- `tools/check_global_cache.py` (timing adapter source): `9890e033bfab6f9bd7a4c3bfdf42cab005fdddb97fa248cbd107d0d82c1a21af`
- `tools/check_pair.py` (tolerance/hash helper): `6f4d72e7da0489b6b3d0b6fd1d5c15d7ba3b7e6179a14ef453b7e8ccb762e646`
