# Radial-L0 exact reuse checker: fresh Luna review

**Result: PASS for the targeted pure-data checker gate.** The helper keeps `check_pair.assess()` as its first gate, including the unchanged `2e-13` scaled field-parity tolerance and unchanged physical metadata checks. Only a parity PASS proceeds to reuse attribution; the helper then re-reads and hashes the candidate `result.json` against the verified candidate receipt before parsing it. The five local unit tests passed. No build, SSH, or science run was performed.

The helper fails closed unless direct response reuse was explicitly requested and used, spool and both fallback flags are explicitly false, the baseline is explicitly not forward-only, symbolic mode matches the requested mode, and the result reports 52 physical bins plus the radial-L0 bypass. For degree `L`, it requires exactly `53*(2L)` nonradial entries/misses and `52*(2L)` hits, 105 bypass calls, and `53*(2L+1)+52` factorization calls. It checks retained bytes ≤ estimated preflight bytes ≤ budget, receiver-basis residual ≤ `1e-12`, matching result/profile bypass counts, and symbolic analyses plus reuses equal factorization calls. Symbolic-on requires a reuse; symbolic-off requires zero reuses.

The five pure-data tests cover the valid L=3 case, counter drift, both-pass bypass accounting, missing/true fallback and spool routing, plus basis and symbolic-closure failures. One small coverage limit: the checker enforces the cache-budget inequality, but these five tests do not inject an over-budget result to assert rejection. The guard itself is present and correctly ordered; this remains a test-coverage note, not a source-review hold.

## Hash bindings

- `tools/check_reuse_radial_exact.py`: `a985177a244eef296ec1b80751880ba4e6b23db6d4c8f0345f1c9fe430420cec`
- `tools/test_check_reuse_radial_exact.py`: `1691d0de6b3a3ecdd851cfa91c41924111604b3b12bd4eff9f3089ba014c7df6`
- `tools/check_pair.py` (unchanged parity implementation): `6f4d72e7da0489b6b3d0b6fd1d5c15d7ba3b7e6179a14ef453b7e8ccb762e646`
- Final-v3 radial-L0 route: `engine_headers/Windowed3DFields.h` `79a79bd14d0a8efe7a120ac9637bd5c3169da603d5b9cccd5e03cae721247dbd`; `physical_driver.cpp` `1834483513308b81346b53844a015c1990efe5b6741a5b5cb46a144fd13c3a11`

## Gate boundary

This PASS supports attribution only when the helper returns PASS on receipt-verified jobs. It does not upgrade the preserved failed L128 receipt or establish scientific correctness or performance. Five local pure-data tests passed; no build, remote access, or science execution was involved.
