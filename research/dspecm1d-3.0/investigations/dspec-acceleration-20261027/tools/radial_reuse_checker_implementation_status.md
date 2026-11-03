# Radial-L0 response reuse checker status

## Change log

- Preserved the reviewed `check_reuse.py` bytes in `preserved-reuse-checker-v1/` before adding the separate `check_reuse_radial_exact.py` helper.
- Reuses the unchanged `check_pair.assess()` parity gate and `2e-13` threshold, then validates the final-v3 route with cache/spool/fallback booleans, forward-only-off baseline, symbolic mode, budget/preflight/storage bounds, and the `1e-12` receiver-basis reconstruction gate.
- Applies the final-v3 radial-L0 arithmetic-preservation accounting: nonradial cache misses/preflight entries are `53*(2L)`, nonradial hits `52*(2L)`, and sparse factorisations `53*(2L+1)+52`. Requires the source-reported radial bypass flag and 105 bypass calls (53 baseline plus 52 kernel), with matching profiler counter and complete symbolic analysis/reuse closure.
- Added five focused pure-data tests for the accepted formulas, count drift, bypass counts, fail-closed route/fallback flags, basis, storage, and symbolic closure.

## Validation

- `python3 -m py_compile check_reuse_radial_exact.py test_check_reuse_radial_exact.py`: passed.
- `python3 -m unittest -v test_check_reuse_radial_exact.py`: 5 tests passed.
- No build, SSH, or scientific execution was performed.
