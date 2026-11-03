# Final Luna High follow-up: pair checker

2026-10-27. Fresh read-only review of the repaired checker plus its committed data-only synthetic suite. No SSH, build, or scientific execution was performed. The earlier `layout-and-checker-fresh-luna.md` HOLD is preserved.

**Verdict: HOLD.** The main repaired checks now cover the previously identified qdot/(2E) recomputation, physical metadata expansion, per-job domain arithmetic, and native radial-side/material identity. However, two calculation-mode fields remain omitted from the physical equivalence gate, and a receipted counterexample passes.

Reviewed source snapshot:

- `tools/check_pair.py`: SHA-256 `975af4d59cc6b427b28170281d1883a8bbf2e2609a53e7437eb5ee17fd944168`
- `tools/test_check_pair.py`: SHA-256 `8a4d623f07c413cc7f700b35568b680a74e998d3cdaa286523bdb1a5f7ba8b48`

The committed synthetic suite passed all 7 tests with `python3 -m unittest -v test_check_pair`. I also created and ran `/tmp/test_check_pair_review.py` as a focused extra case: both valid job receipts were refreshed after changing `pole_only_nonradial_adjoint` to `true` for the reference and `false` for the candidate, with `forward_only_baseline` fixed false. `assess(exact=True)` returned `PASS`, so the test's expected-Fail assertion failed. This confirms that the omission can allow non-equivalent adjoint calculations through even when all compared arrays happen to be byte-identical.

`PHYSICAL` is still missing `pole_only_nonradial_adjoint` and `forward_only_baseline`, both written in `physical_driver.cpp:433`. Add them to `PHYSICAL` and add regression coverage for mismatches in both flags. Then rerun the synthetic suite and obtain a fresh independent review. The checker also leaves byte exactness optional (`--exact`); for these unchanged-arithmetic packing/nodal routes, Tunya invocation must explicitly require exact outputs or encode that policy in the checker.

The previously reviewed candidate source finding remains bounded static PASS only. No implementation source was changed in this follow-up.
