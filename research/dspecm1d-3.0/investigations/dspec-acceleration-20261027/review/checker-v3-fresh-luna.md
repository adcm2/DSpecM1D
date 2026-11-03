# Fresh Luna High review: pair checker v3

2026-10-27. Fresh read-only review of the pair checker and its synthetic regressions, following the prior HOLD reports. No SSH, build, or scientific execution was performed. Earlier review notes remain intact.

**PASS for the current packing/nodal matched-route comparison gate, provided the checker is invoked with `--exact`.** The two formerly omitted solve-mode fields, `pole_only_nonradial_adjoint` and `forward_only_baseline`, are now included in `PHYSICAL` (`tools/check_pair.py:12-23`). The new regression test independently changes each field, refreshes the candidate receipt, and requires rejection (`tools/test_check_pair.py:90-95`). This closes the demonstrated counterexample. The result still retains the `2e-13` comparison gate; `--exact` additionally requires byte identity for every compared output array and sampling table.

Reviewed snapshot hashes:

- `tools/check_pair.py`: `6f4d72e7da0489b6b3d0b6fd1d5c15d7ba3b7e6179a14ef453b7e8ccb762e646`
- `tools/test_check_pair.py`: `933605f851222ead419c89b9d50ef5a023623eb0b86d0840c230e38a8d6d3485`

Independent verification: `python3 -m unittest -v test_check_pair` passed all 8 synthetic tests, including receipt tampering, numerical regression, recomputed qdot/(2E), inconsistent sampled point product, native radial/material mismatch, cleanup qualification, exact-pair acceptance, and both solve-flag mismatches.

This PASS is scoped to current packing/nodal routes that use matching physical and solve methodology, with exact result-array equality requested. It does not pre-approve the later shared-response implementation; that route must carry its separately reviewed configuration and methodology evidence. No numerical evidence or compilation result is implied.
