# Fresh review: global pair checker

**Decision: HOLD.** This is a data-only static review of the checker against the combined-v3 producer schema. The implementation must be corrected and independently re-reviewed before build/science. No implementation files were changed and no science run was performed.

Reviewed SHA-256:

| File | SHA-256 |
|---|---|
| `tools/check_global_pair.py` | `99993f5d9440eb305ca2f8bfd935e923762ef685efc0f3c2f092a5fd473c4c0e` |
| `tools/test_check_global_pair.py` | `f1383cc4efc827b49e0fddbb0ed196faf7bd156cf621a2c6c2eadc85eb95ebf0` |
| `tools/check_pair.py` | `6f4d72e7da0489b6b3d0b6fd1d5c15d7ba3b7e6179a14ef453b7e8ccb762e646` |
| `candidates/combined-v3/physical_driver.cpp` | `357e0038c5e56188d24b47de99697c79fe4a1b3bd7e3b2b340cfb274b928d5cf` |
| `candidates/global-v3/global_kernel_consumer.hpp` | `ab741d8fda380b31c473a75034d100fd314c422f66031703df4b29c1301ec994` |
| `candidates/global-v3/pole_response_spool.hpp` | `39a8e1bd6de5d52257fdfcb365f68fb09c9609439696dbba3e9c81b7d3605ea4` |

## Blocking mismatches

1. `validate_physical_pair` hardcodes degree 3, SEM degree 12, 179 sides, and 53 output frequency bins. The producer cases use SEM degree 6; low-L comparison cases include L=3 and L=12. `result.json` reports 180 native element/node sides (`NE*NN`, including the center), whereas the positive-radius native table contains 179. Physical output `frequency_bins` is 52 (`last+1`); the spool separately contains 53 records because `endBin=last+2` retains the first out-of-band guard. These are distinct counts and the checker currently conflates/rejects them.

2. `validate_global_manifest` and `validate_direct_domain` require the positive-radius CSV row count to equal `result.native_element_node_sides`. The producer writes 179 positive-radius rows but reports the 180-entry `NE*NN` count. Validate the table count against its explicit positive-radius/available-side metadata and separately validate the 180 total-side field.

3. `validate_frequency_records` rejects zero `adjoint_mixing` at every bin. The retained guard is deliberately outside the band: `bandWeight` is zero there and the producer therefore emits exactly zero complex q weight. Continue requiring contiguous identity, finite complex values and bounded receiver-basis residuals; accept/require zero q mixing on the designated guard and retain the nonzero check for in-band bins.

4. Ensure the expected checked matrix of L=3/L=12 full complete angular grids and eight selected global radial sides is explicit. `--probe64` should remain a separate incomplete-angle probe contract, not be mistaken for the production complete-grid comparison. Preserve the existing `2e-13` numerical limit.

## What passed static inspection

The direct binary is angle-major with the selected radial side inside; the global binary is radial-side-major with theta-major angles. The checker maps native `(element,node)` identities and transposes/reorders into that reference layout before channel comparisons. The synthetic permutation/transpose test exercises this mapping. Required global/direct payloads are bound to their completion receipts and the checker rechecks receipt stability around comparison. `validate_energy` recomputes qdot/(2E); global route checks require spool enabled and cache/reuse disabled; factorization attribution counts the producer's 53 solved bins. Keep these checks when correcting schema expectations.

The synthetic suite passes: `python3 -m unittest -v test_check_global_pair.py` — 6 tests. It does not cover the mismatches above or validate a real producer manifest, guard record, 179/180 count split, or the full/probe case matrix, so it does not clear this HOLD.

## Required remediation before re-review

Make output-bin, spool-bin, and guard-bin semantics explicit; validate total native sides separately from positive-radius table rows; use the actual SEM and supported L case matrix; and add synthetic fixtures derived from the combined-v3 field names/counts, including the exact-zero guard q weight. Re-run the unit suite and provide the resulting checker/test hashes and producer manifest schema for a fresh review. No speedup or scientific result is established by this review.
