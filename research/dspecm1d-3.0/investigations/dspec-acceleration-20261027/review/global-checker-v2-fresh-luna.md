# Fresh Luna High review: global checker v2

2026-10-27. Read-only review of the corrected global checker against combined-v4 producer output, its 14 synthetic tests, and `configs/global-run-plan.json`. No SSH, build, or scientific execution was performed. The earlier `global-checker-fresh-luna.md` HOLD remains preserved.

**Decision: HOLD before Tunya use.** The earlier schema mismatches are corrected: the checker accepts the intended L=3/L=12 and SEM-6 comparison, 180 total element/node sides versus 179 positive-radius rows, 52 physical output bins versus 53 spool records with exact-zero guard q at bin 52, and the global eight-side subset/layout. Its direct-to-global path compares trace, spectrum, energy/qdot and kernel channels at the existing `2e-13` limit; receipts, profile identity, spool/cache mode and factorization counts are checked.

Reviewed hashes:

- `tools/check_global_pair.py`: `903a6abdb7fcc661ec41ac2660e7e8b43f18c690d4e6deeaa2ff74a2c302929c`
- `tools/test_check_global_pair.py`: `0b73ea66ae788a185dc08ccb58d27a593249f934d7242bbda0933a3e7c5c8861`

Local verification: `python3 -m unittest -v test_check_global_pair.py` passed all 14 data-only tests. The plan has 12 jobs: six L=3/L=12 solver jobs (direct reference, global w1, global w4) and six checker jobs (each global-to-direct comparison plus the w1-to-w4 parallel comparison at each degree). The configured parallel cases request one and four workers over eight selected sides, with the stated 32 GiB working-memory cap.

Blocking gaps in the planned global-to-global parallel gate:

1. **Compare all unchanged physical outputs.** The `global_reference` branch (`check_global_pair.py:479-502`) requires byte-exact `kernel.bin` and `trace.csv`, but does not compare receipted `spectrum.csv`, `energy_m2_s`, or `energy_adjoint_qdotp`. The combined-v4 driver computes raw spectrum and energy audit before entering the global consumer (`candidates/combined-v4/physical_driver.cpp:536-557,574-578`); changing worker count only changes the consumer. Require exact spectrum and exact energy/qdot scalars for the current w1-vs-w4 route (or at minimum the fixed `2e-13` scalar gate if exact scalar representation is intentionally not required), and report those comparisons.
2. **Prove the requested worker counts were used.** `validate_worker_timings` checks manifest `workers == len(timings.workers)` (`:424-428`), but neither it nor `validate_global_result` matches that count to `result.global_consumer_workers_requested` (`:258-269`). The `--global-reference` path does not require the plan's reference/candidate to have actual counts 1 and 4. Require result requests and actual manifest counts to match the reviewed plan: w1 actual=1, w4 actual=4; for this memory-nonbinding plan, also require actual count to equal `min(requested, selected_side_count)`. Keep any future memory-capped worker policy separately reviewed.
3. **Validate combined-v4 timing rollups.** `validate_worker_timings` accepts the producer's per-worker records and stage sums but does not check their rollups against each other (`:417-440`); it also defaults a missing `worker_elapsed_sum_s` to 0.0. Require the field and verify it equals the sum of `side_elapsed_s`; verify each `worker_stage_elapsed_sums` value matches the corresponding sum over per-worker stage fields within a small roundoff allowance. This protects reported timing evidence without summing overlapping worker durations into wall time.

No claims about numerical results, speedup, compilation, or runtime are established by this review.
