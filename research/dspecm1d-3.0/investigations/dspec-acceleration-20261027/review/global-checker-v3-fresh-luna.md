# Targeted Luna High review: global checker v3

2026-10-27. Bounded independent re-review of the three prior HOLD items against the corrected checker, its tests, and `configs/global-run-v2-plan.json`. Earlier review notes are preserved. No SSH, build, or scientific execution was performed.

**PASS for the reviewed global run-plan gates.** The v3 checker closes each prior issue:

- The global-to-global path now requires byte-exact `kernel.bin`, `trace.csv`, and `spectrum.csv`, plus exact `energy_m2_s` and `energy_adjoint_qdotp` scalars (`tools/check_global_pair.py:302-315`, used at `:561-562`). The direct-reference path retains the `2e-13` numerical gate.
- The worker-1 versus worker-4 gate checks requested counts of 1 and 4, eight selected sides, and actual worker count equal to `min(requested, selected_sides, 48)`, matching the reviewed memory-nonbinding plan (`:272-299`). The plan specifies eight sides, a 32 GiB working cap, and a 64 GiB job cap for all four global jobs; its parallel checks compare `global-w1-l3` to `global-w4-l3` and `global-w1-l12` to `global-w4-l12`.
- The worker timing validator now requires `worker_elapsed_sum_s`, checks its sum against per-worker side durations, and checks each stage aggregate against its per-worker values within a tight roundoff allowance (`:463-508`). The producer schema in combined-v4 emits these fields.

The v2 plan parses as 12 jobs: four direct reference runs, four global runs built from the reviewed `kernel-acceleration-20261027-combined-03` binary path, and four check jobs. Its two global-to-global checks use the intended 1-worker reference and 4-worker candidate for L=3 and L=12.

Reviewed source hashes:

- `tools/check_global_pair.py`: `719f3581c24a1e9ebf4cd609ca6dbba1c67ec9fae09c5a8f3e1d160fe3dad74a`
- `tools/test_check_global_pair.py`: `9c3a7c9719d273ddfae93dae18870426649cd2ef3b61c2017fad2b7f4d041fff`

Independent verification: `python3 -m unittest -v test_check_global_pair.py` passed all 16 data-only tests, including exact output/scalar parity, worker-count rejection, and timing-rollup rejection. This confirms checker behavior only; no numerical comparison or speedup result is established.
