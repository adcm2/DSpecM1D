# Final short Sol review: isolated analytic-pole control and bounded L3 collection

27 October 2026. Read-only static review after `analytic-pole-control-v3-fresh-luna.md` PASS, against `fullkernel-reference-control-math-sol.md`. No implementation edit, build, scientific execution or SSH performed.

**PASS for the isolated control build, extracted receiver gate, supervisor-v4 delta and ten-job low-L3 pre-assessment plan.** This covers bounded collection only. The separate assessment adapter remains pending; no high-degree dispatch or numerical acceptance is released by this review.

Both original and control trees contain the same 12 files. The control manifest verifies 12/12. Their only source difference is the five-line early branch in `Windowed3DFields.h::harmonic`: degree>=1, valid m/n and exact theta=phi=0 return normalized Kronecker delta. The degree-zero path and all other arguments retain the verbatim old Wigner fallback. The old S/T receiver aggregation and scalar-call graph remain unchanged. This implements the proposed isolation experiment without copying the v6 receiver-matrix implementation or modifying frozen/production sources.

The generated gate contains all four independently extracted bodies verbatim, with matching recorded hashes: old harmonic, control harmonic, unchanged control/v4 receiver and v6 receiver. Every file hash in the gate provenance matches. The harness uses independently rounded long-double rho, asserts wider long-double precision, checks finite 3x3 S/T matrices at L1/2/3/12/384/1536, applies eight-ulp bounds to nonzero components and exact zeros to forbidden/zero-phase components, and compares control against both analytic matrices and v6. It correctly avoids unsupported degree-zero receiver calls. The existing receiver-v2 gate remains a separate first job with old-reference structural diagnostics.

The supervisor diff from v3 removes only automatic creation of declared output parent directories. Job directories still exist, output-relative-path/existence checks remain, and each producer now owns its `data` directory lifecycle. This is consistent with DSpec creating its target and YSpec requiring an absent target. Prefix CSV/receipt outputs go directly into the already-created adapter job directory. Logging/draining fixes, identity cleanup, locking, serial jobs, hash-bound evidence and maximum 50-thread userspace watchdog remain unchanged.

The dedicated build root contains only the native v4 analytic-control driver and standalone control gate. Commands retain reviewed compiler flags/dependencies and one direct sequential compiler per variant through the reviewed builder/receipt scripts. The low plan contains exactly ten jobs: two receiver gates, four original DSpec native/profiled L3 jobs, one control L3 job, the baseline prefix adapter and two YSpec L3 jobs. All physical solver degrees are 3, all numerical executables are pinned to CPU 0, numerical library threads are one and cap 50 remains. The plan has no assessment job and no high-degree solver invocation. The adapter/YSpec declared input lists contain literal `profile-04` paths, so this exact plan is bound to a fresh **profile-04** task, even though its command arguments use `{TASK_DIR}`. Do not launch it under another task ID without reviewing a corresponding path revision.

The original v4/v6 C/F comparison remains **FAIL at the unchanged 2e-13 criterion**. This new collection does not redefine that result. A later separately reviewed assessment must compare the complete control/v6 L3 outputs at the same fixed per-channel/column `2e-13`, preserve exact native/profile gates and metadata/energy checks, and save original/control/v6 localization evidence. If control/v6 passes, label it as that distinct comparison; do not label original v4/v6 as passed. Neither cancellation amplification nor a reference-bias explanation is established merely by this source review.

Reviewed digests:

```text
52219c41cd6b666e069567cfc7826e9bbd0587223b237e3bd9e9914c3384de1c  dspec/source-v4-analytic-pole-control/engine_headers/Windowed3DFields.h
cc9ec12748f7527e172be6c25e84f8758217dac2c4c3501e957a98a6d88142c7  dspec/source-v4-analytic-pole-control.sha256
5fb4bcb5247196ea776c362872b24f71da924aa760d93ed45e8851c5828d9837  dspec/validation-v3-control/control_receiver_gate.cpp
a6c7ea515bc144a41a81d6b05d2b30737696abd4902c984d7025b80af1ef6bd4  tools-v4/supervise_campaign.py
b2ab881c9985da2aff15fcab37eb2a5eeec72aa22ee291cc81de93d60f122346  configs-v6/build-control.json
0b80060bdcbe3cb777752f2815e466f3dc2c661d904161251f7f5efd20fc10dc  configs-v6/build-control-plan.json
789abf59c2a45ee1152d4393d9ecca7533311fd6078d76190acb3dfe6532ce90  configs-v6/low-qualification-plan.pre-assessment.json
```

Parent may perform this bounded build/collection within the existing authorization after exact transfer/hash verification, successful new binary receipts, confirmed runtime/CPU/disk readiness and exclusive detached logging. Preserve all earlier task folders and partial failures; use fresh build/task roots and the reviewed new supervisor. Userspace polling is not kernel aggregate enforcement. Stop after the ten low jobs pending assessment review and saved numerical results; physical `1e-3` thresholds and high-degree/full-global claims remain outside this bounded acceptance.
