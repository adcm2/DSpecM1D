# Final short Sol review: receiver gate v2 and profile-02 delta

27 October 2026. Static read-only review after `receiver-gate-v2-fresh-luna.md` PASS. No implementation edits, builds, numerical/scientific execution or SSH performed.

**PASS for the receiver-gate v2 source/math contract, dedicated build configuration and profile-02 plan delta.** This accepts the versioned validation correction justified in `receiver-normalization-math-sol.md`; it does not retroactively change the preserved profile-01 FAIL or establish that v2 has compiled/passed.

The independent oracle computes pi with `acos(-1.0L)`, uses long-double arithmetic and `sqrt`, then rounds rho once to double. A static assertion requires more long-double than double significand bits. Analytic S/T matrices retain the exact derived phase/sign/support pattern. New receiver shape and finiteness are checked; each present component must agree within eight ulps, and forbidden entries and zero phase components must be exactly zero. The new correctness requirement is not weakened.

Old-reference amplitude and old/new differences are explicit diagnostics. Old shape/finiteness, forbidden support at the existing `5e-15*rho` bound, correct nonzero sign and pure-real/pure-imaginary phase are still decisive checks. The separate old harmonic loop covers every valid n,m in [-2,2]; degrees 1,2,3,12,384,1536 and both receiver families remain. No arbitrary old amplitude envelope is introduced. `REFERENCE_DRIFT` versus `DIAGNOSTIC_WITHIN_2E-14` labels classify amplitude diagnostics; the final gate result separately requires all structural/new-oracle/payload checks to pass.

Independently extracted each of the four bound functions from its frozen header using brace matching: v4 harmonic, v4 receiver, v6 receiver and v6 Wigner estimator appear verbatim in the generated C++. Their recorded body hashes match. All five files listed in `provenance.json::v2_files` match their digests, including generator, template and generated source. Payload tests still compare capped/full exact `GSHIndices<All>` counts with the extracted estimator. Timed matrices now contribute their squared norm to a volatile checksum, preventing wholly unused results; these three-sample microtimings remain diagnostics rather than physical-kernel speed evidence.

The dedicated build config compiles only `validation-v2/receiver_gate-v2.cpp` into fresh root `kernel-profile-20261027-receiver-v2`, variant `receiver-gate-v2`, output name `receiver-gate`. That exact output path is used by the new plan. It uses the already reviewed v3 builder/receipt tools with a single supervised build job and cap 50; the six existing solver binaries are not rebuilt or replaced.

Parsed profile-02 plan equality check passed: after only changing literal prior-result paths from `profile-01` to `profile-02` and replacing the receiver executable/input path, the plan is identical to v3 apart from its explanatory revision note. The six solver executables, single-thread environment, serial CPU-0 scheduling, full L3-before-checker-before-64-angle-probe sequence, all physical arguments and strict L3 acceptance checker are unchanged. Checker SHA-256 remains `3af8a351c66aadd1bcc86a2cc04e303f0cdde27a4aba4b634cede871139901a7`.

Reviewed hashes:

```text
c1bfa12905783399af16c158f9b422fc1f4f81de9feff75d90073895ba63d3c9  dspec/validation-v2/receiver_gate-v2.cpp
b46b28ec67c26bdb7532ba3d0db96d880dd40b57e775a4fa8ad64a2b7bb8dac3  configs-v4/build-receiver.json
31204f3e41bb1bc4a333da33e5171a868b38bc760578ce401d7aa0e2f60c1ee2  configs-v4/build-receiver-plan.json
1e9f5440a27bce2a4a6657357128c38794e9efd968716195c0d7794cdc69443d  configs-v4/qualification-and-profile-plan.json
```

Parent may perform the dedicated reviewed gate build, then launch the staged profile-02 plan after successful receipt/current-binary verification and exact transferred source/config/checker hashes. Preserve failed profile-01, its log, all original gate source, frozen implementations and prior builds. Use a fresh task directory and the reviewed detached supervisor with an exclusive log; verify session/resource accounting as planned. The earlier `wait_and_launch.py` is explicitly profile-01/seven-variant-specific and does not cover this new one-gate rebuild/profile-02 handoff; use the reviewed supervisor directly after parent-owned prerequisite verification. No broader launcher edit is required by this review.

Numerical acceptance still requires the actual v2 receiver gate and strict complete L3 gates to pass. This review authorizes no frozen-source change, Wigner-library repair, physical-tolerance relaxation or high-degree complete-global-kernel claim.
