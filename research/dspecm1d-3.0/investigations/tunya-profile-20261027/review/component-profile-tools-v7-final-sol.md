# Final short Sol review: material-aware component analysis v7

27 October 2026. Narrow read-only review after `component-profile-tools-v7-fresh-luna.md` PASS and the preserved v6 deployment HOLD. No SSH, build, scientific solver execution or implementation edit performed. Independently reran the local data-only synthetic suite with bytecode writing disabled and numerical threads set to one: `component profile synthetic checks passed`.

**PASS for the corrected v7 collector within the intended trusted profile-06 → verified snapshot → data-only collector/finalizer workflow.** This resolves the v6 material-comparison reporting blocker. It does not validate unfinished jobs, resolve observed solid-kernel discrepancies or establish global scientific accuracy.

Verified SHA-256:

```text
894a9185530c34ff67c2863673884eb5d82b55ebbde6237e92f7531359a1d13b  tools-v7/collect_component_profile.py
0e347c21adb38034b9f697da1794db6c9cf44d7e43f7dff242196551f24ca195  tools-v7/finalize_profile.py
ddb2d02b3f99d853f93bbb7aa1af25051ab80580f0b152235ded1722b6363a64  tools-v7/test_component_profile.py
```

The collector delta is limited to material/radius-aware descriptive comparisons and supporting radial/shape/metadata checks. It requires the state-listed radial table, 179 sides, 64 angles, matching point counts, expected layouts and exact binary lengths before reshaping as angle × radial side × channel. Matching radial inputs remain separately checked by the unchanged collector workflow.

It compares A/C/F/L/N only on solid rows and DSpec kappa versus the explicitly labelled YSpec A+C+F sum on fluid rows. Inclusive regular solid 4200–6200 km and outer-core 2200–3200 km subsets supplement all-solid/all-fluid results. Selected-point maxima map back to original angle/radial indices with element/node/radius/material and both values; each region retains absolute differences, scale and counts. No unlike-fluid-channel maxima are presented as direct A/C/F agreement. All results remain descriptive with no kernel-accuracy tolerance/pass claim. The synthetic suite exercises the corrected fluid convention while retaining a deliberately visible solid mismatch.

Finalizer is byte-identical to v6. Existing exact native/profile checks, waveform/time/Hann/energy and relative-L2 `1e-3` gates, DSpec q-dot-P/(2E) `1e-10`, finite evidence, pending-job handling and profile exclusive-time closure at `max(1e-8,2e-12*root)` are unchanged. Completion requires all 26 expected successful/input-invariant/cleanup-verified jobs, safe hash-verified snapshot copying and complete passing collector output. The bounded 12-hour waiter executes only data analysis; it never dispatches a scientific job or modifies the queue. Failures/timeouts cannot create a completion receipt.

Fresh Luna reports that the existing immutable partial snapshot returns INCOMPLETE with four pending YSpec L1536 jobs and material-aware L128/L384 summaries. Those saved-artifact checks were not independently repeated here. Pending L1536 data cannot be inferred from smaller degrees. Parent-reported large solid C/F discrepancies remain unresolved; the reporting correction does not turn them into a scientific PASS. Any final summary PASS refers only to the explicitly listed component/operational gates. Contacts-off 64-angle diagnostics remain distinct from a full global physical kernel certificate, degree convergence or representative throughput claim.

The deployment remains intentionally bounded: use the trusted fixed profile-06 task and I10/tools-v7 files, retain the withdrawn I9/v6 copy as DO_NOT_RUN, and verify the live I8 plan and state plan SHA against `1ecfc12aef8c146515082844564f652392e49fb57d868fc91e961e0f5c22c834` before launch. Parent's manual command comparison and the reviewed pinned supervisor supply command/degree provenance; these scripts are not a general validator for crafted input receipts. The collector may consume only trusted snapshots produced by this finalizer, which creates safe relative receipt keys. No general CLI hardening is required for this stated scope.

Within those conditions, parent may launch the one-off detached finalizer on CPU 1 with one analysis thread and an exclusive fresh output/log path, leaving the scientific queue unchanged. Verify the exact deployed hashes and detached session, retain all failed/partial evidence, and report userspace resource-accounting limits honestly. No production/frozen source or scientific threshold is changed by this acceptance.
