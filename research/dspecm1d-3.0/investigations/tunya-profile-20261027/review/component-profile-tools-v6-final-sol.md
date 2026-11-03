# Sol review: component analysis tools v6 — deployment HOLD

27 October 2026. Read-only source review after the fresh Luna intended-workflow PASS. No deployment, scientific execution, builds or SSH performed. Synthetic tests were not rerun because the parent withdrew this version from deployment while the material-aware reporting correction is prepared.

**HOLD deployment of tools-v6 due to a material reporting error.** Preserve this version and review; review the separately versioned correction through fresh Luna then short Sol before deploying the analysis waiter.

`collect_component_profile.py::cross_code_channels` compares DSpec A/C/F/L/N directly against YSpec A/C/F/L/N over **all** sampled radial rows. Under the saved representation contract, fluid DSpec values reside in kappa while YSpec retains A/C/F components. Direct fluid-row A/C/F comparisons therefore include a representation mismatch and can produce misleading channel maxima, even though the section is called descriptive. A descriptive label does not make unlike channels a meaningful comparison.

The narrow correction should compare matching A/C/F/L/N channels on jointly solid rows and DSpec kappa against YSpec A+C+F on jointly fluid rows. Retain absolute errors, scales, corresponding values and physical/point locations. Optional regular solid 4200–6200 km and fluid 2200–3200 km summaries can help localize errors, while broader solid/fluid summaries remain visible. No accuracy gate or existing tolerance should change. Parent reports remaining large solid C/F discrepancies; these remain unresolved scientific evidence and must not be explained away by the fluid correction.

The intended finalizer → trusted immutable snapshot → collector flow is otherwise consistent with bounded data analysis: exact expected job IDs, successful/input-invariant/cleanup-verified receipts, copied log/output hashes, fresh output roots, finite evidence, exact native/profile pairs, waveform and energy gates, q-dot-P and exclusive-timer closure, and pending versus complete status. The finalizer waits at most 12 hours and launches only its data collector after completed campaign state; it does not dispatch science or alter the queue. Completion receipt requires a complete passing summary. These source findings do not independently establish the current remote campaign state or final results.

Keep the fresh Luna limitations scoped to this intended trusted flow: standalone collector receipt paths are trusted because the finalizer creates safe keys; command/degree agreement is supplied by the hash-pinned reviewed supervisor plan plus parent verification, rather than by a general untrusted-input validator. Parent's stated preflight verifies the live I8 plan and state plan SHA against `1ecfc12aef8c146515082844564f652392e49fb57d868fc91e961e0f5c22c834`. No broader CLI hardening is requested by this review. These caveats are not the reason for HOLD; the all-material channel comparison is.

Reviewed withdrawn v6 identities:

```text
802cfeaaa8cbede81d2c0ed93bac1c064d92cea795ffa3430c06997d856c5a41  tools-v6/collect_component_profile.py
0e347c21adb38034b9f697da1794db6c9cf44d7e43f7dff242196551f24ca195  tools-v6/finalize_profile.py
ce33f018fa45f1baf0711e8f9693c610232b3ed7dcd504f4500af060ec183312  tools-v6/test_component_profile.py
```

Do not launch the staged v6 waiter. Preserve it unchanged and review v7's material-aware reporting delta; all scientific thresholds, original failures, incomplete-job handling and no-full-global-kernel limits remain in force.
