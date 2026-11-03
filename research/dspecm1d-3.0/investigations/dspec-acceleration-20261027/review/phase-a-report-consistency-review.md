# Phase A report consistency review

27 October 2026. **Two minor corrections requested; substantive evidence, comparisons and boundaries PASS.** Read-only review of the report source, generator and saved datasets. No scientific source edit, SSH, solver run, report regeneration or compilation was performed. PDF rendering/layout was not reviewed here.

## Concrete corrections

1. `main.tex:154` / generator resolution paragraph: “degree-384 to 1536 energy change is about 2.93%” leaves its denominator unspecified. The saved energies are E384=8.824364556391933e-10 and E1536=8.573420696708428e-10. The usual coarse-reference change `(E1536-E384)/E384` is **-2.843761%**. The current 2.926998% is the absolute difference divided by the *fine* E1536 value. Use “energy decreases by about 2.84% relative to degree 384,” or explicitly identify the fine-reference normalization. This does not change the exploratory/nonconvergence conclusion.

2. `main.tex:65` / generator response-cache paragraph: “observed native peak RSS is 8.736 GiB” corresponds to the profiled v5 observation (8.7363 GiB). The unprofiled native v5 observation is 8.7352 GiB, rounding to **8.735 GiB**. Use “observed profiled peak RSS is 8.736 GiB” or quote the native value as 8.735 GiB. Logical response-cache storage 10.901 GiB remains correct and is clearly distinguished from RSS.

Update the generator and regenerate the TeX together so the corrections survive rebuilding.

## Checked and consistent

- All ten input-file SHA256 values recorded in the comparison JSON match the current referenced files. All 4,026 CSV component rows match their saved inclusive/exclusive/call values. All 119 native/YSpec appendix component rows match the JSON under the report's stated rounding.
- Main native timings agree with saved process-wall telemetry. The original unprofiled 536.442363 s and final native 313.360937 s support the opening 41.6% decrease. The separate instrumented baseline 532.983888 s is correctly labelled in its own table. The prior native/profile equality hashes establish the original trace/spectrum/kernel chain used in the opening claim.
- Isolated v4/v5 repeats are correctly used at degrees 128/384, with about 11.0%/6.5% v5 regressions explicitly retained. The report does not explain them away as overlap or assert an established compiler/allocation/cache cause. The rejected nodal-recovery route and unresolved global route are disclosed; “approved acceleration campaign is complete” therefore describes completed experiments, not universal success or global scientific acceptance.
- DSpec/YSpec baseline fairness is explicit: YSpec waveform is an independent run, its kernel run uses the frozen shared DSpec measurement baseline, those times are not presented as a measured sum, and DSpec includes its baseline. The comparison is on matching 64-angle/179-side diagnostic points, with different configured algorithms/tolerances and substantial local-kernel disagreements disclosed. The report does not treat the roughly 10.65x diagnostic timing ratio at 1536 as an equal-accuracy end-to-end result.
- Waveform agreement, energy agreement and the material-aware C/F/L kernel differences match the prior saved summary. Fluid bulk comparisons and structural zeros are not misrepresented as local-kernel agreement.
- Global full-angular outputs are clearly separated from the native 64-angle diagnostics. Eight/two/one selected sides, the fluid-only 1536 case, source-contact omission, and the absent 179-side global volume are stated. The 151.47-GiB full-volume estimate is consistent with 179*18,929,680*6*8 bytes.
- Generated/Stored process times and memory values match the receipts. The 270.44-GiB shared table is distinct from the 11.293-GiB response spool. The modest 5.0% Stored end-to-end saving and much higher measured peak RSS are fairly described. Thread configuration, 34 observed threads, cap 40 plus follow-up cap 4, and user-space rather than kernel aggregate resource accounting are explicit.
- The manual consumer table correctly includes worker residuals, worker scheduling/wall gaps, consumer residuals and process outside consumer. It is not added to overlapping MATHPROFILE scopes. Nested exclusive closure is distinguished from incomplete fine-grained manual timing coverage. No completed serial-1536 reference or 32x speedup is claimed.
- The frozen original-direct/global 2e-13 gates remain explicitly failed in L/kappa, despite passing internal thread/storage parity. No threshold relaxation, global scientific promotion or numerical convergence conclusion is implied.
- Resolution output is approximate bilinear angular interpolation without radial interpolation, with endpoint exclusions and 8–18% descriptive channel differences. Proposed section sampling (4L+1 around a full circle), native slice depths and the fixed physical configuration match `configs/proposed-figure-case.json`.
- The proposed figure configuration explicitly requires Phase B approval, and no figures/full volume were generated. The report ends by excluding automatic figure generation or further optimization. Deferred receiver lifts and Sobolev/Riesz work remain deferred.

## Review scope and source snapshot

This is a report consistency review only. It carries forward the final evidence review's limits: large binary payloads are remotely hash-bound rather than locally reread; high-degree global accuracy and convergence remain unestablished. A passing review after the two small clarifications must not be read as changing those scientific gates.

- `reports/acceleration-comparison-20261027/main.tex`: `876fbc6482202ff3945e2dd071b48f1f60fb9c9105b99245b01aaaa472a17509`
- `tools/build_phase_a_report.py`: `b3efbbc212714a70b00e83159dded859a90e86d0ae074788528c6749cd2e31f5`
- `analyses/phase-a-final-comparison.json`: `99f72bf301ac6d72c76ff275cccf53b309b9cb178a11d57562f3f1a52acb688d`


Parent resolution after review: the final report explicitly states the energy difference normalization (2.93% relative to degree1536;2.84% relative to384), and labels8.736GiB as profiled-native RSS. It also includes the measured small-degree32-thread slowdowns. No scientific evidence or thresholds changed. Final TeX SHA256: 46bb60fdfac182d8ee322cb0de1fe051dbd532e85678dae8198419f5ee5c5a24.

Final wording clarification: requested high-degree full-volume non-completion distinguished from completed low-degree whole-domain parity tests. Final TeX SHA256 bd7ba762ef377928f4a31012f7b9c3d5fa70c09f683a9e92fe6708365523890d. Numerical tables unchanged.
