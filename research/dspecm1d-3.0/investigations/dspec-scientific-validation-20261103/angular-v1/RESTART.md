# DSpecM1D angular convergence restart

Saved 3 November2026 after the human's brief sleep/save instruction relayed by CCDSM. The human immediately relayed **carry on** through CCDSM; the save hold is lifted. Continue the authorized study. This file is a checkpoint, not a claim that execution is complete.

## Objective and authority

Execute isolated L400/L600/L800 angular convergence on Tunya. Luna high implements, Sol6.1 high independently reviews. Establish the0.1% observed-change target only where supported; no automatic L>800, production integration, optimization or fullH1 reformulation. User emphasizes dependence on frequency filtering AND measurement window. Frozen fixture remains2/4/8/10mHz filtered transverse-velocity traveltime, source29.5km/solid-side receiver3km,60degree separation,1000..1200s window/40s taper,35FFT bins7..41 (32active/3zero-weight). Other illustrative bands are not new jobs or certified.

## Current scientific/build state

- Source frozen in candidate/dspec. Only driver and batched angular header changed from Study2: driver6f30aa9818d82dfaa50d4566a6f1e88444532ec7ccff4d552c1d73a2fd687988; header2b0f36e65fb7d81e41ee686013f91b10b6ff11ebd08cda322002d6179fb5fdbb. Controls unchanged88600b2bcfea127e4ff0ccc759d5e141e805fe6463dfe91885707da62343e743.
- New explicitly bounded angular-v1 control accepts only exact L400/600/800 cases and L800NQ8/10 controls; full78-angle route retained. New candidate32MiB geometry/64MiB angular cache+temporary budget, default16/32unchanged. No solver/angular mathematics changed. Wigner L800 within pinned double limit1827; transient recursion arrays accounted.
- Build PASSED on Tunya, one timer-disabled variant,96.7806s, no compiler stderr,13tracked/12external dependencies unchanged, cleanup verified. Binary SHA d7ffb6cbd3138b922831915906c91fb199e7a17adf103131572318ad1e2631aa.
- Build receipt receipts/angular-build-v1-evidence.json SHA dd5edb758835c57bfeed45df4ced0ef0e9b6318e0a1b2fbc13d88170e4670b2b. Registered slice peak2,667,933,696B, zero kernel max/OOM events. Build unit dspec-science-angular-v1-build.service last verified finished/inactive, MainPID0; historical livePID1283785 is NOT a running-job claim.
- **No angular numerical jobs launched.** Producer/checker stages/controls are proposed, not created yet. Inspect actual remote state before any action; never duplicate existing state/output.
- Exact8-job producer plan configs/angular-run-v1-plan.json is still under review. Checker tools/check_angular_v1.py was corrected for complete-case measurement arguments, higher-cutoff denominators, exact grids/tapers, non-gating99sample diagnostic, fullNQ6/8/10 comparisons without early report exit,32active/3placeholder unit-product bins, independent memory bounds. It is NOT sealed yet; see IMPLEMENTATION_HANDOVER.md and REVIEW_HANDOVER.md for final saved hashes/open findings.

## Completed independent review gates

- Plan/API/feasibility: reviews/angular-plan-feasibility-independent-review-20261103.txt SHA1bc147d24b519ba20db8e129375ee3ec6de2e9fa165dd571051e0e29c688c618.
- Private shared-resource registration: angular-shared-resource-registration-independent-review-20261103.txt SHA9aa919fbb70ebe7935a77b39a117fddb2ee761d252fe647eb2f9107651bd38ad.
- Source/build-config: angular-source-build-v1-independent-review-20261103.txt SHA2d3f1cef310ffe4239838e6cc47046e0a7945e3f9084ad73c00e44d5110259d9.
- Kernel launch hook: angular-kernel-launch-hook-independent-review-20261103.txt SHAb91aa780fa07443510f2e808d17a16901ba8c889b46b65089983487ac372b67e.
- Actual build: angular-build-v1-evidence-independent-review-20261103.txt SHAe0842d70542f6e3240105a8e6f357a2abe86735ff057b70879496c55fc3d2290.
- Numerical producer/checker exact paths/closure/resource-wrapped operators and AFTER-results/report claims remain pending. Do not rerun source/build review for checker-only changes.

## Shared RAM contract: human-authorized liaison completed

Actual chats: CCDSM01a18251-e70e-73a1-9212-9ab80eb289a1; DSpecM3D Implementation01a18879-2a7b-7810-bb21-9d8528e356a2; both host local. Our chat01a19297-d693-76a0-9c5f-43a42427320d.

Resource ADMITTED receipt resource-coordination/deployment_receipt.json SHAdc59cc0f82a84409e5338c2f7c73f37a683a82e4aaf6b060e2c0ce35f38503e1; read resource-coordination/ADMISSION.md before launch. CCDSM owner deployed exactly03421ab7cb49f3e28cfa08c3907ff264cae34d3a4b6f9ea170120b1976989d51 inventory-only extension to the SAME V4 guardian, tested coverage, preserved old state/actions/source and running native job. No second guardian.

- Human ceiling900,000,000,000B ACTIVE RAM; exclude reclaimable Linux/evictableARC using existing valid metric, not occupied RAM. Pause825GB/urgent875GB. Existing owned-work ceiling600GB retained; combined reservations318GB, separate slice caps339GB.
- Our separate dspecscienceangular20261103.slice:40,000,000,000B,12tasks,swap0,CPU800%. One build32GiB OR numerical campaign16GiB OR checker8GiB, exact service dspec-science-angular-v1-{build,run,check}.service.8outer/allinnerlibraries1. All descendants inside service kernel cgroup; no simultaneous M1D services.
- Pilot/M3D retains SEPARATE8GB/12task allowance. Existing CCDSM shared299GB/90task slice and260GB/67task child unchanged. Do not borrow/relabel reservations or change other processes.
- Before EACH service launch: fresh<5s validMONITORING/armed V4 state, no paused/frozen/unknown/conflicting owned units, active guardian exact sourcehash/epoch, exact slice/childcaps, no liveM1Dservice. Conservative forecast active+318GB<825GB. Atomic cross-service flock plus exact service names prevent overlap. Never thaw/kill/freeze a slice or start a second freezing guardian.
- Remote guardian source /space/adcm2/research/runs/CCDSM/extended_elastic_20261031_v1/sp12_20261101_v1/controller/monitor_tunya_global_ram_v4.py; service ccdsm-node-memory-guard-v4.service. State under same SP/resource_coordination_20261101_v1/global_ram_guard_v4/state.json. It restarts a new peak/actions epoch; prior epoch preserved, do not conflate peaks.
- operations/resource_admission.py SHA47311038ba6555e2df4331dc2f8f3f6f9360fb120b90d1db28ecfd694c161339 contains reviewed remote admission/launch functions, readonly inventory import derives root from actual guardian cgroup. Disable bytecode writes. Embed/import into new numerical operators; its resource_launch supports kind run/check. operations/build_stage.py775fc9dc0a6da2169593a1d9343478534860d61ab1a62159ebfa3418f3c54dec is already packaged/frozen; do not mutate or rerun build.

## Paths and next steps

BASE=/space/adcm2/research/runs/dspec-yspec-kernels/2026-10-27-tunya
Existing source stage BASE/scientific-inputs/dspec-science-angular-v1
Existing build /space/adcm2/research/builds/dspec-science-angular-v1
Existing build control BASE/dspec-science-angular-v1-build-control
Proposed validation stage BASE/scientific-inputs/dspec-science-angular-v1-validation
Proposed producer control BASE/dspec-science-angular-v1-control; task dspec-science-angular-v1-run-v1
Proposed checker stage BASE/scientific-inputs/dspec-science-angular-v1-check; control BASE/dspec-science-angular-v1-check-control; task dspec-science-angular-v1-check-v1

Eight jobs: angular-l400-dspec-nq6, angular-l400-yspec-nq6, angular-l600-dspec-nq6, angular-l600-yspec-nq6, angular-l800-dspec-nq6, angular-l800-yspec-nq6, angular-l800-dspec-nq8, angular-l800-dspec-nq10. Same fixed h and168probes,68along/10across; coarse wedges are diagnostic, not full angular bandwidth certification.

Parent now owns NEW numerical staging/finalization operators (not yet written as of this checkpoint); Luna owns exact producer/immutable-terminal-check plans, checker and minimal allowlist archives. Proposed manifests angular-producer-v1-manifest.json and angular-validation-v1-manifest.json. Align with Luna before creating operators; do not overwrite work.

1. Refresh actual shared guard/caps and remote state, preserve existing jobs. Read both agent handovers and current status/hashes.
2. Finish Sol gate for checker/exact8-job plan/reference closure/terminal immutable-state plan. Prepare+review minimal producer/checker archives and resource-wrapped operators.
3. Stage/verify, then launch exactly one kernel-contained16GiB campaign via reviewed helper. All numerical jobs sequential; no build overlap. Monitor guard/cgroup progress. Preserve any failed job; fix minimally with review.
4. Finalize exact completed producer state to an exclusive immutable snapshot; hash all declared data/logs. Run separate8GiB checker only from immutable snapshot, under fresh admission. Retrieve/hashverify raw evidence even if a checker fails.
5. Sol independently reviews actual metrics and acceptance. 0.1% requires two successive increments plus candidate-vs-L800 on taper-window waveform/normalization/P and89far-solid-interior probes+eightdepthstrata (higher-cutoff denominator, exactzero classes/no floors). Source-shell/interfaces/near-contact remain diagnostic. High-L3D action unavailable in current API, declared BEFORE results; no pretend volume norm/derivative validation.
6. Conditional NQ10 angular ladder is only if NQ6 radial gate fails AND NQ8->10 passes: add L400/L600 NQ10 (reuse800NQ10), total10<=12. Requires fresh focused source-whitelist/build/run review; not an automatic current queue extension. No mixed-NQ ladder or L>800.
7. Produce filter/window-specific tables, standard-style report and readable great-circle kernels/difference plots with LINEAR percentage axes and explicit denominators/saturation. No cutoff claim unless evidence passes; report measured limits otherwise. Notify CCDSM/Pilot when our40GB allocation is released at actual completion.

Transport: unchanged investigations/tunya-profile-20261027/ssh_transport.py; use exec_command require_escalated for SSH/SCP, existing /tmp/ccdsm-tunya-20261027/ssh_config. Never change auth/home/security state. No commits/pushes/deletion/production edits. Log every code change immediately.

## Latest execution update (supersedes the initial no-NUM checkpoint)

Stage/recovery now complete. Correct terminal plan is configs/angular-terminal-check-v2-plan.json584f9c94..., validation-v4; oldv1/v2/v3/stoppedattempt preserved. Numerical campaign dspec-science-angular-v1-run-v1 LAUNCHED once inside dspec-science-angular-v1-run.service; inspect receipts/angular-run-v1-launch.json for exact live PID/starttime and currentremote state before any action. DO NOT launch again or rerunstage. Current stage operator e896d21... and finalizer d511331... have focused review c27bf0.../b42e4d...; use finalizer only after all8clean jobs.

## FINAL NUMERICAL / RESOURCE CHECKPOINT — supersedes all launch instructions above

All8 jobs and separate terminal-v2 checker COMPLETED cleanly; DO NOT launch, stage, recover, finalize or check again. Raw producer142 files/checker8 files retrieved and independently hash-verified. AFTER-results review282595ff seals actual metrics: no tested cutoff passes0.1%; YSpec independently shows large degree dependence; NQ6 fails atL800 and NQ8-to10 fails50km0.187392%, so conditional extension trigger FALSE. No further numerical jobs authorized within this stage. Owner release receipt474455aa confirms logical40GB/12tasks RELEASED; all services inactive/tasks0. Preserved empty kernel slice is NOT future resource admission. Final local report/plots generated; five-page PDF e38a87cc rendered all pages and visual QA PASS. Only final Sol report-claim seal, delivery copy/manifest and completion message remain.

## DELIVERED / STOP — final checkpoint

Final Sol report-claim PASS15cfbfed obtained; exactbyte delivery PDF e38a87cc copied to output/pdf/dspec-angular-convergence-20261103.pdf and reports/delivery-manifest.json binds sources/figures/raw metrics/QA/reviews/resource receipts. No remaining stage work. Stop before any further numerical/source/integration work; wait for explicit next-phase authorization and fresh40GB logical admission (old reservation released). No cutoff certified.
