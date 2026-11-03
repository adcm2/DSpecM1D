# Short Sol review: original-route data controls plan

2026-10-27. **PASS for this read-only data orchestration plan.** Static review only; no science, SSH, build, or implementation changes. The wrapper mathematics/data gate is separately approved in `global-original-checker-short-sol.md`.

The plan has one bounded job under its own control root. Its inline program parses and receives the intended four arguments: reviewed wrapper, original reference root, current Estimate candidate root, and new job output directory. It invokes the wrapper for `reference-l3`/`reference-l12` in `accel-reuse-run-01` versus `global-generated-w1-l3`/`global-generated-w1-l12` in `accel-global-benchmark-02`. Parent-reported completed/stable state is not freshly checked remotely here; the wrapper independently requires qualified completed receipts and unchanged artifacts for each comparison. No producer or resume command is present.

Both independent controls run sequentially with `check=False`, intentionally retaining the first FAIL while completing the second control. Reports are `l3-checks.json` and `l12-checks.json`; their exact names and `controls-summary.json` match declared outputs. Summary records both return codes and exits 1 if either is nonzero. The wrapper/summary use exclusive output creation in the new job directory. Missing reports still fail the supervisor's output qualification. This preserves and exposes failures rather than promoting them.

All four imported checker source files are declared inputs and are included by the reviewed staging helper's checker-file selection beneath `global-original-checks-v1`; command, cwd and staged paths agree. Source/data inputs are read-only to this job; saved numerical outputs are not write targets. The existing helper retains package/tool hash checks and detached launch; the supervisor records job/outputs and applies cap 6, 8 GiB and 600-second limits. Numerical thread variables remain forced to one. The stated control/direct/benchmark caps sum to 22 (6+4+12), without claiming aggregate kernel enforcement or freshly measured resource availability.

Plan SHA-256: `ecf5de25aed1b079600004a7ef465bef33383067be0c30a320a1ab0868c25b2e`.

This permits only the requested completed-data controls. The original `2e-13` gate remains fixed; actual cross-route status remains unresolved until the reports are assessed. No science producer, figure, numerical acceptance, speedup or global scientific promotion is authorized or established here.
