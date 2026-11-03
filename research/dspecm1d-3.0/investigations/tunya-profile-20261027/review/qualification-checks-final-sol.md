# Final short Sol review: qualification checker and numerical launch gate

27 October 2026. Reviewed after `qualification-checks-v2-fresh-luna.md` PASS. Static inspection of final checker, relevant actual frozen/portable DSpec and YSpec driver schemas, and `configs-v3/qualification-and-profile-plan.json`. No implementation edits, compilation, scientific execution or SSH performed.

**PASS for the final checker and staged numerical-plan launch contract.** This supersedes the checker-related HOLD in `tunya-launch-sol.md` and the later path/portability supplements. It does not claim the build completed, that any numerical gate passed, or that scientific accuracy/performance has been established.

Reviewed artifacts:

```text
3af8a351c66aadd1bcc86a2cc04e303f0cdde27a4aba4b634cede871139901a7  tools/qualification_checks.py
e9a0e8eb947895c435c87697f9bcafca159b789b8b4614286bae4bb8080b9d7a  configs-v3/qualification-and-profile-plan.json
```

The prior acceptance blockers are closed:

- Both DSpec native/profile pairs require exact kernel bytes and exact numeric trace/spectrum values. YSpec native/profile kernel bytes and stable result metadata also require exact equality.
- V4/v6 kernel differences are scaled independently for A, C, F, L, N and kappa, at `2e-13`. Trace/spectrum comparisons use separate column scales, preventing time/bin magnitudes from masking displacement or complex-spectrum differences.
- DSpec layout parsing matches actual `point-major ... doubles[A,C,F,L,N,kappa]` metadata; YSpec matches its native-endian, angle-major then radial-node then A/C/F/L/N driver layout. Full L3 jobs require seven theta rows, complete target/actual point counts, all 179 radial sides, consistent binary lengths, and finite kernels/tables.
- DSpec native/profile stable result/domain metadata is exact after the documented path/timing/preflight exclusions. Cross-version checks require each listed physical option, including qex/DC/source/receiver and element/route flags, plus domain options to match explicitly. Energy and q-dot-P fields are compared independently rather than under a larger unrelated scale. Each DSpec result must meet `abs(qdotp-2E)/abs(2E)<=1e-10` with positive energy.
- Profile records check actual driver solver/objective identities, finite nonnegative values and exclusive-sum closure against the physical root. Failure emits FAIL and nonzero exit. The report is written exclusively and refuses to overwrite prior evidence.

YSpec kernel-mode baseline trace is an input finite check, not a computed YSpec waveform. Its `trace_rows` summary therefore counts baseline-input rows. Cross-code A/C/F/L/N and fluid-shell kappa comparisons remain descriptive; no global scientific accuracy threshold is inferred from them. The metric named `relative_l2_error` in the per-channel comparison uses maximum absolute channel magnitude as denominator, so interpret it as a scaled error norm rather than the usual relative L2 norm; acceptance is governed by the correctly implemented per-channel maximum-error gate.

The 30-job plan has the correct fail-closed sequence: receiver oracle first, four DSpec plus two YSpec complete L3 kernel jobs, then the checker, and only then L128/384/1536 probes. The reviewed supervisor stops on nonzero exit, missing output, changed inputs or cleanup failure, so a failed receiver/L3/checker cannot fall through to high degree. Every numerical executable is pinned to CPU 0; jobs are serial and numerical thread environments are one, under the campaign-wide maximum of 50. The data-only checker need not have numerical-executable affinity. All high-degree kernel jobs use explicit positive 64-angle diagnostics; separate YSpec waveform jobs are not misrepresented as kernel results.

Concrete launch conditions remain parent-owned: build-03 must finish successfully with a verified receipt and binary/source/dependency hashes; stage the reviewed exact checker and v3 plan at the bound remote paths; confirm CPU 0 and runtime readiness; use the fresh profile task directory and detached reviewed supervisor with preservation of every failed/partial task. Build03's reported live compilation and five observed threads are parent-provided progress, not independently verified here. Userspace polling remains supervision rather than kernel aggregate enforcement.

After those conditions, the reviewed staged numerical plan may launch. Numerical acceptance still requires its receiver and complete L3 gates to return PASS. Inspect the saved new receiver-scope call counts before interpreting component timings, retain serial same-CPU controls for performance claims, and keep the broader angular-convergence/full-global-materialization questions separate from these bounded diagnostics.
