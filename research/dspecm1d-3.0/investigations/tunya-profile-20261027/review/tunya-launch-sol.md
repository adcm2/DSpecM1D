# Short final Sol review: Tunya launch gates

27 October 2026. Read-only review after fresh Luna's v6 tools, receiver-gate and exact-config reviews. Inspected frozen v6 delta, active generated receiver gate, exact build/run JSON, `prepare_configs.py`, `build_receipt.py`, builder, supervisor and the current qualification checker. No SSH, build, scientific execution or source edit was performed.

## Disposition

| Scope | Decision |
|---|---|
| Frozen DSpec v6 estimator and unchanged pole formula | **PASS, static source consistency** |
| Receiver-gate implementation | **PASS for build and execution as a qualification test** |
| Exact build-only plan and reviewed supervision tools | **PASS for build-only launch after parent verifies transferred inputs and host readiness** |
| Entire qualification-and-profile plan | **HOLD before launch: qualification checker still needs its fresh review and the findings below resolved** |
| Numerical equivalence, performance or full high-degree kernel completion | **Not established** |

The build-only PASS is independent of the unfinished checker: `configs/build-plan.json` invokes only `build_receipt.py` and the reviewed sequential builder. It launches no scientific executable. Preserve all existing build roots; its requested output root must be fresh and empty. Parent-owned remote checks must establish exact transfer hashes, source/dependency paths, compiler/runtime availability, write space and the detached supervised session. This review does not assert those remote facts.

## Source and receiver findings

The v6 fix resolves the earlier P1 mismatch without changing the allocation or numerical algorithm: the full-order sampler now estimates with `(L,L)`, while the pole and tensor sites estimate with `(L,min(2,L))`, matching each constructor's actual mMax. The estimate remains payload accounting, not whole-process RSS or a strict allocation ceiling: even a one-angle minimum can exceed the 32 MiB target. The v6 stage map now includes the three receiver scopes. Frozen v5 remains preserved.

The receiver body is unchanged from v5. With `rho=sqrt((2L+1)/(4*pi))`, columns m=(-1,0,+1), its S rows are `(0,rho,0)`, `(-rho,0,rho)`, `(-i*rho,0,-i*rho)`; T rows are zero, `(i*rho,0,i*rho)`, `(-rho,0,rho)`. The independently written analytic matrices encode those six channels directly. The separate harmonic-indicator loop checks every degree-valid n,m pair in [-2,2], including both signs at order 2, against `rho*delta(n,m)`. Present components require eight-ulp agreement, old forbidden components have the `5e-15*rho` bound, and new forbidden components require zero. Degrees 1,2,3,12,384,1536 are covered. The old/new comparisons also impose a scaled error bound. The corrected includes are complete by inspection.

Degree zero is deliberately absent from this standalone receiver test: the old helper's n=±1 requests are outside the degree-zero Wigner range. Radial l=0 must be covered through the actual physical route. Payload checks compare both full and capped order counts; because they reuse the library indexing primitive, they verify caller wiring and byte arithmetic rather than independently proving the indexing library. No numerical PASS exists until the compiled test succeeds. Its three-sample dead-result microtimings are diagnostic only and cannot establish a physical-kernel speedup.

## Build/config/supervision findings

`configs/build.json` has seven sequential direct compiler invocations: DSpec v4/v6 native/profiled, YSpec v6 native/profiled, and receiver gate. Native commands alone define `MATHPROFILE_DISABLED`; DSpec retains `DSFAST_DISABLE_PROFILE`. The commands retain separate PlanetModel dependency roots and required YSpec definitions. No matrix-BLAS backend macro, compiler job-parallelism flag or LTO is enabled. Source/dependency trees are fingerprinted before and after each compile. The receipt requires builder completion and rehashes each reported binary before exclusively creating its result.

The supervisor validates cap 1–50, fixes numerical thread variables to one, holds a campaign-wide nonblocking lock, runs one job at a time, and records command/input/output/log hashes. It tracks `(pid,start-time)` ancestry with subreaper adoption and pidfd signaling, verifies cleanup, and rejects incomplete resume entries rather than overwriting them. This is a periodic userspace watchdog: short-lived excess threads can escape polling, and RSS sums can double-count shared memory. It is not kernel aggregate enforcement.

The generated numerical plan pins jobs to CPU 0, uses explicit non-FTZ environment, runs the receiver gate and all six L3 jobs before the checker, then schedules bounded L128/384/1536 diagnostics. High-degree kernels are explicitly 64-angle probes. Availability of CPU 0 in Tunya's allowed affinity mask remains a parent-owned host check. No high-degree probe result can be represented as a complete global kernel.

## Numerical-plan blockers in current checker

1. **Instrumentation parity is tolerant for DSpec.** In `qualification_checks.py::run`, both DSpec native/profiled pairs pass `exact=False`, so numerical differences up to `2e-13` can pass and the report still sets `native_profile_parity=True`. The prior receiver review requires exact numerical instrumentation parity. Use exact kernel bytes and exact numeric trace/spectrum values (with byte equality recorded where applicable), or obtain an explicit change to that acceptance contract. V4-versus-v6 may retain the reviewed numerical tolerance.
2. **All six DSpec kernel channels share one global scale.** `compare_arrays` compares the interleaved file as one vector. A small A/C/F/L/N/kappa channel can therefore fail its own scale-relative gate while the combined comparison passes. Apply the v4/v6 tolerance independently to each of the six channels and record per-channel errors. The descriptive cross-code summaries do not repair this same-code qualification gap.
3. **The checker review is unfinished.** Fresh Luna's config PASS explicitly excluded it. Complete that independent review and then obtain the short Sol follow-up on the final checker before launching the whole plan. The DSpec energy-adjoint check, complete-file/179-side checks, finite-value checks and exclusive timer closure are useful existing safeguards, but do not remove blockers 1–2. Relevant DSpec energy/result metadata and expected receiver-scope call counts should also be verified/reported, as required by the earlier numerical-gate specification.

The narrowest concrete next action is the reviewed build-only plan under the authorized detached supervisor after remote preflight. Keep the numerical plan on HOLD until the checker findings and its review sequence are complete; then execute the receiver/L3 gates before allowing its later high-degree jobs. Preserve failures and partial outputs, do not edit frozen inputs, and make no numerical or speedup claim from these static passes.
