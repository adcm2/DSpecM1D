# Fresh independent Luna review: DSpec v6 and campaign tools v2

27 October 2026. Read-only review of the frozen DSpec v6 estimator change, pole receiver math and receiver-gate generator, plus `tools/` v2. Preserved `tools-v1/`, frozen source, prior review evidence, and transferred material were not edited. No build, scientific calculation, solver/Earth run, or SSH/remote activity was performed. Two tiny supervisor preflight probes used disposable directories under `/tmp`.

## Decisions

| Area | Decision |
|---|---|
| Frozen v6 Wigner estimator correction | **PASS — static source consistency only.** |
| Pole receiver formula | **PASS — algebraic/source review only.** |
| Generated receiver gate | **HOLD — generated C++ has a syntax error and does not implement the full analytic oracle requested by the prior review.** |
| Campaign tools v2 | **PASS — the five prior tool blockers are addressed by source inspection.** |
| Exact campaign configuration and remote readiness | **HOLD — no exact Tunya plan/build manifest exists for review.** |

## DSpec v6 estimator and receiver

Both 12-entry frozen-source manifests verify against their files: all v6 entries and all transferred v5 entries passed SHA-256 checks. The recorded v5 manifest digest, v6 manifest digest, and source-diff digest agree with `source-v6-provenance.json`.

The v5-to-v6 diff is limited to `engine_headers/Windowed3DFields.h` and `stage-to-math-map.json`. The generic all-orders sampler allocates `AngularWigner(maxAngularDegree,maxAngularDegree,min(2,maxAngularDegree),...)` and now estimates with `angularWignerBytesPerAngle(maxAngularDegree,maxAngularDegree)`. The capped pole and tensor allocation sites use `min(2,maxDegree)` for their estimate. `angularTableBytesPerAngle` forwards the same capped bound at its sole call site. The changed estimates use the existing exact `GSHIndices<All>(L,mMax,n).Size()` payload sum plus offset storage. The three preserved receiver timer names and mappings appear in the v6 JSON map and stage-map note. These checks support the specific memory-estimate correction; they do not establish build correctness or numerical parity.

The v6 `receiverAtPole` body is unchanged from v5. Its normalization is `sqrt((2L+1)/(4*pi))`; at the pole the nonzero harmonics are the matching order indicators for `m=0` and `m=±1`. The spheroidal radial row and both spheroidal/toroidal tangential combinations retain the previous signs and factors of `i`. This is algebraically consistent with the prior source review's six component matrices. No numerical execution was performed, so receiver equivalence and the advertised timing comparison remain unverified.

## Receiver-gate blockers

1. `validation/generate_receiver_gate.py` emits an unterminated include, `#include <iostream` (template near line 90); the generated `validation/receiver_gate.cpp` repeats it at line 17. The generated translation unit is syntactically invalid and cannot be used as the proposed gate until corrected and regenerated.
2. The gate only compares v4 values against v6 values. It does not compare either implementation with the analytic pole matrices, inspect forbidden entries against the prior bounds, or check pole normalization independently. Thus a shared defect could pass the v4/v6 comparison. The previous review's analytic and forbidden-entry checks remain outstanding.
3. The gate intentionally omits degree zero and directs radial degree zero to a future full physical gate. That boundary is not covered by this standalone receiver gate; preserve the previous review's restriction that radial l=0 uses its separate path.

## Campaign tools v2

Static review confirms the five blockers from `campaign-tools-luna.md` have been addressed: resume requires a task directory and valid schema-1 state with matching task ID and plan hash; cleanup logs pidfd signal outcomes, retries/escalates, reaps adopted children, and verifies identities gone or records `unclean`; the supervisor enables a Linux child subreaper and surfaces `/proc` inspection errors; README detached logs use a task-specific path with exclusive creation and append; the builder rejects wrappers, parallel-job flags and LTO and runs one direct compiler invocation per variant sequentially. The supervisor caps `thread_cap` at 50, forces common numerical thread environment variables to one, holds a campaign-wide nonblocking lock, and runs jobs sequentially. These are userspace observations and safeguards, not kernel aggregate enforcement.

Two safe `/tmp` probes confirmed that an existing task with missing state fails with `cannot resume: task state is missing`, and a configured cap of 51 fails with `thread_cap must be between 1 and 50`. An empty-job plan was rejected before campaign activity. I did not run cleanup/descendant probes in this review, nor test a detached session on Tunya.

The plan, exact compiler commands, complete source/dependency fingerprints, output declarations, resource thresholds, and remote detach policy still need review once written. The campaign config review is outstanding; the static tools PASS is not launch authorization. Preserve the watchdog caveat that short-lived threads can escape polling and RSS can double-count shared pages.

## Disposition

Do not accept the generated receiver gate until its C++ include is fixed and the analytic-oracle checks are added. The estimator and pole-formula findings are source-level passes only; the numerical gates listed by the earlier review and a compiled gate remain outstanding. Tools v2 resolves the listed design blockers, while the exact Tunya build/run plan and host/session evidence remain pre-launch gates. No implementation or source files were changed in this review.
