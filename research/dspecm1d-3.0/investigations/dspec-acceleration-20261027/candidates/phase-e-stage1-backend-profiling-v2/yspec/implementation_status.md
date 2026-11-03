# Phase E Stage 1 YSpec source candidate

2026-10-31: This fresh source tree was extracted from `git archive` of clean pinned commit `280eab0200495391539c0f97e62f978083c307bc`. The exact archive SHA-256 is `6cb81f5fe1c71a799def806f3107f27630f5a4419f38936fdcb7cfa6272ca11d`, recorded along with the pinned source-tree and library hashes in `phase-e-yspec-provenance.json`. The live YSpec checkout was not edited.

2026-10-31: Plan review approved a Stage 1-only instrumentation design: fresh timer-disabled and timer-enabled DSO variants from this same source, coarse function/degree/frequency scopes, thread-local recorder with one merge per existing frequency work item after all scopes close, one main-thread merge after join, and explicit DSO reset/finalize/write boundary. The YSpec API still computes its configured velocity parameter vector and the driver selects `vs`; formulas and scheduling remain unchanged. No instrumentation code, build, or numerical run has started.

2026-10-31 Phase E code change 1: added `YSpec/PhaseEProfile.hpp`, a fixed scope/work-item API with separate forward/kernel coverage. This defines instrumentation interfaces only; no numerical expressions or scheduling change. Implementation and driver wiring follow.

2026-10-31 Phase E code change 2: added `WorkItem::Complete` so the recorder can distinguish completed from skipped exact-zero bins and incomplete/failed items. It only records metadata after work and does not alter the API's operations.

2026-10-31 Phase E code change 3: implemented `PhaseEProfile.cpp` with fixed-array thread-local inclusive/exclusive scope stats and one aggregate merge per completed existing work item. It tracks attempted/completed forward and kernel bins plus exact-zero skipped bins. `FinalizeJson()` is intended to run after the outer API returns, when workers have joined. Instrumentation remains unbuilt and independent review is pending.

2026-10-31 Phase E code change 4: added the recorder translation unit to the fresh YSpec `YSpec` library source list in `CMakeLists.txt`; the disabled/enabled builds will use the same source list. No compilation has started.

2026-10-31 Phase E code change 5: added explicit `Scope::Stop()` so profiling scopes can bracket setup and transform calls without extending to their enclosing function return. Timer-disabled implementation remains a no-op.

2026-10-31 Phase E code change 6: instrumented YSpec seismogram synthesis at the existing frequency and degree work boundaries, depth-solver setup, and post-frequency transform. Each work item marks its exact frequency bin complete after its child scopes close; exception paths leave it incomplete. Arithmetic and OpenMP schedule are unchanged. No build has started.

2026-10-31 Phase E code change 7: instrumented sampled-kernel frequency work with exact-bin/zero-skip coverage, R/T/S radial solve, degree-batched strain-table creation, Wigner/Legendre table construction, angular evaluation, material contraction, ordered block reduction, and output transform scopes. All scopes wrap existing calls/loops; accumulation and traversal order are unchanged. No build has started.

2026-10-31 Phase E code change 9: instrumented the existing YSpec R/T/S radial solver calls and corrected the recorder’s OpenMP-master ownership boundary. Each frequency work item snapshots the master thread’s preexisting TLS counters and merges only its own delta, preserving caller/setup scopes for the final post-join caller merge. This prevents classifying setup as worker work or clearing the active outer API frame. No formulas or scheduling changed; build/run pending.

2026-10-31 Phase E code change 10: added driver-facing setup, measurement, and output-I/O scope IDs, plus actual YSpec returned parameter-name capture. The matched driver resets/finalizes the YSpec recorder only for the explicitly enabled YSpec-profile job and writes solver-option/layer/depth-group metadata beside the worker/caller scope totals. Native jobs use the same calls against the no-timer DSO. No numerical expressions or scheduling changed; independent source review pending.

2026-10-31 Phase E code change 11: compiled YSpec `Scope` and `WorkItem` instrumentation out inline in the timer-disabled build, and added a master-thread work-item time/count bridge to the enabled recorder so worker-zero frequency tasks are explicitly identified as overlapping the caller API scope. Scope formulas and scheduling are unchanged.

2026-10-31 Phase E code change 12: refined the master-thread bridge to sum the exact direct frequency-root scope inclusive deltas seen by the caller frame, rather than measuring the WorkItem object lifetime. This aligns the reported bridge with the caller child-time accounting while leaving caller wall and worker sums separate. Timer-only metadata; build/run pending review.

2026-10-31 Phase E code change 13: added fixed-size OpenMP thread-slot coverage to the existing WorkItem merge: each slot reports its attempted/completed forward/kernel bins, item counts, and exact-zero kernel bins. The slot is keyed by the existing OpenMP thread id and merged under the existing once-per-WorkItem lock; no per-scope event lock or dynamic per-event map was introduced.
