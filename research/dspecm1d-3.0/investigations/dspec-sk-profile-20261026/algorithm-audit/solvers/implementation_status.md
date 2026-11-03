# Solver prototype implementation status

This directory contains diagnostic-only solver alternatives for the approved
SK profiling campaign. The accepted production checkout remains read-only.

## Change log

- 2026-10-26: Read `algorithm-audit/theory-before-implementation.md`,
  `algorithm-audit/theory-solver-before-implementation.md`, and
  `algorithm-audit/protocol.json` before implementation. The theory supports
  equilibrated general-band LU, NaturalOrdering SparseLU with exact symbolic
  pattern reuse, and exact element-interior static condensation. No compilation
  is authorized until the parent assigns one of the two available compiler
  slots; numerical validation and benchmark evidence are therefore pending.
- 2026-10-26: Began a header-only backend layer matching the existing fast
  prototype's matrix/RHS interface (`compute`/`factorize`/`solve`/
  `solveAdjoint`). No claim of successful compilation or numerical equivalence
  is made yet.
- 2026-10-26: Added `solver_backends.hpp` with (a) the existing row-max
  symmetric equilibration, standard pivoted `ZGBTRF`/`ZGBTRS` band LU using
  `ldab=2*kl+ku+1` and shared forward/adjoint factors, (b) NaturalOrdering
  Eigen SparseLU with exact CSC-pattern checks before numeric-only reuse, and
  (c) a pivoted element-interior Schur elimination/recovery path with local
  rank, reciprocal-condition, and local-solve residual fallback gates. Added
  templated LTG-based SEM partition builders for radial, spheroidal, and
  toroidal families; endpoint DOFs remain traces. LAPACK symbols were found
  in the installed OpenBLAS library (`zgbtrf_`, `zgbtrs_`; package
  `libopenblas-dev 0.3.26+ds-1ubuntu0.1`, non-64-bit symbol ABI). This is
  static ABI inspection only; compilation, backend solve checks, true SEM
  partition checks, original-operator residuals, field comparison, and
  benchmarks remain pending the parent's compiler/runtime slots.
- 2026-10-26: Added separate accumulated stage timers for band scaling/packing/
  factorization/substitution, NaturalOrdering scaling/analysis/factorization/
  substitution, and condensation scaling/local factors/Schur assembly/Schur
  factorization/substitution/recovery. Condensation rejects cross-element
  interior couplings before elimination and reuses Schur symbolic analysis only
  after exact Schur CSC-pattern equality. This remains uncompiled and unreviewed.
- 2026-10-26: Added original-operator relative-residual and baseline-relative
  field-difference helpers. Their denominators use the actual matrix/solution/
  RHS or reference-field norms and handle the exact all-zero case explicitly;
  no `max(1,norm)` floor masks small physical fields. Timings and helpers have
  not yet been exercised.
- 2026-10-26: Added `sem_cases.hpp` to form the actual radial, spheroidal, and
  toroidal `H=K-z^2 M` matrices, source-treatment RHS, and element partition
  directly from SEM accessors. It checks dimensions and finite matrix/RHS
  values. Receiver RHS construction, end-to-end field/kernel comparisons,
  compilation, and real-case execution remain outstanding.
- 2026-10-26: Added `sem_solver_benchmark.cpp`, a bounded driver over the actual
  profile source/model, displacement-map adjoint weights, modal receiver RHS,
  and SEM source-treatment loads. It compares COLAMD reference, NaturalOrdering,
  general-band, and condensation forward/adjoint fields; emits original-system
  residuals, baseline-relative field errors, bandwidth/nonzeros, dimensions,
  stage timing fields, and process peak RSS. Input `f2` is extended only as
  needed for the selected first 1/3/12/48 raw bins. This driver has not been
  compiled or run; compiler slot is currently assigned but the parent has
  requested no numerical timings until its YSpec window closes.
- 2026-10-26: Parent granted a compiler slot. The first compile surfaced a
  constness/adjoint-view issue in NaturalOrdering SparseLU and dependent Eigen
  casts in the SEM case builder; after fixing those and the condensed Schur
  solver's mutable substitution state, the combined benchmark driver compiled
  successfully with the audit flags plus `-lopenblas`. Compilation alone does
  not validate numerical behavior; only bounded actual-SEM smoke cases are
  next, with no headline timings until independent process runs are approved.
- 2026-10-26: The first actual SEM smoke exposed an interior-recovery defect
  in static condensation: the already solved local base `Aii^{-1} bi` was
  incorrectly passed through the local inverse a second time. Recovery now
  subtracts only the solved trace correction for both forward and Hermitian
  adjoint solves. The smoke also showed that condensation timing initialization
  erased the measured full-system scaling stage; that timer is now preserved.
  These are unverified fixes pending a rebuild and repeated residual checks.
- 2026-10-26: Rebuilt the benchmark after the recovery/timer fixes (success).
  Re-ran the actual homogeneous-shell SEM case at LMAX=3 and one raw bin with
  one numerical thread and bounded memory/CPU. All 28 backend/family/degree
  records completed. Natural LU, band LU, and static condensation all had
  original-system residuals at reference levels: about 1e-20 for radial and
  spheroidal and 1e-17 for toroidal. Condensation's largest full-vector
  reference difference was 1.35e-11; band's was 1.45e-11. Both forward and
  adjoint solved fields were exercised; the condensation scale timer is now
  nonzero. This is a bounded correctness smoke only: frozen sampled-kernel and
  pointwise gates have not yet been evaluated, and no performance conclusion
  follows. Logs are preserved in `smoke-l3-b1-fixed/`.
- 2026-10-26: Extended the one-bin correctness smoke to LMAX=6. All 52
  backend/family/degree records completed. Across these records, the largest
  original-system residuals were 3.18e-17 forward and 5.34e-21 adjoint for
  condensation; its largest full-vector relative difference was 1.35e-11
  forward and 1.23e-11 adjoint. This still does not evaluate the stricter
  sampled kernel/pointwise gates. Logs are in `smoke-l6-b1.stdout.log` and
  `smoke-l6-b1.stderr.log`; timings are preliminary.
- 2026-10-26: Ran additional LMAX=3/bin=1 actual cases on the homogeneous
  ocean/seafloor deck and `prem.200.noatten.txt`. Both completed all backend
  families, including toroidal source solves in the solid source region.
  Static condensation passed its local gates and original-system residuals on
  both interface-bearing models; PREM radial/spheroidal residuals were below
  2.3e-20 and toroidal below 3.0e-17. Full-vector errors remained below
  8.5e-12. These are source/interface smoke cases, not the frozen sampled
  kernel/pointwise acceptance checks. Logs are `smoke-ocean-l3-b1.*.log` and
  `smoke-prem-l3-b1.*.log`.
- 2026-10-26: Completed a bounded PREM all-bin run through LMAX=12 and raw
  bin 48 (one thread, 3 GiB address-space limit, 120 s CPU limit). It emitted
  4,800 successful records: 1,200 each for COLAMD reference, NaturalOrdering,
  band LU, and static condensation; no backend failure/fallback was hidden.
  Every raw bin 1..48 and all degrees 1..12 were present. Maximum full-system
  residuals were 4.72e-17 forward / 1.93e-20 adjoint for condensation and
  4.97e-17 / 2.11e-20 for band LU. Maximum full-vector relative differences
  were 1.46e-11 / 6.35e-12 for condensation and 1.60e-11 / 9.51e-12 for band
  LU. This is an actual PREM source/receiver RHS solver packet, with
  preliminary in-process timings only; it does not establish the stricter
  sampled kernel/pointwise gates or independent-process headline timings.
  Logs: `prem-l12-b48.stdout.log` and `prem-l12-b48.stderr.log`.
- 2026-10-26: Updated the actual SEM benchmark driver to retain NaturalOrdering
  and condensation factors by family/degree across raw bins, use numeric-only
  refactorization when exact sparse patterns match, and emit explicit
  `symbolic_reused` evidence. This makes reuse observable in actual frequency
  sequences instead of only exposing an unexercised class method. Rebuild and
  sequence validation are pending.
- 2026-10-26: Rebuilt the updated driver successfully and ran a three-bin
  PREM sequence at LMAX=3. All 84 records passed with no failures. For all
  seven family/degree streams, NaturalOrdering and condensation reported
  `symbolic_reused=false` at bin 1 and `true` at bins 2 and 3. Maximum
  original-system residuals across nonreference backends were 4.45e-17
  forward and 5.50e-21 adjoint. Logs: `smoke-prem-l3-b3-reuse.stdout.log` and
  `smoke-prem-l3-b3-reuse.stderr.log`.
- 2026-10-26: Added `batch_response_cache.hpp` as a stack-scoped, explicitly
  bounded response cache for one objective batch. Its key includes SEM object
  identity, family, degree, and exact frequency bit patterns; each entry checks
  exact forward/unit-basis RHS signatures, retains only forward and solved
  unit-adjoint response blocks plus bounded signatures, and supports thin-QR
  reconstruction of weighted adjoint RHSs with an actual relative residual
  gate. Key/signature, capacity, rank, and span mismatches raise explicit
  exceptions for caller-controlled direct-solve fallback. API alignment and
  compilation/actual batch checks remain pending with the fast integrator.
- 2026-10-26: Added a small API smoke source for `BatchResponseCache` covering
  first compute, exact-key hit without callback, thin-basis weighted response,
  exact force-signature mismatch rejection, and out-of-span rejection. It is
  implementation validation only, not a scientific solver result; compile/run
  is pending the shared compiler semaphore.
- 2026-10-26: Kept the smoke fixture's complex mixing initialization explicit
  by assigning each Eigen matrix entry, avoiding ambiguous initializer-list
  parsing before the first compiler run.
- 2026-10-26: The first cache smoke compile found GCC's nested-type default
  argument restriction for `Limits{}`; replaced that default argument with a
  delegating no-argument constructor. Rebuild and execution remain pending.
- 2026-10-26: Recompiled the cache smoke through the shared two-slot semaphore
  successfully and ran it successfully. The checks cover exact cache hit,
  computed weighted response against a known thin mix, force-signature change,
  and out-of-span rejection. This validates the cache API mechanics only; the
  fast SEM batch integration and actual objective gates remain pending.
- 2026-10-26: Raised the explicit entry cap to 2,048 so a single bounded
  48-bin, LMAX=12 objective batch can represent all radial and spheroidal/
  toroidal family-degree-frequency keys without silent eviction. The 512 MiB
  byte cap remains an independent explicit stop/fallback bound.
- 2026-10-26: Added cache hit/miss counters and the maximum successful
  weighted-RHS basis reconstruction residual so the fast driver can report
  actual batch reuse and reconstruction quality. Parent's frozen production
  fixture review then supplied the decisive factor-backend result: the shared
  accepted EquilibratedSparseLU passes all 336 field and 240 contact groups
  per bridge, but band LU fails 252 field/106 contact groups (maximum scaled
  forward error 1.35e-8, contact L2 5.38e-9), NaturalOrdering fails 246/106
  (maximum field 9.21e-9, contact 6.78e-9), and condensation fails 258/104
  (maximum adjoint-strain error 1.17e-8, contact 8.36e-10). These alternatives
  remain diagnostic prototypes and are not accepted replacements, regardless
  of their small original-system residuals. No tolerance changes or extended
  LU refinement will be pursued.
- 2026-10-26: Extended the cache API smoke to check one miss, one hit, and the
  maximum successful basis-reconstruction counter. Revalidation is included
  with the fast agent's current shared-semaphore compile of the integrated
  batch path.
- 2026-10-26: Independently recompiled and ran the cache API smoke after adding
  the counters; both succeeded. The fast agent is concurrently compiling the
  integrated batch path via the shared semaphore.
- 2026-10-26: Fast v3 integration now exercises the cache in a bounded LMAX=3,
  one-bin, three-objective tensor-real case with contacts disabled for the
  basis-product comparison. Saved result reports 7 unique entries, 7 misses,
  14 hits, and maximum basis reconstruction 3.51e-16; same-engine field and
  objective errors are zero, kernel L2 1.14e-14 and scaled pointwise
  4.68e-13. This validates volume-only batch arithmetic on that copied-engine
  case, while batch contact-bearing behavior remains intentionally opt-in and
  unaccepted; a two-objective packet is still requested.
- 2026-10-26: Final disposition after the parent's independent frozen-field
  comparison: the equilibrated COLAMD `EquilibratedSparseLU` backend is the
  only accepted factor route. The band, NaturalOrdering, and static-condensed
  solver prototypes remain implemented and reproduce the actual SEM equations
  with very small algebraic residuals, but all three fail the frozen field and
  contact gates; they are diagnostic-only and must not replace production LU.
  No tolerance changes or extended iterative-refinement attempt is authorized.
  The response cache is a separate explicit opt-in: its object lifetime is one
  N-objective batch; keys use the SEM object's identity, family, degree, and
  exact real/imaginary frequency bits; stored force and unit-basis RHS
  signatures must match bit-for-bit. The cache retains only solved forward and
  minimal unit-adjoint responses plus signatures, with caps of 2,048 entries
  and 512 MiB. Missing keys compute both response blocks from one factor;
  hits reuse them. A weighted adjoint is formed by pivoted thin QR of the
  supplied unit basis, and its un-floored relative RHS reconstruction must be
  <=1e-12. Key/signature/rank/span/capacity errors are explicit exceptions;
  the caller must use the direct weighted-adjoint path or report failure.
  Batch output passed one same-engine L3/F1/3-objective volume-only smoke, but
  batch contacts remain disabled there because the tiny F contact term is
  sensitive to basis-product roundoff. This does not authorize a contact-safe
  replacement or production default.

- 2026-10-26: Added `evidence_manifest.json` with actual SEM packet paths,
  reusable-pattern evidence, fast batch-cache smoke metrics, and the parent
  frozen-production factor-gate disposition. Validated the manifest as JSON;
  no further numerical runs were made.
