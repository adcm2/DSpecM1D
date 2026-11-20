# Implementation Status

## 2026-07-29 — MeshModel gravity integration correction

- Added focused analytic regression coverage for constant density, linearly
  varying density, and a duplicated two-region density interface.
- Corrected each gravity quadrature contribution to evaluate the existing
  layer-density interpolant at the mapped quadrature radius, splitting only at
  density-spline knots inside the current integration interval.
- Preserved cumulative centre-to-surface storage and exact copying at duplicated
  element and material boundaries, without adding a new model-interface
  requirement.
- Verification completed:
  - all 3 focused analytic gravity tests passed;
  - all 53 registered project tests passed;
  - `git diff --check` passed.

## 2026-07-21 — Corrected upstream Interpolation pin

- Updated the project build and installed package configuration to pin
  `https://github.com/da380/Interpolation.git` at immutable merge commit
  `4508a7b3d2205f5647e87cdbaa4496cc0283c389`.
- This upstream commit contains the corrected cubic-spline implementation merged
  in Interpolation pull request 4 while preserving the existing include path,
  C++ namespace, constructor interface, and CMake target.
- Verification completed from a clean `/tmp/dspec-cubic-spline-fix-TqjqBV`
  build tree using the public HTTPS dependency URLs:
  - CMake fetched Interpolation at the exact pinned commit.
  - The full project build passed.
  - All 50 registered tests passed: 46 unit/component tests, 3 smoke tests,
    and 1 migration reference test.

## 2026-07-15 — GitHub Actions dependency checkout fix

- Removed shallow cloning from commit-pinned `FetchContent` dependencies in the
  project build and installed package configuration.
- Retained shallow cloning for Eigen because it is pinned by the `3.4.0` tag.
- Reason: a shallow clone of GSHTrans could no longer check out the older pinned
  commit after its default branch advanced, causing `build-test-smoke` to fail
  during CMake configuration.
- Verification completed from a clean `/tmp/dspecm1d-ci-fix` build tree:
  - CI-equivalent CMake configuration passed with freshly cloned dependencies.
  - Full build passed.
  - All 46 unit/component tests passed.
  - All 3 smoke tests passed.
  - The `website` documentation target passed.

## 2026-10-01 — Unsupported historical benchmark reports

Copied benchmark_report.pdf and pkp_diffraction_report.pdf from Documents/PhD_Codes into benchmarks/unsupported/reports without altering or removing the originals. Added provenance.json with SHA-256 hashes and README.md describing historical scope, limitations and a future reproducible QSSP benchmark plan. Verified both copies byte-for-byte. No numerical source, build configuration, supported benchmark targets or website content changed; no simulations rerun.

## 2026-10-07 — Cowling Stage 1: two-field LTG declaration

- Added the `SEM::ltgSC` declaration for the compact `(U, V)` Cowling map; no
  existing map or behavior changed.
- Implemented the compact map as `2*element*(NN-1) + 2*node + field + offset`;
  at the first node above a fluid-solid boundary, U subtracts one from the
  shared offset to remain continuous, while V keeps the extra interface DOF.
  The existing `ltgS` implementation was left unchanged.
- Added independently numbered quotient-topology regressions for the existing
  solid-fluid-solid and all-solid fixtures. They check the DOF formula,
  contiguous coverage without gaps, ordinary endpoint sharing, U continuity,
  and V continuity only across same-material interfaces.
- Added a PREM topology case that first confirms at least one solid-solid
  material boundary is present, then checks endpoint sharing across it.
- Verification: the focused debug `SEMComponentTests.CowlingMap*` filter passed
  all 3 tests with `OMP_NUM_THREADS=1`. Running the full debug component binary
  separately reaches the existing `ltgT` assertion in
  `ReceiverVectorsHaveExpectedDimensions`; no toroidal code was changed.
- Parent verification: all 52 configured Release unit/component tests passed.
  The Debug receiver assertion was also reproduced using the original HEAD's
  headers and test source; it is a pre-existing limitation.

## 2026-10-07 — Cowling Stage 2: stored spheroidal matrices

- Extended the existing spheroidal quadrature loops to accumulate compact Cowling
  stiffness, attenuation, and inertia bases alongside the unchanged full-gravity
  triplets. Cowling retains the elastic and background-gravity terms while
  omitting potential couplings, potential stiffness, the exterior potential
  boundary term, and the `4*pi*G*rho^2` self-gravity contribution.
- Added the `hSC(l)`, `hSCa(l)`, and `pSC(l)` accessors with the existing
  `k^2` decomposition. The compact matrices use the Stage 1 `ltgSC` layout.
- Added all-solid and solid-fluid-solid matrix checks and a direct weak-form
  oracle for an all-solid interior U diagonal entry. The oracle evaluates
  quadrature/derivative contributions directly and checks separately that the
  full matrix adds the expected perturbed self-gravity term.
- Verification is in progress; final commands and results will be recorded
  after the focused checks complete.
- Test adjustment: the public mesh accessor is const and does not expose `GLL()`
  as a const method, so the independent oracle constructs the same GLL rule from
  the requested node count directly. The initial focused build caught this
  test compile error; production code was unaffected.
- Focused test refinement: center-node inertia is exactly zero because the
  spherical mass weight contains `r^2`; tests now require nonnegative diagonal
  mass, positive entries away from the center, and the expected center zero.
  Attenuation and elastic matrices have equal sparsity on the chosen fixture, so
  the test checks equal pattern size and that both remain sparser than full
  gravity matrices instead of assuming attenuation has fewer entries.
- Oracle indexing correction after parent observation: the direct test now
  evaluates the fixed test basis derivative at each quadrature node explicitly,
  matching the weak-form integral orientation without reusing production's
  `mat_d` indexing. This supersedes the initial transposed test helper.
- Focused Debug verification: all 6 `SEMComponentTests.Cowling*` tests passed
  with `OMP_NUM_THREADS=1`, including the corrected weak-form oracle and both
  all-solid and layered matrix checks.
- Full Release verification: the configured build completed and `ctest` passed
  all 55 tests (55/55).

- Structural matrix measurements at l=2 on the controlled fixtures:
  - all-solid: full 39 DOFs, 453 nonzeros, lower/upper bandwidth 11/11;
    Cowling 26 DOFs, 244 nonzeros, bandwidth 7/7.
  - solid-fluid-solid: full 59 DOFs, 685 nonzeros, lower/upper bandwidth
    12/12; Cowling 40 DOFs, 370 nonzeros, bandwidth 8/8.
- Parent regression against original HEAD headers: full `hS`, `hSa`, `pS`,
  source vectors/coefficients and receiver interpolation dumps were byte-identical
  at l=1,2,7 on PREM with and without attenuation material coefficients.

## 2026-10-07 — Cowling Stage 2 review finding

- Fresh independent Stage 2 review returned `CHANGES REQUIRED` for test coverage
  only. It found that mass and attenuation were checked structurally, while the
  direct weak-form oracle covered only a diagonal stiffness coefficient.
- Scope of correction is limited to the existing focused test: add direct
  quadrature coefficient checks for one mass entry, one nonzero attenuation
  entry, and one off-diagonal U-V stiffness entry. No production changes or
  optional review suggestions are included.
- Added the three requested direct coefficient oracles to the existing
  all-solid weak-form test. Mass checks U and V entries including the `k^2=6`
  angular factor at l=2; attenuation checks a nonzero U-U quadrature integral;
  off-diagonal U-V stiffness checks the l=1 coefficient including retained
  `rho*g*r` coupling. Production code is unchanged.
- Verification after the review correction: all 6 focused Debug Cowling tests
  passed, the Release SEM component target rebuilt, all 55 Release tests passed,
  and `git diff --check` passed.

## 2026-10-07 — Cowling Stage 3: source and receiver entry points

- Extended the existing spheroidal source and base receiver method declarations
  with an optional Cowling layout selector, defaulting to the unchanged full
  `(U, V, P)` layout.
- The existing four-column radial source assembly now writes its U/V entries
  through either `ltgS` or `ltgSC`; the angular source coefficient methods are
  reused unchanged. The existing full receiver interpolation similarly switches
  only its radial bounds and U/V indices, leaving interpolation and angular
  weighting unchanged.
- The source vector dimension is computed locally in the extended method, so
  choosing the compact layout does not store formulation-dependent size state.
- Added a regression that compares every mapped U/V source row exactly, checks
  full and compact vector dimensions, compares every corresponding receiver
  column, and verifies displacement evaluation agrees when the two layouts
  contain the same physical U/V solution.
- This work is on `dev/cowling` against immutable baseline
  `889452b6790d2c60afacc156b4c7512990758aee`; checkpoint commit/push handling is
  reserved for the parent agent after the Stage 1–3 review checkpoint.
- The first focused build identified that the SEM-local complex alias is private;
  the test now uses `std::complex<double>` directly. The next compile check also
  confirmed the radial mesh exposes its surface through the final element's
  upper radius, which the synthetic receiver solution now uses. Verification
  then completed: all 7 `SEMComponentTests.Cowling*` tests passed in Debug, and
  the full Release build and all 56 Release tests passed. `git diff --check`
  passed. The separate Debug `ReceiverVectorsHaveExpectedDimensions` assertion
  remains the previously reproduced baseline `ltgT` range assertion.
- The Stage 3-only diff against the Stage 2 snapshot is saved at
  `/tmp/dspecm1d-cowling-ZwV33K/stage3-only.diff` (116 insertions, 12 deletions
  across `SEM.h`, `SEMForceSpheroidal.h`, `SEMReceivers.h`, the SEM component
  tests, and this status file).
- Parent original-HEAD regression after Stage 3: full matrices, source vectors
  and receiver interpolation dumps at l=1,2,7 remain byte-identical on both
  PREM fixtures. Full single-SEM and multi-SEM spheroidal spectrum dumps are
  also byte-identical for l=1..4, 2..8 mHz, two receivers, attenuation off/on,
  and one OpenMP thread. These checks cover full-path preservation only.
- Independent Luna-high gates completed: Stage 1 passed with non-blocking
  notes; Stage 2 passed a fresh rereview after the required test-only additions;
  Stage 3 passed. Integration Review A must inspect the complete accumulated
  diff against the immutable baseline before the first remote checkpoint.

## 2026-10-07 — Cowling Stage 4: start-index / truncation support

- Added an optional Cowling layout selector to each `allIndicesSph` overload.
  All overloads retain the existing `startElementSph` physical selection and
  `nskip` cadence, mapping the selected element's first U DOF through `ltgSC`
  only when requested; existing calls continue to map through `ltgS`.
- Added explicit `<vector>` and `StartElement.h` includes in the SEM component
  test source to support direct truncation regressions.
- Added checks for all three `allIndicesSph` overloads, both formulations,
  source-constrained and unconstrained starts, `nskip=1` and cadence reuse,
  plus the source-element-zero edge case. The fixture asserts that its chosen
  finite positive frequencies produce a nonzero selected start element.
- The first run showed degree 4 did not truncate on this mesh; increased the
  test degree further to 100 after that guard still observed no truncation.
- Verification: all 8 focused Debug `SEMComponentTests.Cowling*` tests passed
  with `OMP_NUM_THREADS=1`; the full Release build completed and all 57 Release
  CTest tests passed. `git diff --check` passed. Focused and full Release logs
  and the Stage 4-only diff are saved under
  `/tmp/dspecm1d-cowling-ZwV33K/stage4-{debug-focused,release-ctest}.log` and
  `/tmp/dspecm1d-cowling-ZwV33K/stage4-only.diff`.

## 2026-10-07 — Cowling Stage 5: standalone single-SEM solve

- Added an opt-in `cowling` boolean to the low-level `SpectraRunContext` plus
  single-SEM solver overload. It defaults to false; the existing public and
  legacy entry points therefore retain the full-gravity calculation.
- In the existing spheroidal loop, the selected formulation now supplies its
  own matrices, global and receiver bounds, source vector, and frequency
  truncation indices. Both formulations continue through the same attenuation,
  sparse solve, receiver accumulation, and output code.
- Added a PREM single-SEM regression that exercises both layouts across the
  selected frequency band, checks default output equality, finite/nonzero
  results, a nonzero Cowling start index, and the attenuation-enabled Cowling
  path. It also samples full/Cowling relative differences in increasing
  frequency order and reports their endpoint trend without imposing a test
  acceptance threshold.
- The first end-to-end Debug run exposed that Cowling has two displacement
  fields but retains four source right-hand sides. The receiver product now
  uses all four force columns, as in the existing full-gravity solve.
- The PREM comparison band is 5–80 mHz at l=1–2. A separate one-degree l=100
  Cowling solve checks a nonzero start element and finite, nonzero receiver
  output after reduction.
- Initial sampled differences decreased overall from the first to last bin but
  fluctuated between bins; the test accepts only this fixture-specific endpoint
  relationship and does not require monotonic differences.
- Fresh Stage 5 review required the test to enforce the observed endpoint
  trend. Added one fixture-specific assertion that the final sampled relative
  difference is below the initial one, after the existing two-sample guard;
  intermediate fluctuations remain allowed. No production code changed.
- After the review correction, the focused Debug solver test passed with one
  OpenMP thread, and the matching Release CTest passed (1/1).
- Focused Debug `PreferredSolverApiTests.StandaloneCowlingSpheroidalSolveOnPrem`
  passed with one OpenMP thread. It exercised 15 PREM elements with 3 GLL nodes
  per element, a 5–80 mHz requested band (sampled bins 6.25–75 mHz), l=1–2,
  attenuation off and on, and an actual l=100 reduced solve. Default context,
  explicit full-gravity, and existing API outputs were identical on this run.
- Relative full/Cowling spectrum differences by sampled frequency were:
  6.25: 2.08132e-5; 12.5: 1.33465e-5; 18.75: 8.95589e-6; 25: 5.74779e-6;
  31.25: 4.80171e-6; 37.5: 6.83836e-6; 43.75: 7.97056e-6; 50: 6.29792e-6;
  56.25: 3.70785e-6; 62.5: 5.74482e-6; 68.75: 9.42314e-6; 75: 1.00119e-5.
  These coarse-mesh, low-degree measurements fluctuate between bins; they are
  diagnostic evidence for this fixture, not broad waveform validation.
- Full Release build passed and all 58 CTest tests passed with one OpenMP
  thread; `git diff --check` passed. Stage 5 independent review is pending.
- Original-HEAD regression probes are byte-identical for full matrices, source
  vectors and receiver vectors at l=1,2,7, and for full single-/multi-SEM
  spheroidal spectra over l=1–4 and 2–8 mHz with attenuation off and on. The
  evidence is in `/tmp/dspecm1d-cowling-ZwV33K/stage5-baseline-regression.log`.

## 2026-10-07 — Cowling Stage 6: configurable cutoff setup

- Added an optional final ordered parameter value, `cowling_frequency_mhz`,
  defaulting to zero (disabled). The parser accepts comments/blanks around the
  optional value, rejects malformed values and trailing tokens, and validates
  that the cutoff is finite and nonnegative. Added getter/setter forwarding
  through `InputParametersNew` for programmatic configuration.
- Added a `FreqFull::timeNorm()` accessor so cutoff comparisons can convert
  the frequency helper's dimensionless cycles to physical mHz using its own
  normalization. The cutoff is documented as applying to the single-SEM path
  at this checkpoint.
- Updated the single-SEM spheroidal solve to route bins below the configured
  cutoff through full gravity and bins at or above it through Cowling. Both
  regions use fresh solver state; each frequency region restarts its own
  `nskip` counter and computes its own start-index cadence anchors. The
  standalone all-Cowling selector remains available.
- Documented the optional physical-mHz value and its single-SEM-only scope in
  the parameter-file guide. Added parser coverage for the disabled default,
  valid optional values, malformed/trailing input, and setter validation.
- Added single-SEM cutoff regressions for disabled and above-range exact
  preservation, all-Cowling selection below range, and a mixed band whose exact
  boundary and adjacent bins are checked against standalone formulations.
  These run with attenuation both off and on. Added a high-degree mixed-band
  case for `nskip=1,2,3` plus a cadence longer than the band. It verifies
  nonzero truncation on the full-gravity side, each layout's region-local
  endpoint anchor, and finite/nonzero mixed output.
- Adjusted the high-degree test cutoff to begin near the low end of the band
  after the first run showed the Cowling-only high-frequency segment did not
  reach a truncated start element at the original switch bin.
- The high-degree test uses a one-bin full-gravity region so that its cadence
  restarts at the cutoff for `nskip=2,3`; this fixture's resulting global start
  index also happens to match the old whole-band anchor, so the test directly
  checks the region-local start mapping rather than claiming a differing index.
- Added one cadence larger than the entire band. In that case the whole-band
  anchor is reused throughout, while the one-bin full region must recompute its
  nonzero start at the switch; the test compares those indices explicitly.
- For this long-cadence case, the mixed solver's full-gravity boundary bin is
  also compared exactly with an independent full solve using `nskip=1`.
- Verification on the configured build trees: Debug parser and solver API
  tests passed 15/15, Release parser and solver API tests passed 15/15, and the
  full Release CTest suite passed 62/62 with `OMP_NUM_THREADS=1`. `git diff
  --check` passed. Independent Stage 6 review is pending.

## 2026-10-07 — Cowling Stages 4–6 review checkpoint preparation

- Fresh independent Luna-high gates passed: Stage 4 `PASS`; Stage 5
  `PASS WITH NON-BLOCKING NOTES` after its required endpoint-trend test
  correction and fresh rereview; Stage 6 `PASS WITH NON-BLOCKING NOTES`.
  Non-blocking notes were reported without automatic code changes.
- Parent original-baseline regression after Stage 6 remains byte-identical for
  full matrices, source vectors and receiver weights at l=1,2,7, and full
  single-/multi-SEM spheroidal spectra at l=1–4 over 2–8 mHz, attenuation
  off/on, with one OpenMP thread. Evidence:
  `/tmp/dspecm1d-cowling-ZwV33K/stage6-baseline-regression.log`.
- All 62 Release tests and the 15 focused Debug parser/solver tests passed.
  The known original-baseline full-Debug receiver `ltgT` assertion is unchanged.
- Production SEM headers from Stages 1–3 are unchanged from the approved
  checkpoint. The approved direct `<vector>` test include was added in Stage 4.
- Review B must inspect the complete feature against fixed baseline
  `889452b6790d2c60afacc156b4c7512990758aee`, including earlier stages again.
  The previous checkpoint remains `5bc5ed840b6466b9c7af60bb5cda95ef41027e78`;
  Stages 4–6 remain uncommitted on `dev/cowling`. Stage 7 is not authorized
  before the next human review.

## 2026-10-07 — Low-level Cowling selector naming

- Renamed only the `SparseFSpec::spectra(SpectraRunContext, SEM, bool)`
  parameter from `cowling` to `forceCowling` and updated its local uses. The
  API comment states that true forces all spheroidal frequencies through
  Cowling, while false follows the configured cutoff, including zero-disabled.
- Added a test comment naming `forceCowling=true/false` and clarifying the
  existing explicit-selector behavior;
  no generic Cowling result variables, layout selectors, or semantics changed.
- Verification: rebuilt the Debug `dspecm1d_solver_api_tests` target with `-j2`;
  all 8 `PreferredSolverApiTests` passed with `OMP_NUM_THREADS=1`. Build and test
  logs are saved at `/tmp/dspecm1d-cowling-ZwV33K/rename-debug-build.log` and
  `/tmp/dspecm1d-cowling-ZwV33K/rename-debug-preferred-tests.log`.

## 2026-10-07 — Cowling Stage 7 implementation

- Extended only the multi-SEM spheroidal solve to split each existing chunk at
  the configured physical-frequency cutoff. Each nonempty region selects its
  full-gravity or Cowling matrices, maps, source and receiver operators while
  retaining the original chunk assignment and output-column positions.
- Each region now has a local SparseLU instance and restarts its existing
  chunk-derived `nskip`/truncation cadence. A zero cutoff retains one
  full-gravity region; empty chunks are skipped by the region bounds.
- Updated `docs/site/parameter-files.html` to state that the configured
  physical-mHz cutoff applies in both single-SEM and multi-SEM paths.
- Added focused solver API coverage for cutoff disabled/above/below/crossing,
  equality and adjacent bins, attenuation on/off, single- versus multi-SEM
  agreement on a shared mesh, and a three-chunk mixed-region case. Fixture
  inspection found two valid bins in its 2–8 mHz
  requested band, so its minimum-bin assertion now reflects the actual grid.
  The longer three-chunk fixture uses a 60-minute window and asserts its
  computed chunk count and derived `nskip > 1`; the test prints both counts.
- The new focused test passes in Debug. Its one-chunk 2–8 mHz cases use the
  same capped 0.05 mesh step in both solver paths and agree within `1e-10` for
  disabled, all-full, all-Cowling, and mixed cutoffs with attenuation both off
  and on. The longer mixed fixture reports 157 bins, 3 chunks, and derived
  `nskip=7`; the same-mesh single-/multi-SEM outputs agree within `1e-10`.
- All 9 `PreferredSolverApiTests` pass in Debug with `OMP_NUM_THREADS=1`.
  The full Release build succeeded with `-j2`, and all 63 Release CTest tests
  pass with `OMP_NUM_THREADS=1`. The root agent's original-baseline regression
  is separate and has not yet been reported. Stage 8 has not begun.
- Saved logs: `/tmp/dspecm1d-cowling-ZwV33K/stage7-debug-focused-tests.log`
  and `/tmp/dspecm1d-cowling-ZwV33K/stage7-release-ctest.log`. The Stage 7-only
  patch against the approved start snapshot is
  `/tmp/dspecm1d-cowling-ZwV33K/stage7-only.diff` (244 insertions, 42
  deletions across the implementation, focused test, guide, and this status
  entry).

## 2026-10-07 — Stage 7 parent regression and stopping point

- Original-baseline full single-SEM and multi-SEM spheroidal spectrum outputs
  remain byte-identical at l=1–4 over 2–8 mHz, two receivers, attenuation
  off/on, with one OpenMP thread. Evidence:
  `/tmp/dspecm1d-cowling-ZwV33K/stage7-baseline-regression.log`.
- The latest human instruction changes the stopping point to after Stage 7
  verification and its independent Luna-high review. Stages 8–9 have not
  begun. No checkpoint commit or push is authorized for this work.

## 2026-10-07 — Cowling Stage 8 validation harness

- Added an explicitly disabled focused PREM validation case to the existing
  solver API tests. The opt-in case compares full and Cowling receiver
  components at nearest available bins to 5, 10, 20, 40, and 80 mHz for
  degrees 1, 20, and 100; emits complex, amplitude, and phase differences;
  checks a hard cutoff against standalone results on its lower and upper
  adjacent bins; and compares selected l=1 and l=100 responses on two mesh
  resolutions. The existing short-window test helper now accepts optional
  window length, GLL order, and mesh step for this campaign.
- This is diagnostic evidence only: the case sets no universal acceptance or
  monotonicity threshold. The test is disabled in ordinary CTest runs and is
  invoked explicitly with `--gtest_also_run_disabled_tests`.
- Reproducible source is `tests/test_solver_api.cpp`; the measured CSV, exact
  input/model provenance, run command, and interpretation will be recorded in
  `validation/cowling/` after execution.
- After parent review, the diagnostic was tightened to use l=1, 20, and 100 on
  the campaign mesh, test an explicit disabled GTest entry point, assert mixed
  cutoff routing at the lower/equal/upper bins, and report switch excess beside
  the natural full/full adjacent-bin change. The resolution table now records
  receiver component and actual sampled frequencies plus each mesh's element
  and formulation DOF counts; the refinement comparison is l=1/100 at 5/80
  mHz. Results remain pending.
- Initial execution passed in 0.261 s but showed that an absolute `1e-30`
  response floor suppressed phase reporting for the SI-scale output, while a
  `0.025` mesh step was not adequate to call the 80 mHz row resolved. Changed
  the diagnostics to a `1e-12` per-component row-maximum floor and tightened
  the campaign mesh to `nq=5, maxstep=0.01` with a selected-resolution pair at
  `0.01` and `0.005`; final evidence is pending.
- A second run passed in 0.949 s with 22 MB maximum RSS. The 0.01/0.005 mesh
  pair changed 80 mHz components by 3.6–8.8% at l=1 and 11–44% at l=100,
  although the relative full/Cowling differences were similar on both meshes.
  To test further, raised the campaign mesh to `nq=5, maxstep=0.005` and the
  selected resolution pair to `0.005`/`0.0025`; the final resolution result is
  pending.
- The `0.005`/`0.0025` run passed in 1.76 s at 47 MB maximum RSS. It used 207
  and 406 elements (2490/4878 full-gravity DOFs; 1661/3253 Cowling DOFs).
  At 79.6875 mHz, the maximum selected response change across refinement was
  0.190% for l=1 and 0.785% for l=100. Parent review requested the same
  adjacent-bin cutoff check at l=20 and 100 as well as l=1; the opt-in case now
  reuses all three campaign spectra for that check. Final rerun pending.
- Final opt-in run passed in 1.885 s at 50.3 MB maximum RSS; the 9 ordinary
  `PreferredSolverApiTests.*` cases also passed after the helper changes. At
  the computed bins 4.6875, 9.375, 20.3125, 40.625, and 79.6875 mHz, maximum
  component-wise relative complex errors across l=1/20/100 were respectively
  2.47e-2, 1.12e-2, 1.13e-2, 1.01e-2, and 3.59e-4. The l=100 series is
  nonmonotone at 20–40 mHz; the endpoints show stronger high-frequency
  agreement. At the 20.3125 mHz cutoff, maximum switch excess over components
  0–2 was 1.19e-5, 0.353%, and 0.748% for l=1/20/100, alongside natural
  adjacent-bin changes up to 19.9%, 20.8%, and 33.8%. Exact mixed/full/Cowling
  routing assertions passed for bins 18.75, 20.3125, and 21.875 mHz.
- Saved reproducible inputs, raw stdout, separated CSV tables, mesh/DOF records,
  hashes, commands, interpretation, MINEOS-reference limitation, and damping
  caveats under `validation/cowling/`. The independent parent full-gravity
  baseline regression passed byte-identically with output SHA-256
  `ec838d8173387a019879565053f52fa7d7dcbe9cc449c8bae81e8370317e9c2b`.
  The focused ordinary solver API filter passed 9/9 and `git diff --check`
  passed; its output is also saved in the validation folder. Independent Stage
  8 review is now pending; Stage 9 has not started.

## 2026-10-07 — Cowling Stage 8 independent review

- Fresh separate Luna-high review returned `PASS WITH NON-BLOCKING NOTES`,
  with no blocking findings. The reviewer independently reran the campaign,
  reproduced its CSV data byte-for-byte, checked provenance hashes, and passed
  the nine focused solver API tests. Report: `validation/cowling/stage8_review.md`.
- Notes concern diagnostic code size and the bounded resolution/cutoff/MINEOS
  evidence. No automatic changes were made in response to these notes.
- The human instruction to finish this phase authorizes Stage 9 and the final
  complete-baseline Sol-high review after its own independent Luna-high gate.

## 2026-10-07 — Cowling Stage 9 performance diagnostic (in progress)

- Added a disabled opt-in PREM performance diagnostic to the existing solver
  API test fixture. It reports full/Cowling DOFs, the matched physical
  frequency, truncated matrix sizes and bandwidth/nonzero counts, repeated
  factorization and four-RHS solve medians, shared SEM mesh/both-layout build
  time, and prepared-SEM whole-spectrum medians. Symbolic analysis is included
  in factorization; ordinary CTest does not run the diagnostic. The initial
  run selected 18.75 mHz (nearest to 20 mHz), with 2490/1661 full/Cowling DOFs
  and 2430/1621 reduced matrix rows. Final measurements follow after rebuilding
  this construction-timing addition.
- Results, reproducibility commands, build/compiler provenance, and the
  simplification inspection are pending execution and documentation.
- Corrected benchmark lookup to use the frequency-band-relative index returned
  by `allIndicesSph`, and added a start-element/layout assertion against the
  production helper before measuring; the initial run is superseded.
- Removed one trailing space flagged by `git diff --check` in the new diagnostic.
- Added the actual frequency-bin count and physical endpoint output so the
  padded solver band is explicit; the preceding source hash/logs are superseded.
- Final source hash is recorded in `validation/cowling/README.md`; the final
  Release diagnostic passed with 13 bins over 6.25–81.25 mHz. At the matched
  18.75 mHz row, full/Cowling DOFs are 2490/1661 (1.4991 ratio); production
  truncation leaves 2406/1605 rows from starts 84/56. Stored nonzeros are
  35224/19213 and lower/upper bandwidths 15/15 and 10/10. Median factorization
  including symbolic analysis was 2.138/0.609 ms; four-RHS solve was
  0.341/0.143 ms; prepared-SEM 13-bin spectrum was 41.645/16.112 ms. The
  shared SEM build took 8.777 ms and constructs both layouts. These timings
  are bounded to this warm-cache single-thread local case; all RHS residual
  and finite checks passed.
- Inspected the accumulated production feature against the Stage 9
  simplification checklist. Shared full/Cowling traversal and selector-based
  layout access remain concise and no production simplification was justified;
  Stage 9 changed no production source. The Stage 8 CSV section remained
  byte-identical (89 lines); the 9 focused solver API tests and
  `git diff --check` passed. Evidence and reproducibility details are in
  `validation/cowling/README.md`, `stage9_performance.log`, and
  `stage8_after_stage9.log`. Stage 9 independent review is pending.

## 2026-10-07 — Parent checks before Stage 9 review

- Rebuilt current Release and Debug targets from this repository using cached
  dependencies and `-j2`. All 63 active Release tests passed (two opt-in
  diagnostics are disabled in ordinary CTest), plus nine Debug solver API and
  eight Debug Cowling SEM component tests. `git diff --check` passed.
- Saved logs and scope/provenance in `validation/cowling/final_checks.md`.
  Production source remains unchanged during this phase; the refreshed
  original-baseline full-gravity probe remains byte-identical.

## 2026-10-07 — Stage 9 independent review and final checkpoint

- Fresh separate Luna-high Stage 9 review returned `PASS WITH NON-BLOCKING
  NOTES`, with no blocking findings. It reran the benchmark and nine focused
  tests, verified structural counts/provenance and the preserved Stage 8 CSV,
  and assessed the cumulative production diff for simplification. Report:
  `validation/cowling/stage9_review.md`. No automatic code changes followed.
- Stage 7's prior Luna-high gate and the fresh Stage 8–9 gates have passed.
  The final fresh GPT-5.6 Sol-high reviewer will independently reassess the
  complete feature against immutable baseline
  `889452b6790d2c60afacc156b4c7512990758aee`, including the earlier code again.
  Its complete report will be `validation/cowling/final_review.md`.
- Review statistics and patches include uncommitted implementation/evidence
  files. The raw baseline diff also contains four historical benchmark files
  already committed at phase entry; they were not modified in Stages 8–9.
- After the final review report, stop for human review without code changes,
  commit, push, PR, merge, or further development.

## 2026-10-07 — Cowling remote checkpoint hygiene

- Human authorization covers packaging, a checkpoint commit and a normal push
  only. Numerical implementation and reviewed tests remain unchanged; no new
  validation or review-finding corrections are being performed.
- Added narrowly scoped ignore rules for local Cowling logs, timing output,
  derived result CSVs, temporary diffs and the three machine-specific archived
  input templates. These inputs are recreated by the retained test harness;
  intentional fixtures/reference data elsewhere remain visible to Git.
- Retained the existing baseline-regression probe source by copying it verbatim
  from the ignored local audit directory to `validation/cowling/full_regression.cpp`.
  Reports and methodology remain source-controlled; generated artifacts remain
  locally available. No tracked generated files require untracking.
- Documentation packaging notes and reproduction instructions distinguish
  local outputs from committed sources/reports. The historical benchmark PDFs
  and their accompanying files already in branch history are unchanged.
- Updated only packaging documentation: a portable Release build recipe,
  explicit local-only output/input references, and a reproduction recipe for
  the preserved baseline probe. Historical scientific results and review
  verdicts remain unchanged. No command in those validation recipes was run.

## 2026-10-08 — Approved final Cowling PR cleanup (in progress)

- Restored one chunk-local spheroidal accumulator per degree/chunk in
  `FullSpecMultiSem.h`, with one protected global block addition after both
  full/Cowling regions. Frequency indices are assigned consecutively by the
  unchanged chunk builder; empty chunks skip the final addition. Region split,
  factorization/truncation cadence, weak forms and scaling remain unchanged.
- Corrected three public cutoff accessor comments to cover single- and
  multi-SEM. No API signature or configuration behavior changed.
- Added one focused multi-degree chunk-accumulation test using the existing
  parameter helper and single-SEM reference. It covers three chunks, full-only,
  Cowling-only and a crossing cutoff, attenuation off/on, finite output and
  untouched columns outside the band; execution with one/four threads follows.
- Retained a small standalone validation driver in
  `validation/cowling/multi_sem_cleanup.cpp`, using the existing input generator
  to compare unchanged pre-cleanup/current headers for six routing/attenuation
  cases and time one fixed four-thread mixed-cutoff multi-SEM case. Generated
  outputs stay local; this adds no production or benchmark framework API.
- Removed only the four feature-added `benchmarks/unsupported/` paths from
  Git tracking with `git rm --cached`; all local bytes are preserved. Exact
  path exclusions in `.git/info/exclude` keep those local copies out of status.
  These paths did not exist in the immutable baseline. No history was rewritten.
- Validation completed: 64 active Release tests; 10 focused solver tests in
  Release/Debug with one/four threads. Serial pre-cleanup/current outputs for
  six routing/attenuation cases and the original full-gravity baseline probe
  are byte-identical; parallel differences are at most about `3.1e-16` relative.
- The isolated four-thread whole multi-SEM benchmark measured 0.465767 s before
  and 0.470966 s after (median of five warm calls). This is 1.1% slower in the
  fixed case, with no demonstrated speed-up. Locks reduce to one per nonempty
  degree/chunk. Methodology/results and reproduction driver are retained in
  `validation/cowling/cleanup_checks.md` and `multi_sem_cleanup.cpp`.
- Fresh independent Luna-high cleanup review returned `PASS`, checked all eight
  requested areas and reran the focused tests in Release/Debug with one/four
  threads. No code corrections followed. Final fresh Sol-high interaction,
  OpenMP, evidence and complete-PR scope review is next; commit/push pending.
- Fresh independent Sol-6.1-high review returned `PASS WITH NON-BLOCKING NOTES`,
  with no required corrections. It reassessed the complete immutable-baseline
  feature and cleanup, independently passed solver/component tests in both
  builds with one/four threads, and verified snapshots, numerical comparisons,
  performance wording and historical-file preservation. No code changed after
  either review. Only these gate records were updated before the authorized
  checkpoint commit/push and PR metadata update; CI must pass, then stop.
- Validation, benchmark-file untracking, fresh Luna-high/Sol-high review and
  the explicitly authorized normal commit/push are pending. Do not merge.

## 2026-11-18 — Stage 1 diffracted-wave slowness-window pilot (Luna High implementer)

- This Stage 1 was performed in a new Tunya worktree only:
  `/space/adcm2/research/projects/DSpecM1D-diffracted-wave-window`, branch
  `feature/diffracted-wave-slowness-window`, pinned at
  `b7d87b8dba1817aaa82cf30fac6c5980d228a7ba` (`Generate readable reports for
  paired YSpecpp experiments (#25)`). No commit or push was made.
- The protected canonical checkout remained unchanged before and after: HEAD
  `a778b13b22dca9709e1f96604b62925707654e5d`, branch
  `feature/radial-sensitivity-kernels-k1`, with the pre-existing modifications
  in `DSpecM1D/src/SensitivityKernels.h`,
  `DSpecM1D/src/Windowed3DFields.h`, `DSpecM1D/src/Windowed3DKernels.h`,
  `implementation_status.md`, `tests/test_sensitivity_kernels.cpp`, and
  untracked `validation/phase0/`. The pinned object remained unchanged.
- Only the Stage 1 whitelist was written on the new branch:
  `validation/diffracted_wave/{README.md,run_pilot.py,check_outputs.py,
  toroidal.params,spheroidal.params}`. No production source was modified.
  Final file SHA-256 values were: `README.md`
  `3b879e9d10e4c6364bb5391850ebcb216efa5c2c60e80dd1ce80c16e6d4466be`,
  `run_pilot.py`
  `646372e9717b5a0eded94f00e0b8d5f66fe6b5638aa2a6c67f7a576222ad1fa5`,
  `check_outputs.py`
  `c709ebb3579349cf2e92270a742ace3668c25d49cb2af568ded3db69e12a8aec`,
  `toroidal.params`
  `fb69961cf35d552e0a584e1e0e1f7fd2d63f2751a0718d1eb3f28428434e9765`, and
  `spheroidal.params`
  `6197d5ad1dbfe3805c3c980015f4289225d956d549f30e190a8b7bb982d56436`.
- The final input is an explicitly unconverged elastic PREM proxy: clean
  isotropic `data/models/prem.200.isotropic`, attenuation off, full gravity,
  equatorial source `(depth 300 km, lat 0, lon 0)`, receivers at depth 3 km
  and `(lat,lon)=(0,95),(0,130),(0,165)`, target period about 10.5 s via the
  80–113 mHz band, 10-minute/2-second record, and filter corners
  80–85–108–113 mHz. Toroidal uses `l=284..292`; spheroidal uses `l=145..153`.
  These are smoke settings only; no convergence or waveform-arrival claim is
  made.
- Offline Release configure/build succeeded with GCC 15.2.0, CMake 3.25.1,
  Ninja, and build parallelism 2. The dedicated build was
  `/space/adcm2/research/builds/dspecm1d-diffracted-wave-window/20261118-stage1`;
  the executable SHA-256 was
  `fba72734701d0cdc6c519a6569c7a73fd2a7f6a2e6090243cf36e04f95b8f786`.
  Exact pinned offline dependency sources were used. The first built-in smoke
  invocation exposed only the GCC runtime loader mismatch; rerunning with
  `LD_LIBRARY_PATH=/opt/gcc-15.2.0/lib64` passed
  `smoke_generate_seismogram` (1/1). Solver API tests were not built because a
  pinned GoogleTest source was not present on Tunya; the existing built-in
  smoke was the applicable diagnostic.
- Final pilot receipts are preserved without overwriting earlier failed
  attempts under
  `/space/adcm2/research/runs/dspecm1d-diffracted-wave-window/20261118-stage1-r4`.
  Summary SHA-256 is
  `fa937e5adaebd927112b2285c2ba42928f3487d518aa2cbba03dc1bed0419404`. Both
  cases passed finite/nonzero checks with 302 rows, 10 columns, 2-second
  spacing, and final sample 602 s (the pinned writer emits one sample before
  stopping):
  - toroidal: output SHA-256
    `28fb8a5d664313ac6827d041cd649e3d9591c992d7a4900b6f4438702052ed5d`,
    58,847 bytes, 0.14 s wall, 66,384 kB maximum RSS; vertical peaks were
    exactly zero and horizontal peaks per receiver were
    `(9.5049e-12,2.6e-15)`, `(6.5727e-12,1.64e-14)`,
    `(4.8294e-12,7.42e-14)`.
  - spheroidal: output SHA-256
    `7e20c3e6364cc3c386c11695bac110139c30f2185fe039440977d376f4a3ca63`,
    59,273 bytes, 0.95 s wall, 83,664 kB maximum RSS; component peaks per
    receiver `(6.1508e-12,4.42e-14,8.6526e-12)`,
    `(6.1142e-12,3.71e-14,5.8042e-12)`, and
    `(6.6108e-12,8.04e-14,4.7196e-12)`.
  Per-case resource receipts use `/usr/bin/busybox time -v` because Tunya has
  no `/usr/bin/time`; the receipt files record user/system time, maximum RSS,
  and exit status 0.
- `git diff --check` passed on the new worktree. Stage 1 is complete for the
  requested smoke scope and remains unconverged; no Phase 2 or production
  change was started.

## 2026-11-18 — Stage 1 solver-API evidence and pilot-safeguard closure

- The solver-API evidence gap was closed in a separate external build tree,
  preserving the existing CLI build and r4 evidence. Clean local GoogleTest
  source was verified at `v1.14.0`, commit
  `f8d7d77c06936315286eb55f8de22cd23c188571`, with no worktree changes. Its
  streamed archive SHA-256 was
  `e5f89b43d9cc0916723695c8503430fdd5cc25e8bd50881dbd242cef468a0533`.
  Tunya received it only at
  `/space/adcm2/research/dependencies/dspecm1d-diffracted-wave-window/20261118-stage1-googletest`;
  extracted `CMakeLists.txt` SHA-256 is
  `e9934d724b827530afbc2227a56e7459e608a49947ff0c4c802345eaa3ac7efb`.
- The dedicated solver-API build was
  `/space/adcm2/research/builds/dspecm1d-diffracted-wave-window/20261118-stage1-solver-api`.
  It used the pinned DSpecM1D dependency sources, GCC 15.2.0, Release mode,
  and at most two build workers. The initial CMake GoogleTest discovery hit
  the known GCC runtime loader mismatch; rebuilding the same target with
  `LD_LIBRARY_PATH=/opt/gcc-15.2.0/lib64` completed successfully. Solver test
  executable SHA-256:
  `41494979d490dbd7e41fbae70dfd561be2e037302b304e82675fa29d55909f88`.
- The ordinary focused command used one worker for OpenMP and BLAS, plus
  `timeout --kill-after=5s 600s` outside BusyBox timing. It ran
  `--gtest_filter=PreferredSolverApiTests.*`: 15 tests passed; the only two
  omitted tests were the source-marked disabled
  `DISABLED_CowlingValidationCampaignOnPrem` and
  `DISABLED_CowlingPerformanceOnPrem`. Receipts are preserved under
  `/space/adcm2/research/runs/dspecm1d-diffracted-wave-window/20261118-stage1-solver-api`.
  Focused-test log SHA-256 is
  `1da2b14d4922c5c83e59221819800ea2e0979bf5cc9b2fe93735aada0e513fbb`;
  test resource receipt SHA-256 is
  `6236ebbbf3d72a79c37a4c03b5732fc929d450781f1b41bda6a65c5ff90acac0`;
  maximum RSS was 54,016 kB and exit status was 0.
- The pilot driver was hardened without production changes: it now requires a
  nonexistent run directory (`exist_ok=False`), pins
  `OMP_NUM_THREADS=1`, `OMP_THREAD_LIMIT=1`, `OPENBLAS_NUM_THREADS=1`, and
  `MKL_NUM_THREADS=1`, and wraps each solver invocation in GNU
  `timeout --kill-after=5s 300s` outside `/usr/bin/busybox time -v`. Launcher
  stdout was redirected to a sibling log so it could not pre-create the run
  directory. Final driver SHA-256 is
  `cf994e87411e5ed4bd9692049da1e27babee8da7858b713e73fc1b2e90b20290` and
  final README SHA-256 is
  `89a4644cf6003acb777552691e408983e4e9bfb7b63b89cba982ed9f4fbfc6d1`.
- Both proxies were rerun once with the safeguarded driver in fresh
  `/space/adcm2/research/runs/dspecm1d-diffracted-wave-window/20261118-stage1-r5`.
  Summary SHA-256 is
  `7b3b0ae94e76e7a3ac6c068d06eaec3c68e3dbbd5d5af2ac8de2dc185b0c0ace`.
  Both checks passed with 302 rows, 10 columns, and final sample 602 s.
  Toroidal output SHA-256 remains
  `28fb8a5d664313ac6827d041cd649e3d9591c992d7a4900b6f4438702052ed5d`,
  with 67,920 kB maximum RSS; spheroidal output SHA-256 remains
  `7e20c3e6364cc3c386c11695bac110139c30f2185fe039440977d376f4a3ca63`,
  with 83,648 kB maximum RSS. Both resource receipts report exit status 0.
- No production source, canonical checkout, commit, push, or next stage was
  touched. The final new-worktree diff remains limited to
  `implementation_status.md` and `validation/diffracted_wave/`; the canonical
  checkout and pinned source commit remain unchanged. Stage 1 is frozen for
  independent review.

## 2026-11-18 — Stage 1 timeout advisory closure

- Applied the sole requested post-review code delta in
  `validation/diffracted_wave/run_pilot.py`: Python `subprocess.run` grace is
  now `timeout=315`, while the GNU wrapper remains
  `/usr/bin/timeout --kill-after=5s 300s` (TERM at 300 s, KILL after 5 s).
  The resource command is accurately recorded as BusyBox `time -v` wrapping
  that GNU timeout command.
- Ran both unchanged proxies once in fresh
  `/space/adcm2/research/runs/dspecm1d-diffracted-wave-window/20261118-stage1-r6`;
  r5 was preserved. Both checks passed. The r6 toroidal and spheroidal raw
  outputs are byte-identical to r5 (`cmp` passed for both), with unchanged
  SHA-256 values `28fb8a5d664313ac6827d041cd649e3d9591c992d7a4900b6f4438702052ed5d`
  and `7e20c3e6364cc3c386c11695bac110139c30f2185fe039440977d376f4a3ca63`.
  R6 summary/launcher log SHA-256 is
  `0e5466a5f8a5bb5c214ad48f3ada5e2752f4b5a43ec68790fc73205fe9cc930e`.
- Final driver SHA-256 is
  `e766d21ddddb3ebfbadb0ddee37dbfd5788126dac3dae943e3f911d48433521a`.
  R6 resource receipts report exit status 0, 69,296 kB maximum RSS for
  toroidal and 83,696 kB for spheroidal. No other source, input, dependency,
  build, or output was changed; no commit or push was made.

### Independent Stage 1 closeout

The independent reviewer returned PASS on the Stage 1 workflow and PASS on
the timeout correction. The final review verified the six scoped files,
the pinned base and protected canonical branch, the CLI smoke (1/1), the
enabled solver API tests (15/15), and the final r6 proxy evidence. The
reviewed r6 summary SHA-256 is
`0e5466a5f8a5bb5c214ad48f3ada5e2752f4b5a43ec68790fc73205fe9cc930e`.
Both r6 waveforms match the r5 bytes. The pilot qualifies execution and
finite, mode-appropriate output only. The next stage requires user approval
and establishes numerical reference convergence and measurement accuracy.

### 2026-11-19 — Opt-in post-undamping SH solve pruning

The 60 s and 70 s Stage 2 proxy used 495 physical frequency bins and degrees
1–420. Offline removal of pre-undamping pairs below the 26 s/degree cutoff
retained 80,543 of 207,900 pairs. Independent review reproduced the counts and
receiver scoring. Across both time-end ramps and FFT lengths, the zero-guard
mask gave worst errors of 0.001495 s/degree and 0.0003933 per degree against the
unwindowed reference, within the 0.003 and 0.004 limits. Both the primary groups
and the 145–165 degree stress group passed.

The implementation on `feature/sh-post-undamping-guard` adds
`generate_seismogram <parameter-file> --experimental-sh-window` and an explicit
`SparseFSpec::postUndampingSH` API. The command accepts toroidal mode 2 and uses
the validated frequency sampling setting qex=4. It solves only the retained
high-frequency suffix of each degree/chunk, uses the original radial truncation
and symbolic-analysis cadence, removes exponential damping per degree on the
half-step synthesis grid, and then applies the four-corner cosine taper with
corners 0, 0.01, 23.5, 26 s/degree. A cosine ramp covers the final 2.5% of the
time record. The experimental `_sh_window_t.out` has SI Z/N/E samples and is
written before an additional analysis bandpass. The ordinary command retains
its existing filter and output contract.

Gate before timing: verify cutoff equality, empty-window and wrong-mode cases,
ordinary output parity, actual solve counts, and agreement of native traces
with the independent offline mask for both 60–70 s proxy models. This evidence
applies to the archived long-period SH grid and proxy source/receiver setup.

Independent runtime gate PASS for the scoped experiment. The native executable
solved exactly 80,543 pairs in each model. The single-worker test traces matched
the independently masked archived north component within 4e-15 peak-relative.
Re-picked PREM and model-anomaly maxima were 0.001257 s/degree and 0.0001605
per degree across the tested periods and distance groups. The ordinary tool
produced identical output bytes on the regression fixture. Exact equality at
degree 292 and 31.25 mHz, empty-window output, and spheroidal rejection passed.
The Tunya native receipts were archived for review.

Existing regression suites also passed against the exact source snapshot.
The solver API suite passed all 15 active tests and the output writer suite
passed all 10. Two existing opt-in Cowling diagnostics were disabled.

The scoped performance gate passed with a persistent picker on Tunya. With
eight OpenMP workers, one BLAS worker, and three measured repeats after a
warm-up, full workflow medians including reading, filtering and both period
picks were 3.410 to 3.041 s on PREM and 3.420 to 3.056 s on basal-Vs.
Every measured repeat met the scientific gates. Median radial time fell
from about 3.19 to 2.07 s, and the post-undamping synthesis took about
0.75 s. The benchmark with a fresh Python picker for each trace gave mixed
end-to-end times due to picker-process startup variability. Both benchmark
designs and run receipts are preserved under the Tunya run directory.

The accepted scope is the long-period SH proxy at 60 and 70 s. The
experimental path and its receipts are on the separate feature worktree.

### 2026-11-19 — Collaborator stage 1 and Wigner high-degree support

The GSHTrans pin moves to 991dd1eaba2d0c47e4200356a0329c51dfe5fb5e.
The other dependency pins stay fixed. A focused regression exercises the
production Wigner order shape through degree 3600 and compares double against
long double at three distances with both signs of the supported orders.

Sol planned the scope. GPT-6-Luna High implemented it. An independent review
and the complete available regression suite support the engineering checkpoint.
The suite passed 85 active tests; two existing opt-in diagnostics are disabled.
The old double pin rejected degree 4096 at its 1827 guard. On the candidate,
21,150 angular values through degree 4096 agreed with the old long-double
calculation to 2.7e-12 peak-relative. A degree-3600 waveform smoke at 248–250 mHz,
19.5 km source depth and 141 receivers produced finite traces. Full-gravity and
Cowling proxy smokes were also exercised.

The collaborator’s exact decks, input settings, reference traces and adaptive picker
were unavailable. The full-band collaborator reference and inversion-scale
accuracy and timing remain pending. Subsequent engineering work must label its
stand-in cases and keep the Cowling approximation separate from taper error.

### 2026-11-19 — Collaborator stage 2 and component-specific taper

Sol replanned following stage 1. GPT-6-Luna High added opt-in mixed-family
transverse Sdiff for mode 4 and spheroidal vertical Pdiff for mode 3.
`--experimental-sdiff-window` and `--experimental-pdiff-window` take four
explicit slowness corners and an optional `--pre-solve-guard-s-per-deg`.
The default in each new route retains the full requested solve grid. The
explicit guard expands both outer corners for choosing raw solve pairs and
preserves the original chunks, radial truncation and Cowling boundaries.
The cosine taper acts after inverse transforming, removing exponential
damping, and transforming back to real frequencies. The mixed route
sums both families in the physical great-circle T projection.

An independent per-degree reconstruction on oblique paths agreed with the
native output to 4.6e-13 peak-relative over 16 tested combinations. Tests
covered full gravity, all-Cowling and a mixed cutoff, with no guard, a
two-sided guard and an empty guard. Both exact endpoints and four-worker
contractions passed. Legacy and original SH output were byte-identical on
the compatibility fixture. The full suite passed 85 active tests and narrow
Cowling smoke runs returned finite nonzero T and Z at degrees 3596–3600 for
141 receivers. Independent code review supports the scoped engineering gate.

The comparison with collaborator decks and picker remains pending. These
proxy checks establish numerical agreement with the chosen masked
calculation. They do not validate any guard, taper or Cowling cutoff for
the collaborator's measured phase slowness or decay.

### 2026-11-19 — Collaborator stage 3, raw output and in-memory trials

Sol replanned after the stage 2 review. Two GPT-6-Luna High agents added
independent raw output and in-memory model changes. A separate GPT-6-Luna
High agent reviewed the final engineering scope.

`generate_seismogram <parameter-file> --raw-time --raw-spectrum` writes
both outputs from one standard solve. Either flag can be used individually.
The unfiltered time output applies exponential undamping once and writes
samples in seconds up to the requested record duration. The complex
output preserves the damped one-sided grid through Nyquist, including
unsolved zeros, with an axis in Hz and amplitudes in the selected SI
observable times seconds. Output columns are grouped by receiver in Z, N,
E order. The console reports epsilon in inverse seconds.

`EarthModels::IsotropicModel` takes radius, density, P and S velocities,
both Q arrays and the reference period as in-memory input. It applies the
existing layered spline conversion with independent values on the two sides
of repeated radii. The solver overloads accept one candidate at a time for
standard spectra or the phase-window route, using the prepared source,
receivers, frequency grid and Cowling cutoff. Each trial builds its own
model-dependent SEM state and leaves the prepared context usable after a
rejected model.

All 87 active tests passed, with two existing Cowling diagnostics disabled.
The exported raw spectrum agreed with the independently scaled stage 2
solver spectra to 1.6e-16 peak-relative error or better. Its NumPy inverse
transform reproduced displacement, velocity and acceleration time traces
to 3.1e-15 peak-relative error or better. A fractional sample interval
aligned exactly with the final record time retained the endpoint. The independent PREM and
perturbed-deck experiment found exact agreement for raw spectra and
windowed T and Z. The perturbation changed Vs, an internal interface and
the attenuation reference period. Four invalid trial models were rejected
on each route. Stage 2 regression retained legacy filtered and original SH
file identity and passed the independent reconstruction and solve-count
checks for both mixed-family T and spheroidal Z.

The engineering gate used pinned dependencies in an isolated local build.
The collaborator's matching model decks, traces and adaptive picker are
needed for the 0.001–0.002 s/deg acceptance gate and a fair end-to-end
runtime comparison. This stage provides the C++ trial API and retains
the existing spline interpolation. Model-specific SEM and factorization
work remains in each trial.

### 2026-11-19 — Deep-source benchmark follow-up, Stage 1

The spheroidal phase-window path in `FullSpecMultiSem.h` now evaluates the
receiver contribution once before extracting Z or T. The existing signs,
component indices and output factor are retained. This is the only code
hunk in this stage.

Sol specified the implementation and review gates. One GPT-6-Luna High agent
implemented the hunk and ran the frozen controls. A separate GPT-6-Luna High
agent independently reviewed the source, replayed the reduced controls and
audited the full native waveforms, solve counts and benchmark receipts.

On Tunya, the matched 605 km, 65 mHz, degree-375 PREM controls all matched
byte for byte on the 141-receiver quarter-second waveform grid. Pdiff
retained 199125 pairs with all solves and 127308 with the zero guard.
Mixed Sdiff retained 398250 and 254616 pairs respectively. The independent
reduced cases cover one and four workers and a nonzero spheroidal T
contribution. Serial original SH and raw controls matched exactly.
Both source snapshots passed all 87 active tests.

Two alternating four-worker all-solves Pdiff timing repeats averaged
401.215 s for the baseline and 97.743 s for the candidate. User CPU averaged
1531.61 and 322.665 s. Synthesis averaged about 20.5 s in both builds.
The supplemental four-worker raw comparison differed at floating-point
scale, with relative L2 error 1.787e-15 versus 1.634e-15 between repeats of
the unchanged baseline, exact axes and no peak-index shifts. Independent
review accepted raw compatibility and maintained the exact window gate.

The timing and waveform receipts were archived for review.
The matched full-band Cowling test and picker/taper investigations have
separate approval gates.
