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
