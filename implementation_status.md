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
