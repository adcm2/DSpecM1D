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
