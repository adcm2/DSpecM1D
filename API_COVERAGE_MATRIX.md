# DSpecM1D API Coverage Matrix

Last updated: 2026-04-23 (Europe/London)

This file records how the current release-facing `DSpecM1D` surface is covered after the API hardening work.

## Coverage Labels

- `Direct modular tests`: Covered by the small unit/component suites under `tests/`.
- `Indirect modular tests`: Not called directly by name in a test, but exercised through a directly tested preferred wrapper.
- `Paper examples`: Covered mainly by `ex1` through `ex7` rather than by a narrow modular test.
- `Smoke`: Covered mainly by `t1` or the smoke build checks.
- `Legacy`: Kept for compatibility and not a primary hardening target for this pass.

## Preferred Release-Facing API

| Surface | Coverage | Notes |
| --- | --- | --- |
| `InputParametersNew` construction and path resolution | Direct modular tests | Covered by `InputParametersNewTests.ResolvesAbsoluteModelPathAndBuildsContext`. |
| `InputParametersNew::normFactor()` | Direct modular tests | Covered for displacement, velocity, and acceleration output modes. |
| `InputParametersNew::setNq()` / `setNskip()` | Direct modular tests | Covered by the clamping test. |
| `InputParametersNew` bundled context accessors (`inputParameters`, `cmt`, `srInfo`, `earthModel`, `freqFull`, `tref`) | Direct modular tests | Covered by the `InputParametersNew` and `PreferredSolverApi` suites. |
| `SPARSESPEC::SparseFSpec::spectra(InputParametersNew &)` | Direct modular tests | Covered by `PreferredSolverApiTests.PreferredSolverOverloadsReturnStableShapes`. |
| `SPARSESPEC::SparseFSpec::spectra(InputParametersNew &, Full1D::SEM &)` | Direct modular tests | Covered directly in the preferred solver API suite. |
| `SPARSESPEC::SparseFSpec::spectra(Full1D::SEM &, InputParametersNew &)` | Direct modular tests | Covered directly and checked for consistency with the reused-SEM overload above. |
| `SPARSESPEC::SparseFSpec::spectra(const SpectraRunContext &, Full1D::SEM &)` | Indirect modular tests | The preferred SEM-reuse overloads collapse into this path, so it is exercised indirectly. |
| `SPARSESPEC::SpectraRunContext` constructor and accessors | Direct modular tests | Covered by `PreferredSolverApiTests.SpectraRunContextExposesWorkflowObjects`. |
| `Full1D::SEM(const InputParametersNew &)` | Direct modular tests | Covered by `PreferredSolverApiTests.SemConstructorFromInputParametersNewUsesSettings`. |
| `DSpecM::FilterOptions` defaults | Direct modular tests | Covered by `SignalFilteringTests.FilterOptionsExposeExpectedDefaults`. |
| `DSpecM::applyFilter(const Eigen::MatrixXcd &, ...)` | Direct modular tests | Covered for valid input, zero input, invalid width, multiple passes, and `enforceRealSignal`. |
| `DSpecM::applyFilter(const Eigen::MatrixXd &, ...)` | Direct modular tests | Covered for valid input, invalid width, dimension preservation, and finite output. |
| `DSpecM::loadYSpecTimeSeries(...)` | Direct modular tests | Covered for truncation, padding, and zero-column requests. |
| `DSpecM::loadSpecnmTimeSeries(...)` | Direct modular tests | Covered for parsing, padding, and zero-column requests. |
| `DSpecM::loadMineosTimeSeries(...)` | Direct modular tests | Covered for scaling, truncation to shortest component, and zero-column requests. |
| `DSpecM::loadReferenceTimeSeries(...)` | Direct modular tests | Covered directly by the bundled reference-series test. |
| `DSpecM::writeFrequencyComparison(...)` low-level two-way overload | Direct modular tests | Covered for truncation, layout, and invalid output path handling. |
| `DSpecM::writeFrequencyComparison(...)` low-level three-way overload | Direct modular tests | Covered for three-way layout and truncation behavior. |
| `DSpecM::writeFrequencyComparison(...)` `InputParametersNew` overloads | Direct modular tests | Covered for both two-way and three-way convenience paths. |
| `DSpecM::writeTimeComparison(...)` low-level overload | Direct modular tests | Covered for output-window truncation, shortest-input truncation, and invalid output path handling. |
| `DSpecM::writeTimeComparison(...)` `InputParametersNew` overload | Direct modular tests | Covered directly by the convenience writer test. |

## Helper Utilities

| Surface | Coverage | Notes |
| --- | --- | --- |
| `InputParser` helper functions (`get_next_value_line`, `read_required_scalar`, `read_required_lat_lon`) | Direct modular tests | Covered directly by `InputParserTests`. |
| `InputParameters` legacy parser | Direct modular tests | Covered for valid parsing and invalid latitude rejection. |
| `SPARSESPEC::resolveModeFlags(...)` | Direct modular tests | Covered by `SpecHelpersTests.ResolveModeFlagsHonorsModeTypeAndAngularRange`. |
| `SPARSESPEC::outputFactor(...)` | Direct modular tests | Covered for displacement, velocity, and acceleration modes. |
| `SPARSESPEC::SpecConstants` | Direct modular tests | Covered for the reference-period derived constants. |
| `SPARSESPEC::factorizeOrCompute(...)` | Direct modular tests | Covered with a fake solver to validate call selection. |

## SEM Low-Level API

| Surface | Coverage | Notes |
| --- | --- | --- |
| `SEM` legacy constructor `(model, maxstep, NQ, lmax)` | Direct modular tests | Used by the SEM component suite and the preferred SEM-construction comparison test. |
| `SEM::ltgS`, `ltgR` | Direct modular tests | Covered for monotonic and contiguous ordering invariants. |
| `SEM::ltgT` | Direct modular tests | Covered indirectly through toroidal receiver-base shape checks. |
| `SEM::receiverElements(...)` / `sourceElement(...)` | Direct modular tests | Covered for in-mesh lookup behavior. |
| `SEM::mesh()` / `meshModel()` | Direct modular tests | `mesh()` is covered directly in SEM tests and `meshModel()` now has a small direct density/gravity contract check in the SEM component suite. |
| `SEM::hS`, `pS`, `hTk`, `pTk`, `hR`, `pR` | Direct modular tests | Covered for square shape compatibility and structural non-emptiness. |
| `SEM::rvZR`, `rvRedZR`, `rvBaseFull`, `rvBaseFullT` | Direct modular tests | Covered for shape and non-zero support. |
| `SEM::rvFull`, `rvFullT`, `rvBaseZ`, `rvValZ`, `rvBaseTheta`, `rvValTheta`, `rvBaseThetaT`, `rvValThetaT`, `rvValPhi`, `rvBasePhiT`, `rvValPhiT`, `rvThetaT`, `rvPhiT` | Paper examples / indirect solver coverage | These remain lower-level helpers used inside solver/example workflows and are not individually targeted by narrow tests yet. |
| `SEM::calculateForce*` families | Paper examples / indirect solver coverage | These are exercised through the solver and paper examples, but not individually unit-tested. |

## Legacy Solver Surface

| Surface | Coverage | Notes |
| --- | --- | --- |
| `SparseFSpec::spectra(FreqFull &, Full1D::SEM &, model1d &, EarthquakeCMT &, InputParameters &, int)` | Direct modular tests | Now covered directly by the small legacy single-SEM regression test in `test_solver_api.cpp`. |
| `SparseFSpec::spectra(FreqFull &, model1d &, EarthquakeCMT &, InputParameters &, int, SRInfo &, double)` | Direct modular tests | Now covered directly by the small legacy multi-SEM regression test in `test_solver_api.cpp`. |

## Example And Smoke Coverage

| Workflow | Coverage role | Notes |
| --- | --- | --- |
| `t1` | Smoke | Protects the self-contained end-to-end synthesis path without external comparison data. |
| `ex1` | Paper example + modern public-path example | Exercises the clearest preferred release-facing workflow. |
| `ex2`, `ex3`, `ex7` | Paper examples | Protect legacy/paper comparison workflows, including STF-specific or multi-reference output paths. |
| `ex4` | Paper example | Protects the record-section style workflow. |
| `ex5` | Paper example | Protects the convergence-study workflow. |
| `ex6` | Paper example | Protects the low-level SEM/matrix workflow discussed in the paper. |

## Residual Risks And Intentional Gaps

- The remaining uncovered low-level SEM helper methods are intentionally not unit-tested individually because they are tightly coupled to the solver internals and are already exercised through the solver and paper examples.
- The legacy `SparseFSpec` overloads now have small direct regression coverage, but the paper examples remain the authoritative scientific check for those workflows.
- `meshModel()` is publicly reachable through `SEM` and now has a small direct contract check, but its richer behavior is still mainly relevant to specialist low-level workflows such as `ex6`.
- The paper examples remain the authoritative protection for scientific workflow fidelity; the modular tests are designed to catch smaller regressions closer to their source, not to replace the paper-reproduction checks.

## Docs Integration

- The documentation site now exposes the testing material as multiple linked pages rather than one long page:
  - testing overview
  - smoke and integration checks
  - parser and workflow-context tests
  - reference I/O and output-writer tests
  - filtering, solver API, and SEM/component tests
  - testing coverage summary
