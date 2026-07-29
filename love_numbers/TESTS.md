# Love-number test inventory

This inventory is based on `ctest --show-only=json-v1` from a configured
paper-validation build. It contains 58 normal tests and 11
paper-validation tests: 69 Love-number tests in total. Core DSpecM1D tests
are excluded. Classifications and descriptions were checked against the
registered commands and the corresponding test bodies.

## Radial state, material model, and gravity

- **normal** `RadialStateTests.PolynomialOrderUsesOneMoreGllNode` — `tests/test_radial_state.cpp` — Checks that polynomial degree `p` creates exactly `p+1` GLL nodes per element.
- **normal** `RadialStateTests.MeshEndsAtSurfaceWithoutExteriorBall` — `tests/test_radial_state.cpp` — Checks that the mesh contains no exterior ball and ends at the model outer radius.
- **normal** `RadialStateTests.SampledModelValuesAreFinite` — `tests/test_radial_state.cpp` — Samples density, TI moduli, gravity, and the density derivative at every element node and requires finite values.
- **normal** `RadialStateTests.HomogeneousSphereGravityIsAnalyticAtEveryNode` — `tests/test_radial_state.cpp` — Compares every nodal gravity value with the homogeneous-sphere analytic solution, including the finite zero at the centre.
- **normal** `RadialStateTests.ControlledSolidFluidSolidGravityIsAnalytic` — `tests/test_radial_state.cpp` — Compares nodal gravity with the piecewise analytic solid--fluid--solid density integral.
- **normal** `RadialStateTests.MinimumPolynomialOrderExactlyIntegratesCubicDensityGravity` — `tests/test_radial_state.cpp` — Checks roundoff-level analytic surface gravity for cubic density at the minimum polynomial order `p=3`.
- **normal** `RadialStateTests.ControlledCentralFluidGravityIsAnalytic` — `tests/test_radial_state.cpp` — Checks the analytic gravity of a constant-density central fluid and solid exterior.
- **normal** `RadialStateTests.GravityIsContinuousAcrossEveryBoundary` — `tests/test_radial_state.cpp` — Requires identical gravity values on both copies of every element and material boundary.
- **normal** `RadialStateTests.DensityDerivativeUsesElementLayerSpline` — `tests/test_radial_state.cpp` — Compares each sampled derivative with the derivative of that element's density spline.
- **normal** `RadialStateTests.MaterialInterfaceUsesOneSidedDensityDerivatives` — `tests/test_radial_state.cpp` — Checks that duplicated-radius interface nodes use the lower and upper layer derivatives independently.
- **normal** `RadialStateTests.RejectsSurfaceFluid` — `tests/test_radial_state.cpp` — Requires a clear exception when the outermost layer is fluid.
- **normal** `RadialStateTests.RejectsInvalidConfiguration` — `tests/test_radial_state.cpp` — Rejects negative maximum degree and non-positive radial step.

## DOF mapping

- **normal** `LoveDofMapTests.DegreeZeroInterleavesEveryNodeAndMatchesRadialSemMap` — `tests/test_love_dof_map.cpp` — Checks continuous interleaved `U,P` indexing and exact agreement with the existing radial SEM map.
- **normal** `LoveDofMapTests.DegreeZeroRetainsCentreFields` — `tests/test_love_dof_map.cpp` — Verifies independent `U(0)` and `P(0)` coefficients.
- **normal** `LoveDofMapTests.AllSolidDegreeTwoInterleavesEveryField` — `tests/test_love_dof_map.cpp` — Checks all-solid `U,V,P` interleaving, dimension, and exact agreement with the existing spheroidal SEM map.
- **normal** `LoveDofMapTests.SolidCentreRetainsThreeIndependentFields` — `tests/test_love_dof_map.cpp` — Verifies independent `U(0)`, `V(0)`, and `P(0)` coefficients.
- **normal** `LoveDofMapTests.PotentialIsSharedAcrossEveryElementBoundary` — `tests/test_love_dof_map.cpp` — Checks continuity of the single potential coefficient across solid, fluid, and material boundaries.
- **normal** `LoveDofMapTests.DisplacementsAreSharedAcrossSolidBoundaries` — `tests/test_love_dof_map.cpp` — Checks that adjacent solid elements share their `U` and `V` boundary coefficients.
- **normal** `LoveDofMapTests.DisplacementsAreAbsentFromFluidElements` — `tests/test_love_dof_map.cpp` — Requires `U` and `V` to return `nullopt` at every fluid node.
- **normal** `LoveDofMapTests.InterfacesKeepDisplacementsOnlyOnSolidSide` — `tests/test_love_dof_map.cpp` — Exercises both interface orientations and checks that displacement exists only on the solid side while `P` is shared.
- **normal** `LoveDofMapTests.InternalFluidDimensionsMatchTopologyFormulas` — `tests/test_love_dof_map.cpp` — Compares explicit map sizes with the radial-node and separate-solid-run dimension formulas.
- **normal** `LoveDofMapTests.DegreeOneOnlyOmitsSurfacePotential` — `tests/test_love_dof_map.cpp` — Compares degree one and degree two and requires the sole missing coefficient to be surface `P`.
- **normal** `LoveDofMapTests.DegreeOneRetainsSurfaceDisplacements` — `tests/test_love_dof_map.cpp` — Verifies that degree-one surface `U` and `V` remain present.
- **normal** `LoveDofMapTests.RejectsNegativeDegree` — `tests/test_love_dof_map.cpp` — Requires a negative degree to throw `std::invalid_argument`.

## All-solid operator

- **normal** `StaticOperatorTests.MatchesIsotropicReference` — `tests/test_static_operator.cpp` — Compares degree 2 and 10 isotropic matrices with `Full1D::SEM::hS(l)`, including dimensions and roundoff symmetry.
- **normal** `StaticOperatorTests.MatchesTransverselyIsotropicReference` — `tests/test_static_operator.cpp` — Repeats the degree 2 and 10 `hS(l)` matrix comparison for a genuinely TI constant-density model.
- **normal** `StaticOperatorTests.RetainsCentreFieldsAndMatchesSemOrdering` — `tests/test_static_operator.cpp` — Checks the three centre coefficients and exact all-solid SEM ordering.
- **normal** `StaticOperatorTests.RejectsNegativeDegree` — `tests/test_static_operator.cpp` — Requires the private static operator to reject a negative degree.

## Fluid operator and interfaces

- **normal** `StaticFluidOperatorTests.InterfaceSignsUseFluidSideDensity` — `tests/test_static_fluids.cpp` — Isolates both solid--fluid interface orientations and verifies sign, fluid-side density, shared `P`, and solid-side `U`.
- **normal** `StaticFluidOperatorTests.NonconstantFluidBlockUsesPositiveDensityGradient` — `tests/test_static_fluids.cpp` — Compares the fluid potential block with an independent quadrature assembly containing `+r^2 rho'/g`.
- **normal** `StaticFluidOperatorTests.ConstantDensityCentralFluidIsFinite` — `tests/test_static_fluids.cpp` — Checks zero stratification, a finite central-fluid matrix, and safe handling of `r=g=0`.
- **normal** `StaticFluidOperatorTests.InternalFluidMatricesAreCompleteAndSymmetric` — `tests/test_static_fluids.cpp` — Requires internal-fluid and no-ocean PREM matrices to have the map dimension, roundoff symmetry, and no empty row or column.

## Forcing, solve, residuals, and reciprocity

- **normal** `StaticForcingTests.SurfaceLoadsHaveExpectedEntriesAndSigns` — `tests/test_static_solve.cpp` — Checks the displacement and gravitational surface load entries after the single `KX=-L` negation.
- **normal** `StaticForcingTests.AllSolidTideMatchesIndependentAssembly` — `tests/test_static_solve.cpp` — Compares the all-solid tidal column with a test-local GLL assembly at degrees 2 and 10.
- **normal** `StaticForcingTests.InternalFluidTideMatchesBothInterfaceOrientations` — `tests/test_static_solve.cpp` — Independently assembles solid, fluid-gradient, and both oriented interface tidal terms.
- **normal** `StaticSolveTests.OneFactorizationSolvesAllColumnsAndExtractsResults` — `tests/test_static_solve.cpp` — Checks three-column residuals, combined-load linearity, dimensional extraction, reciprocity, and finite results for isotropic, TI, fluid, and PREM models.
- **normal** `StaticSolveTests.RejectsNegativeDegree` — `tests/test_static_solve.cpp` — Requires the private solve function to reject a negative degree.

## Degree zero

- **normal** `StaticRadialOperatorTests.MatchesRadialSemForSupportedModels` — `tests/test_static_radial.cpp` — Compares degree-zero matrices with `Full1D::SEM::hR()` for constant-density isotropic and TI models, with symmetry and ordering checks.
- **normal** `StaticRadialOperatorTests.RetainsCentreFields` — `tests/test_static_radial.cpp` — Verifies independent radial `U(0)` and `P(0)` coefficients.
- **normal** `StaticRadialForcingTests.SurfaceLoadsAndZeroTideAreExact` — `tests/test_static_radial.cpp` — Checks the two surface loads and requires an identically zero degree-zero tidal column.
- **normal** `StaticRadialSolveTests.SolvesSupportedModels` — `tests/test_static_radial.cpp` — Checks degree-zero residuals, reciprocity, finite model coverage, and exact zero `h_t,k_t`.

## Degree one

- **normal** `StaticDegreeOneOperatorTests.AllSolidMatricesMatchConstrainedSem` — `tests/test_static_degree_one.cpp` — Compares isotropic and TI matrices with the principal submatrix of `hS(1)` after removing surface `P`, and evaluates the full translation null vector.
- **normal** `StaticDegreeOneOperatorTests.RetainsCentreAndSurfaceFields` — `tests/test_static_degree_one.cpp` — Checks independent centre `U,V,P`, retained surface `U,V`, and absent surface `P`.
- **normal** `StaticDegreeOneOperatorTests.InternalFluidMatricesHaveExpectedSizeAndSymmetry` — `tests/test_static_degree_one.cpp` — Checks map dimension and roundoff symmetry for an internal fluid and no-ocean PREM.
- **normal** `StaticDegreeOneForcingTests.SurfaceLoadsAndIndependentTidesAreCorrect` — `tests/test_static_degree_one.cpp` — Verifies the displacement load, exactly zero gravitational column, and independent all-solid and fluid tidal assemblies.
- **normal** `StaticDegreeOneSolveTests.SolvesSupportedModels` — `tests/test_static_degree_one.cpp` — Checks constrained residuals, finite remaining responses, and exact zero potential-related outputs.

## Public API

- **normal** `LoveNumbersTests.LoadSumsCombineGeneralizedResponses` — `tests/test_love_numbers.cpp` — Checks that `h_load()` and `k_load()` return the sums of their generalized components.
- **normal** `LoveNumbersTests.PublicHeaderCalculatesOrderedDegreeRange` — `tests/test_love_numbers.cpp` — Uses only the public header and checks an ordered result for every degree from zero through the configured maximum.
- **normal** `LoveNumbersTests.MaximumDegreeZeroIsValid` — `tests/test_love_numbers.cpp` — Checks that a degree-zero-only public calculation succeeds and returns one result.
- **normal** `LoveNumbersTests.SupportedModelsReturnFiniteValues` — `tests/test_love_numbers.cpp` — Runs the public API on isotropic, TI, internal-fluid, and no-ocean PREM models and checks special-degree zeros and finite values.
- **normal** `LoveNumbersTests.RejectsInvalidConfigurationAndSurfaceFluid` — `tests/test_love_numbers.cpp` — Checks public rejection of invalid configuration values and a surface fluid.
- **normal** `LoveNumbersTests.EnforcesMinimumPolynomialOrder` — `tests/test_love_numbers.cpp` — Requires the public calculation to reject polynomial orders below three and accept `p=3`.
- **normal** `LoveNumbersPublicPrivateTests.PublicResultsMatchPrivateSolves` — `tests/test_public_private.cpp` — Compares every public component directly with private solves at selected degrees including 0, 1, 2, and 10.

## CLI

- **normal** `LoveNumbersCliTests.ExplicitOptionsWriteFileAndMatchPublicResults` — `tests/test_cli.cpp` — Runs explicit optional arguments, parses the file metadata and nine columns, and compares every value with `calculate`.
- **normal** `LoveNumbersCliTests.DefaultsWriteStandardOutput` — `tests/test_cli.cpp` — Runs the default polynomial order and radial step through standard output and checks ordered rows.
- **normal** `LoveNumbersCliTests.InvalidArgumentsAndMissingModelFail` — `tests/test_cli.cpp` — Requires non-zero status and clear diagnostics for malformed arguments and a missing model.
- **normal** `LoveNumbersCliTests.RejectsUnderintegratedPolynomialOrders` — `tests/test_cli.cpp` — Requires non-zero status and the `p >= 3` diagnostic for `p=1` and `p=2`.
- **normal** `LoveNumbersCliTests.SurfaceFluidIsRejectedClearly` — `tests/test_cli.cpp` — Runs the CLI on a surface-fluid model and checks the explicit no-ocean error.

## gia3D validation

- **paper-validation** `LoveNumbersGia3DReference.ReportsDirectComparison` — `validation/gia3d/compare_reference.py` — Runs the pinned original-truncation PREM driver, reports three independent solutions and stock combined-load agreement, and requires finite residuals and degree-one zeros.
- **paper-validation** `LoveNumbersGia3DReference.ReportsControlledSphereComparison` — `validation/compare_controlled_sphere.py` — Compares matched homogeneous-sphere meshes, dimensions, gravity, residuals, reciprocity, and eight response components.
- **paper-validation** `LoveNumbersGia3DReference.ReportsSolidFluidSolidComparison` — `validation/compare_controlled_sphere.py` — Repeats the matched comparison for two interface orientations and a non-zero fluid density derivative, now through `N=256`.
- **paper-validation** `LoveNumbersGia3DReference.ReportsCentralFluidComparison` — `validation/compare_controlled_sphere.py` — Checks a central fluid for finite matrices and responses, zero stratification, matching topology, and no invalid floating-point flags.
- **paper-validation** `LoveNumbersGia3DReference.ReportsIsotropicPremComparison` — `validation/compare_isotropic_prem.py` — Performs the analytic-gia3D, sampled-gia3D, and sampled-DSpecM1D PREM comparison across mesh and knot refinements.

## ELLN validation

- **paper-validation** `LoveNumbersELLNReference.ValidatesConversion` — `validation/elln/validate_elln.py` — Algebraically round-trips the dimensional-to-conventional `h_l,k_l` conversions and their signs.
- **paper-validation** `LoveNumbersELLNReference.AuditsOfficialExample` — `validation/elln/validate_elln.py` — Verifies the official archive, source and model checksums, runs the non-GUI example when possible, and otherwise reports the publisher table.
- **paper-validation** `LoveNumbersELLNReference.AuditsDegreeOneConvention` — `validation/elln/validate_elln.py` — Traces ELLN's deformed-centre subtraction and demonstrates why it is not the `P(a)=0` generic conversion.
- **paper-validation** `LoveNumbersELLNReference.ReportsTiPremComparison` — `validation/elln/compare_ti_prem.py` — Validates generated TI decks and reports degree-10 Table S7 differences, mesh/knot convergence, residuals, and the differing gravitational constants.
- **paper-validation** `LoveNumbersELLNReference.ReportsGravitationalConstantSensitivity` — `validation/elln/compare_gravitational_constants.py` — Compares production and ELLN `G`, reports sensitivity and explained discrepancy, and audits exported `rho,A,C,F,L,N` against the official model.

## Other gated diagnostics

- **paper-validation** `LoveNumbersReferenceValidation.ReportsPremConvergence` — `tests/test_reference_validation.cpp` — Reports selected-degree convergence against the published `love.dat` subset without treating the sampled-model comparison as a regression oracle.
