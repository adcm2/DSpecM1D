#include <gtest/gtest.h>
#include <cmath>
#include <DSpecM1D/src/SpecHelpers.h>

TEST(SpecHelpersTests, ResolveModeFlagsHonorsModeTypeAndAngularRange) {
  const auto allModes = SPARSESPEC::resolveModeFlags(4, 0, 6);
  EXPECT_TRUE(allModes.inc_rad);
  EXPECT_TRUE(allModes.inc_tor);
  EXPECT_TRUE(allModes.inc_sph);

  const auto torOnly = SPARSESPEC::resolveModeFlags(2, 1, 6);
  EXPECT_FALSE(torOnly.inc_rad);
  EXPECT_TRUE(torOnly.inc_tor);
  EXPECT_FALSE(torOnly.inc_sph);

  const auto radialSuppressed = SPARSESPEC::resolveModeFlags(1, 1, 6);
  EXPECT_FALSE(radialSuppressed.inc_rad);
}

TEST(SpecHelpersTests, OutputFactorMatchesRequestedQuantity) {
  const std::complex<double> myi{0.0, 1.0};
  const std::complex<double> w{2.0, 0.0};

  EXPECT_EQ(SPARSESPEC::outputFactor(1, w, myi), std::complex<double>(1.0, 0.0));
  EXPECT_EQ(SPARSESPEC::outputFactor(2, w, myi), std::complex<double>(0.0, 2.0));
  EXPECT_EQ(SPARSESPEC::outputFactor(0, w, myi), std::complex<double>(0.0, -0.5));
}

TEST(SpecHelpersTests, SpecConstantsUsesReferencePeriod) {
  SPARSESPEC::SpecConstants constants(0.1, 10.0);
  EXPECT_NEAR(constants.w0, (2.0 * M_PI) / 10.0, 1e-12);
  EXPECT_NEAR(constants.twodivpi, 2.0 / M_PI, 1e-12);
}

TEST(SpecHelpersTests, FactorizeOrComputeChoosesExpectedPath) {
  struct FakeSolver {
    int computeCalls = 0;
    int factorizeCalls = 0;

    void compute(int) { ++computeCalls; }
    void factorize(int) { ++factorizeCalls; }
  };

  FakeSolver solver;
  SPARSESPEC::factorizeOrCompute(solver, 42, 0, 4);
  SPARSESPEC::factorizeOrCompute(solver, 42, 1, 4);

  EXPECT_EQ(solver.computeCalls, 1);
  EXPECT_EQ(solver.factorizeCalls, 1);
}
