#include <gtest/gtest.h>
#include <cmath>
#include <DSpecM1D/src/SpecHelpers.h>
#include <DSpecM1D/src/StartElement.h>

namespace {

struct FakeStartQuadrature {
  double W(int) const { return 1.0; }
};

struct FakeStartMesh {
  int NE() const { return 4; }
  int NN() const { return 2; }
  double EW(int) const { return 100000.0; }
  double NodeRadius(int idxe, int idxq) const {
    return 1000.0 * static_cast<double>(idxe + 1) + idxq;
  }
  FakeStartQuadrature GLL() const { return {}; }
};

struct FakeStartModel {
  double PSlow(int, int) const { return 0.0; }
  double SSlow(int, int) const { return 0.0; }
};

struct FakeStartSem {
  FakeStartMesh mesh() const { return {}; }
  FakeStartModel meshModel() const { return {}; }
  int el() const { return 0; }
};

}   // namespace

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

TEST(SpecHelpersTests, StartElementCadenceTraversesHighToLow) {
  const std::vector<double> frequencies{1.0, 2.0, 3.0, 4.0,
                                        5.0, 6.0, 7.0};
  std::vector<double> recomputed;

  const auto starts = SpectralTools::detail::allIndicesImpl(
      frequencies, 1, static_cast<int>(frequencies.size()), 3,
      [&](double frequency) {
        recomputed.push_back(frequency);
        return static_cast<int>(frequency * 10.0);
      });

  EXPECT_EQ(recomputed, (std::vector<double>{7.0, 4.0}));
  EXPECT_EQ(starts, (std::vector<int>{40, 40, 40, 70, 70, 70}));
}

TEST(SpecHelpersTests, SourceConstrainedStartRetainsSourceElement) {
  FakeStartSem sem;

  EXPECT_EQ(SpectralTools::startElementTor(sem, 1, 1.0), 3);
  EXPECT_EQ(SpectralTools::startElementTor(sem, 1, 1.0, 2), 1);
}
