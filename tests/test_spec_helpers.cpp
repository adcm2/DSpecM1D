#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <DSpecM1D/src/SpecHelpers.h>

namespace {

std::complex<double> independentAttenOracle(double omega, double epsilon,
                                            double w0) {
  const double pi = std::acos(-1.0);
  return {(2.0 / pi) * std::log(omega / w0) +
              (1.0 / pi) * std::log1p((epsilon / omega) * (epsilon / omega)),
          1.0 - (2.0 / pi) * std::atan(epsilon / omega)};
}

void expectComplexNear(std::complex<double> actual,
                       std::complex<double> expected, double absTol,
                       double relTol) {
  EXPECT_NEAR(std::abs(actual - expected), 0.0,
              absTol + relTol * std::abs(expected));
}

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

TEST(SpecHelpersTests, ComplexAttenuationMatchesIndependentRealOracle) {
  constexpr double w0 = 0.37;
  for (const auto &[omega, epsilon] :
       {std::pair{0.2, 0.001}, std::pair{1.7, 0.13},
        std::pair{5.0, 0.8}}) {
    const auto z = std::complex<double>(omega, -epsilon);
    const auto actual = SPARSESPEC::attenFactor(z, w0);
    const auto expected = independentAttenOracle(omega, epsilon, w0);
    expectComplexNear(actual, expected, 1e-12, 1e-12);
    EXPECT_GT(actual.real() -
                  (2.0 / std::acos(-1.0)) * std::log(omega / w0),
              0.0);
    EXPECT_LT(actual.imag(), 1.0);
  }
}

TEST(SpecHelpersTests, ComplexAttenuationHasPositiveAxisLimitAndReality) {
  constexpr double w0 = 0.61;
  constexpr double omega = 1.4;
  const auto nearAxis = SPARSESPEC::attenFactor(
      std::complex<double>(omega, -1e-15), w0);
  const std::complex<double> axisLimit{
      (2.0 / std::acos(-1.0)) * std::log(omega / w0), 1.0};
  expectComplexNear(nearAxis, axisLimit, 1e-12, 1e-12);

  const std::complex<double> z{2.3, -0.17};
  const auto positiveFrequency = SPARSESPEC::attenFactor(z, w0);
  const auto negativeFrequency = SPARSESPEC::attenFactor(-std::conj(z), w0);
  expectComplexNear(negativeFrequency, std::conj(positiveFrequency), 1e-12,
                    1e-12);
}

TEST(SpecHelpersTests, ComplexAttenuationIsUnitScaleInvariant) {
  const std::complex<double> zSi{1.25, -0.09};
  constexpr double w0Si = 0.42;
  constexpr double unitsPerSi = 3600.0;
  const auto si = SPARSESPEC::attenFactor(zSi, w0Si);
  const auto nondimensional = SPARSESPEC::attenFactor(
      zSi / unitsPerSi, w0Si / unitsPerSi);
  expectComplexNear(nondimensional, si, 1e-12, 1e-12);
}

TEST(SpecHelpersTests, RealAttenuationOverloadRemainsSourceCompatible) {
  const std::complex<double> myi{0.0, 1.0};
  constexpr double omega = 1.9;
  constexpr double w0 = 0.55;
  constexpr double twodivpi = 2.0 / 3.14159265358979323846;
  const auto legacy =
      SPARSESPEC::attenFactor(omega, w0, twodivpi, myi);
  const auto expected =
      myi + twodivpi * std::log(omega / w0);
  EXPECT_EQ(legacy, expected);
}

TEST(SpecHelpersTests, ComplexAttenuationRejectsInvalidDomainAndReference) {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  EXPECT_THROW(SPARSESPEC::attenFactor(std::complex<double>{0.0, 0.0}, 1.0),
               std::invalid_argument);
  EXPECT_THROW(SPARSESPEC::attenFactor(std::complex<double>{nan, -0.1}, 1.0),
               std::invalid_argument);
  EXPECT_THROW(SPARSESPEC::attenFactor(std::complex<double>{1.0, -inf}, 1.0),
               std::invalid_argument);
  for (const double w0 : {0.0, -1.0, nan, inf})
    EXPECT_THROW(SPARSESPEC::attenFactor(std::complex<double>{1.0, -0.1}, w0),
                 std::invalid_argument);
  for (const double tref : {0.0, -1.0, nan, inf,
                            std::numeric_limits<double>::min()})
    EXPECT_THROW(SPARSESPEC::validateAttenuationReferencePeriod(tref),
                 std::invalid_argument);
  EXPECT_NO_THROW(SPARSESPEC::validateAttenuationReferencePeriod(3600.0));
}
