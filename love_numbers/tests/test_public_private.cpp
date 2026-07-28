#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>

// Keep the public implementation and private inline solver in one translation
// unit because EarthMesh defines non-inline functions in its public header.
#include "../src/LoveNumbers.cpp"

namespace {

using DSpecM1D::LoveNumbers::Config;
using DSpecM1D::LoveNumbers::DegreeResult;
using DSpecM1D::LoveNumbers::calculate;
using DSpecM1D::LoveNumbers::detail::RadialState;
using DSpecM1D::LoveNumbers::detail::solveStaticDegree;

double relativeResultError(const DegreeResult &actual,
                           const DegreeResult &expected) {
  const std::array<double, 6> actual_values{
      actual.h_u,   actual.k_u, actual.h_phi,
      actual.k_phi, actual.h_t, actual.k_t,
  };
  const std::array<double, 6> expected_values{
      expected.h_u,   expected.k_u, expected.h_phi,
      expected.k_phi, expected.h_t, expected.k_t,
  };
  double difference = 0.0;
  double scale = 1.0;
  for (std::size_t index = 0; index < actual_values.size(); ++index) {
    difference =
        std::max(difference,
                 std::abs(actual_values[index] - expected_values[index]));
    scale = std::max(scale, std::abs(expected_values[index]));
  }
  return difference / scale;
}

TEST(LoveNumbersPublicPrivateTests, PublicResultsMatchPrivateSolves) {
  const EarthModels::ModelInput<double> model(
      DSPECM1D_LOVE_NUMBERS_TI_MODEL);
  const Config config{
      .maximum_degree = 10,
      .polynomial_order = 2,
      .maximum_radial_step = 0.4,
  };
  const std::vector<DegreeResult> public_results =
      calculate(model, config);
  RadialState state(model, config);

  ASSERT_EQ(public_results.size(), 11);
  for (const int degree : {0, 1, 2, 10}) {
    const DegreeResult private_result =
        solveStaticDegree(model, state, degree);
    const double error =
        relativeResultError(public_results[degree], private_result);

    std::cout << std::setprecision(17) << "l=" << degree
              << " public_private_error=" << error << '\n';
    EXPECT_EQ(public_results[degree].degree, degree);
    EXPECT_LE(error, 1.0e-13);
  }
}

}   // namespace
