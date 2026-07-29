#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <DSpecM1D/LoveNumbers>

namespace {

using DSpecM1D::LoveNumbers::Config;
using DSpecM1D::LoveNumbers::DegreeResult;
using DSpecM1D::LoveNumbers::calculate;

struct ReferenceResult {
  int degree;
  std::array<double, 6> values;
};

std::vector<ReferenceResult> readReference() {
  std::ifstream input(DSPECM1D_LOVE_NUMBERS_REFERENCE);
  std::string header;
  std::getline(input, header);

  std::vector<ReferenceResult> reference;
  ReferenceResult result;
  while (input >> result.degree >> result.values[0] >> result.values[1] >>
         result.values[2] >> result.values[3] >> result.values[4] >>
         result.values[5]) {
    reference.push_back(result);
  }
  return reference;
}

std::array<double, 6> values(const DegreeResult &result) {
  return {
      result.h_u,   result.k_u, result.h_phi,
      result.k_phi, result.h_t, result.k_t,
  };
}

double relativeError(double calculated, double reference) {
  const double difference = calculated - reference;
  if (reference != 0.0) {
    return std::abs(difference / reference);
  }
  return difference == 0.0
             ? 0.0
             : std::numeric_limits<double>::infinity();
}

TEST(LoveNumbersReferenceValidation, ReportsPremConvergence) {
  const EarthModels::ModelInput<double> model(
      DSPECM1D_LOVE_NUMBERS_REFERENCE_MODEL);
  const std::vector<ReferenceResult> reference = readReference();
  const std::array<const char *, 6> names{
      "h_u", "k_u", "h_phi", "k_phi", "h_t", "k_t",
  };

  ASSERT_FALSE(reference.empty());
  EXPECT_DOUBLE_EQ(model.OuterRadius() * model.LengthNorm(), 6368000.0);

  for (const double maximum_radial_step :
       {0.02, 0.01, 0.005, 0.0025, 0.00125}) {
    const Config config{
        .maximum_degree = reference.back().degree,
        .polynomial_order = 6,
        .maximum_radial_step = maximum_radial_step,
    };
    const std::vector<DegreeResult> calculated =
        calculate(model, config);

    for (const ReferenceResult &expected : reference) {
      const DegreeResult &actual = calculated[expected.degree];
      const auto actual_values = values(actual);
      ASSERT_EQ(actual.degree, expected.degree);

      for (std::size_t component = 0; component < names.size();
           ++component) {
        const double difference =
            actual_values[component] - expected.values[component];
        const double relative_error =
            relativeError(actual_values[component],
                          expected.values[component]);
        EXPECT_TRUE(std::isfinite(actual_values[component]));
        std::cout << std::setprecision(17)
                  << "step=" << maximum_radial_step
                  << " l=" << expected.degree
                  << " component=" << names[component]
                  << " calculated=" << actual_values[component]
                  << " reference=" << expected.values[component]
                  << " signed_difference=" << difference
                  << " relative_error=" << relative_error << '\n';
      }
    }
  }
}

}   // namespace
