#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <DSpecM1D/LoveNumbers>

namespace {

using DSpecM1D::LoveNumbers::Config;
using DSpecM1D::LoveNumbers::DegreeResult;
using DSpecM1D::LoveNumbers::calculate;

Config smallConfig(int maximum_degree) {
  return Config{
      .maximum_degree = maximum_degree,
      .polynomial_order = 2,
      .maximum_radial_step = 0.4,
  };
}

TEST(LoveNumbersTests, LoadSumsCombineGeneralizedResponses) {
  const DegreeResult result{
      .degree = 2,
      .h_u = 1.25,
      .k_u = -0.5,
      .h_phi = 0.75,
      .k_phi = 0.125,
      .h_t = 0.0,
      .k_t = 0.0,
  };

  EXPECT_DOUBLE_EQ(result.h_load(), 2.0);
  EXPECT_DOUBLE_EQ(result.k_load(), -0.375);
}

TEST(LoveNumbersTests, PublicHeaderCalculatesOrderedDegreeRange) {
  const EarthModels::ModelInput<double> model(
      DSPECM1D_LOVE_NUMBERS_ISOTROPIC_MODEL);
  const std::vector<DegreeResult> results =
      calculate(model, smallConfig(2));

  ASSERT_EQ(results.size(), 3);
  for (int degree = 0; degree <= 2; ++degree) {
    EXPECT_EQ(results[degree].degree, degree);
  }

  EXPECT_EQ(results[0].h_t, 0.0);
  EXPECT_EQ(results[0].k_t, 0.0);
  EXPECT_EQ(results[1].h_phi, 0.0);
  EXPECT_EQ(results[1].k_u, 0.0);
  EXPECT_EQ(results[1].k_phi, 0.0);
  EXPECT_EQ(results[1].k_t, 0.0);
}

TEST(LoveNumbersTests, MaximumDegreeZeroIsValid) {
  const EarthModels::ModelInput<double> model(
      DSPECM1D_LOVE_NUMBERS_ISOTROPIC_MODEL);
  const std::vector<DegreeResult> results =
      calculate(model, smallConfig(0));

  ASSERT_EQ(results.size(), 1);
  EXPECT_EQ(results.front().degree, 0);
}

TEST(LoveNumbersTests, SupportedModelsReturnFiniteValues) {
  const std::vector<std::pair<std::string, std::string>> models{
      {"isotropic", DSPECM1D_LOVE_NUMBERS_ISOTROPIC_MODEL},
      {"TI", DSPECM1D_LOVE_NUMBERS_TI_MODEL},
      {"S-F-S", DSPECM1D_LOVE_NUMBERS_SFS_MODEL},
      {"PREM no-ocean", DSPECM1D_LOVE_NUMBERS_SOLID_SURFACE_MODEL},
  };

  for (const auto &[name, path] : models) {
    const EarthModels::ModelInput<double> model(path);
    Config config = smallConfig(2);
    if (name == "PREM no-ocean") {
      config.maximum_radial_step = 0.2;
    }
    const std::vector<DegreeResult> results = calculate(model, config);

    ASSERT_EQ(results.size(), 3);
    for (const DegreeResult &result : results) {
      const std::array<double, 6> values{
          result.h_u,   result.k_u, result.h_phi,
          result.k_phi, result.h_t, result.k_t,
      };
      EXPECT_TRUE(std::all_of(
          values.begin(), values.end(),
          [](double value) { return std::isfinite(value); }));
      EXPECT_DOUBLE_EQ(result.h_load(), result.h_u + result.h_phi);
      EXPECT_DOUBLE_EQ(result.k_load(), result.k_u + result.k_phi);
    }
  }
}

TEST(LoveNumbersTests, RejectsInvalidConfigurationAndSurfaceFluid) {
  const EarthModels::ModelInput<double> model(
      DSPECM1D_LOVE_NUMBERS_ISOTROPIC_MODEL);

  EXPECT_THROW(calculate(
                   model,
                   Config{.maximum_degree = -1,
                          .polynomial_order = 2,
                          .maximum_radial_step = 0.4}),
               std::invalid_argument);
  EXPECT_THROW(calculate(
                   model,
                   Config{.maximum_degree = 2,
                          .polynomial_order = 0,
                          .maximum_radial_step = 0.4}),
               std::invalid_argument);
  EXPECT_THROW(calculate(
                   model,
                   Config{.maximum_degree = 2,
                          .polynomial_order = 2,
                          .maximum_radial_step = 0.0}),
               std::invalid_argument);

  const EarthModels::ModelInput<double> ocean_model(
      DSPECM1D_LOVE_NUMBERS_OCEAN_MODEL);
  EXPECT_THROW(calculate(ocean_model, smallConfig(2)),
               std::invalid_argument);
}

}   // namespace
