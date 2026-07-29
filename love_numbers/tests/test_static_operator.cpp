#include <gtest/gtest.h>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <string>

#include <DSpecM1D/ModelInput>
#include <DSpecM1D/src/SEM/SEM.h>

#include "LoveDofMap.h"
#include "RadialState.h"
#include "StaticOperator.h"

namespace {

using DSpecM1D::LoveNumbers::Config;
using DSpecM1D::LoveNumbers::detail::Field;
using DSpecM1D::LoveNumbers::detail::LoveDofMap;
using DSpecM1D::LoveNumbers::detail::RadialState;
using DSpecM1D::LoveNumbers::detail::assembleStaticOperator;

constexpr double comparison_tolerance = 1.0e-12;

Config testConfig() {
  return Config{
      .maximum_degree = 10,
      .polynomial_order = 4,
      .maximum_radial_step = 0.25,
  };
}

void compareWithReference(const std::string &path,
                          const std::string &model_name) {
  const EarthModels::ModelInput<double> model(path);
  const Config config = testConfig();
  RadialState state(model, config);
  const Full1D::SEM sem(model, config.maximum_radial_step,
                        config.polynomial_order + 1,
                        config.maximum_degree);

  for (const int degree : {2, 10}) {
    const Eigen::SparseMatrix<double> matrix =
        assembleStaticOperator(state, degree);
    const Eigen::SparseMatrix<double> reference = sem.hS(degree);
    const Eigen::SparseMatrix<double> difference = matrix - reference;
    const Eigen::SparseMatrix<double> asymmetry =
        matrix - Eigen::SparseMatrix<double>(matrix.transpose());
    const double relative_difference =
        difference.norm() / std::max(1.0, reference.norm());
    const double symmetry_error =
        asymmetry.norm() / std::max(1.0, matrix.norm());

    std::cout << std::setprecision(17) << model_name << " l=" << degree
              << " relative_difference=" << relative_difference
              << " symmetry_error=" << symmetry_error << '\n';

    EXPECT_EQ(matrix.rows(), reference.rows());
    EXPECT_EQ(matrix.cols(), reference.cols());
    EXPECT_LE(relative_difference, comparison_tolerance);
    EXPECT_LE(symmetry_error, comparison_tolerance);
  }
}

TEST(StaticOperatorTests, MatchesIsotropicReference) {
  compareWithReference(DSPECM1D_LOVE_NUMBERS_ISOTROPIC_MODEL, "isotropic");
}

TEST(StaticOperatorTests, MatchesTransverselyIsotropicReference) {
  compareWithReference(DSPECM1D_LOVE_NUMBERS_TI_MODEL, "TI");
}

TEST(StaticOperatorTests, RetainsCentreFieldsAndMatchesSemOrdering) {
  const EarthModels::ModelInput<double> model(
      DSPECM1D_LOVE_NUMBERS_TI_MODEL);
  const Config config = testConfig();
  RadialState state(model, config);
  const Full1D::SEM sem(model, config.maximum_radial_step,
                        config.polynomial_order + 1,
                        config.maximum_degree);
  const LoveDofMap map(state.mesh(), 2);

  ASSERT_TRUE(map.index(Field::U, 0, 0));
  ASSERT_TRUE(map.index(Field::V, 0, 0));
  ASSERT_TRUE(map.index(Field::P, 0, 0));
  EXPECT_NE(*map.index(Field::U, 0, 0), *map.index(Field::V, 0, 0));
  EXPECT_NE(*map.index(Field::U, 0, 0), *map.index(Field::P, 0, 0));
  EXPECT_NE(*map.index(Field::V, 0, 0), *map.index(Field::P, 0, 0));

  for (int element = 0; element < state.mesh().NE(); ++element) {
    for (int node = 0; node < state.mesh().NN(); ++node) {
      EXPECT_EQ(*map.index(Field::U, element, node),
                sem.ltgS(0, element, node));
      EXPECT_EQ(*map.index(Field::V, element, node),
                sem.ltgS(1, element, node));
      EXPECT_EQ(*map.index(Field::P, element, node),
                sem.ltgS(2, element, node));
    }
  }
}

TEST(StaticOperatorTests, RejectsNegativeDegree) {
  const EarthModels::ModelInput<double> model(
      DSPECM1D_LOVE_NUMBERS_ISOTROPIC_MODEL);
  const Config config = testConfig();
  RadialState state(model, config);

  EXPECT_THROW(assembleStaticOperator(state, -1), std::invalid_argument);
}

}   // namespace
