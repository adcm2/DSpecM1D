#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/SparseLU>

#include <DSpecM1D/src/SEM/SEM.h>

#include "LoveDofMap.h"
#include "RadialState.h"
#include "StaticForcing.h"
#include "StaticOperator.h"
#include "StaticSolve.h"

namespace {

using DSpecM1D::LoveNumbers::Config;
using DSpecM1D::LoveNumbers::DegreeResult;
using DSpecM1D::LoveNumbers::detail::Field;
using DSpecM1D::LoveNumbers::detail::LoveDofMap;
using DSpecM1D::LoveNumbers::detail::RadialState;
using DSpecM1D::LoveNumbers::detail::assembleStaticOperator;
using DSpecM1D::LoveNumbers::detail::assembleStaticRightHandSides;
using DSpecM1D::LoveNumbers::detail::solveStaticDegree;

constexpr double matrix_tolerance = 1.0e-12;
constexpr double solve_tolerance = 1.0e-11;

Config smallModelConfig() {
  return Config{
      .maximum_degree = 2,
      .polynomial_order = 4,
      .maximum_radial_step = 0.4,
  };
}

Config premConfig() {
  return Config{
      .maximum_degree = 2,
      .polynomial_order = 3,
      .maximum_radial_step = 0.2,
  };
}

void compareWithRadialReference(const std::string &path,
                                const std::string &name,
                                const Config &config) {
  const EarthModels::ModelInput<double> model(path);
  RadialState state(model, config);
  const LoveDofMap dofs(state.mesh(), 0);
  const Full1D::SEM sem(model, config.maximum_radial_step,
                        config.polynomial_order + 1,
                        config.maximum_degree);
  const Eigen::SparseMatrix<double> matrix =
      assembleStaticOperator(state, dofs, 0);
  const Eigen::SparseMatrix<double> reference = sem.hR();
  const Eigen::SparseMatrix<double> difference = matrix - reference;
  const Eigen::SparseMatrix<double> asymmetry =
      matrix - Eigen::SparseMatrix<double>(matrix.transpose());
  const double relative_difference =
      difference.norm() / std::max(1.0, reference.norm());
  const double symmetry_error =
      asymmetry.norm() / std::max(1.0, matrix.norm());

  std::cout << std::setprecision(17) << name
            << " radial_matrix_error=" << relative_difference
            << " symmetry_error=" << symmetry_error << '\n';

  EXPECT_EQ(matrix.rows(), reference.rows());
  EXPECT_EQ(matrix.cols(), reference.cols());
  EXPECT_EQ(matrix.rows(), dofs.size());
  EXPECT_LE(relative_difference, matrix_tolerance);
  EXPECT_LE(symmetry_error, matrix_tolerance);

  for (int element = 0; element < state.mesh().NE(); ++element) {
    for (int node = 0; node < state.mesh().NN(); ++node) {
      EXPECT_EQ(*dofs.index(Field::U, element, node),
                sem.ltgR(0, element, node));
      EXPECT_EQ(*dofs.index(Field::P, element, node),
                sem.ltgR(1, element, node));
    }
  }
}

TEST(StaticRadialOperatorTests, MatchesRadialSemForSupportedModels) {
  compareWithRadialReference(DSPECM1D_LOVE_NUMBERS_ISOTROPIC_MODEL,
                             "isotropic", smallModelConfig());
  compareWithRadialReference(DSPECM1D_LOVE_NUMBERS_TI_MODEL, "TI",
                             smallModelConfig());
}

TEST(StaticRadialOperatorTests, RetainsCentreFields) {
  const EarthModels::ModelInput<double> model(
      DSPECM1D_LOVE_NUMBERS_SFS_MODEL);
  RadialState state(model, smallModelConfig());
  const LoveDofMap dofs(state.mesh(), 0);

  ASSERT_TRUE(dofs.index(Field::U, 0, 0));
  ASSERT_TRUE(dofs.index(Field::P, 0, 0));
  EXPECT_EQ(*dofs.index(Field::U, 0, 0), 0);
  EXPECT_EQ(*dofs.index(Field::P, 0, 0), 1);
}

TEST(StaticRadialForcingTests, SurfaceLoadsAndZeroTideAreExact) {
  const EarthModels::ModelInput<double> model(
      DSPECM1D_LOVE_NUMBERS_SFS_MODEL);
  RadialState state(model, smallModelConfig());
  const LoveDofMap dofs(state.mesh(), 0);
  const Eigen::MatrixXd right_hand_sides =
      assembleStaticRightHandSides(state, dofs, 0);
  const int surface_element = state.mesh().NE() - 1;
  const int surface_node = state.mesh().NN() - 1;
  const int surface_u =
      *dofs.index(Field::U, surface_element, surface_node);
  const int surface_p =
      *dofs.index(Field::P, surface_element, surface_node);
  const double radius = state.surfaceRadius();

  EXPECT_DOUBLE_EQ(right_hand_sides(surface_u, 0),
                   -radius * radius * state.surfaceGravity());
  EXPECT_DOUBLE_EQ(right_hand_sides(surface_p, 1),
                   -radius * radius);
  EXPECT_EQ(right_hand_sides.col(2).norm(), 0.0);

  for (int row = 0; row < dofs.size(); ++row) {
    if (row != surface_u) {
      EXPECT_EQ(right_hand_sides(row, 0), 0.0);
    }
    if (row != surface_p) {
      EXPECT_EQ(right_hand_sides(row, 1), 0.0);
    }
  }
}

void checkRadialSolve(const std::string &path, const std::string &name,
                      const Config &config) {
  const EarthModels::ModelInput<double> model(path);
  RadialState state(model, config);
  const LoveDofMap dofs(state.mesh(), 0);
  const Eigen::SparseMatrix<double> matrix =
      assembleStaticOperator(state, dofs, 0);
  const Eigen::MatrixXd right_hand_sides =
      assembleStaticRightHandSides(state, dofs, 0);

  Eigen::SparseLU<Eigen::SparseMatrix<double>> factorization;
  factorization.compute(matrix);
  ASSERT_EQ(factorization.info(), Eigen::Success);
  const Eigen::MatrixXd solutions =
      factorization.solve(right_hand_sides);
  ASSERT_EQ(factorization.info(), Eigen::Success);

  const double displacement_residual =
      (matrix * solutions.col(0) - right_hand_sides.col(0)).norm() /
      std::max(1.0, right_hand_sides.col(0).norm());
  const double gravitational_residual =
      (matrix * solutions.col(1) - right_hand_sides.col(1)).norm() /
      std::max(1.0, right_hand_sides.col(1).norm());
  const DegreeResult result = solveStaticDegree(model, state, 0);
  const double gravity_si =
      state.surfaceGravity() * model.AccelerationNorm();
  const double reciprocity_left = gravity_si * result.h_phi;
  const double reciprocity_scale =
      std::max({1.0, std::abs(reciprocity_left),
                std::abs(result.k_u)});
  const double reciprocity_error =
      std::abs(reciprocity_left - result.k_u) / reciprocity_scale;

  const std::array<double, 6> values{
      result.h_u, result.k_u, result.h_phi,
      result.k_phi, result.h_t, result.k_t,
  };
  EXPECT_TRUE(std::all_of(values.begin(), values.end(),
                          [](double value) {
                            return std::isfinite(value);
                          }));
  EXPECT_EQ(result.degree, 0);
  EXPECT_LE(displacement_residual, solve_tolerance);
  EXPECT_LE(gravitational_residual, solve_tolerance);
  EXPECT_LE(reciprocity_error, solve_tolerance);
  EXPECT_EQ(right_hand_sides.col(2).norm(), 0.0);
  EXPECT_EQ(solutions.col(2).norm(), 0.0);
  EXPECT_EQ(result.h_t, 0.0);
  EXPECT_EQ(result.k_t, 0.0);

  std::cout << std::setprecision(17) << name
            << " radial_displacement_residual="
            << displacement_residual
            << " radial_gravitational_residual="
            << gravitational_residual
            << " radial_reciprocity_error=" << reciprocity_error
            << " tidal_rhs_norm=" << right_hand_sides.col(2).norm()
            << '\n';
}

TEST(StaticRadialSolveTests, SolvesSupportedModels) {
  const std::vector<std::pair<std::string, std::string>> models{
      {"isotropic", DSPECM1D_LOVE_NUMBERS_ISOTROPIC_MODEL},
      {"TI", DSPECM1D_LOVE_NUMBERS_TI_MODEL},
      {"S-F-S", DSPECM1D_LOVE_NUMBERS_SFS_MODEL},
      {"PREM no-ocean", DSPECM1D_LOVE_NUMBERS_SOLID_SURFACE_MODEL},
  };

  for (const auto &[name, path] : models) {
    checkRadialSolve(path, name,
                     name == "PREM no-ocean" ? premConfig()
                                              : smallModelConfig());
  }
}

}   // namespace
