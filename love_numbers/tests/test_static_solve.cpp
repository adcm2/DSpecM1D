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

constexpr double solve_tolerance = 1.0e-11;

Config smallModelConfig() {
  return Config{
      .maximum_degree = 10,
      .polynomial_order = 4,
      .maximum_radial_step = 0.4,
  };
}

double relativeError(const Eigen::VectorXd &actual,
                     const Eigen::VectorXd &expected) {
  return (actual - expected).norm() / std::max(1.0, expected.norm());
}

Eigen::VectorXd independentTidalRightHandSide(
    RadialState &state, const LoveDofMap &dofs, int degree,
    int *solid_to_fluid_interfaces = nullptr,
    int *fluid_to_solid_interfaces = nullptr) {
  const auto &mesh = state.mesh();
  const auto &model = state.meshModel();
  const auto &quadrature = state.quadrature();
  const double surface_radius = state.surfaceRadius();
  const double angular_degree = degree * (degree + 1);
  Eigen::VectorXd load_functional =
      Eigen::VectorXd::Zero(dofs.size());

  for (int element = 0; element < mesh.NE(); ++element) {
    const double half_width = 0.5 * mesh.EW(element);
    for (int point = 0; point < mesh.NN(); ++point) {
      const double radius = mesh.NodeRadius(element, point);
      const double weight = half_width * quadrature.W(point);
      const double tide_potential =
          std::pow(radius / surface_radius, degree);

      if (mesh.IsSolid(element)) {
        const double value =
            weight * model.Density(element, point) * radius *
            tide_potential;
        const int u = *dofs.index(Field::U, element, point);
        const int v = *dofs.index(Field::V, element, point);
        load_functional(u) += degree * value;
        load_functional(v) += angular_degree * value;
      } else {
        double density_gradient_coefficient = 0.0;
        if (radius != 0.0) {
          density_gradient_coefficient =
              radius * radius * state.densityDerivative(element, point) /
              state.gravity(element, point);
        }
        const int p = *dofs.index(Field::P, element, point);
        load_functional(p) +=
            weight * density_gradient_coefficient * tide_potential;
      }
    }
  }

  for (int lower = 0; lower + 1 < mesh.NE(); ++lower) {
    const int upper = lower + 1;
    int orientation = 0;
    int solid_element = 0;
    int solid_node = 0;
    int fluid_element = 0;
    int fluid_node = 0;

    if (mesh.IsSolid(lower) && mesh.IsFluid(upper)) {
      orientation = -1;
      solid_element = lower;
      solid_node = mesh.NN() - 1;
      fluid_element = upper;
      fluid_node = 0;
      if (solid_to_fluid_interfaces != nullptr) {
        ++*solid_to_fluid_interfaces;
      }
    } else if (mesh.IsFluid(lower) && mesh.IsSolid(upper)) {
      orientation = 1;
      solid_element = upper;
      solid_node = 0;
      fluid_element = lower;
      fluid_node = mesh.NN() - 1;
      if (fluid_to_solid_interfaces != nullptr) {
        ++*fluid_to_solid_interfaces;
      }
    } else {
      continue;
    }

    const double radius = mesh.NodeRadius(solid_element, solid_node);
    const double density = model.Density(fluid_element, fluid_node);
    const double tide_potential =
        std::pow(radius / surface_radius, degree);
    const int u = *dofs.index(Field::U, solid_element, solid_node);
    load_functional(u) +=
        orientation * density * radius * radius * tide_potential;
  }

  return -load_functional;
}

void expectNear(double actual, double expected) {
  EXPECT_NEAR(actual, expected,
              solve_tolerance * std::max(1.0, std::abs(expected)));
}

TEST(StaticForcingTests, SurfaceLoadsHaveExpectedEntriesAndSigns) {
  const EarthModels::ModelInput<double> model(
      DSPECM1D_LOVE_NUMBERS_ISOTROPIC_MODEL);
  RadialState state(model, smallModelConfig());
  const LoveDofMap dofs(state.mesh(), 2);
  const Eigen::MatrixXd right_hand_sides =
      assembleStaticRightHandSides(state, dofs, 2);
  const int surface_element = state.mesh().NE() - 1;
  const int surface_node = state.mesh().NN() - 1;
  const int surface_u =
      *dofs.index(Field::U, surface_element, surface_node);
  const int surface_p =
      *dofs.index(Field::P, surface_element, surface_node);
  const double radius = state.surfaceRadius();

  EXPECT_EQ(right_hand_sides.cols(), 3);
  EXPECT_DOUBLE_EQ(right_hand_sides(surface_u, 0),
                   -radius * radius * state.surfaceGravity());
  EXPECT_DOUBLE_EQ(right_hand_sides(surface_p, 1),
                   -radius * radius);
  for (int row = 0; row < dofs.size(); ++row) {
    if (row != surface_u) {
      EXPECT_EQ(right_hand_sides(row, 0), 0.0);
    }
    if (row != surface_p) {
      EXPECT_EQ(right_hand_sides(row, 1), 0.0);
    }
  }
}

TEST(StaticForcingTests, AllSolidTideMatchesIndependentAssembly) {
  const EarthModels::ModelInput<double> model(
      DSPECM1D_LOVE_NUMBERS_TI_MODEL);
  RadialState state(model, smallModelConfig());

  for (const int degree : {2, 10}) {
    const LoveDofMap dofs(state.mesh(), degree);
    const Eigen::MatrixXd right_hand_sides =
        assembleStaticRightHandSides(state, dofs, degree);
    const Eigen::VectorXd expected =
        independentTidalRightHandSide(state, dofs, degree);
    const double error =
        relativeError(right_hand_sides.col(2), expected);

    std::cout << std::setprecision(17) << "all-solid l=" << degree
              << " tidal_rhs_error=" << error << '\n';
    EXPECT_LE(error, solve_tolerance);
  }
}

TEST(StaticForcingTests,
     InternalFluidTideMatchesBothInterfaceOrientations) {
  const EarthModels::ModelInput<double> model(
      DSPECM1D_LOVE_NUMBERS_SFS_MODEL);
  RadialState state(model, smallModelConfig());

  for (const int degree : {2, 10}) {
    const LoveDofMap dofs(state.mesh(), degree);
    const Eigen::MatrixXd right_hand_sides =
        assembleStaticRightHandSides(state, dofs, degree);
    int solid_to_fluid_interfaces = 0;
    int fluid_to_solid_interfaces = 0;
    const Eigen::VectorXd expected = independentTidalRightHandSide(
        state, dofs, degree, &solid_to_fluid_interfaces,
        &fluid_to_solid_interfaces);
    const double error =
        relativeError(right_hand_sides.col(2), expected);

    std::cout << std::setprecision(17) << "S-F-S l=" << degree
              << " tidal_rhs_error=" << error << '\n';
    EXPECT_EQ(solid_to_fluid_interfaces, 1);
    EXPECT_EQ(fluid_to_solid_interfaces, 1);
    EXPECT_LE(error, solve_tolerance);
  }
}

void checkSolve(const std::string &path, const std::string &name,
                const Config &config, int degree) {
  const EarthModels::ModelInput<double> model(path);
  RadialState state(model, config);
  const LoveDofMap dofs(state.mesh(), degree);
  const Eigen::SparseMatrix<double> matrix =
      assembleStaticOperator(state, dofs, degree);
  const Eigen::MatrixXd right_hand_sides =
      assembleStaticRightHandSides(state, dofs, degree);

  Eigen::SparseLU<Eigen::SparseMatrix<double>> factorization;
  factorization.compute(matrix);
  ASSERT_EQ(factorization.info(), Eigen::Success);
  const Eigen::MatrixXd solutions =
      factorization.solve(right_hand_sides);
  ASSERT_EQ(factorization.info(), Eigen::Success);

  const double residual =
      (matrix * solutions - right_hand_sides).norm() /
      std::max(1.0, right_hand_sides.norm());
  const Eigen::VectorXd combined_solution =
      factorization.solve(right_hand_sides.col(0) +
                          right_hand_sides.col(1));
  ASSERT_EQ(factorization.info(), Eigen::Success);
  const double linearity_error =
      (combined_solution - solutions.col(0) - solutions.col(1)).norm() /
      std::max(1.0,
               (solutions.col(0) + solutions.col(1)).norm());

  const int surface_element = state.mesh().NE() - 1;
  const int surface_node = state.mesh().NN() - 1;
  const int surface_u =
      *dofs.index(Field::U, surface_element, surface_node);
  const int surface_p =
      *dofs.index(Field::P, surface_element, surface_node);
  const double potential_scale =
      model.LengthNorm() /
      (model.DensityNorm() * model.TimeNorm() * model.TimeNorm());
  const DegreeResult result = solveStaticDegree(model, state, degree);

  expectNear(result.h_u,
             solutions(surface_u, 0) / model.DensityNorm());
  expectNear(result.k_u, solutions(surface_p, 0) * potential_scale);
  expectNear(result.h_phi,
             solutions(surface_u, 1) / model.DensityNorm());
  expectNear(result.k_phi, solutions(surface_p, 1) * potential_scale);
  expectNear(result.h_t,
             solutions(surface_u, 2) * model.TimeNorm() *
                 model.TimeNorm() / model.LengthNorm());
  expectNear(result.k_t, solutions(surface_p, 2));

  const double gravity_si =
      state.surfaceGravity() * model.AccelerationNorm();
  const double reciprocity_left = gravity_si * result.h_phi;
  const double reciprocity_scale =
      std::max(std::abs(reciprocity_left), std::abs(result.k_u));
  ASSERT_GT(reciprocity_scale, 0.0);
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
  EXPECT_LE(residual, solve_tolerance);
  EXPECT_LE(linearity_error, solve_tolerance);
  EXPECT_LE(reciprocity_error, solve_tolerance);

  std::cout << std::setprecision(17) << name << " l=" << degree
            << " solve_residual=" << residual
            << " linearity_error=" << linearity_error
            << " reciprocity_error=" << reciprocity_error << '\n';
}

TEST(StaticSolveTests, OneFactorizationSolvesAllColumnsAndExtractsResults) {
  const std::vector<std::pair<std::string, std::string>> models{
      {"isotropic", DSPECM1D_LOVE_NUMBERS_ISOTROPIC_MODEL},
      {"TI", DSPECM1D_LOVE_NUMBERS_TI_MODEL},
      {"S-F-S", DSPECM1D_LOVE_NUMBERS_SFS_MODEL},
      {"PREM no-ocean", DSPECM1D_LOVE_NUMBERS_SOLID_SURFACE_MODEL},
  };

  for (const auto &[name, path] : models) {
    const Config config =
        name == "PREM no-ocean"
            ? Config{.maximum_degree = 10,
                     .polynomial_order = 3,
                     .maximum_radial_step = 0.2}
            : smallModelConfig();
    for (const int degree : {2, 10}) {
      checkSolve(path, name, config, degree);
    }
  }
}

TEST(StaticSolveTests, RejectsNegativeDegree) {
  const EarthModels::ModelInput<double> model(
      DSPECM1D_LOVE_NUMBERS_ISOTROPIC_MODEL);
  RadialState state(model, smallModelConfig());
  const LoveDofMap degree_zero(state.mesh(), 0);

  EXPECT_THROW(assembleStaticRightHandSides(state, degree_zero, -1),
               std::invalid_argument);
  EXPECT_THROW(solveStaticDegree(model, state, -1),
               std::invalid_argument);
}

}   // namespace
