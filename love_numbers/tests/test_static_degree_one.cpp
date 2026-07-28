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

#include <DSpecM1D/ModelInput>
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
constexpr double translation_tolerance = 2.0e-7;

Config smallModelConfig() {
  return Config{
      .maximum_degree = 2,
      .polynomial_order = 4,
      .maximum_radial_step = 0.1,
  };
}

Config premConfig() {
  return Config{
      .maximum_degree = 2,
      .polynomial_order = 2,
      .maximum_radial_step = 0.2,
  };
}

Eigen::SparseMatrix<double> removeRowAndColumn(
    const Eigen::SparseMatrix<double> &matrix, int removed) {
  std::vector<Eigen::Triplet<double>> entries;
  for (int outer = 0; outer < matrix.outerSize(); ++outer) {
    for (Eigen::SparseMatrix<double>::InnerIterator entry(matrix, outer);
         entry; ++entry) {
      if (entry.row() == removed || entry.col() == removed) {
        continue;
      }
      const int row = entry.row() - (entry.row() > removed ? 1 : 0);
      const int column =
          entry.col() - (entry.col() > removed ? 1 : 0);
      entries.emplace_back(row, column, entry.value());
    }
  }

  Eigen::SparseMatrix<double> result(matrix.rows() - 1,
                                     matrix.cols() - 1);
  result.setFromTriplets(entries.begin(), entries.end());
  return result;
}

double relativeError(const Eigen::VectorXd &actual,
                     const Eigen::VectorXd &expected) {
  return (actual - expected).norm() / std::max(1.0, expected.norm());
}

Eigen::VectorXd independentTidalRightHandSide(
    RadialState &state, const LoveDofMap &dofs,
    int *solid_to_fluid_interfaces = nullptr,
    int *fluid_to_solid_interfaces = nullptr) {
  const auto &mesh = state.mesh();
  const auto &model = state.meshModel();
  const auto &quadrature = state.quadrature();
  const double surface_radius = state.surfaceRadius();
  constexpr double angular_degree = 2.0;
  Eigen::VectorXd load_functional =
      Eigen::VectorXd::Zero(dofs.size());

  for (int element = 0; element < mesh.NE(); ++element) {
    const double half_width = 0.5 * mesh.EW(element);
    for (int point = 0; point < mesh.NN(); ++point) {
      const double radius = mesh.NodeRadius(element, point);
      const double weight = half_width * quadrature.W(point);
      const double tide_potential = radius / surface_radius;

      if (mesh.IsSolid(element)) {
        const double coefficient =
            weight * model.Density(element, point) * radius *
            tide_potential;
        const int u = *dofs.index(Field::U, element, point);
        const int v = *dofs.index(Field::V, element, point);
        load_functional(u) += coefficient;
        load_functional(v) += angular_degree * coefficient;
        continue;
      }

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
    const double tide_potential = radius / surface_radius;
    const int u = *dofs.index(Field::U, solid_element, solid_node);
    load_functional(u) +=
        orientation * density * radius * radius * tide_potential;
  }

  return -load_functional;
}

void compareDegreeOneMatrix(const std::string &path,
                            const std::string &name) {
  const EarthModels::ModelInput<double> model(path);
  const Config config = smallModelConfig();
  RadialState state(model, config);
  const LoveDofMap dofs(state.mesh(), 1);
  const Full1D::SEM sem(model, config.maximum_radial_step,
                        config.polynomial_order + 1,
                        config.maximum_degree);
  const Eigen::SparseMatrix<double> unconstrained = sem.hS(1);
  const int surface_p =
      sem.ltgS(2, sem.mesh().NE() - 1, sem.mesh().NN() - 1);
  const Eigen::SparseMatrix<double> reference =
      removeRowAndColumn(unconstrained, surface_p);
  const Eigen::SparseMatrix<double> matrix =
      assembleStaticOperator(state, dofs, 1);
  const Eigen::SparseMatrix<double> difference = matrix - reference;
  const Eigen::SparseMatrix<double> asymmetry =
      matrix - Eigen::SparseMatrix<double>(matrix.transpose());
  const double matrix_error =
      difference.norm() / std::max(1.0, reference.norm());
  const double symmetry_error =
      asymmetry.norm() / std::max(1.0, matrix.norm());

  Eigen::VectorXd translation =
      Eigen::VectorXd::Zero(unconstrained.rows());
  for (int element = 0; element < sem.mesh().NE(); ++element) {
    for (int node = 0; node < sem.mesh().NN(); ++node) {
      translation(sem.ltgS(0, element, node)) = 1.0;
      translation(sem.ltgS(1, element, node)) = 1.0;
      translation(sem.ltgS(2, element, node)) =
          -sem.meshModel().Gravity(element, node);
    }
  }
  const double translation_residual =
      (unconstrained * translation).norm() /
      std::max(1.0, unconstrained.norm() * translation.norm());

  std::cout << std::setprecision(17) << name
            << " l=1 matrix_error=" << matrix_error
            << " translation_residual=" << translation_residual
            << " symmetry_error=" << symmetry_error << '\n';

  EXPECT_EQ(matrix.rows(), dofs.size());
  EXPECT_EQ(matrix.rows(), reference.rows());
  EXPECT_EQ(matrix.cols(), reference.cols());
  EXPECT_LE(matrix_error, matrix_tolerance);
  EXPECT_LE(translation_residual, translation_tolerance);
  EXPECT_LE(symmetry_error, matrix_tolerance);
}

TEST(StaticDegreeOneOperatorTests,
     AllSolidMatricesMatchConstrainedSem) {
  compareDegreeOneMatrix(DSPECM1D_LOVE_NUMBERS_ISOTROPIC_MODEL,
                         "isotropic");
  compareDegreeOneMatrix(DSPECM1D_LOVE_NUMBERS_TI_MODEL, "TI");
}

TEST(StaticDegreeOneOperatorTests, RetainsCentreAndSurfaceFields) {
  const EarthModels::ModelInput<double> model(
      DSPECM1D_LOVE_NUMBERS_TI_MODEL);
  RadialState state(model, smallModelConfig());
  const LoveDofMap dofs(state.mesh(), 1);
  const int surface_element = state.mesh().NE() - 1;
  const int surface_node = state.mesh().NN() - 1;

  ASSERT_TRUE(dofs.index(Field::U, 0, 0));
  ASSERT_TRUE(dofs.index(Field::V, 0, 0));
  ASSERT_TRUE(dofs.index(Field::P, 0, 0));
  EXPECT_NE(*dofs.index(Field::U, 0, 0),
            *dofs.index(Field::V, 0, 0));
  EXPECT_NE(*dofs.index(Field::U, 0, 0),
            *dofs.index(Field::P, 0, 0));
  EXPECT_NE(*dofs.index(Field::V, 0, 0),
            *dofs.index(Field::P, 0, 0));
  EXPECT_TRUE(dofs.index(Field::U, surface_element, surface_node));
  EXPECT_TRUE(dofs.index(Field::V, surface_element, surface_node));
  EXPECT_FALSE(dofs.index(Field::P, surface_element, surface_node));
}

TEST(StaticDegreeOneOperatorTests,
     InternalFluidMatricesHaveExpectedSizeAndSymmetry) {
  const std::vector<std::pair<std::string, std::string>> models{
      {"S-F-S", DSPECM1D_LOVE_NUMBERS_SFS_MODEL},
      {"PREM no-ocean", DSPECM1D_LOVE_NUMBERS_SOLID_SURFACE_MODEL},
  };

  for (const auto &[name, path] : models) {
    const EarthModels::ModelInput<double> model(path);
    RadialState state(
        model, name == "PREM no-ocean" ? premConfig()
                                        : smallModelConfig());
    const LoveDofMap dofs(state.mesh(), 1);
    const Eigen::SparseMatrix<double> matrix =
        assembleStaticOperator(state, dofs, 1);
    const Eigen::SparseMatrix<double> asymmetry =
        matrix - Eigen::SparseMatrix<double>(matrix.transpose());
    const double symmetry_error =
        asymmetry.norm() / std::max(1.0, matrix.norm());

    std::cout << std::setprecision(17) << name
              << " l=1 symmetry_error=" << symmetry_error << '\n';
    EXPECT_EQ(matrix.rows(), dofs.size());
    EXPECT_EQ(matrix.cols(), dofs.size());
    EXPECT_LE(symmetry_error, matrix_tolerance);
  }
}

TEST(StaticDegreeOneForcingTests,
     SurfaceLoadsAndIndependentTidesAreCorrect) {
  const std::vector<std::pair<std::string, std::string>> models{
      {"TI", DSPECM1D_LOVE_NUMBERS_TI_MODEL},
      {"S-F-S", DSPECM1D_LOVE_NUMBERS_SFS_MODEL},
  };

  for (const auto &[name, path] : models) {
    const EarthModels::ModelInput<double> model(path);
    RadialState state(model, smallModelConfig());
    const LoveDofMap dofs(state.mesh(), 1);
    const Eigen::MatrixXd right_hand_sides =
        assembleStaticRightHandSides(state, dofs, 1);
    const int surface_element = state.mesh().NE() - 1;
    const int surface_node = state.mesh().NN() - 1;
    const int surface_u =
        *dofs.index(Field::U, surface_element, surface_node);
    const double radius = state.surfaceRadius();
    int solid_to_fluid_interfaces = 0;
    int fluid_to_solid_interfaces = 0;
    const Eigen::VectorXd expected_tide =
        independentTidalRightHandSide(
            state, dofs, &solid_to_fluid_interfaces,
            &fluid_to_solid_interfaces);
    const double tidal_error =
        relativeError(right_hand_sides.col(2), expected_tide);

    std::cout << std::setprecision(17) << name
              << " l=1 tidal_rhs_error=" << tidal_error << '\n';
    EXPECT_DOUBLE_EQ(
        right_hand_sides(surface_u, 0),
        -radius * radius * state.surfaceGravity());
    EXPECT_EQ(right_hand_sides.col(1).norm(), 0.0);
    EXPECT_LE(tidal_error, solve_tolerance);
    if (name == "S-F-S") {
      EXPECT_EQ(solid_to_fluid_interfaces, 1);
      EXPECT_EQ(fluid_to_solid_interfaces, 1);
    }
  }
}

void checkDegreeOneSolve(const std::string &path,
                         const std::string &name,
                         const Config &config) {
  const EarthModels::ModelInput<double> model(path);
  RadialState state(model, config);
  const LoveDofMap dofs(state.mesh(), 1);
  const Eigen::SparseMatrix<double> matrix =
      assembleStaticOperator(state, dofs, 1);
  const Eigen::MatrixXd right_hand_sides =
      assembleStaticRightHandSides(state, dofs, 1);

  Eigen::SparseLU<Eigen::SparseMatrix<double>> factorization;
  factorization.compute(matrix);
  ASSERT_EQ(factorization.info(), Eigen::Success);
  const Eigen::MatrixXd solutions =
      factorization.solve(right_hand_sides);
  ASSERT_EQ(factorization.info(), Eigen::Success);
  const double residual =
      (matrix * solutions - right_hand_sides).norm() /
      std::max(1.0, right_hand_sides.norm());
  const DegreeResult result = solveStaticDegree(model, state, 1);
  const std::array<double, 2> remaining_responses{
      result.h_u, result.h_t,
  };

  std::cout << std::setprecision(17) << name
            << " l=1 solve_residual=" << residual
            << " h_phi_abs=" << std::abs(result.h_phi)
            << " k_u=" << result.k_u
            << " k_phi=" << result.k_phi
            << " k_t=" << result.k_t << '\n';

  EXPECT_EQ(result.degree, 1);
  EXPECT_TRUE(std::all_of(
      remaining_responses.begin(), remaining_responses.end(),
      [](double value) { return std::isfinite(value); }));
  EXPECT_LE(residual, solve_tolerance);
  EXPECT_EQ(right_hand_sides.col(1).norm(), 0.0);
  EXPECT_EQ(solutions.col(1).norm(), 0.0);
  EXPECT_LE(std::abs(result.h_phi), 1.0e-14);
  EXPECT_EQ(result.k_u, 0.0);
  EXPECT_EQ(result.k_phi, 0.0);
  EXPECT_EQ(result.k_t, 0.0);
}

TEST(StaticDegreeOneSolveTests, SolvesSupportedModels) {
  const std::vector<std::pair<std::string, std::string>> models{
      {"isotropic", DSPECM1D_LOVE_NUMBERS_ISOTROPIC_MODEL},
      {"TI", DSPECM1D_LOVE_NUMBERS_TI_MODEL},
      {"S-F-S", DSPECM1D_LOVE_NUMBERS_SFS_MODEL},
      {"PREM no-ocean", DSPECM1D_LOVE_NUMBERS_SOLID_SURFACE_MODEL},
  };

  for (const auto &[name, path] : models) {
    checkDegreeOneSolve(
        path, name,
        name == "PREM no-ocean" ? premConfig() : smallModelConfig());
  }
}

}   // namespace
