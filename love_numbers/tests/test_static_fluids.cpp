#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <Interpolation/Lagrange>

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

Config smallModelConfig() {
  return Config{
      .maximum_degree = 10,
      .polynomial_order = 4,
      .maximum_radial_step = 0.4,
  };
}

double relativeSymmetryError(const Eigen::SparseMatrix<double> &matrix) {
  const Eigen::SparseMatrix<double> difference =
      matrix - Eigen::SparseMatrix<double>(matrix.transpose());
  return difference.norm() / std::max(1.0, matrix.norm());
}

double potentialBlockError(RadialState &state, int degree,
                           const Eigen::SparseMatrix<double> &matrix) {
  const auto &mesh = state.mesh();
  const auto &quadrature = state.quadrature();
  const auto lagrange = Interpolation::LagrangePolynomial(
      quadrature.Points().begin(), quadrature.Points().end());
  const LoveDofMap map(mesh, degree);
  const double angular_degree = degree * (degree + 1);
  const double potential_scale =
      1.0 / (4.0 * 3.14159265358979323846 *
             state.dimensionlessGravitationalConstant());
  Eigen::MatrixXd expected =
      Eigen::MatrixXd::Zero(map.size(), map.size());

  for (int element = 0; element < mesh.NE(); ++element) {
    const double half_width = 0.5 * mesh.EW(element);
    const double derivative_scale = 1.0 / half_width;

    for (int point = 0; point < mesh.NN(); ++point) {
      const double radius = mesh.NodeRadius(element, point);
      const double weight = half_width * quadrature.W(point);
      double density_gradient_coefficient = 0.0;
      if (mesh.IsFluid(element) && radius != 0.0) {
        density_gradient_coefficient =
            radius * radius * state.densityDerivative(element, point) /
            state.gravity(element, point);
      }

      for (int row_node = 0; row_node < mesh.NN(); ++row_node) {
        const int row = *map.index(Field::P, element, row_node);
        const double row_value = row_node == point ? 1.0 : 0.0;
        const double row_derivative =
            derivative_scale *
            lagrange.Derivative(row_node, quadrature.X(point));

        for (int column_node = 0; column_node < mesh.NN();
             ++column_node) {
          const int column =
              *map.index(Field::P, element, column_node);
          const double column_value =
              column_node == point ? 1.0 : 0.0;
          const double column_derivative =
              derivative_scale *
              lagrange.Derivative(column_node, quadrature.X(point));
          double integrand =
              potential_scale *
              (radius * radius * row_derivative * column_derivative +
               angular_degree * row_value * column_value);
          integrand += density_gradient_coefficient * row_value *
                       column_value;
          expected(row, column) += weight * integrand;
        }
      }
    }
  }

  const int surface_element = mesh.NE() - 1;
  const int surface_node = mesh.NN() - 1;
  const int surface_p =
      *map.index(Field::P, surface_element, surface_node);
  expected(surface_p, surface_p) +=
      state.surfaceRadius() * (degree + 1.0) * potential_scale;

  std::vector<int> potential_indices;
  std::vector<bool> seen(map.size(), false);
  for (int element = 0; element < mesh.NE(); ++element) {
    for (int node = 0; node < mesh.NN(); ++node) {
      const int index = *map.index(Field::P, element, node);
      if (!seen[index]) {
        seen[index] = true;
        potential_indices.push_back(index);
      }
    }
  }

  double difference_squared = 0.0;
  double expected_squared = 0.0;
  for (const int row : potential_indices) {
    for (const int column : potential_indices) {
      const double difference =
          matrix.coeff(row, column) - expected(row, column);
      difference_squared += difference * difference;
      expected_squared += expected(row, column) * expected(row, column);
    }
  }
  return std::sqrt(difference_squared) /
         std::max(1.0, std::sqrt(expected_squared));
}

double densityGradientContributionNorm(RadialState &state) {
  const auto &mesh = state.mesh();
  const auto &quadrature = state.quadrature();
  double squared_norm = 0.0;

  for (int element = 0; element < mesh.NE(); ++element) {
    if (!mesh.IsFluid(element)) {
      continue;
    }
    for (int point = 0; point < mesh.NN(); ++point) {
      const double radius = mesh.NodeRadius(element, point);
      if (radius == 0.0) {
        continue;
      }
      const double contribution =
          0.5 * mesh.EW(element) * quadrature.W(point) * radius *
          radius * state.densityDerivative(element, point) /
          state.gravity(element, point);
      squared_norm += contribution * contribution;
    }
  }
  return std::sqrt(squared_norm);
}

double solidVolumeUpEntry(RadialState &state, int element, int node) {
  const auto &mesh = state.mesh();
  const auto &model = state.meshModel();
  const auto &quadrature = state.quadrature();
  const auto lagrange = Interpolation::LagrangePolynomial(
      quadrature.Points().begin(), quadrature.Points().end());
  const double radius = mesh.NodeRadius(element, node);
  const double derivative =
      2.0 / mesh.EW(element) *
      lagrange.Derivative(node, quadrature.X(node));
  return 0.5 * mesh.EW(element) * quadrature.W(node) *
         model.Density(element, node) * radius * radius * derivative;
}

void expectNoEmptyRowsOrColumns(
    const Eigen::SparseMatrix<double> &matrix) {
  std::vector<bool> nonempty_rows(matrix.rows(), false);
  std::vector<bool> nonempty_columns(matrix.cols(), false);
  for (int column = 0; column < matrix.outerSize(); ++column) {
    for (Eigen::SparseMatrix<double>::InnerIterator entry(matrix, column);
         entry; ++entry) {
      if (entry.value() != 0.0) {
        nonempty_rows[entry.row()] = true;
        nonempty_columns[entry.col()] = true;
      }
    }
  }
  EXPECT_TRUE(
      std::all_of(nonempty_rows.begin(), nonempty_rows.end(),
                  [](bool nonempty) { return nonempty; }));
  EXPECT_TRUE(
      std::all_of(nonempty_columns.begin(), nonempty_columns.end(),
                  [](bool nonempty) { return nonempty; }));
}

TEST(StaticFluidOperatorTests, InterfaceSignsUseFluidSideDensity) {
  const EarthModels::ModelInput<double> model(
      DSPECM1D_LOVE_NUMBERS_SFS_MODEL);
  RadialState state(model, smallModelConfig());
  const auto &mesh = state.mesh();
  const auto &mesh_model = state.meshModel();
  const LoveDofMap map(mesh, 2);
  const Eigen::SparseMatrix<double> matrix =
      assembleStaticOperator(state, 2);
  int interfaces_checked = 0;

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
    } else if (mesh.IsFluid(lower) && mesh.IsSolid(upper)) {
      orientation = 1;
      solid_element = upper;
      solid_node = 0;
      fluid_element = lower;
      fluid_node = mesh.NN() - 1;
    } else {
      continue;
    }

    const int u = *map.index(Field::U, solid_element, solid_node);
    const int p = *map.index(Field::P, lower, mesh.NN() - 1);
    const double radius = mesh.NodeRadius(solid_element, solid_node);
    const double fluid_density =
        mesh_model.Density(fluid_element, fluid_node);
    const double solid_density =
        mesh_model.Density(solid_element, solid_node);
    const double interface_entry =
        matrix.coeff(u, p) -
        solidVolumeUpEntry(state, solid_element, solid_node);
    const double sign_ratio =
        interface_entry / (fluid_density * radius * radius);
    const double recovered_density =
        interface_entry / (orientation * radius * radius);

    std::cout << std::setprecision(17)
              << (orientation < 0 ? "solid-fluid" : "fluid-solid")
              << " interface_sign_ratio=" << sign_ratio
              << " recovered_fluid_density=" << recovered_density << '\n';

    EXPECT_NEAR(sign_ratio, static_cast<double>(orientation),
                comparison_tolerance);
    EXPECT_NEAR(recovered_density, fluid_density,
                comparison_tolerance * std::max(1.0, fluid_density));
    EXPECT_GT(std::abs(recovered_density - solid_density),
              0.5 * std::abs(fluid_density - solid_density));
    ++interfaces_checked;
  }

  EXPECT_EQ(interfaces_checked, 2);
}

TEST(StaticFluidOperatorTests,
     NonconstantFluidBlockUsesPositiveDensityGradient) {
  const EarthModels::ModelInput<double> model(
      DSPECM1D_LOVE_NUMBERS_SFS_MODEL);
  RadialState state(model, smallModelConfig());
  const Eigen::SparseMatrix<double> matrix =
      assembleStaticOperator(state, 2);
  const double error = potentialBlockError(state, 2, matrix);
  const double density_gradient_norm =
      densityGradientContributionNorm(state);

  std::cout << std::setprecision(17)
            << "S-F-S density_gradient_norm=" << density_gradient_norm
            << " fluid_block_error=" << error << '\n';
  EXPECT_GT(density_gradient_norm, 100.0 * comparison_tolerance);
  EXPECT_LE(error, comparison_tolerance);
}

TEST(StaticFluidOperatorTests, ConstantDensityCentralFluidIsFinite) {
  const EarthModels::ModelInput<double> model(
      DSPECM1D_LOVE_NUMBERS_CENTRAL_FLUID_MODEL);
  RadialState state(model, smallModelConfig());
  const auto &mesh = state.mesh();

  ASSERT_TRUE(mesh.IsFluid(0));
  EXPECT_EQ(mesh.NodeRadius(0, 0), 0.0);
  EXPECT_EQ(state.gravity(0, 0), 0.0);
  for (int element = 0; element < mesh.NE(); ++element) {
    if (!mesh.IsFluid(element)) {
      continue;
    }
    for (int node = 0; node < mesh.NN(); ++node) {
      EXPECT_NEAR(state.densityDerivative(element, node), 0.0,
                  comparison_tolerance);
    }
  }

  const Eigen::SparseMatrix<double> matrix =
      assembleStaticOperator(state, 2);
  for (int column = 0; column < matrix.outerSize(); ++column) {
    for (Eigen::SparseMatrix<double>::InnerIterator entry(matrix, column);
         entry; ++entry) {
      EXPECT_TRUE(std::isfinite(entry.value()));
    }
  }

  const double error = potentialBlockError(state, 2, matrix);
  std::cout << std::setprecision(17)
            << "central constant-fluid block_error=" << error << '\n';
  EXPECT_LE(error, comparison_tolerance);
}

TEST(StaticFluidOperatorTests, InternalFluidMatricesAreCompleteAndSymmetric) {
  const std::vector<std::pair<std::string, std::string>> models{
      {"S-F-S", DSPECM1D_LOVE_NUMBERS_SFS_MODEL},
      {"central-fluid", DSPECM1D_LOVE_NUMBERS_CENTRAL_FLUID_MODEL},
      {"PREM no-ocean", DSPECM1D_LOVE_NUMBERS_SOLID_SURFACE_MODEL},
  };

  for (const auto &[name, path] : models) {
    const EarthModels::ModelInput<double> model(path);
    const Config config =
        name == "PREM no-ocean"
            ? Config{.maximum_degree = 2,
                     .polynomial_order = 3,
                     .maximum_radial_step = 0.2}
            : smallModelConfig();
    RadialState state(model, config);
    const LoveDofMap map(state.mesh(), 2);
    const Eigen::SparseMatrix<double> matrix =
        assembleStaticOperator(state, 2);
    const double symmetry_error = relativeSymmetryError(matrix);

    std::cout << std::setprecision(17) << name
              << " symmetry_error=" << symmetry_error << '\n';
    EXPECT_EQ(matrix.rows(), map.size());
    EXPECT_EQ(matrix.cols(), map.size());
    EXPECT_LE(symmetry_error, comparison_tolerance);
    expectNoEmptyRowsOrColumns(matrix);
  }
}

}   // namespace
