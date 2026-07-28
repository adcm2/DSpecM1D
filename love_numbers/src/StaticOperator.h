#ifndef DSPECM1D_LOVE_NUMBERS_STATIC_OPERATOR_H
#define DSPECM1D_LOVE_NUMBERS_STATIC_OPERATOR_H

#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

#include <Eigen/Sparse>
#include <Interpolation/Lagrange>

#include "LoveDofMap.h"
#include "RadialState.h"

namespace DSpecM1D::LoveNumbers::detail {

inline Eigen::SparseMatrix<double>
assembleStaticRadialOperator(RadialState &state,
                             const LoveDofMap &dof_map) {
  const auto &mesh = state.mesh();
  const auto &model = state.meshModel();
  const auto &quadrature = state.quadrature();
  const auto lagrange = Interpolation::LagrangePolynomial(
      quadrature.Points().begin(), quadrature.Points().end());
  const double gravitational_constant =
      state.dimensionlessGravitationalConstant();
  const double potential_scale =
      1.0 / (4.0 * 3.14159265358979323846 * gravitational_constant);

  using FormEntry = std::pair<int, double>;
  using LocalForm = std::vector<FormEntry>;
  std::vector<Eigen::Triplet<double>> triplets;

  auto add_outer_product =
      [&triplets](double coefficient, const LocalForm &left,
                  const LocalForm &right) {
        for (const auto &[row, left_value] : left) {
          for (const auto &[column, right_value] : right) {
            const double value = coefficient * left_value * right_value;
            if (value != 0.0) {
              triplets.emplace_back(row, column, value);
            }
          }
        }
      };

  for (int element = 0; element < mesh.NE(); ++element) {
    const double element_width = mesh.EW(element);
    const double derivative_scale = 2.0 / element_width;

    for (int point = 0; point < mesh.NN(); ++point) {
      const double radius = mesh.NodeRadius(element, point);
      const double integration_weight =
          0.5 * element_width * quadrature.W(point);
      const double density = model.Density(element, point);
      const double a = model.A(element, point);
      const double c = model.C(element, point);
      const double f = model.F(element, point);
      const double n = model.N(element, point);
      const double gravity = state.gravity(element, point);

      LocalForm displacement;
      LocalForm displacement_derivative;
      LocalForm potential_derivative;

      for (int node = 0; node < mesh.NN(); ++node) {
        const int u = *dof_map.index(Field::U, element, node);
        const int p = *dof_map.index(Field::P, element, node);
        const double value = node == point ? 1.0 : 0.0;
        const double derivative =
            derivative_scale *
            lagrange.Derivative(node, quadrature.X(point));
        displacement.emplace_back(u, value);
        displacement_derivative.emplace_back(u, derivative);
        potential_derivative.emplace_back(p, derivative);
      }

      add_outer_product(integration_weight * c * radius * radius,
                        displacement_derivative,
                        displacement_derivative);
      add_outer_product(integration_weight * 2.0 * f * radius,
                        displacement_derivative, displacement);
      add_outer_product(integration_weight * 2.0 * f * radius,
                        displacement, displacement_derivative);
      add_outer_product(
          integration_weight * 4.0 *
              (density *
                   (3.14159265358979323846 * gravitational_constant *
                        density * radius -
                    gravity) *
                   radius +
               a - n),
          displacement, displacement);
      add_outer_product(integration_weight * density * radius * radius,
                        displacement, potential_derivative);
      add_outer_product(integration_weight * density * radius * radius,
                        potential_derivative, displacement);
      add_outer_product(
          integration_weight * potential_scale * radius * radius,
          potential_derivative, potential_derivative);
    }
  }

  const int surface_element = mesh.NE() - 1;
  const int surface_node = mesh.NN() - 1;
  const int surface_p =
      *dof_map.index(Field::P, surface_element, surface_node);
  triplets.emplace_back(surface_p, surface_p,
                        state.surfaceRadius() * potential_scale);

  Eigen::SparseMatrix<double> matrix(dof_map.size(), dof_map.size());
  matrix.setFromTriplets(triplets.begin(), triplets.end());
  matrix.makeCompressed();
  return matrix;
}

inline Eigen::SparseMatrix<double>
assembleStaticOperator(RadialState &state, const LoveDofMap &dof_map,
                       int degree) {
  if (degree == 0) {
    return assembleStaticRadialOperator(state, dof_map);
  }
  if (degree < 0) {
    throw std::invalid_argument(
        "Static operator requires a non-negative degree.");
  }

  const auto &mesh = state.mesh();
  const auto &model = state.meshModel();
  const auto &quadrature = state.quadrature();
  const auto lagrange = Interpolation::LagrangePolynomial(
      quadrature.Points().begin(), quadrature.Points().end());
  const double angular_degree = degree * (degree + 1);
  const double gravitational_constant =
      state.dimensionlessGravitationalConstant();
  const double potential_scale =
      1.0 / (4.0 * 3.14159265358979323846 * gravitational_constant);

  using FormEntry = std::pair<int, double>;
  using LocalForm = std::vector<FormEntry>;
  std::vector<Eigen::Triplet<double>> triplets;

  auto add_outer_product =
      [&triplets](double coefficient, const LocalForm &left,
                  const LocalForm &right) {
        for (const auto &[row, left_value] : left) {
          for (const auto &[column, right_value] : right) {
            const double value = coefficient * left_value * right_value;
            if (value != 0.0) {
              triplets.emplace_back(row, column, value);
            }
          }
        }
      };

  for (int element = 0; element < mesh.NE(); ++element) {
    const double element_width = mesh.EW(element);
    const double derivative_scale = 2.0 / element_width;

    for (int point = 0; point < mesh.NN(); ++point) {
      const double radius = mesh.NodeRadius(element, point);
      const double integration_weight =
          0.5 * element_width * quadrature.W(point);

      LocalForm radial_derivative;
      LocalForm shear;
      LocalForm trace;
      LocalForm radial_displacement;
      LocalForm tangential_displacement;
      LocalForm potential;
      LocalForm potential_derivative;

      for (int node = 0; node < mesh.NN(); ++node) {
        const double value = node == point ? 1.0 : 0.0;
        const double derivative =
            derivative_scale *
            lagrange.Derivative(node, quadrature.X(point));

        if (const auto p = dof_map.index(Field::P, element, node)) {
          potential.emplace_back(*p, value);
          potential_derivative.emplace_back(*p, derivative);
        }

        if (mesh.IsFluid(element)) {
          continue;
        }

        const int u = *dof_map.index(Field::U, element, node);
        const int v = *dof_map.index(Field::V, element, node);
        radial_derivative.emplace_back(u, radius * derivative);
        shear.emplace_back(u, value);
        shear.emplace_back(v, radius * derivative - value);
        trace.emplace_back(u, 2.0 * value);
        trace.emplace_back(v, -angular_degree * value);
        radial_displacement.emplace_back(u, value);
        tangential_displacement.emplace_back(v, value);
      }

      add_outer_product(
          integration_weight * potential_scale * radius * radius,
          potential_derivative, potential_derivative);
      add_outer_product(
          integration_weight * potential_scale * angular_degree, potential,
          potential);

      if (mesh.IsFluid(element)) {
        double density_gradient_coefficient = 0.0;
        if (radius != 0.0) {
          density_gradient_coefficient =
              radius * radius * state.densityDerivative(element, point) /
              state.gravity(element, point);
          if (!std::isfinite(density_gradient_coefficient)) {
            throw std::runtime_error(
                "Fluid density-derivative coefficient is non-finite away "
                "from r = 0.");
          }
        }
        add_outer_product(
            integration_weight * density_gradient_coefficient, potential,
            potential);
        continue;
      }

      const double density = model.Density(element, point);
      const double gravity = state.gravity(element, point);
      const double a = model.A(element, point);
      const double c = model.C(element, point);
      const double f = model.F(element, point);
      const double l = model.L(element, point);
      const double n = model.N(element, point);

      add_outer_product(integration_weight * c, radial_derivative,
                        radial_derivative);
      add_outer_product(integration_weight * angular_degree * l, shear,
                        shear);
      add_outer_product(integration_weight * (a - n), trace, trace);
      add_outer_product(integration_weight * f, radial_derivative, trace);
      add_outer_product(integration_weight * f, trace, radial_derivative);
      add_outer_product(
          integration_weight * angular_degree * (angular_degree - 2.0) * n,
          tangential_displacement, tangential_displacement);

      add_outer_product(
          integration_weight * 4.0 * density *
              (3.14159265358979323846 * gravitational_constant * density *
                   radius -
               gravity) *
              radius,
          radial_displacement, radial_displacement);
      add_outer_product(
          integration_weight * angular_degree * density * gravity * radius,
          radial_displacement, tangential_displacement);
      add_outer_product(
          integration_weight * angular_degree * density * gravity * radius,
          tangential_displacement, radial_displacement);

      add_outer_product(integration_weight * density * radius * radius,
                        radial_displacement, potential_derivative);
      add_outer_product(integration_weight * density * radius * radius,
                        potential_derivative, radial_displacement);
      add_outer_product(
          integration_weight * angular_degree * density * radius,
          tangential_displacement, potential);
      add_outer_product(
          integration_weight * angular_degree * density * radius, potential,
          tangential_displacement);
    }
  }

  for (int lower_element = 0; lower_element + 1 < mesh.NE();
       ++lower_element) {
    const int upper_element = lower_element + 1;
    int orientation = 0;
    int solid_element = 0;
    int solid_node = 0;
    int fluid_element = 0;
    int fluid_node = 0;

    if (mesh.IsSolid(lower_element) && mesh.IsFluid(upper_element)) {
      orientation = -1;
      solid_element = lower_element;
      solid_node = mesh.NN() - 1;
      fluid_element = upper_element;
      fluid_node = 0;
    } else if (mesh.IsFluid(lower_element) &&
               mesh.IsSolid(upper_element)) {
      orientation = 1;
      solid_element = upper_element;
      solid_node = 0;
      fluid_element = lower_element;
      fluid_node = mesh.NN() - 1;
    } else {
      continue;
    }

    const int u =
        *dof_map.index(Field::U, solid_element, solid_node);
    const int p =
        *dof_map.index(Field::P, lower_element, mesh.NN() - 1);
    const double radius = mesh.NodeRadius(solid_element, solid_node);
    const double density = model.Density(fluid_element, fluid_node);
    const double gravity = state.gravity(fluid_element, fluid_node);
    const double coefficient =
        orientation * density * radius * radius;

    triplets.emplace_back(u, u, coefficient * gravity);
    triplets.emplace_back(u, p, coefficient);
    triplets.emplace_back(p, u, coefficient);
  }

  if (degree >= 2) {
    const int surface_element = mesh.NE() - 1;
    const int surface_node = mesh.NN() - 1;
    const int surface_p =
        *dof_map.index(Field::P, surface_element, surface_node);
    const double surface_radius = state.surfaceRadius();
    triplets.emplace_back(
        surface_p, surface_p,
        surface_radius * (degree + 1.0) * potential_scale);
  }

  Eigen::SparseMatrix<double> matrix(dof_map.size(), dof_map.size());
  matrix.setFromTriplets(triplets.begin(), triplets.end());
  matrix.makeCompressed();
  return matrix;
}

inline Eigen::SparseMatrix<double>
assembleStaticOperator(RadialState &state, int degree) {
  const LoveDofMap dof_map(state.mesh(), degree);
  return assembleStaticOperator(state, dof_map, degree);
}

}   // namespace DSpecM1D::LoveNumbers::detail

#endif   // DSPECM1D_LOVE_NUMBERS_STATIC_OPERATOR_H
