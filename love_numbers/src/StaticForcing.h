#ifndef DSPECM1D_LOVE_NUMBERS_STATIC_FORCING_H
#define DSPECM1D_LOVE_NUMBERS_STATIC_FORCING_H

#include <cmath>
#include <stdexcept>

#include <Eigen/Dense>

#include "LoveDofMap.h"
#include "RadialState.h"

namespace DSpecM1D::LoveNumbers::detail {

inline Eigen::MatrixXd
assembleStaticRightHandSides(RadialState &state, const LoveDofMap &dofs,
                             int degree) {
  if (degree < 0) {
    throw std::invalid_argument(
        "Static forcing requires a non-negative degree.");
  }

  constexpr int displacement_column = 0;
  constexpr int gravitational_column = 1;
  constexpr int tide_column = 2;

  const auto &mesh = state.mesh();
  const auto &model = state.meshModel();
  const auto &quadrature = state.quadrature();
  const double surface_radius = state.surfaceRadius();
  Eigen::MatrixXd load_functionals =
      Eigen::MatrixXd::Zero(dofs.size(), 3);

  const int surface_element = mesh.NE() - 1;
  const int surface_node = mesh.NN() - 1;
  const int surface_u =
      *dofs.index(Field::U, surface_element, surface_node);
  load_functionals(surface_u, displacement_column) +=
      surface_radius * surface_radius * state.surfaceGravity();
  if (const auto surface_p =
          dofs.index(Field::P, surface_element, surface_node)) {
    load_functionals(*surface_p, gravitational_column) +=
        surface_radius * surface_radius;
  }

  if (degree == 0) {
    return -load_functionals;
  }

  const double angular_degree = degree * (degree + 1);
  for (int element = 0; element < mesh.NE(); ++element) {
    const double half_width = 0.5 * mesh.EW(element);
    for (int point = 0; point < mesh.NN(); ++point) {
      const double radius = mesh.NodeRadius(element, point);
      const double weight = half_width * quadrature.W(point);
      const double tide_potential =
          std::pow(radius / surface_radius, degree);

      if (mesh.IsSolid(element)) {
        const double coefficient =
            weight * model.Density(element, point) * radius *
            tide_potential;
        const int u = *dofs.index(Field::U, element, point);
        const int v = *dofs.index(Field::V, element, point);
        load_functionals(u, tide_column) += degree * coefficient;
        load_functionals(v, tide_column) +=
            angular_degree * coefficient;
        continue;
      }

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
      const int p = *dofs.index(Field::P, element, point);
      load_functionals(p, tide_column) +=
          weight * density_gradient_coefficient * tide_potential;
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

    const double radius = mesh.NodeRadius(solid_element, solid_node);
    const double density = model.Density(fluid_element, fluid_node);
    const double tide_potential =
        std::pow(radius / surface_radius, degree);
    const int u = *dofs.index(Field::U, solid_element, solid_node);
    load_functionals(u, tide_column) +=
        orientation * density * radius * radius * tide_potential;
  }

  return -load_functionals;
}

}   // namespace DSpecM1D::LoveNumbers::detail

#endif   // DSPECM1D_LOVE_NUMBERS_STATIC_FORCING_H
