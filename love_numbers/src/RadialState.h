#ifndef DSPECM1D_LOVE_NUMBERS_RADIAL_STATE_H
#define DSPECM1D_LOVE_NUMBERS_RADIAL_STATE_H

#include <stdexcept>
#include <vector>

#include <DSpecM1D/EarthMesh>
#include <DSpecM1D/LoveNumbers>
#include <DSpecM1D/src/model_info/MeshModel.h>
#include <GaussQuad/All>

namespace DSpecM1D::LoveNumbers::detail {

class RadialState {
public:
  RadialState(const EarthModels::ModelInput<double> &model,
              const Config &config)
      : mesh_(model, validatedNodeCount(model, config), 1.0,
              config.maximum_radial_step, false),
        mesh_model_(mesh_, model),
        density_derivatives_(
            mesh_.NE(), std::vector<double>(mesh_.NN(), 0.0)),
        gravity_(mesh_.NE(), std::vector<double>(mesh_.NN(), 0.0)),
        dimensionless_gravitational_constant_(
            6.67230e-11 / model.GravitationalConstant()) {
    for (int element = 0; element < mesh_.NE(); ++element) {
      const int layer = mesh_.LayerNumber(element);
      const auto density = model.Density(layer);
      for (int node = 0; node < mesh_.NN(); ++node) {
        const double radius = mesh_.NodeRadius(element, node);
        density_derivatives_[element][node] = density.Derivative(radius);
      }
    }
    sampleGravity(model);
  }

  const EarthMesh::RadialMesh &mesh() const noexcept { return mesh_; }

  const MeshModel &meshModel() const noexcept { return mesh_model_; }

  const GaussQuad::Quadrature1D<double> &quadrature() noexcept {
    return mesh_.GLL();
  }

  double dimensionlessGravitationalConstant() const noexcept {
    return dimensionless_gravitational_constant_;
  }

  double densityDerivative(int element, int node) const {
    return density_derivatives_[element][node];
  }

  double gravity(int element, int node) const {
    return gravity_[element][node];
  }

  double surfaceRadius() const {
    return mesh_.NodeRadius(mesh_.NE() - 1, mesh_.NN() - 1);
  }

  double surfaceGravity() const {
    return gravity(mesh_.NE() - 1, mesh_.NN() - 1);
  }

private:
  static double integrateDensityMoment(
      const EarthModels::ModelInput<double>::InterpA &density,
      double lower_radius, double upper_radius,
      const GaussQuad::Quadrature1D<double> &quadrature) {
    const double half_width = 0.5 * (upper_radius - lower_radius);
    const double centre = 0.5 * (upper_radius + lower_radius);
    double integral = 0.0;
    for (int point = 0; point < quadrature.N(); ++point) {
      const double radius =
          centre + half_width * quadrature.X(point);
      integral += quadrature.W(point) * density(radius) * radius * radius;
    }
    return half_width * integral;
  }

  void sampleGravity(const EarthModels::ModelInput<double> &model) {
    // Three Gauss--Legendre points integrate degree-five polynomials exactly.
    const auto quadrature =
        GaussQuad::GaussLegendreQuadrature1D<double>(3);
    constexpr double pi = 3.14159265358979323846;
    double enclosed_density_moment = 0.0;

    for (int element = 0; element < mesh_.NE(); ++element) {
      const int layer = mesh_.LayerNumber(element);
      const auto density = model.Density(layer);
      const auto knots = model.LayerRadii(layer);
      if (element > 0) {
        gravity_[element][0] =
            gravity_[element - 1][mesh_.NN() - 1];
      }

      for (int node = 1; node < mesh_.NN(); ++node) {
        double lower_radius = mesh_.NodeRadius(element, node - 1);
        const double upper_radius = mesh_.NodeRadius(element, node);
        for (const double knot : knots) {
          if (knot > lower_radius && knot < upper_radius) {
            enclosed_density_moment += integrateDensityMoment(
                density, lower_radius, knot, quadrature);
            lower_radius = knot;
          }
        }
        enclosed_density_moment += integrateDensityMoment(
            density, lower_radius, upper_radius, quadrature);
        gravity_[element][node] =
            4.0 * pi * dimensionless_gravitational_constant_ *
            enclosed_density_moment / (upper_radius * upper_radius);
      }
    }
  }

  static int validatedNodeCount(
      const EarthModels::ModelInput<double> &model, const Config &config) {
    if (config.maximum_degree < 0) {
      throw std::invalid_argument("maximum_degree must be non-negative.");
    }
    if (config.polynomial_order < 1) {
      throw std::invalid_argument("polynomial_order must be at least 1.");
    }
    if (!(config.maximum_radial_step > 0.0)) {
      throw std::invalid_argument("maximum_radial_step must be positive.");
    }
    if (model.IsFluid(model.NumberOfLayers() - 1)) {
      throw std::invalid_argument(
          "Surface fluids and oceans are not supported.");
    }

    return config.polynomial_order + 1;
  }

  EarthMesh::RadialMesh mesh_;
  MeshModel mesh_model_;
  std::vector<std::vector<double>> density_derivatives_;
  std::vector<std::vector<double>> gravity_;
  double dimensionless_gravitational_constant_;
};

}   // namespace DSpecM1D::LoveNumbers::detail

#endif   // DSPECM1D_LOVE_NUMBERS_RADIAL_STATE_H
