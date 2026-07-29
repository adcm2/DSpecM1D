#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "RadialState.h"

namespace {

using DSpecM1D::LoveNumbers::Config;
using DSpecM1D::LoveNumbers::detail::RadialState;

constexpr double gravity_tolerance = 2.0e-13;
constexpr double pi = 3.14159265358979323846;

Config gravityConfig() {
  return Config{
      .maximum_degree = 2,
      .polynomial_order = 3,
      .maximum_radial_step = 0.1,
  };
}

struct GravityErrors {
  double maximum_absolute_si = 0.0;
  double maximum_scaled = 0.0;
};

template <typename EnclosedDensityMoment>
GravityErrors measureGravityErrors(
    const EarthModels::ModelInput<double> &model,
    const RadialState &state,
    const EnclosedDensityMoment &enclosed_density_moment) {
  GravityErrors errors;
  const auto &mesh = state.mesh();
  const double gravitational_constant =
      state.dimensionlessGravitationalConstant();

  for (int element = 0; element < mesh.NE(); ++element) {
    for (int node = 0; node < mesh.NN(); ++node) {
      const double radius = mesh.NodeRadius(element, node);
      const double expected =
          radius == 0.0
              ? 0.0
              : 4.0 * pi * gravitational_constant *
                    enclosed_density_moment(radius) / (radius * radius);
      const double scale = std::max(1.0, std::abs(expected));
      const double error =
          std::abs(state.gravity(element, node) - expected);
      errors.maximum_absolute_si =
          std::max(errors.maximum_absolute_si,
                   error * model.AccelerationNorm());
      errors.maximum_scaled =
          std::max(errors.maximum_scaled, error / scale);

      EXPECT_TRUE(std::isfinite(state.gravity(element, node)));
      EXPECT_LE(error / scale, gravity_tolerance);
    }
  }
  return errors;
}

double exactSurfaceGravity(
    const EarthModels::ModelInput<double> &model,
    const RadialState &state) {
  const auto quadrature =
      GaussQuad::GaussLegendreQuadrature1D<double>(3);
  double density_moment = 0.0;
  for (int layer = 0; layer < model.NumberOfLayers(); ++layer) {
    const auto density = model.Density(layer);
    const auto radii = model.LayerRadii(layer);
    for (std::size_t interval = 1; interval < radii.size(); ++interval) {
      const double half_width =
          0.5 * (radii[interval] - radii[interval - 1]);
      const double centre =
          0.5 * (radii[interval] + radii[interval - 1]);
      for (int point = 0; point < quadrature.N(); ++point) {
        const double radius =
            centre + half_width * quadrature.X(point);
        density_moment +=
            half_width * quadrature.W(point) *
            density(radius) * radius * radius;
      }
    }
  }

  const double radius = state.surfaceRadius();
  return 4.0 * pi * state.dimensionlessGravitationalConstant() *
         density_moment / (radius * radius);
}

TEST(RadialStateTests, PolynomialOrderUsesOneMoreGllNode) {
  const EarthModels::ModelInput<double> model{
      DSPECM1D_LOVE_NUMBERS_SOLID_SURFACE_MODEL};
  const Config config{
      .maximum_degree = 2,
      .polynomial_order = 4,
      .maximum_radial_step = 0.2,
  };

  const RadialState state(model, config);

  EXPECT_EQ(state.mesh().NN(), config.polynomial_order + 1);
}

TEST(RadialStateTests, MeshEndsAtSurfaceWithoutExteriorBall) {
  const EarthModels::ModelInput<double> model{
      DSPECM1D_LOVE_NUMBERS_SOLID_SURFACE_MODEL};
  const Config config{
      .maximum_degree = 2,
      .polynomial_order = 3,
      .maximum_radial_step = 0.2,
  };

  const RadialState state(model, config);
  const auto &mesh = state.mesh();

  EXPECT_EQ(mesh.OuterPlanetaryElement(), mesh.NE() - 1);
  EXPECT_TRUE(mesh.HasFluid());
  EXPECT_DOUBLE_EQ(mesh.OR(), model.OuterRadius());
  EXPECT_DOUBLE_EQ(state.surfaceRadius(), model.OuterRadius());
}

TEST(RadialStateTests, SampledModelValuesAreFinite) {
  const EarthModels::ModelInput<double> model{
      DSPECM1D_LOVE_NUMBERS_SOLID_SURFACE_MODEL};
  const Config config{
      .maximum_degree = 2,
      .polynomial_order = 3,
      .maximum_radial_step = 0.2,
  };

  const RadialState state(model, config);
  const auto &mesh = state.mesh();
  const auto &mesh_model = state.meshModel();

  for (int element = 0; element < mesh.NE(); ++element) {
    for (int node = 0; node < mesh.NN(); ++node) {
      SCOPED_TRACE(testing::Message()
                   << "element " << element << ", node " << node);
      EXPECT_TRUE(std::isfinite(mesh_model.Density(element, node)));
      EXPECT_TRUE(std::isfinite(mesh_model.A(element, node)));
      EXPECT_TRUE(std::isfinite(mesh_model.C(element, node)));
      EXPECT_TRUE(std::isfinite(mesh_model.F(element, node)));
      EXPECT_TRUE(std::isfinite(mesh_model.L(element, node)));
      EXPECT_TRUE(std::isfinite(mesh_model.N(element, node)));
      EXPECT_TRUE(std::isfinite(state.gravity(element, node)));
    }
  }
  EXPECT_TRUE(std::isfinite(state.surfaceGravity()));
}

TEST(RadialStateTests, HomogeneousSphereGravityIsAnalyticAtEveryNode) {
  const EarthModels::ModelInput<double> model{
      DSPECM1D_LOVE_NUMBERS_HOMOGENEOUS_MODEL};
  const RadialState state(model, gravityConfig());
  const double density = model.Density(0)(0.0);
  const auto errors = measureGravityErrors(
      model, state, [density](double radius) {
        return density * radius * radius * radius / 3.0;
      });

  std::cout << std::setprecision(17)
            << "homogeneous_gravity_max_abs_si="
            << errors.maximum_absolute_si
            << " max_scaled=" << errors.maximum_scaled << '\n';
  EXPECT_EQ(state.gravity(0, 0), 0.0);
}

TEST(RadialStateTests, ControlledSolidFluidSolidGravityIsAnalytic) {
  const EarthModels::ModelInput<double> model{
      DSPECM1D_LOVE_NUMBERS_CONTROLLED_SFS_MODEL};
  const RadialState state(model, gravityConfig());
  const double first_interface = model.UpperRadius(0);
  const double second_interface = model.UpperRadius(1);
  const double inner_density = model.Density(0)(0.0);
  const double fluid_lower_density =
      model.Density(1)(first_interface);
  const double fluid_upper_density =
      model.Density(1)(second_interface);
  const double fluid_slope =
      (fluid_upper_density - fluid_lower_density) /
      (second_interface - first_interface);
  const double outer_density =
      model.Density(2)(second_interface);

  const auto enclosed_density_moment =
      [=](double radius) {
        const double inner_moment =
            inner_density * std::pow(first_interface, 3) / 3.0;
        if (radius <= first_interface) {
          return inner_density * std::pow(radius, 3) / 3.0;
        }

        const double capped_radius =
            std::min(radius, second_interface);
        const double cubic_difference =
            std::pow(capped_radius, 3) -
            std::pow(first_interface, 3);
        const double fluid_moment =
            fluid_lower_density * cubic_difference / 3.0 +
            fluid_slope *
                ((std::pow(capped_radius, 4) -
                  std::pow(first_interface, 4)) /
                     4.0 -
                 first_interface * cubic_difference / 3.0);
        if (radius <= second_interface) {
          return inner_moment + fluid_moment;
        }
        return inner_moment + fluid_moment +
               outer_density *
                   (std::pow(radius, 3) -
                    std::pow(second_interface, 3)) /
                   3.0;
      };
  const auto errors =
      measureGravityErrors(model, state, enclosed_density_moment);
  const double radius = state.surfaceRadius();
  const double expected_surface =
      4.0 * pi * state.dimensionlessGravitationalConstant() *
      enclosed_density_moment(radius) / (radius * radius);

  std::cout << std::setprecision(17)
            << "controlled_SFS_gravity_max_abs_si="
            << errors.maximum_absolute_si
            << " max_scaled=" << errors.maximum_scaled << '\n';
  EXPECT_NEAR(state.surfaceGravity(), expected_surface,
              gravity_tolerance * std::max(1.0,
                                           std::abs(expected_surface)));
}

TEST(RadialStateTests,
     MinimumPolynomialOrderExactlyIntegratesCubicDensityGravity) {
  const EarthModels::ModelInput<double> model{
      DSPECM1D_LOVE_NUMBERS_SFS_MODEL};
  const Config config{
      .maximum_degree = 2,
      .polynomial_order = 3,
      .maximum_radial_step = 0.125,
  };
  const RadialState state(model, config);
  const double expected = exactSurfaceGravity(model, state);
  const double error =
      std::abs(state.surfaceGravity() - expected) /
      std::max(1.0, std::abs(expected));

  std::cout << std::setprecision(17)
            << "cubic_density_p3_surface_gravity_error="
            << error << '\n';
  EXPECT_LE(error, 32.0 * std::numeric_limits<double>::epsilon());
}

TEST(RadialStateTests, ControlledCentralFluidGravityIsAnalytic) {
  const EarthModels::ModelInput<double> model{
      DSPECM1D_LOVE_NUMBERS_CONTROLLED_CENTRAL_FLUID_MODEL};
  const RadialState state(model, gravityConfig());
  const double interface_radius = model.UpperRadius(0);
  const double fluid_density = model.Density(0)(0.0);
  const double solid_density = model.Density(1)(interface_radius);
  const auto enclosed_density_moment =
      [=](double radius) {
        if (radius <= interface_radius) {
          return fluid_density * std::pow(radius, 3) / 3.0;
        }
        return fluid_density * std::pow(interface_radius, 3) / 3.0 +
               solid_density *
                   (std::pow(radius, 3) -
                    std::pow(interface_radius, 3)) /
                   3.0;
      };
  const auto errors =
      measureGravityErrors(model, state, enclosed_density_moment);
  const double radius = state.surfaceRadius();
  const double expected_surface =
      4.0 * pi * state.dimensionlessGravitationalConstant() *
      enclosed_density_moment(radius) / (radius * radius);

  std::cout << std::setprecision(17)
            << "controlled_central_gravity_max_abs_si="
            << errors.maximum_absolute_si
            << " max_scaled=" << errors.maximum_scaled << '\n';
  EXPECT_EQ(state.gravity(0, 0), 0.0);
  EXPECT_TRUE(std::isfinite(state.gravity(0, 0)));
  EXPECT_NEAR(state.surfaceGravity(), expected_surface,
              gravity_tolerance * std::max(1.0,
                                           std::abs(expected_surface)));
}

TEST(RadialStateTests, GravityIsContinuousAcrossEveryBoundary) {
  for (const char *path :
       {DSPECM1D_LOVE_NUMBERS_CONTROLLED_SFS_MODEL,
        DSPECM1D_LOVE_NUMBERS_CONTROLLED_CENTRAL_FLUID_MODEL}) {
    const EarthModels::ModelInput<double> model(path);
    const RadialState state(model, gravityConfig());
    const auto &mesh = state.mesh();
    int material_boundaries = 0;

    for (int element = 0; element < mesh.NE(); ++element) {
      for (int node = 0; node < mesh.NN(); ++node) {
        EXPECT_DOUBLE_EQ(
            state.gravity(element, node),
            state.meshModel().Gravity(element, node));
        EXPECT_TRUE(std::isfinite(state.gravity(element, node)));
      }
    }
    EXPECT_DOUBLE_EQ(state.gravity(0, 0), 0.0);
    EXPECT_DOUBLE_EQ(
        state.surfaceGravity(),
        state.meshModel().Gravity(mesh.NE() - 1, mesh.NN() - 1));

    for (int lower = 0; lower + 1 < mesh.NE(); ++lower) {
      const int upper = lower + 1;
      EXPECT_DOUBLE_EQ(state.gravity(lower, mesh.NN() - 1),
                       state.gravity(upper, 0));
      if (mesh.LayerNumber(lower) != mesh.LayerNumber(upper)) {
        ++material_boundaries;
      }
    }
    EXPECT_GT(material_boundaries, 0);
  }
}

TEST(RadialStateTests, DensityDerivativeUsesElementLayerSpline) {
  const EarthModels::ModelInput<double> model{
      DSPECM1D_LOVE_NUMBERS_SOLID_SURFACE_MODEL};
  const Config config{
      .maximum_degree = 2,
      .polynomial_order = 3,
      .maximum_radial_step = 0.2,
  };

  const RadialState state(model, config);
  const auto &mesh = state.mesh();

  for (int element = 0; element < mesh.NE(); ++element) {
    const int layer = mesh.LayerNumber(element);
    const auto density = model.Density(layer);
    for (int node = 0; node < mesh.NN(); ++node) {
      const double radius = mesh.NodeRadius(element, node);
      EXPECT_DOUBLE_EQ(state.densityDerivative(element, node),
                       density.Derivative(radius));
    }
  }
}

TEST(RadialStateTests, MaterialInterfaceUsesOneSidedDensityDerivatives) {
  const EarthModels::ModelInput<double> model{
      DSPECM1D_LOVE_NUMBERS_SOLID_SURFACE_MODEL};
  const Config config{
      .maximum_degree = 2,
      .polynomial_order = 3,
      .maximum_radial_step = 0.2,
  };

  const RadialState state(model, config);
  const auto &mesh = state.mesh();
  bool checked_interface = false;

  for (int lower_element = 0; lower_element + 1 < mesh.NE();
       ++lower_element) {
    const int upper_element = lower_element + 1;
    const int lower_layer = mesh.LayerNumber(lower_element);
    const int upper_layer = mesh.LayerNumber(upper_element);
    if (lower_layer == upper_layer) {
      continue;
    }

    const int lower_node = mesh.NN() - 1;
    const double radius = mesh.NodeRadius(lower_element, lower_node);
    const double upper_radius = mesh.NodeRadius(upper_element, 0);
    const double lower_derivative =
        model.Density(lower_layer).Derivative(radius);
    const double upper_derivative =
        model.Density(upper_layer).Derivative(upper_radius);
    if (lower_derivative == upper_derivative) {
      continue;
    }

    EXPECT_DOUBLE_EQ(radius, upper_radius);
    EXPECT_DOUBLE_EQ(
        state.densityDerivative(lower_element, lower_node),
        lower_derivative);
    EXPECT_DOUBLE_EQ(state.densityDerivative(upper_element, 0),
                     upper_derivative);
    checked_interface = true;
    break;
  }

  EXPECT_TRUE(checked_interface);
}

TEST(RadialStateTests, RejectsSurfaceFluid) {
  const EarthModels::ModelInput<double> model{
      DSPECM1D_LOVE_NUMBERS_OCEAN_MODEL};
  const Config config{
      .maximum_degree = 2,
      .polynomial_order = 3,
      .maximum_radial_step = 0.2,
  };

  try {
    static_cast<void>(RadialState(model, config));
    FAIL() << "RadialState accepted an outer fluid layer";
  } catch (const std::invalid_argument &error) {
    EXPECT_STREQ(error.what(),
                 "Surface fluids and oceans are not supported.");
  } catch (...) {
    FAIL() << "RadialState threw an unexpected exception type";
  }
}

TEST(RadialStateTests, RejectsInvalidConfiguration) {
  const EarthModels::ModelInput<double> model{
      DSPECM1D_LOVE_NUMBERS_SOLID_SURFACE_MODEL};

  EXPECT_THROW(
      RadialState(model,
                  Config{.maximum_degree = -1,
                         .polynomial_order = 3,
                         .maximum_radial_step = 0.2}),
      std::invalid_argument);
  EXPECT_THROW(
      RadialState(model,
                  Config{.maximum_degree = 2,
                         .polynomial_order = 3,
                         .maximum_radial_step = 0.0}),
      std::invalid_argument);
}

}   // namespace
