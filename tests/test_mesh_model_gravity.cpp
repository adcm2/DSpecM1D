#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include <DSpecM1D/EarthMesh>
#include <DSpecM1D/src/model_info/MeshModel.h>

namespace {

constexpr double kGravitationalConstant = 6.67230e-11;
constexpr double kPi = 3.14159265358979323846;

class AnalyticGravityModel {
public:
  using Field = std::function<double(double)>;

  AnalyticGravityModel(std::vector<double> layerBounds,
                       std::vector<Field> densities,
                       std::vector<std::vector<double>> densityKnots)
      : m_layerBounds(std::move(layerBounds)),
        m_densities(std::move(densities)),
        m_densityKnots(std::move(densityKnots)) {}

  double LengthNorm() const { return 1.0; }
  double MassNorm() const { return 1.0; }
  double TimeNorm() const { return 1.0; }
  double GravitationalConstant() const { return 1.0; }

  int NumberOfLayers() const {
    return static_cast<int>(m_densities.size());
  }
  double LowerRadius(int layer) const { return m_layerBounds[layer]; }
  double UpperRadius(int layer) const { return m_layerBounds[layer + 1]; }
  double OuterRadius() const { return m_layerBounds.back(); }
  bool IsSolid(int) const { return true; }
  bool IsFluid(int) const { return false; }

  const std::vector<double> &LayerRadii(int layer) const {
    return m_densityKnots[layer];
  }
  Field Density(int layer) const { return m_densities[layer]; }

  Field VPV(int) const { return constantField(2.0); }
  Field VSV(int) const { return constantField(1.0); }
  Field VPH(int) const { return constantField(2.0); }
  Field VSH(int) const { return constantField(1.0); }
  Field VP(int) const { return constantField(2.0); }
  Field VS(int) const { return constantField(1.0); }
  Field QKappa(int) const { return constantField(1000.0); }
  Field QMu(int) const { return constantField(1000.0); }
  Field Eta(int) const { return constantField(1.0); }
  Field A(int) const { return constantField(4.0); }
  Field N(int) const { return constantField(1.0); }
  Field L(int) const { return constantField(1.0); }
  Field C(int) const { return constantField(4.0); }
  Field F(int) const { return constantField(2.0); }
  Field Kappa(int) const { return constantField(2.0); }
  Field Mu(int) const { return constantField(1.0); }

private:
  static Field constantField(double value) {
    return [value](double) { return value; };
  }

  std::vector<double> m_layerBounds;
  std::vector<Field> m_densities;
  std::vector<std::vector<double>> m_densityKnots;
};

double gravityTolerance(double expected) {
  return 5e-13 * std::max(std::abs(expected), kGravitationalConstant);
}

std::string preciseString(double value) {
  std::ostringstream stream;
  stream << std::setprecision(17) << value;
  return stream.str();
}

}  // namespace

TEST(MeshModelGravityTests, ConstantDensityMatchesAnalyticProfile) {
  constexpr double density = 3.25;
  AnalyticGravityModel model(
      {0.0, 1.0}, {[=](double) { return density; }}, {{0.0, 1.0}});
  EarthMesh::RadialMesh mesh(model, 5, model.OuterRadius(), 1.0, false);
  MeshModel meshModel(mesh, model);

  EXPECT_DOUBLE_EQ(meshModel.Gravity(0, 0), 0.0);
  for (int element = 0; element < mesh.NE(); ++element) {
    for (int node = 0; node < mesh.NN(); ++node) {
      const double radius = mesh.NodeRadius(element, node);
      const double expected =
          4.0 * kPi * kGravitationalConstant * density * radius / 3.0;
      const double actual = meshModel.Gravity(element, node);
      EXPECT_TRUE(std::isfinite(actual));
      EXPECT_NEAR(actual, expected, gravityTolerance(expected));
    }
  }

  const double surfaceExpected =
      4.0 * kPi * kGravitationalConstant * density / 3.0;
  EXPECT_NEAR(meshModel.Gravity(mesh.NE() - 1, mesh.NN() - 1),
              surfaceExpected, gravityTolerance(surfaceExpected));
}

TEST(MeshModelGravityTests, LinearDensityMatchesAnalyticProfile) {
  constexpr double density0 = 2.0;
  constexpr double alpha = 3.0;
  AnalyticGravityModel model(
      {0.0, 1.0},
      {[=](double radius) { return density0 + alpha * radius; }},
      {{0.0, 0.37, 0.79, 1.0}});
  EarthMesh::RadialMesh mesh(model, 5, model.OuterRadius(), 1.0, false);
  MeshModel meshModel(mesh, model);

  double maxAbsoluteError = 0.0;
  for (int element = 0; element < mesh.NE(); ++element) {
    for (int node = 0; node < mesh.NN(); ++node) {
      const double radius = mesh.NodeRadius(element, node);
      const double expected =
          4.0 * kPi * kGravitationalConstant *
          (density0 * radius / 3.0 + alpha * radius * radius / 4.0);
      const double actual = meshModel.Gravity(element, node);
      maxAbsoluteError =
          std::max(maxAbsoluteError, std::abs(actual - expected));
      EXPECT_NEAR(actual, expected, gravityTolerance(expected));
    }
  }

  const double surfaceExpected =
      4.0 * kPi * kGravitationalConstant *
      (density0 / 3.0 + alpha / 4.0);
  const double surfaceError =
      std::abs(meshModel.Gravity(mesh.NE() - 1, mesh.NN() - 1) -
               surfaceExpected);
  testing::Test::RecordProperty("max_absolute_error",
                                preciseString(maxAbsoluteError));
  testing::Test::RecordProperty("surface_absolute_error",
                                preciseString(surfaceError));
  EXPECT_NEAR(meshModel.Gravity(mesh.NE() - 1, mesh.NN() - 1),
              surfaceExpected, gravityTolerance(surfaceExpected));
}

TEST(MeshModelGravityTests, DuplicatedInterfacePreservesEnclosedMass) {
  constexpr double interfaceRadius = 0.5;
  constexpr double innerDensity = 2.0;
  constexpr double outerDensity = 5.0;
  AnalyticGravityModel model(
      {0.0, interfaceRadius, 1.0},
      {[=](double) { return innerDensity; },
       [=](double) { return outerDensity; }},
      {{0.0, interfaceRadius}, {interfaceRadius, 1.0}});
  EarthMesh::RadialMesh mesh(model, 5, model.OuterRadius(), 0.5, false);
  MeshModel meshModel(mesh, model);

  ASSERT_EQ(mesh.NE(), 2);
  ASSERT_DOUBLE_EQ(mesh.NodeRadius(0, mesh.NN() - 1),
                   mesh.NodeRadius(1, 0));
  const double lowerInterfaceGravity =
      meshModel.Gravity(0, mesh.NN() - 1);
  const double upperInterfaceGravity = meshModel.Gravity(1, 0);
  EXPECT_NEAR(lowerInterfaceGravity, upperInterfaceGravity,
              gravityTolerance(lowerInterfaceGravity));
  EXPECT_DOUBLE_EQ(lowerInterfaceGravity, upperInterfaceGravity);

  const double enclosedMassIntegral =
      (innerDensity * std::pow(interfaceRadius, 3.0) +
       outerDensity * (1.0 - std::pow(interfaceRadius, 3.0))) /
      3.0;
  const double surfaceExpected =
      4.0 * kPi * kGravitationalConstant * enclosedMassIntegral;
  EXPECT_NEAR(meshModel.Gravity(1, mesh.NN() - 1), surfaceExpected,
              gravityTolerance(surfaceExpected));
}
