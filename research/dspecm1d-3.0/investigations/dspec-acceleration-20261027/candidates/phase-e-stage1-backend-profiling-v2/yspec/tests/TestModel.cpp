#include <gtest/gtest.h>

#include <YSpec/Model.hpp>
#include <cmath>
#include <numbers>

namespace {

using namespace YSpec;

const auto kPrem = std::string(YSPEC_EXAMPLES_DIR) + "/prem.200";
const auto kPremFluidInnerCore =
    std::string(YSPEC_EXAMPLES_DIR) + "/prem.200.fic";
const auto kPremNoOcean = std::string(YSPEC_EXAMPLES_DIR) + "/prem.200.no";

TEST(Model, FindsTheLayersAndRegionsOfPrem) {
  const auto model = Model::ReadDeck(kPrem);
  EXPECT_EQ(model.NumberOfLayers(), 13);
  const auto regions = model.Regions();
  ASSERT_EQ(regions.size(), 4u);
  EXPECT_EQ(regions[0].phase, Phase::Solid);  // inner core
  EXPECT_EQ(regions[1].phase, Phase::Fluid);  // outer core
  EXPECT_EQ(regions[2].phase, Phase::Solid);  // mantle and crust
  EXPECT_EQ(regions[3].phase, Phase::Fluid);  // ocean
  EXPECT_EQ(regions[3].first, 12);
  EXPECT_DOUBLE_EQ(model.SurfaceRadius(), 1.0);
  const auto scales = model.GetScales();
  EXPECT_NEAR(regions[0].upper * scales.Length(), 1221.5e3, 1);
  EXPECT_NEAR(regions[1].upper * scales.Length(), 3480.0e3, 1);
  EXPECT_NEAR(regions[2].upper * scales.Length(), 6368.0e3, 1);
}

TEST(Model, SeesAFluidInnerCoreAndAMissingOcean) {
  const auto fic = Model::ReadDeck(kPremFluidInnerCore);
  ASSERT_EQ(fic.Regions().size(), 3u);
  EXPECT_EQ(fic.Regions()[0].phase, Phase::Fluid);
  const auto noOcean = Model::ReadDeck(kPremNoOcean);
  EXPECT_EQ(noOcean.Regions().back().phase, Phase::Solid);
}

TEST(Model, ReproducesTheKnotValues) {
  const auto model = Model::ReadDeck(kPrem);
  const auto scales = model.GetScales();
  // The first knot of the deck: r = 0, rho = 13088.5, vpv = 11262.2,
  // vsv = 3667.8, Qkappa = 1327.7, Qmu = 84.6.
  const auto m = model(0, 0.0);
  EXPECT_NEAR(m.rho * scales.Density(), 13088.5, 1e-9);
  const auto vs = std::sqrt(m.L / m.rho) * scales.Velocity();
  const auto vp = std::sqrt(m.C / m.rho) * scales.Velocity();
  EXPECT_NEAR(vs, 3667.8, 1e-9);
  EXPECT_NEAR(vp, 11262.2, 1e-9);
  EXPECT_NEAR(m.qMu, 1 / 84.6, 1e-15);
  EXPECT_NEAR(m.qKappa, 1 / 1327.7, 1e-15);
  // In the fluid outer core the shear moduli and their attenuation vanish.
  const auto c = model(1, 0.4);
  EXPECT_EQ(c.L, 0);
  EXPECT_EQ(c.mu, 0);
  EXPECT_EQ(c.qMu, 0);
  EXPECT_NEAR(c.qA, c.qKappa, 1e-15);
}

TEST(Model, GivesTheEarthsMassAndSurfaceGravity) {
  const auto model = Model::ReadDeck(kPrem);
  const auto scales = model.GetScales();
  const auto surface = model.SurfaceRadius();
  const auto top = model.NumberOfLayers() - 1;
  const auto g = model(top, surface).g * scales.Acceleration();
  EXPECT_NEAR(g, 9.815, 0.005);
  // g = G M / a^2, so M follows from g.
  const auto a = surface * scales.Length();
  const auto mass = g * a * a / GravitationalConstant;
  EXPECT_NEAR(mass / 5.972e24, 1.0, 2e-3);
  // Gravity is continuous across discontinuities and vanishes at the centre.
  EXPECT_EQ(model(0, 0.0).g, 0.0);
  for (auto i = Int{0}; i + 1 < model.NumberOfLayers(); ++i) {
    const auto r = model.GetLayer(i).Upper();
    EXPECT_NEAR(model(i, r).g, model(i + 1, r).g, 1e-14);
  }
}

TEST(Model, FindsTheLayerOfARadius) {
  const auto model = Model::ReadDeck(kPrem);
  const auto icb = model.Regions()[0].upper;
  EXPECT_EQ(model.LayerIndex(icb, Side::Below), 0);
  EXPECT_EQ(model.LayerIndex(icb, Side::Above), 1);
  EXPECT_EQ(model.LayerIndex(0.1), 0);
  EXPECT_EQ(model.LayerIndex(model.SurfaceRadius()), 12);
  EXPECT_EQ(model.LayerIndex(model.SurfaceRadius(), Side::Above), 12);
  EXPECT_THROW(model.LayerIndex(1.1), std::out_of_range);
}

}  // namespace
