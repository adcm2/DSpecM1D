// A homogeneous solid sphere, whose toroidal modes are known in closed form:
// a check that needs no other code, and of a solid region reaching the centre,
// which the Fortran YSpec never meets.

#include <gtest/gtest.h>

#include <YSpec/RadialSolver.hpp>
#include <cmath>
#include <functional>
#include <numbers>
#include <string>

namespace {

using namespace YSpec;

constexpr auto kRadius = 6371e3;  // m
constexpr auto kVs = 4500.0;      // m/s

Model HomogeneousSphere() {
  // Six columns: r, rho, vp, vs, Qkappa, Qmu.
  const auto text = std::string(R"(homogeneous solid sphere
  0 1.0 1
  3 0 0
        0. 5000. 8000. 4500. 1000. 500.
  3185500. 5000. 8000. 4500. 1000. 500.
  6371000. 5000. 8000. 4500. 1000. 500.
)");
  return Model::ReadDeckText(text);
}

// The first root above x0 of (l - 1) j_l(x) - x j_{l+1}(x), the condition
// that the traction of a toroidal mode vanishes at the surface.
double ToroidalRoot(int l, double x0) {
  auto g = [l](double x) {
    return (l - 1) * std::sph_bessel(l, x) - x * std::sph_bessel(l + 1, x);
  };
  auto a = x0;
  auto b = a + 0.01;
  while (g(a) * g(b) > 0) {
    a = b;
    b += 0.01;
  }
  for (auto i = 0; i < 100; ++i) {
    const auto c = (a + b) / 2;
    (g(a) * g(c) <= 0 ? b : a) = c;
  }
  return (a + b) / 2;
}

double Peak(const std::function<double(double)>& response, double f0, double f1,
            int n) {
  auto best = 0;
  auto values = std::vector<double>(n + 1);
  for (auto i = 0; i <= n; ++i) {
    values[i] = response(f0 + (f1 - f0) * i / n);
    if (values[i] > values[best]) best = i;
  }
  const auto a = values[best - 1], b = values[best], c = values[best + 1];
  const auto shift = 0.5 * (a - c) / (a - 2 * b + c);
  return f0 + (f1 - f0) * (best + shift) / n;
}

TEST(HomogeneousSphere, HasItsToroidalModesWhereTheBesselFunctionsSay) {
  const auto model = HomogeneousSphere();
  ASSERT_EQ(model.Regions().size(), 1u);
  const auto scales = model.GetScales();
  const auto solver = RadialSolver(model, 1 - 100e3 / kRadius,
                                   1 - 10e3 / kRadius, Gravity::None, false);
  const auto toMHz = 1000 / scales.Time() / (2 * std::numbers::pi);
  for (auto l : {1, 2, 5}) {
    for (auto overtone : {0, 1}) {
      // l = 1 has a rigid rotation at zero frequency; skip it.
      const auto x0 = overtone == 0 ? 0.5 : ToroidalRoot(l, 0.5) + 0.1;
      const auto x = ToroidalRoot(l, x0);
      const auto expected = x * kVs / (2 * std::numbers::pi * kRadius) * 1e3;
      const auto omega = [&](double mHz) {
        return Complex(mHz / toMHz, -1e-4 * expected / toMHz);
      };
      const auto f = Peak(
          [&](double mHz) {
            return std::abs(solver.Toroidal(l, omega(mHz))(0, 0).W[0]);
          },
          0.97 * expected, 1.03 * expected, 300);
      EXPECT_NEAR(f / expected, 1, 1e-5)
          << "l = " << l << ", overtone " << overtone;
    }
  }
}

TEST(HomogeneousSphere, StartsLowDegreesAtTheCentre) {
  // At low degree nothing decays enough above the centre, so the toroidal
  // integration starts there.
  const auto model = HomogeneousSphere();
  const auto solver = RadialSolver(model, 0.9, 0.99, Gravity::Full);
  const auto omega = Complex(5.0, -0.1);
  EXPECT_EQ(solver.StartPoint(1, omega.real(), solver.Source(), 0).r, 0.0);
  for (auto l : {1, 2, 3}) {
    const auto t = solver.Toroidal(l, omega)(0, 0).W[0];
    EXPECT_TRUE(std::isfinite(std::abs(t)) && std::abs(t) > 0);
    const auto u = solver.Spheroidal(l, omega)(0, 0).U[0];
    EXPECT_TRUE(std::isfinite(std::abs(u)) && std::abs(u) > 0);
  }
}

TEST(HomogeneousSphere, AllowsASourceInTheInnerCoreOfPrem) {
  const auto model =
      Model::ReadDeck(std::string(YSPEC_EXAMPLES_DIR) + "/prem.200");
  const auto solver = RadialSolver(model, 0.1, 0.15, Gravity::Full);
  for (auto l : {1, 2, 10}) {
    const auto t = solver.Toroidal(l, Complex(10.0, -0.5))(0, 0);
    EXPECT_TRUE(std::isfinite(std::abs(t.W[1])));
  }
}

}  // namespace
