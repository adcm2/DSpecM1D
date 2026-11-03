// The responses peak at the normal-mode frequencies of PREM: a check of the
// physics that is independent of the Fortran YSpec. The reference values are
// the well-known periods of the gravest modes of PREM (Dziewonski & Anderson
// 1981): 0S2 53.9 min, 0T2 43.9 min and 0S0 20.5 min.

#include <gtest/gtest.h>

#include <YSpec/RadialSolver.hpp>
#include <cmath>
#include <functional>
#include <numbers>
#include <string>

namespace {

using namespace YSpec;

// The frequency (mHz) in [f0, f1] at which |response(f)| is largest, refined
// by a parabola through the largest sample and its neighbours.
double Peak(const std::function<double(double)>& response, double f0, double f1,
            int n) {
  auto best = 0;
  auto values = std::vector<double>(n + 1);
  for (auto i = 0; i <= n; ++i) {
    values[i] = response(f0 + (f1 - f0) * i / n);
    if (values[i] > values[best]) best = i;
  }
  if (best == 0 || best == n) return f0 + (f1 - f0) * best / n;
  const auto a = values[best - 1], b = values[best], c = values[best + 1];
  const auto shift = 0.5 * (a - c) / (a - 2 * b + c);
  return f0 + (f1 - f0) * (best + shift) / n;
}

class NormalModes : public testing::Test {
 protected:
  Model model = Model::ReadDeck(std::string(YSPEC_EXAMPLES_DIR) + "/prem.200");
  double surface = model.SurfaceRadius();
  double scale = 1000 / model.GetScales().Time();  // mHz per unit frequency

  // A small damping, much less than the separation of the modes.
  Complex Omega(double mHz) const {
    return {2 * std::numbers::pi * mHz / scale,
            -2 * std::numbers::pi * 2e-4 / scale};
  }

  RadialSolver Solver(Gravity gravity) const {
    return {model, surface - 0.005, surface - 0.001, gravity};
  }
};

TEST_F(NormalModes, FindsTheFootballMode) {
  const auto solver = Solver(Gravity::Full);
  const auto f = Peak(
      [&](double mHz) {
        return std::abs(solver.Spheroidal(2, Omega(mHz))(0, 0).U[0]);
      },
      0.30, 0.32, 100);
  EXPECT_NEAR(f, 1000 / (53.9 * 60), 0.0005) << "0S2 at " << f << " mHz";
}

TEST_F(NormalModes, NeedsSelfGravitationForTheFootballMode) {
  // The Cowling approximation is poor at so low a degree: it moves 0S2 by
  // far more than the tolerance above.
  const auto solver = Solver(Gravity::Cowling);
  const auto f = Peak(
      [&](double mHz) {
        return std::abs(solver.Spheroidal(2, Omega(mHz))(0, 0).U[0]);
      },
      0.25, 0.45, 400);
  EXPECT_GT(std::abs(f - 1000 / (53.9 * 60)), 0.005)
      << "0S2 at " << f << " mHz";
}

TEST_F(NormalModes, FindsTheGravestToroidalMode) {
  const auto solver = Solver(Gravity::Full);
  const auto f = Peak(
      [&](double mHz) {
        return std::abs(solver.Toroidal(2, Omega(mHz))(0, 0).W[1]);
      },
      0.37, 0.39, 100);
  EXPECT_NEAR(f, 1000 / (43.9 * 60), 0.0005) << "0T2 at " << f << " mHz";
}

TEST_F(NormalModes, FindsTheBreathingMode) {
  const auto solver = Solver(Gravity::Full);
  const auto f = Peak(
      [&](double mHz) {
        return std::abs(solver.Radial(Omega(mHz))(0, 0).U[0]);
      },
      0.80, 0.83, 150);
  EXPECT_NEAR(f, 1000 / (20.5 * 60), 0.002) << "0S0 at " << f << " mHz";
}

}  // namespace
