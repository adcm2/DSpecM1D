#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <numbers>

#include "Bessel.hpp"
#include "RungeKutta.hpp"

namespace {

using namespace YSpec;
using namespace std::complex_literals;

TEST(BesselRatio, MatchesTheStandardLibraryForRealArguments) {
  for (auto l : {0, 1, 2, 5, 20, 100}) {
    for (auto x : {1e-3, 0.1, 1.0, 3.7, 12.5, 40.0}) {
      const auto jl = std::sph_bessel(l, x);
      if (std::abs(jl) < 1e-200) continue;
      const auto expected = x * std::sph_bessel(l + 1, x) / jl;
      if (!std::isfinite(expected)) continue;  // the reference underflowed
      const auto q = BesselRatio(l, x * x);
      EXPECT_NEAR(q, expected, 1e-11 * std::max(1.0, std::abs(expected)))
          << "l = " << l << ", x = " << x;
    }
  }
}

TEST(BesselRatio, SatisfiesItsRecurrenceForComplexArguments) {
  for (auto l : {0, 3, 50, 600}) {
    for (auto z : {0.5 + 0.1i, -20.0 + 3.0i, 400.0 - 50.0i, -1e4 + 1.0i}) {
      const auto q = BesselRatio(l, z);
      const auto next = BesselRatio(l + 1, z);
      EXPECT_LT(std::abs(q - z / (2.0 * l + 3.0 - next)), 1e-12 * std::abs(q))
          << "l = " << l << ", z = " << z;
    }
  }
}

TEST(BesselRatio, ApproachesItsSmallArgumentLimit) {
  for (auto l : {0, 2, 10}) {
    const auto z = 1e-8 + 1e-9i;
    // The first two terms of the continued fraction.
    const auto expected = z / (2.0 * l + 3 - z / (2.0 * l + 5));
    EXPECT_LT(std::abs(BesselRatio(l, z) - expected), 1e-24);
  }
}

// y'' = -k^2 y as a first-order system, with a complex wavenumber so that
// the solution both oscillates and grows.
struct Oscillator {
  Complex k;
  Vector<2> operator()(double, const Vector<2>& y) const {
    return {y[1], -k * k * y[0]};
  }
  Vector<2> Exact(double r) const {
    return {std::cos(k * r), -k * std::sin(k * r)};
  }
};

template <class Tableau>
double FixedStepError(int steps) {
  const auto system = Oscillator{3.0 + 0.2i};
  auto y = ScaledVector<2>{{1.0, 0.0}};
  auto h = 2.0 / steps;
  RungeKutta<Tableau>({.adaptive = false}).Integrate(system, 0.0, 2.0, y, h);
  const auto exact = system.Exact(2.0);
  const auto scale = std::exp(y.logScale);
  return std::abs(y.value[0] * scale - exact[0]) / std::abs(exact[0]);
}

TEST(RungeKutta, ConvergesAtFifthOrderWithFixedSteps) {
  for (auto [coarse, fine] : {std::pair{FixedStepError<DormandPrince54>(40),
                                        FixedStepError<DormandPrince54>(80)},
                              std::pair{FixedStepError<CashKarp54>(40),
                                        FixedStepError<CashKarp54>(80)}}) {
    const auto order = std::log2(coarse / fine);
    EXPECT_GT(order, 4.7);
    EXPECT_LT(order, 5.5);
  }
}

TEST(RungeKutta, MeetsItsToleranceAdaptively) {
  const auto system = Oscillator{8.0 - 0.5i};
  for (auto tolerance : {1e-6, 1e-9}) {
    auto y = ScaledVector<2>{{1.0, 0.0}};
    auto h = 0.1;
    RungeKutta<>({.tolerance = tolerance}).Integrate(system, 0.0, 3.0, y, h);
    const auto exact = system.Exact(3.0);
    const auto scale = std::exp(y.logScale);
    const auto error = Norm(Vector<2>{y.value[0] * scale - exact[0],
                                      y.value[1] * scale - exact[1]});
    EXPECT_LT(error / Norm(exact), 100 * tolerance);
  }
}

TEST(RungeKutta, CarriesGrowthBeyondTheRangeOfADouble) {
  // y' = 2000 y over [0, 1] grows by e^2000; integrated backwards it decays.
  const auto grow = [](double, const Vector<1>& y) {
    return Vector<1>{2000.0 * y[0]};
  };
  auto y = ScaledVector<1>{{1.0}};
  auto h = 1e-3;
  RungeKutta<>({.tolerance = 1e-10}).Integrate(grow, 0.0, 1.0, y, h);
  EXPECT_NEAR(y.logScale + std::log(std::abs(y.value[0])), 2000.0, 1e-6);
  RungeKutta<>({.tolerance = 1e-10}).Integrate(grow, 1.0, 0.0, y, h);
  EXPECT_NEAR(y.logScale + std::log(std::abs(y.value[0])), 0.0, 1e-6);
}

TEST(RungeKutta, IntegratesSeveralSolutionsTogether) {
  // y' = lambda y with very different growth for three solutions of one
  // "system" (the rate depends on the component, so they differ), one of
  // them zero; each must come out as if integrated alone.
  const auto rates = Vector<2>{300.0 + 2.0i, -40.0 + 15.0i};
  const auto single = [&](double, const Vector<2>& y) {
    return Vector<2>{rates[0] * y[0], rates[1] * y[1]};
  };
  const auto together = [&](double r, const Vectors<2, 3>& y) {
    auto d = Vectors<2, 3>{};
    for (auto k = 0; k < 3; ++k) d[k] = single(r, y[k]);
    return d;
  };
  auto ys = std::array<ScaledVector<2>, 3>{ScaledVector<2>{{1.0, 0.0}},
                                           ScaledVector<2>{{1e-3, 1.0}},
                                           ScaledVector<2>{{0.0, 0.0}}};
  auto h = 1e-3;
  const auto rk = RungeKutta<>({.tolerance = 1e-10});
  rk.IntegrateTogether(together, 0.0, 1.0, ys, h);
  for (auto k = 0; k < 2; ++k) {
    auto alone =
        ScaledVector<2>{k == 0 ? Vector<2>{1.0, 0.0} : Vector<2>{1e-3, 1.0}};
    auto step = 1e-3;
    rk.Integrate(single, 0.0, 1.0, alone, step);
    for (auto i = 0; i < 2; ++i) {
      if (alone.value[i] == 0.0) {
        EXPECT_EQ(ys[k].value[i], 0.0);
        continue;
      }
      // Compare in logarithms: the first grows by e^300.
      const auto a = std::log(ys[k].value[i]) + ys[k].logScale;
      const auto b = std::log(alone.value[i]) + alone.logScale;
      EXPECT_LT(std::abs(a - b), 1e-6) << k << " " << i;
    }
  }
  EXPECT_EQ(Norm(ys[2].value), 0.0);
  EXPECT_EQ(ys[2].logScale, 0.0);
}

}  // namespace
