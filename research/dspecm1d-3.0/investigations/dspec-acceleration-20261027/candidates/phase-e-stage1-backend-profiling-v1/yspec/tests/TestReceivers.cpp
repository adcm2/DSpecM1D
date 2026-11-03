#include <gtest/gtest.h>

#include <YSpec/Location.hpp>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <vector>

#include "AngularFunctions.hpp"
#include "Geometry.hpp"

namespace {

using namespace YSpec;
constexpr auto pi = std::numbers::pi;

TEST(Geometry, PlacesAReceiverEastOfTheSourceOnTheEquator) {
  const auto g = Geometry({0, 30}, {0, 0});
  EXPECT_NEAR(g.distance, pi / 6, 1e-12);
  // The receiver lies due east (azimuth 90 degrees), so phi = 180 - 90.
  EXPECT_NEAR(g.phi, pi / 2, 1e-12);
  // The source lies due west of the receiver (azimuth 270 degrees).
  EXPECT_NEAR(g.rotation, pi / 2, 1e-12);
}

TEST(Geometry, PlacesAReceiverNorthOfTheSource) {
  const auto g = Geometry({40, 10}, {10, 10});
  EXPECT_NEAR(g.phi, pi, 1e-12);       // azimuth 0
  EXPECT_NEAR(g.rotation, pi, 1e-12);  // the source is due south
  // Geocentric colatitudes make the distance a little less than 30 degrees.
  EXPECT_NEAR(g.distance * 180 / pi, 29.9, 0.1);
}

TEST(Geometry, HandlesCoincidentAndAntipodalPoints) {
  const auto same = Geometry({12, 34}, {12, 34});
  EXPECT_EQ(same.distance, 0.0);
  // As for a receiver approaching the source from any direction.
  EXPECT_NEAR(same.rotation, same.phi, 1e-15);
  const auto near = Geometry({12.000001, 34.000001}, {12, 34});
  EXPECT_NEAR(near.rotation, near.phi, 1e-8);
  const auto antipode = Geometry({0, 180}, {0, 0});
  EXPECT_NEAR(antipode.distance, pi, 1e-12);
  EXPECT_TRUE(std::isfinite(antipode.phi));
  EXPECT_TRUE(std::isfinite(antipode.rotation));
}

TEST(GreatCircle, PlacesPointsAtTheDistancesAndAzimuthAsked) {
  // Measured as the seismograms measure them, from a source at the start.
  const auto source = Location{35, -120};
  const auto azimuth = 50.0;
  const auto points = GreatCircle(source, azimuth, 0, 180, 37);
  ASSERT_EQ(points.size(), 37u);
  EXPECT_EQ(points.front().latitude, source.latitude);
  EXPECT_EQ(points.front().longitude, source.longitude);
  for (auto i = std::size_t{1}; i + 1 < points.size(); ++i) {
    const auto g = Geometry(points[i], source);
    EXPECT_NEAR(g.distance * 180 / pi, 5.0 * static_cast<double>(i), 1e-9);
    // phi is pi less the azimuth of the receiver from the source.
    EXPECT_NEAR(g.phi, pi - azimuth * pi / 180, 1e-9);
  }
  EXPECT_NEAR(Geometry(points.back(), source).distance, pi, 1e-9);
}

TEST(GreatCircle, JoinsTwoPointsWithEqualSteps) {
  const auto from = Location{-30, 10};
  const auto to = Location{50, 140};
  const auto points = GreatCircle(from, to, 9);
  ASSERT_EQ(points.size(), 9u);
  EXPECT_EQ(points.front().latitude, from.latitude);
  EXPECT_EQ(points.back().longitude, to.longitude);
  const auto total = Geometry(to, from).distance;
  for (auto i = std::size_t{0}; i < points.size(); ++i) {
    // Each point lies on the great circle, at its share of the distance.
    const auto along = Geometry(points[i], from).distance;
    EXPECT_NEAR(along, total * static_cast<double>(i) / 8, 1e-9);
    EXPECT_NEAR(along + Geometry(to, points[i]).distance, total, 1e-9);
  }
  EXPECT_EQ(GreatCircle(from, to, 1).size(), 1u);
}

TEST(GreatCircle, RefusesWhatIsUndefined) {
  EXPECT_THROW((void)GreatCircle({0, 0}, {0, 180}, 5), std::invalid_argument);
  EXPECT_THROW((void)GreatCircle({0, 0}, {10, 10}, 0), std::invalid_argument);
  EXPECT_THROW((void)GreatCircle({0, 0}, 0, 0, 10, 0), std::invalid_argument);
}

TEST(AngularFunctions, AgreeWithTheStandardLibrary) {
  const auto lMax = 60;
  const auto theta = std::vector{0.3, 1.2, 2.0, 3.0};
  const auto angular = AngularFunctions(lMax, theta);
  const auto h = 1e-6;
  for (auto i = 0; i < std::ssize(theta); ++i) {
    const auto t = theta[i];
    for (auto l = 0; l <= lMax; ++l) {
      for (auto m = -std::min(l, 2); m <= std::min(l, 2); ++m) {
        const auto sign = (m < 0 && m % 2 != 0) ? -1.0 : 1.0;
        auto X = [&](double angle) {
          return sign * std::sph_legendre(l, std::abs(m), angle);
        };
        const auto& a = angular(i, l, m);
        const auto scale = std::sqrt((2.0 * l + 1) / (4 * pi));
        EXPECT_NEAR(a.value, X(t), 1e-12 * scale) << l << " " << m;
        const auto dx = (X(t + h) - X(t - h)) / (2 * h);
        EXPECT_NEAR(a.dTheta, dx, 1e-6 * scale * (l + 1)) << l << " " << m;
        EXPECT_NEAR(a.mOverSin, m * X(t) / std::sin(t), 1e-11 * scale * (l + 1))
            << l << " " << m;
      }
    }
  }
}

TEST(AngularFunctions, ReachBeyondTheLimitOfADoubleTable) {
  // Above GSHTrans's limit for double the table is built in long double.
  const auto lMax = 2500;
  const auto theta = std::vector{0.4, 1.9};
  const auto angular = AngularFunctions(lMax, theta);
  for (auto i : {0, 1}) {
    for (auto l : {10, 1000, 2000, 2500}) {
      for (auto m : {-2, 0, 1, 2}) {
        const auto sign = (m < 0 && m % 2 != 0) ? -1.0 : 1.0;
        const auto expected =
            sign * std::sph_legendre(l, std::abs(m), theta[i]);
        EXPECT_NEAR(angular(i, l, m).value, expected,
                    1e-10 * std::sqrt((2.0 * l + 1) / (4 * pi)))
            << l << " " << m;
      }
    }
  }
}

TEST(AngularFunctions, AreRegularAtThePoles) {
  const auto angular = AngularFunctions(10, {0.0, pi});
  for (auto i : {0, 1}) {
    for (auto l = 1; l <= 10; ++l) {
      // Only m = 0 survives at a pole, and m = +-1 in the derivatives.
      EXPECT_NEAR(angular(i, l, 2).value, 0, 1e-13);
      EXPECT_NEAR(angular(i, l, 2).mOverSin, 0, 1e-13);
      EXPECT_TRUE(std::isfinite(angular(i, l, 1).mOverSin));
      // At the north pole m X / sin(theta) -> -dX/dtheta for m = 1 ... in
      // magnitude both equal sqrt((2l+1)/4pi) sqrt(l(l+1))/2.
      const auto expected =
          std::sqrt((2.0 * l + 1) / (4 * pi)) * std::sqrt(l * (l + 1.0)) / 2;
      EXPECT_NEAR(std::abs(angular(i, l, 1).mOverSin), expected, 1e-12);
      EXPECT_NEAR(std::abs(angular(i, l, 1).dTheta), expected, 1e-12);
    }
  }
}

TEST(SpinFunctions, SatisfyTheirRecurrences) {
  // d/dtheta X^N = (a_N X^{N-1} - b_N X^{N+1}) / 2 and
  // (m - N cos theta) / sin theta X^N = -(a_N X^{N-1} + b_N X^{N+1}) / 2,
  // with a_N = sqrt((l+N)(l-N+1)) and b_N = sqrt((l-N)(l+N+1)).
  const auto theta = 0.7;
  const auto h = 1e-6;
  const auto x = SpinFunctions(12, {theta - h, theta, theta + h});
  for (auto l = 2; l <= 12; ++l) {
    for (auto m = -2; m <= 2; ++m) {
      for (auto n = -1; n <= 1; ++n) {
        const auto& at = x(1, l, m);
        const auto a = std::sqrt((l + n) * (l - n + 1.0));
        const auto b = std::sqrt((l - n) * (l + n + 1.0));
        const auto derivative =
            (x(2, l, m)[n + 2] - x(0, l, m)[n + 2]) / (2 * h);
        EXPECT_NEAR(derivative, (a * at[n + 1] - b * at[n + 3]) / 2, 1e-7)
            << l << " " << m << " " << n;
        EXPECT_NEAR((m - n * std::cos(theta)) / std::sin(theta) * at[n + 2],
                    -(a * at[n + 1] + b * at[n + 3]) / 2, 1e-12)
            << l << " " << m << " " << n;
      }
    }
  }
}

TEST(SpinFunctions, PairUnderIntegrationOverTheSphere) {
  // With Y^N_lm = X^N_lm(theta) exp(i m phi), the bilinear integral
  // int Y^N_lm Y^-N_l'm' dOmega is (-1)^(N+m) when l' = l and m' = -m, and
  // zero for l' != l.
  const auto count = 2001;
  auto theta = std::vector<double>{};
  for (auto i = 0; i < count; ++i) theta.push_back(pi * i / (count - 1));
  const auto x = SpinFunctions(8, theta);
  for (auto l = 2; l <= 8; ++l) {
    for (auto lp : {l, l - 1}) {
      for (auto m = -2; m <= 2; ++m) {
        for (auto n = -2; n <= 2; ++n) {
          // Simpson's rule in theta; the integral in phi is 2 pi.
          auto sum = 0.0;
          for (auto i = 0; i < count; ++i) {
            const auto weight = (i == 0 || i == count - 1) ? 1.0
                                : (i % 2 == 1)             ? 4.0
                                                           : 2.0;
            sum += weight * x(i, l, m)[n + 2] * x(i, lp, -m)[-n + 2] *
                   std::sin(theta[i]);
          }
          const auto integral = 2 * pi * sum * (pi / (count - 1)) / 3;
          const auto sign = ((n + m) % 2 == 0) ? 1.0 : -1.0;
          EXPECT_NEAR(integral, lp == l ? sign : 0.0, 1e-8)
              << l << " " << lp << " " << m << " " << n;
        }
      }
    }
  }
}

}  // namespace
