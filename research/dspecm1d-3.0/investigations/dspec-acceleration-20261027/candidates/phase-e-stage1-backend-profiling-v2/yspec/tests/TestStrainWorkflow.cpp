#include <gtest/gtest.h>

#include <YSpec/Synthetics.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <string>
#include <vector>

// The full workflow: strains computed through the whole band, checked against
// the gradient of the displacement and against reciprocity. Slow; the quick
// tests of the strain sums on fields given by hand are in TestStrain.cpp.

namespace {

using namespace YSpec;
using S = StrainComponent;

const auto kPrem = std::string(YSPEC_EXAMPLES_DIR) + "/prem.200";
constexpr auto kDegree = std::numbers::pi / 180;

Calculation Small() {
  auto c = Calculation{};
  c.degrees = {0, 40};
  c.band = {0.2e-3, 5e-3};
  c.sampling = {3600, 60};
  // Tightly, and from deep enough that where the integrations start does not
  // matter: the two sides of reciprocity start them at different depths.
  c.solver.tolerance = 1e-10;
  c.solver.startDecay = 25;
  return c;
}

double Peak(std::span<const double> trace) {
  return std::ranges::max(
      trace | std::views::transform([](double x) { return std::abs(x); }));
}

// The geographic latitude of a geocentric one, both in degrees, as the
// geometry converts them.
double Geographic(double geocentric) {
  return std::atan(std::tan(geocentric * kDegree) / 0.993277) / kDegree;
}

TEST(Strain, IsTheGradientOfTheDisplacement) {
  const auto model = Model::ReadDeck(kPrem);
  const auto calculation = Small();
  const auto sources = std::vector<Source>{
      {{10, 20},
       30e3,
       MomentTensor{.rr = 1e18, .rt = 2e18, .rp = -1e18, .tp = 5e17},
       Force{.r = 1e12, .t = -2e12, .p = 3e12}}};
  // Around a point at geocentric colatitude theta and longitude phi, in the
  // mantle: displacement there and at points a little away in each
  // direction, and strain there.
  const auto latitude = 40.0;  // geocentric
  const auto longitude = 75.0;
  const auto depth = 200e3;
  const auto radius = model.SurfaceRadius() * model.GetScales().Length();
  const auto r = radius - depth;
  const auto h = 1e-4;    // radians
  const auto dr = 500.0;  // m
  auto at = [&](double dLatitude, double dLongitude, double dDepth) {
    return Receiver{{Geographic(latitude + dLatitude / kDegree),
                     longitude + dLongitude / kDegree},
                    depth + dDepth};
  };
  const auto receivers = std::vector<Receiver>{
      at(0, 0, 0),  at(h, 0, 0),  at(-h, 0, 0), at(0, h, 0),
      at(0, -h, 0), at(0, 0, dr), at(0, 0, -dr)};
  const auto u = ComputeSeismograms(model, calculation, sources, receivers);
  const auto strain = ComputeStrains(model, calculation, sources,
                                     std::span{receivers}.first(1));

  const auto theta = (90 - latitude) * kDegree;
  // The displacement in (r, theta, phi) at receiver k and sample i.
  auto ur = [&](int k, std::size_t i) { return u(0, k, Vertical, i); };
  auto ut = [&](int k, std::size_t i) { return -u(0, k, North, i); };
  auto up = [&](int k, std::size_t i) { return u(0, k, East, i); };
  // Derivatives: theta increases southwards, and depth downwards.
  auto dTheta = [&](auto f, std::size_t i) {
    return (f(2, i) - f(1, i)) / (2 * h);
  };
  auto dPhi = [&](auto f, std::size_t i) {
    return (f(3, i) - f(4, i)) / (2 * h);
  };
  auto dR = [&](auto f, std::size_t i) {
    return (f(6, i) - f(5, i)) / (2 * dr);
  };

  auto peak = 0.0;
  for (auto c : {S::rr, S::tt, S::pp, S::rt, S::rp, S::tp}) {
    peak = std::max(peak, Peak(strain.Trace(0, 0, c)));
  }
  ASSERT_GT(peak, 0);
  const auto cot = std::cos(theta) / std::sin(theta);
  for (auto i = std::size_t{0}; i < strain.Samples(); i += 3) {
    const auto expected = std::array<double, 6>{
        dR(ur, i),
        (dTheta(ut, i) + ur(0, i)) / r,
        (dPhi(up, i) / std::sin(theta) + ut(0, i) * cot + ur(0, i)) / r,
        (dR(ut, i) - ut(0, i) / r + dTheta(ur, i) / r) / 2,
        (dR(up, i) - up(0, i) / r + dPhi(ur, i) / (r * std::sin(theta))) / 2,
        (dTheta(up, i) - up(0, i) * cot + dPhi(ut, i) / std::sin(theta)) /
            (2 * r)};
    for (auto c : {S::rr, S::tt, S::pp, S::rt, S::rp, S::tp}) {
      EXPECT_NEAR(strain(0, 0, c, i), expected[static_cast<int>(c)],
                  1e-4 * peak)
          << "component " << static_cast<int>(c) << ", sample " << i;
    }
  }
}

TEST(Strain, SatisfiesReciprocity) {
  // The displacement at a receiver from a moment tensor at a source is the
  // moment tensor contracted with the strain at the source from a unit force
  // at the receiver, along the component recorded.
  const auto model = Model::ReadDeck(kPrem);
  const auto calculation = Small();
  const auto moment = MomentTensor{
      .rr = 1e18, .rt = 2e18, .rp = -1e18, .tt = -3e17, .tp = 5e17, .pp = 7e17};
  const auto source = Source{{10, 20}, 30e3, moment};
  for (auto depth : {3e3, 1000e3}) {
    const auto receiver = Receiver{{-25, 80}, depth};
    const auto u = ComputeSeismograms(model, calculation, std::span{&source, 1},
                                      std::span{&receiver, 1});
    // Unit forces, scaled up to keep the numbers ordinary.
    const auto unit = 1e15;
    const auto forces =
        std::vector<Source>{{receiver.location, depth, {}, Force{.r = unit}},
                            {receiver.location, depth, {}, Force{.t = -unit}},
                            {receiver.location, depth, {}, Force{.p = unit}}};
    const auto atSource = Receiver{source.location, source.depth};
    const auto strain =
        ComputeStrains(model, calculation, forces, std::span{&atSource, 1});
    for (auto component : {Vertical, North, East}) {
      // The force along this component is the source of that index.
      const auto k = static_cast<std::size_t>(component);
      const auto peak = Peak(u.Trace(0, 0, component));
      ASSERT_GT(peak, 0);
      for (auto i = std::size_t{0}; i < u.Samples(); ++i) {
        auto e = [&](S c) { return strain(k, 0, c, i) / unit; };
        const auto contracted =
            moment.rr * e(S::rr) + moment.tt * e(S::tt) + moment.pp * e(S::pp) +
            2 * (moment.rt * e(S::rt) + moment.rp * e(S::rp) +
                 moment.tp * e(S::tp));
        EXPECT_NEAR(contracted, u(0, 0, component, i), 1e-7 * peak)
            << "depth " << depth << ", component " << k << ", sample " << i;
      }
    }
  }
}

}  // namespace
