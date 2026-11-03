#include <gtest/gtest.h>

#include <YSpec/Receiver.hpp>
#include <YSpec/Source.hpp>
#include <array>
#include <cmath>
#include <complex>
#include <vector>

#include "Coefficients.hpp"
#include "Kinds.hpp"

// The strain sums on fields given by hand, where the strain is known in
// closed form: rigid motions, whose strain is zero whatever the angular
// functions evaluate to, and a uniform expansion, whose strain is isotropic.
// The full-workflow checks, against the gradient of the displacement and
// against reciprocity, are in TestStrainWorkflow.cpp, which is slow.

namespace {

using namespace YSpec;

const auto kSource = std::vector<Source>{{{10, 20}, 30e3, MomentTensor{}}};
const auto kReceivers = std::vector<Receiver>{
    {{40, 75}, 0}, {{-10, 130}, 0}, {{70, -20}, 0}, {{-55, 200}, 0}};

// An arbitrary coefficient at order m.
Complex Value(int m) { return {0.3 + 0.1 * m, -0.2 + 0.05 * m}; }

// The strain sums at one degree for each receiver, from coefficients given by
// hand at radius r.
std::vector<Strain::Sums> StrainSums(const Coefficients& coefficients, int l,
                                     double r) {
  const auto views = Views<SpinFunctions>(kSource, kReceivers, 2);
  auto sums = std::vector<Strain::Sums>(kReceivers.size());
  for (auto k = std::size_t{0}; k < kReceivers.size(); ++k) {
    Strain::Add(sums[k], coefficients, l, r, false, views, 0, k);
  }
  return sums;
}

// The size of the displacement the same coefficients give, which sets the
// scale that the strain is compared against.
double DisplacementScale(const Coefficients& coefficients, int l) {
  const auto views = Views<AngularFunctions>(kSource, kReceivers, 2);
  auto scale = 0.0;
  for (auto k = std::size_t{0}; k < kReceivers.size(); ++k) {
    auto sums = Displacement::Sums{};
    Displacement::Add(sums, coefficients, l, 1.0, false, views, 0, k);
    for (const auto& value : sums) scale = std::max(scale, std::abs(value));
  }
  return scale;
}

TEST(Strain, VanishesForARigidTranslation) {
  // u = grad(c r Y_1m) is a constant vector: U = V = c at l = 1, with zero
  // radial derivatives.
  const auto r = 0.9;
  auto c = Coefficients{};
  for (auto m : OrderArray<Complex>::Orders(1)) {
    c.u[m] = Value(m);
    c.v[m] = Value(m);
  }
  const auto scale = DisplacementScale(c, 1);
  ASSERT_GT(scale, 0);
  for (const auto& sums : StrainSums(c, 1, r)) {
    for (const auto& component : sums) {
      EXPECT_LT(std::abs(component), 1e-13 * scale / r);
    }
  }
}

TEST(Strain, VanishesForARigidRotation) {
  // u = Omega x x is toroidal at l = 1 with W proportional to r, so
  // W' = W / r.
  const auto r = 0.7;
  auto c = Coefficients{};
  for (auto m : OrderArray<Complex>::Orders(1)) {
    c.w[m] = Value(m);
    c.dw[m] = Value(m) / r;
  }
  const auto scale = DisplacementScale(c, 1);
  ASSERT_GT(scale, 0);
  for (const auto& sums : StrainSums(c, 1, r)) {
    for (const auto& component : sums) {
      EXPECT_LT(std::abs(component), 1e-13 * scale / r);
    }
  }
}

TEST(Strain, IsIsotropicForAUniformExpansion) {
  // u = x is radial at l = 0, U = r, and its strain is the identity: after
  // the turn to the receiver's frame, equal diagonal components and zero
  // shears, whatever X_00 is.
  const auto r = 0.8;
  auto c = Coefficients{};
  c.u[0] = r;
  c.du[0] = 1;
  const auto views = Views<SpinFunctions>(kSource, kReceivers, 2);
  const auto sums = StrainSums(c, 0, r);
  for (auto k = std::size_t{0}; k < kReceivers.size(); ++k) {
    const auto rotation = views.Of(0, k).geometry.rotation;
    const auto strain = Strain::ToReceiver(sums[k], rotation);
    const auto rr = strain[0];
    ASSERT_GT(std::abs(rr), 0);
    EXPECT_LT(std::abs(strain[1] - rr), 1e-14 * std::abs(rr));  // south-south
    EXPECT_LT(std::abs(strain[2] - rr), 1e-14 * std::abs(rr));  // east-east
    for (auto shear : {strain[3], strain[4], strain[5]}) {
      EXPECT_LT(std::abs(shear), 1e-14 * std::abs(rr));
    }
  }
}

TEST(Strain, MatchesTheDisplacementMachineryForAUniformExpansion) {
  // For u = x the radial displacement is r X_00 and the rr strain is X_00:
  // the two tables of angular functions must agree.
  const auto r = 0.8;
  auto c = Coefficients{};
  c.u[0] = r;
  c.du[0] = 1;
  const auto strains = StrainSums(c, 0, r);
  const auto views = Views<AngularFunctions>(kSource, kReceivers, 2);
  for (auto k = std::size_t{0}; k < kReceivers.size(); ++k) {
    auto u = Displacement::Sums{};
    Displacement::Add(u, c, 0, r, false, views, 0, k);
    EXPECT_LT(std::abs(strains[k][0] - u[0] / r), 1e-14 * std::abs(u[0] / r));
  }
}

TEST(Strain, TurnsTheTensorToTheReceiversFrame) {
  // Spin components built from a tensor given by hand, turned by a known
  // angle: the result is R t R^T with R the rotation that turns
  // (theta, phi) to (south, east).
  const auto rr = Complex{1.1, -0.3};
  const auto tt = Complex{0.4, 0.7};
  const auto pp = Complex{-0.8, 0.2};
  const auto rt = Complex{0.6, -0.5};
  const auto rp = Complex{-0.2, 0.9};
  const auto tp = Complex{0.3, 0.1};
  const auto spin = Strain::Sums{rr,
                                 tt + pp,
                                 rt - TimesI(rp),
                                 rt + TimesI(rp),
                                 (tt - pp) - 2.0 * TimesI(tp),
                                 (tt - pp) + 2.0 * TimesI(tp)};
  const auto rotation = 0.6;
  const auto strain = Strain::ToReceiver(spin, rotation);
  const auto c = std::cos(rotation);
  const auto s = std::sin(rotation);
  const auto expected =
      std::array<Complex, 6>{rr,
                             c * c * tt - 2 * c * s * tp + s * s * pp,
                             s * s * tt + 2 * c * s * tp + c * c * pp,
                             c * rt - s * rp,
                             s * rt + c * rp,
                             c * s * (tt - pp) + (c * c - s * s) * tp};
  for (auto i = std::size_t{0}; i < expected.size(); ++i) {
    EXPECT_LT(std::abs(strain[i] - expected[i]), 1e-14) << "component " << i;
  }
}

}  // namespace
