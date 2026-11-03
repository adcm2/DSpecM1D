#include <gtest/gtest.h>

#include <YSpec/Kernels.hpp>
#include <YSpec/Model.hpp>
#include <array>
#include <cmath>
#include <string>
#include <vector>

#include "Coefficients.hpp"
#include "KernelParts.hpp"
#include "KernelProducts.hpp"
#include "Spectra.hpp"

// The pieces of the kernel computation on inputs given by hand: the mesh in
// radius, the products of the forward and adjoint fields, the accumulation
// over frequency, and the weights of the sum over frequency. The full
// workflow, against finite differences, is in TestKernelsWorkflow.cpp (slow)
// and, on a few frequencies, in TestKernels.cpp.

namespace {

using namespace YSpec;

const auto kPrem = std::string(YSPEC_EXAMPLES_DIR) + "/prem.200";

// Coefficients with arbitrary values in every field at the orders of degree
// l, salted by a seed so that two fields differ.
Coefficients Arbitrary(int l, double seed) {
  auto c = Coefficients{};
  auto fill = [&, slot = 0](OrderArray<Complex>& field) mutable {
    ++slot;
    for (auto m : OrderArray<Complex>::Orders(l)) {
      field[m] = {std::sin(seed + slot + 0.7 * m),
                  std::cos(2 * seed - slot + 0.3 * m)};
    }
  };
  fill(c.u);
  fill(c.v);
  fill(c.w);
  fill(c.du);
  fill(c.dv);
  fill(c.dw);
  fill(c.phi);
  fill(c.dphi);
  return c;
}

std::array<Complex, 10> Members(const Products& p) {
  return {p.A,        p.C,        p.F,      p.L,
          p.N,        p.kinetic,  p.radial, p.radialWithPotential,
          p.coupling, p.potential};
}

TEST(Mesh, IntegratesExactly) {
  const auto model = Model::ReadDeck(kPrem);
  auto spans = 0.0;
  for (const auto& layer : model.Layers()) {
    spans += layer.Upper() - layer.Lower();
  }
  const auto surface = model.SurfaceRadius();
  for (auto breaks : {std::vector<double>{},
                      std::vector<double>{0.97 * surface, 0.5 * surface}}) {
    const auto mesh = Mesh(model, {.order = 5, .maxElement = 200e3}, breaks);
    auto total = 0.0;
    auto f = std::vector<double>{};
    for (auto i = std::size_t{0}; i < mesh.Size(); ++i) {
      total += mesh.Weight(i);
      const auto r = mesh.Radius(i);
      f.push_back(r * r);
    }
    EXPECT_NEAR(total, spans, 1e-12 * surface);
    // The integral of r^2 from each node to the surface, exactly.
    const auto above = mesh.IntegralAbove(f);
    for (auto i = std::size_t{0}; i < mesh.Size(); ++i) {
      const auto r = mesh.Radius(i);
      const auto exact = (surface * surface * surface - r * r * r) / 3;
      EXPECT_NEAR(above[i], exact, 1e-11 * exact) << "node " << i;
    }
  }
}

TEST(Products, VanishForARigidTranslationOfTheForwardField) {
  // u = grad(c r Y_1m) is a constant vector: no strain, so no elastic
  // products with any adjoint field, whatever the normalisations.
  const auto r = 0.9;
  auto a = Coefficients{};
  for (auto m : OrderArray<Complex>::Orders(1)) {
    a.u[m] = {0.3 + 0.1 * m, -0.2 + 0.05 * m};
    a.v[m] = a.u[m];
  }
  const auto b = Arbitrary(1, 1.0);
  auto p = Products{};
  Add(p, a, b, 1, r, true, true);
  for (auto product : {p.A, p.C, p.F, p.L, p.N}) {
    EXPECT_LT(std::abs(product), 1e-14);
  }
}

TEST(Products, VanishForARigidRotationOfTheForwardField) {
  // u = Omega x x is toroidal at l = 1 with W proportional to r.
  const auto r = 0.7;
  auto a = Coefficients{};
  for (auto m : OrderArray<Complex>::Orders(1)) {
    a.w[m] = {0.3 + 0.1 * m, -0.2 + 0.05 * m};
    a.dw[m] = a.w[m] / r;
  }
  const auto b = Arbitrary(1, 2.0);
  auto p = Products{};
  Add(p, a, b, 1, r, true, true);
  for (auto product : {p.A, p.C, p.F, p.L, p.N}) {
    EXPECT_LT(std::abs(product), 1e-14);
  }
}

TEST(Products, AreSymmetricInTheTwoFields) {
  // Every product pairs the fields through Pair, which is symmetric: swapping
  // the forward and adjoint fields changes nothing, which is reciprocity at
  // the level of the products.
  const auto r = 0.75;
  const auto a = Arbitrary(2, 3.0);
  const auto b = Arbitrary(2, 4.0);
  auto ab = Products{};
  auto ba = Products{};
  Add(ab, a, b, 2, r, true, true);
  Add(ba, b, a, 2, r, true, true);
  const auto those = Members(ba);
  const auto these = Members(ab);
  for (auto i = std::size_t{0}; i < these.size(); ++i) {
    EXPECT_LT(std::abs(these[i] - those[i]), 1e-14 * std::abs(these[i]))
        << "product " << i;
  }
}

TEST(Accumulate, FollowsItsContract) {
  const auto model = Model::ReadDeck(kPrem);
  const auto x = 0.9;  // in the mantle
  const auto layer = model.LayerIndex(x);
  const auto m = model(layer, x);
  const auto a = Arbitrary(2, 5.0);
  const auto b = Arbitrary(2, 6.0);
  auto p = Products{};
  Add(p, a, b, 2, x, true, true);
  const auto weight = 0.37;
  const auto omega2 = Complex{1.44, 0.1};
  auto at = [](RadialParameter q) { return static_cast<std::size_t>(q); };

  // Without dispersion the kernels of the quality factors vanish and each
  // modulus takes minus its product.
  auto sum = RadialAccumulator{};
  Accumulate(sum, p, m, x, weight, Complex{}, omega2, false);
  for (auto q :
       {RadialParameter::QA, RadialParameter::QF, RadialParameter::QMu}) {
    EXPECT_EQ(sum[at(q)], 0);
  }
  EXPECT_NEAR(sum[at(RadialParameter::A)], -weight * p.A.real(),
              1e-14 * std::abs(p.A));
  EXPECT_NEAR(sum[at(RadialParameter::L)], -weight * p.L.real(),
              1e-14 * std::abs(p.L));
  // Without gravity the density's kernel is the kinetic term less the
  // potential's, and the angular integrals for the adjoint potential stay
  // zero.
  EXPECT_NEAR(sum[at(RadialParameter::Rho)],
              weight * (omega2 * p.kinetic - p.potential).real(),
              1e-14 * std::abs(omega2 * p.kinetic));
  EXPECT_EQ(sum[RadialMomentSlot], 0);
  EXPECT_EQ(sum[FluxSlot], 0);

  // With gravity, what the density's kernel loses is exactly what the two
  // angular integrals record: g times the flux, and 4 pi G times M_rr.
  auto withGravity = RadialAccumulator{};
  Accumulate(withGravity, p, m, x, weight, Complex{}, omega2, true);
  const auto lost =
      sum[at(RadialParameter::Rho)] - withGravity[at(RadialParameter::Rho)];
  EXPECT_NEAR(lost,
              (m.g / m.rho) * withGravity[FluxSlot] +
                  4 * PiG * withGravity[RadialMomentSlot],
              1e-13 * std::abs(lost));
}

TEST(AdjointSpectrum, WeighsEachFrequencyBySpacing) {
  // The weights of the sum over frequency that replaces the integral: 2 df
  // everywhere but df at the Nyquist frequency and nothing above it.
  const auto model = Model::ReadDeck(kPrem);
  auto calculation = Calculation{};
  calculation.degrees = {0, 2};
  calculation.band = {0.5e-3, 2.5e-3};
  calculation.sampling = {1200, 60};
  const auto& scales = model.GetScales();
  const auto grid = MakeGrid(calculation, scales);
  const auto samples = SampleCount(grid, calculation.sampling);
  auto measurement = Measurement{};
  measurement.receiver = {{0, 0}, 0};
  measurement.direction = Along(Vertical);
  measurement.adjointSource.assign(samples, 1.0);
  const auto spectrum = MakeAdjointSpectrum(model, calculation, measurement);
  ASSERT_EQ(spectrum.weight.size(), grid.Count());
  const auto nyquist = grid.Samples() / 2;
  for (auto k = std::size_t{0}; k < grid.Count(); ++k) {
    const auto index = grid.First() + k;
    const auto expected = index > nyquist    ? 0.0
                          : index == nyquist ? grid.Spacing()
                                             : 2 * grid.Spacing();
    EXPECT_EQ(spectrum.weight[k], expected) << "frequency " << k;
  }
}

}  // namespace
