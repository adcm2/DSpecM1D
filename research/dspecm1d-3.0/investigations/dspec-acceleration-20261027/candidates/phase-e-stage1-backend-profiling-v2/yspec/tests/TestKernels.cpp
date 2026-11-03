#include <gtest/gtest.h>

#include <YSpec/Kernels.hpp>
#include <array>
#include <cmath>
#include <vector>

#include "KernelCheck.hpp"

namespace {

using namespace YSpec;
using P = RadialParameter;
using KernelTest::Check;
using KernelTest::ExpectAgreement;
using KernelTest::kPrem;
using KernelTest::TheReceiver;
using KernelTest::TheSource;

// A few frequencies and low degrees: the identity between the kernels and the
// finite differences holds on any grid, so the quick check runs on a tiny
// one. The full band is the workflow test in TestKernelsWorkflow.cpp, which
// is slow.
Calculation Tiny() {
  auto c = Calculation{};
  c.degrees = {0, 8};
  c.band = {0.5e-3, 2.5e-3};
  c.sampling = {1200, 60};
  c.quantity = Quantity::Velocity;
  c.solver.tolerance = 1e-10;
  c.solver.startDecay = 25;
  return c;
}

const auto kWindow = Window{.start = 100, .end = 1100, .taper = 100};

TEST(RadialKernels, PredictFiniteDifferencesOnAFewFrequencies) {
  ExpectAgreement(Check(Gravity::Full, Along(Vertical), Tiny(), kWindow, 5e-5),
                  {
                      {"rho, mantle", {Property::Rho}, 5800e3, 300e3},
                      {"A, mantle", {Property::A}, 5800e3, 300e3},
                      {"C, mantle", {Property::C}, 5800e3, 300e3},
                      {"F, mantle", {Property::F}, 5800e3, 300e3},
                      {"L, mantle", {Property::L}, 5800e3, 300e3},
                      {"N, mantle", {Property::N}, 5800e3, 300e3},
                      {"qMu, mantle", {Property::QMu}, 5800e3, 300e3},
                      // The point parts at the source.
                      {"C, crust", {Property::C}, 6340e3, 10e3},
                      // The fluid, and the inner core through gravity alone.
                      {"kappa, outer core",
                       {Property::A, Property::C, Property::F},
                       2500e3,
                       500e3},
                      {"rho, inner core", {Property::Rho}, 600e3, 300e3},
                  });
}

TEST(RadialKernels, PredictFiniteDifferencesOnAFewFrequenciesWithoutGravity) {
  ExpectAgreement(
      Check(Gravity::None,
            TransverseDirection(TheSource().location, TheReceiver().location),
            Tiny(), kWindow),
      {
          {"rho, mantle", {Property::Rho}, 5800e3, 300e3},
          {"N, mantle", {Property::N}, 5800e3, 300e3},
      });
}

TEST(RadialKernels,
     PredictFiniteDifferencesOnAFewFrequenciesInTheCowlingApproximation) {
  ExpectAgreement(
      Check(Gravity::Cowling,
            RadialDirection(TheSource().location, TheReceiver().location),
            Tiny(), kWindow),
      {
          {"rho, mantle", {Property::Rho}, 5800e3, 300e3},
          {"rho, outer core", {Property::Rho}, 2500e3, 500e3},
      });
}

// The kernels of each set of parameters predict, for a perturbation of those
// parameters, what the model's own predict for the same perturbation: here
// with made-up kernels, at radii in the anisotropic upper mantle, the
// isotropic lower mantle and the fluid outer core.
TEST(RadialKernels, ReparameteriseConsistently) {
  const auto model = Model::ReadDeck(kPrem);
  const auto& scales = model.GetScales();
  auto radii = std::vector<double>{};
  for (auto r = 3000e3; r < 6300e3; r += 150e3) radii.push_back(r);
  auto kernels =
      RadialKernels(radii, std::vector<double>(radii.size(), 1e17), 6000e3);
  constexpr auto count = static_cast<std::size_t>(P::Count);
  for (auto p = std::size_t{0}; p < count; ++p) {
    const auto parameter = static_cast<P>(p);
    for (auto i = std::size_t{0}; i < radii.size(); ++i) {
      kernels[parameter][i] = std::sin(1.0 + p + 0.7 * i) * 1e-20;
    }
    kernels.AtSource(parameter) = std::cos(2.0 + p) * 1e-11;
  }
  // Small relative perturbations of each derived parameter at a radius.
  auto relative = [](std::size_t j, double r) {
    return 1e-3 * std::sin(0.3 * j + r / 1e6);
  };
  // The change they predict, through the model's own parameters and through
  // the set's.
  auto check = [&](ParameterSet set, auto toLove) {
    const auto named = Reparameterise(kernels, model, set, Scaling::Relative);
    auto direct = 0.0;
    auto derived = 0.0;
    auto add = [&](double r, double weight, auto kernel, auto derivedKernel) {
      const auto x = r / scales.Length();
      const auto m = model(model.LayerIndex(x), x);
      const auto love =
          toLove(m, [&](std::size_t j) { return relative(j, r); });
      for (auto p = std::size_t{0}; p < count; ++p) {
        direct += weight * kernel(static_cast<P>(p)) * love[p];
      }
      for (auto j = std::size_t{0}; j < named.names.size(); ++j) {
        derived += weight * derivedKernel(j) * relative(j, r);
      }
    };
    for (auto i = std::size_t{0}; i < radii.size(); ++i) {
      add(
          radii[i], kernels.Weights()[i], [&](P p) { return kernels[p][i]; },
          [&](std::size_t j) { return named.values[j][i]; });
    }
    add(
        kernels.SourceRadius(), 1.0, [&](P p) { return kernels.AtSource(p); },
        [&](std::size_t j) { return named.atSource[j]; });
    EXPECT_NEAR(derived, direct, 1e-10 * std::abs(direct));
  };
  const auto density = scales.Density();
  const auto modulus = scales.Modulus();
  // The changes of rho, A, C, F, L, N, qA, qF and qMu, in SI units.
  using Love = std::array<double, count>;
  // Isotropic changes of kappa and mu, and of qkappa and qmu, by Minos's
  // effective quality factors.
  auto isotropic = [&](const Material& m, double dRho, double dKappa,
                       double dMu, double dqKappa, double dqMu) {
    const auto x = 4 * m.mu / (3 * m.kappa + 4 * m.mu);
    return Love{dRho,
                dKappa + 4 * dMu / 3,
                dKappa + 4 * dMu / 3,
                dKappa - 2 * dMu / 3,
                dMu,
                dMu,
                (1 - x) * dqKappa + x * dqMu,
                ((1 - x) * dqKappa - x * dqMu / 2) / (1 - 1.5 * x),
                dqMu};
  };
  check(ParameterSet::Love, [&](const Material& m, auto d) {
    return Love{d(0) * m.rho * density,
                d(1) * m.A * modulus,
                d(2) * m.C * modulus,
                d(3) * m.F * modulus,
                d(4) * m.L * modulus,
                d(5) * m.N * modulus,
                d(6) * m.qA,
                d(7) * m.qF,
                d(8) * m.qMu};
  });
  check(ParameterSet::Moduli, [&](const Material& m, auto d) {
    return isotropic(m, d(0) * m.rho * density, d(1) * m.kappa * modulus,
                     d(2) * m.mu * modulus, d(3) * m.qKappa, d(4) * m.qMu);
  });
  check(ParameterSet::Velocities, [&](const Material& m, auto d) {
    // rho vp^2 = kappa + 4 mu / 3 and rho vs^2 = mu, to first order.
    const auto dMu = m.mu * (d(0) + 2 * d(2));
    const auto dP = (m.kappa + 4 * m.mu / 3) * (d(0) + 2 * d(1));
    return isotropic(m, d(0) * m.rho * density, (dP - 4 * dMu / 3) * modulus,
                     dMu * modulus, d(3) * m.qKappa, d(4) * m.qMu);
  });
  check(ParameterSet::Anisotropic, [&](const Material& m, auto d) {
    // rho, vpv, vph, vsv, vsh, eta.
    const auto dA = m.A * (d(0) + 2 * d(2));
    const auto dC = m.C * (d(0) + 2 * d(1));
    const auto dL = m.L * (d(0) + 2 * d(3));
    const auto dN = m.N * (d(0) + 2 * d(4));
    const auto eta = m.F / (m.A - 2 * m.L);
    const auto dF = eta * (dA - 2 * dL) + d(5) * eta * (m.A - 2 * m.L);
    auto love = isotropic(m, d(0) * m.rho * density, 0, 0, d(6) * m.qKappa,
                          d(7) * m.qMu);
    love[1] = dA * modulus;
    love[2] = dC * modulus;
    love[3] = dF * modulus;
    love[4] = dL * modulus;
    love[5] = dN * modulus;
    return love;
  });
}

}  // namespace
