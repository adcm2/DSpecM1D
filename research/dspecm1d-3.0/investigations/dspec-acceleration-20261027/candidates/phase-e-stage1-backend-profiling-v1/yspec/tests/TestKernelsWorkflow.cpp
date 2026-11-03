#include <gtest/gtest.h>

#include <YSpec/Kernels.hpp>

#include "KernelCheck.hpp"

// The full workflow: kernels of a traveltime over the whole band, checked
// against finite differences in every parameter and region. Slow; the quick
// counterpart on a few frequencies is in TestKernels.cpp.

namespace {

using namespace YSpec;
using KernelTest::Check;
using KernelTest::ExpectAgreement;
using KernelTest::TheReceiver;
using KernelTest::TheSource;

Calculation Small() {
  auto c = Calculation{};
  c.degrees = {0, 30};
  c.band = {0.2e-3, 4e-3};
  c.sampling = {3600, 60};
  c.taper = Taper{0.3e-3, 0.6e-3, 3e-3, 4e-3};
  c.quantity = Quantity::Velocity;
  c.solver.tolerance = 1e-10;
  c.solver.startDecay = 25;
  return c;
}

const auto kWindow = Window{.start = 700, .end = 1500, .taper = 200};

TEST(RadialKernels, PredictFiniteDifferences) {
  ExpectAgreement(Check(Gravity::Full, Along(Vertical), Small(), kWindow),
                  {
                      {"rho, mantle", {Property::Rho}, 5800e3, 300e3},
                      {"A, mantle", {Property::A}, 5800e3, 300e3},
                      {"C, mantle", {Property::C}, 5800e3, 300e3},
                      {"F, mantle", {Property::F}, 5800e3, 300e3},
                      {"L, mantle", {Property::L}, 5800e3, 300e3},
                      {"N, mantle", {Property::N}, 5800e3, 300e3},
                      {"qA, mantle", {Property::QA}, 5800e3, 300e3},
                      {"qF, mantle", {Property::QF}, 5800e3, 300e3},
                      {"qMu, mantle", {Property::QMu}, 5800e3, 300e3},
                      // Around the source, where the kernels of C, F and L have
                      // a point part, and clear of the ocean, where only
                      // isotropic perturbations have meaning.
                      {"rho, crust", {Property::Rho}, 6340e3, 10e3},
                      {"C, crust", {Property::C}, 6340e3, 10e3},
                      {"F, crust", {Property::F}, 6340e3, 10e3},
                      {"L, crust", {Property::L}, 6340e3, 10e3},
                      {"rho, outer core", {Property::Rho}, 2500e3, 500e3},
                      {"kappa, outer core",
                       {Property::A, Property::C, Property::F},
                       2500e3,
                       500e3},
                      // Reached only through the gravity above it.
                      {"rho, inner core", {Property::Rho}, 600e3, 300e3},
                  });
}

TEST(RadialKernels, PredictFiniteDifferencesInTheCowlingApproximation) {
  ExpectAgreement(
      Check(Gravity::Cowling,
            RadialDirection(TheSource().location, TheReceiver().location),
            Small(), kWindow),
      {
          {"rho, mantle", {Property::Rho}, 5800e3, 300e3},
          {"L, mantle", {Property::L}, 5800e3, 300e3},
          {"N, mantle", {Property::N}, 5800e3, 300e3},
          {"rho, outer core", {Property::Rho}, 2500e3, 500e3},
          {"rho, inner core", {Property::Rho}, 600e3, 300e3},
      });
}

TEST(RadialKernels, PredictFiniteDifferencesWithoutGravity) {
  ExpectAgreement(
      Check(Gravity::None,
            TransverseDirection(TheSource().location, TheReceiver().location),
            Small(), kWindow),
      {
          {"rho, mantle", {Property::Rho}, 5800e3, 300e3},
          {"L, crust", {Property::L}, 6340e3, 10e3},
          {"N, mantle", {Property::N}, 5800e3, 300e3},
          {"rho, outer core", {Property::Rho}, 2500e3, 500e3},
      });
}

}  // namespace
