// The responses of the radial solvers to unit jumps, compared with those of
// the Fortran YSpec (tests/data/solver, made by
// tests/fortran/make_solver_reference.sh with the Fortran's step density
// raised so that it is accurate to about 1e-8).

#include <gtest/gtest.h>

#include <YSpec/RadialSolver.hpp>
#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <map>
#include <numbers>
#include <ranges>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace {

using namespace YSpec;

struct Reference {
  std::string model;
  bool attenuation = true;
  Gravity gravity = Gravity::Full;
  double sourceDepth = 0;    // km
  double receiverDepth = 0;  // km
  double length = 0;         // min
  // Keyed by kind ('T', 'S', 'R'), degree and frequency (mHz); the values
  // in the order the driver prints them.
  std::map<std::tuple<char, int, double>, std::vector<Complex>> values;
};

Reference Read(const std::string& name) {
  auto file = std::ifstream(std::string(YSPEC_TEST_DATA_DIR) + "/solver/" +
                            name + ".txt");
  if (!file) throw std::runtime_error("cannot open reference " + name);
  auto reference = Reference{};
  auto line = std::string{};
  while (std::getline(file, line)) {
    auto in = std::istringstream(line);
    auto tag = std::string{};
    in >> tag;
    if (tag == "#") {
      auto key = std::string{};
      in >> key;
      if (key == "model") in >> reference.model;
      if (key == "attenuation") {
        auto a = 0;
        in >> a;
        reference.attenuation = a == 1;
      }
      if (key == "gravity") {
        auto g = 0;
        in >> g;
        reference.gravity = g == 0   ? Gravity::None
                            : g == 1 ? Gravity::Cowling
                                     : Gravity::Full;
      }
      if (key == "source_depth_km") in >> reference.sourceDepth;
      if (key == "receiver_depth_km") in >> reference.receiverDepth;
      if (key == "time_series_min") in >> reference.length;
      continue;
    }
    auto l = 0;
    auto f = 0.0;
    in >> l >> f;
    if (tag == "S") {
      auto j = 0;
      in >> j;
    }
    auto& values = reference.values[{tag[0], l, f}];
    for (auto re = 0.0, im = 0.0; in >> re >> im;) values.emplace_back(re, im);
  }
  return reference;
}

class SolverAgainstFortran : public testing::TestWithParam<std::string> {};

TEST_P(SolverAgainstFortran, AgreesWithTheFortran) {
  const auto reference = Read(GetParam());
  const auto model =
      Model::ReadDeck(std::string(YSPEC_EXAMPLES_DIR) + "/" + reference.model);
  const auto& scales = model.GetScales();
  const auto surface = model.SurfaceRadius();
  const auto rs = surface - reference.sourceDepth * 1e3 / scales.Length();
  const auto rr = surface - reference.receiverDepth * 1e3 / scales.Length();
  // The imaginary part of the frequency, as in the original YSpec.
  const auto epsilon = 5 / (reference.length * 60 / scales.Time());
  auto options = SolverOptions{};
  options.tolerance = 1e-10;
  const auto solver = RadialSolver(model, rs, rr, reference.gravity,
                                   reference.attenuation, options);
  // The Fortran prints its radially scaled variables: r U, zeta r V, zeta r W
  // and r phi.
  const auto r = solver.Receiver().r;

  for (const auto& [key, expected] : reference.values) {
    const auto [kind, l, f] = key;
    const auto omega =
        Complex(2 * std::numbers::pi * f / 1000 * scales.Time(), -epsilon);
    const auto zeta = std::sqrt(l * (l + 1.0));
    auto computed = std::vector<Complex>{};
    if (kind == 'R') {
      for (auto u : solver.Radial(omega)(0, 0).U) computed.push_back(r * u);
    } else if (kind == 'T') {
      for (auto w : solver.Toroidal(l, omega)(0, 0).W) {
        computed.push_back(zeta * r * w);
      }
    } else {
      const auto y = solver.Spheroidal(l, omega)(0, 0);
      // For each jump, U, V and the potential.
      for (auto j : {0, 1, 2, 3}) {
        computed.push_back(r * y.U[j]);
        computed.push_back(zeta * r * y.V[j]);
        computed.push_back(r * y.Phi[j]);
      }
    }
    ASSERT_EQ(computed.size(), expected.size());
    auto scale = 0.0;
    for (const auto& v : expected) scale = std::max(scale, std::abs(v));
    for (auto i = std::size_t{0}; i < expected.size(); ++i) {
      EXPECT_LT(std::abs(computed[i] - expected[i]), 1e-5 * scale)
          << kind << " l = " << l << " f = " << f << " mHz, entry " << i << ": "
          << computed[i] << " against " << expected[i];
    }
  }
}

INSTANTIATE_TEST_SUITE_P(Cases, SolverAgainstFortran,
                         testing::Values("full", "full_below", "cowling",
                                         "cowling_below", "none", "none_below",
                                         "fic_full", "fic_cowling", "fic_none",
                                         "no_ocean_full", "no_ocean_none"));

TEST(Solver, GivesForManyRadiiWhatItGivesForEach) {
  // Sources and receivers at several depths, above and below one another,
  // solved together and one pair at a time.
  const auto model =
      Model::ReadDeck(std::string(YSPEC_EXAMPLES_DIR) + "/prem.200");
  const auto length = model.GetScales().Length();
  const auto surface = model.SurfaceRadius();
  auto radii = [&](std::initializer_list<double> depths) {
    auto r = std::vector<double>{};
    for (auto d : depths) r.push_back(surface - d * 1e3 / length);
    return r;
  };
  const auto sources = radii({29.5, 300, 600});
  const auto receivers = radii({3, 50, 300, 450, 700});
  auto options = SolverOptions{};
  options.tolerance = 1e-10;
  const auto together =
      RadialSolver(model, sources, receivers, Gravity::Full, true, options);
  auto close = [](Complex a, Complex b, double scale) {
    return std::abs(a - b) <= 1e-6 * scale;
  };
  for (auto [l, omega] : {std::pair{0, Complex(10.0, -0.5)},
                          {2, Complex(5.0, -0.5)},
                          {40, Complex(30.0, -1.0)}}) {
    for (auto s : std::views::iota(std::size_t{0}, sources.size())) {
      for (auto k : std::views::iota(std::size_t{0}, receivers.size())) {
        const auto alone = RadialSolver(model, sources[s], receivers[k],
                                        Gravity::Full, true, options);
        if (l == 0) {
          const auto a = together.Radial(omega)(s, k).U;
          const auto b = alone.Radial(omega)(0, 0).U;
          const auto scale = std::max(std::abs(b[0]), std::abs(b[1]));
          for (auto j : {0, 1}) EXPECT_TRUE(close(a[j], b[j], scale));
          continue;
        }
        const auto ta = together.Toroidal(l, omega)(s, k).W;
        const auto tb = alone.Toroidal(l, omega)(0, 0).W;
        const auto tScale = std::max(std::abs(tb[0]), std::abs(tb[1]));
        for (auto j : {0, 1}) {
          EXPECT_TRUE(close(ta[j], tb[j], tScale))
              << "toroidal l " << l << " source " << s << " receiver " << k;
        }
        const auto sa = together.Spheroidal(l, omega)(s, k);
        const auto sb = alone.Spheroidal(l, omega)(0, 0);
        auto sScale = 0.0;
        for (auto j : {0, 1, 2, 3}) {
          sScale = std::max({sScale, std::abs(sb.U[j]), std::abs(sb.V[j])});
        }
        for (auto j : {0, 1, 2, 3}) {
          EXPECT_TRUE(close(sa.U[j], sb.U[j], sScale) &&
                      close(sa.V[j], sb.V[j], sScale))
              << "spheroidal l " << l << " source " << s << " receiver " << k
              << " jump " << j;
        }
      }
    }
  }
}

TEST(Solver, RefusesASourceInTheOcean) {
  const auto model =
      Model::ReadDeck(std::string(YSPEC_EXAMPLES_DIR) + "/prem.200");
  const auto surface = model.SurfaceRadius();
  const auto oceanFloor = model.Regions().back().lower;
  EXPECT_THROW(
      RadialSolver(model, (surface + oceanFloor) / 2, surface, Gravity::Full),
      std::invalid_argument);
}

TEST(Solver, PutsASourceOnTheCoreMantleBoundaryInTheMantle) {
  const auto model =
      Model::ReadDeck(std::string(YSPEC_EXAMPLES_DIR) + "/prem.200");
  const auto cmb = model.Regions()[1].upper;
  const auto solver =
      RadialSolver(model, cmb, model.Regions()[2].upper, Gravity::Full);
  EXPECT_EQ(solver.Source().layer, model.Regions()[2].first);
  const auto y = solver.Spheroidal(10, Complex(20.0, -1.0))(0, 0);
  EXPECT_TRUE(std::isfinite(std::abs(y.U[0])));
}

TEST(Solver, MatchesTheExteriorPotentialAtAFreeSurface) {
  // Below a free surface phi' = -(l+1) phi / r - 4 pi G rho U, which joins
  // the potential to its exterior solution: the flux variable vanishes there,
  // which is what makes the potential correction at the surface the original
  // YSpec's.
  const auto model =
      Model::ReadDeck(std::string(YSPEC_EXAMPLES_DIR) + "/prem.200.no");
  const auto surface = model.SurfaceRadius();
  const auto solver =
      RadialSolver(model, surface - 0.005, surface, Gravity::Full);
  const auto rho = model(solver.Receiver().layer, surface).rho;
  for (auto l : {2, 20}) {
    const auto y = solver.Spheroidal(l, Complex(20.0, -0.5))(0, 0);
    for (auto j : {0, 1, 2, 3}) {
      const auto exterior =
          -(l + 1.0) * y.Phi[j] / surface - 4 * PiG * rho * y.U[j];
      EXPECT_LT(std::abs(y.DPhi[j] - exterior), 1e-9 * std::abs(y.DPhi[j]))
          << l << " " << j;
    }
  }
}

TEST(Solver, RefusesOptionsOutOfRange) {
  const auto model =
      Model::ReadDeck(std::string(YSPEC_EXAMPLES_DIR) + "/prem.200");
  auto options = SolverOptions{};
  options.tolerance = 0;
  EXPECT_THROW(RadialSolver(model, 0.99, 0.995, Gravity::Full, true, options),
               std::invalid_argument);
}

TEST(Solver, StartsDeeperForLowerDegreesAndHigherFrequencies) {
  const auto model =
      Model::ReadDeck(std::string(YSPEC_EXAMPLES_DIR) + "/prem.200");
  const auto surface = model.SurfaceRadius();
  const auto solver =
      RadialSolver(model, surface - 0.005, surface - 0.001, Gravity::Full);
  const auto omega = 2 * std::numbers::pi * 0.02 * model.GetScales().Time();
  const auto top = solver.Source();
  const auto high = solver.StartPoint(600, omega, top, 0);
  const auto low = solver.StartPoint(100, omega, top, 0);
  EXPECT_GT(high.r, low.r);
  EXPECT_GT(high.r, 0.5);
  // Below about l = 10 at 20 mHz nothing is evanescent enough.
  EXPECT_EQ(solver.StartPoint(2, omega, top, 0).r, 0.0);
}

}  // namespace
