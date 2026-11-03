#include <gtest/gtest.h>

#include <YSpec/Model.hpp>
#include <YSpec/RadialSolver.hpp>
#include <algorithm>
#include <cmath>
#include <complex>
#include <map>
#include <numbers>
#include <ranges>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "Equations.hpp"

// Responses on the far side of boundaries between solid and fluid, against
// direct shooting: the solutions themselves, rather than their minors, are
// integrated from the centre and from the surface with a fixed-step
// Runge-Kutta method, joined at each boundary by the plain conditions on
// vectors, and combined at the source by solving the jump conditions. At a
// low degree and frequency the solutions stay independent, so shooting is
// accurate, and it shares nothing with the solver beyond the coefficients of
// the equations and the starting solutions.

namespace {

using namespace YSpec;

const auto kPrem = std::string(YSPEC_EXAMPLES_DIR) + "/prem.200";

using Column = std::vector<Complex>;
using Matrix = std::vector<Column>;  // by rows

// The matrix [[T, K], [S, -T^T]] of a system.
Matrix Hamiltonian(const Matrix& T, const Column& K, const Matrix& S) {
  const auto h = T.size();
  auto a = Matrix(2 * h, Column(2 * h));
  for (auto i = std::size_t{0}; i < h; ++i) {
    a[i][h + i] = K[i];
    for (auto j = std::size_t{0}; j < h; ++j) {
      a[i][j] = T[i][j];
      a[h + i][j] = S[i][j];
      a[h + i][h + j] = -T[j][i];
    }
  }
  return a;
}

Matrix Of(const Block2& b) { return Hamiltonian({{b.t}}, {b.k}, {{b.s}}); }

Matrix Of(const Block4& b) {
  return Hamiltonian({{b.t11, b.t12}, {b.t21, b.t22}}, {b.k11, b.k22},
                     {{b.s11, b.s12}, {b.s12, b.s22}});
}

Matrix Of(const Block6& b) {
  return Hamiltonian(
      {{b.t11, b.t12, 0.0}, {b.t21, b.t22, 0.0}, {b.t31, 0.0, b.t33}},
      {b.k11, b.k22, b.k33},
      {{b.s11, b.s12, b.s13}, {b.s12, b.s22, b.s23}, {b.s13, b.s23, 0.0}});
}

Column Times(const Matrix& a, const Column& y) {
  auto result = Column(a.size());
  for (auto i = std::size_t{0}; i < a.size(); ++i) {
    for (auto j = std::size_t{0}; j < y.size(); ++j) {
      result[i] += a[i][j] * y[j];
    }
  }
  return result;
}

Column Plus(const Column& y, Complex a, const Column& x) {
  auto result = y;
  for (auto i = std::size_t{0}; i < y.size(); ++i) result[i] += a * x[i];
  return result;
}

// Solves a x = b by Gaussian elimination with partial pivoting.
Column Solve(Matrix a, Column b) {
  const auto n = b.size();
  for (auto k = std::size_t{0}; k < n; ++k) {
    auto pivot = k;
    for (auto i = k + 1; i < n; ++i) {
      if (std::abs(a[i][k]) > std::abs(a[pivot][k])) pivot = i;
    }
    std::swap(a[k], a[pivot]);
    std::swap(b[k], b[pivot]);
    for (auto i = k + 1; i < n; ++i) {
      const auto factor = a[i][k] / a[k][k];
      for (auto j = k; j < n; ++j) a[i][j] -= factor * a[k][j];
      b[i] -= factor * b[k];
    }
  }
  auto x = Column(n);
  for (auto k = n; k-- > 0;) {
    auto sum = b[k];
    for (auto j = k + 1; j < n; ++j) sum -= a[k][j] * x[j];
    x[k] = sum / a[k][k];
  }
  return x;
}

enum class Kind { Spheroidal, Radial };

// A set of solutions carried through the model, with the solutions at the
// receivers passed kept in terms of the coefficients of the current set.
class Shooter {
 public:
  Shooter(const Model& model, const Context& c, Gravity gravity, Kind kind)
      : model_{model}, c_{c}, gravity_{gravity}, kind_{kind} {}

  // Starts with the solutions regular at the centre, a little off it.
  void StartBelow(double r) {
    const auto m = c_(0, r);
    solutions_.clear();
    if (kind_ == Kind::Radial) {
      const auto y = RadialStart(c_, m, gravity_ != Gravity::None);
      solutions_.push_back({y.begin(), y.end()});
    } else if (gravity_ == Gravity::Full) {
      for (const auto& y : SolidStartSolutions6(c_, m)) {
        solutions_.push_back({y.begin(), y.end()});
      }
    } else {
      for (const auto& y : SolidStartSolutions4(c_, m)) {
        solutions_.push_back({y.begin(), y.end()});
      }
    }
    inFluid_ = false;
    here_ = {0, r};
  }

  // Starts with the solutions that leave the surface free: unit
  // displacements (and potential) with no traction (and flux).
  void StartAbove() {
    const auto top = model_.NumberOfLayers() - 1;
    inFluid_ = Fluid(top);
    const auto size = Size(inFluid_);
    solutions_.clear();
    for (auto i = std::size_t{0}; i < size / 2; ++i) {
      auto y = Column(size);
      y[i] = 1;
      solutions_.push_back(y);
    }
    here_ = {top, model_.SurfaceRadius()};
  }

  // Integrates to a point, keeping the solutions at the stops, which lie
  // between in the order of travel.
  void WalkTo(Point to,
              const std::vector<std::pair<std::size_t, Point>>& stops) {
    for (const auto& [slot, point] : stops) {
      Walk(point);
      kept_[slot] = {inFluid_, solutions_};
    }
    Walk(to);
  }

  const std::vector<Column>& Solutions() const { return solutions_; }

  // The response at a kept stop for given coefficients of the current set.
  std::pair<bool, Column> Response(std::size_t slot,
                                   const Column& coefficients) const {
    const auto& [fluid, columns] = kept_.at(slot);
    auto y = Column(columns.front().size());
    for (auto i = std::size_t{0}; i < columns.size(); ++i) {
      y = Plus(y, coefficients[i], columns[i]);
    }
    return {fluid, y};
  }

 private:
  const Model& model_;
  const Context& c_;
  Gravity gravity_;
  Kind kind_;
  std::vector<Column> solutions_;
  bool inFluid_ = false;
  Point here_{0, 0};
  std::map<std::size_t, std::pair<bool, std::vector<Column>>> kept_;

  bool Fluid(Int layer) const {
    return kind_ == Kind::Spheroidal && model_.GetLayer(layer).IsFluid();
  }
  std::size_t Size(bool fluid) const {
    if (kind_ == Kind::Radial) return 2;
    const auto solid = gravity_ == Gravity::Full ? std::size_t{6} : 4;
    return fluid ? solid - 2 : solid;
  }

  Matrix At(Int layer, double r) const {
    const auto m = c_(layer, r);
    if (kind_ == Kind::Radial) {
      return Of(RadialBlock(c_, m, gravity_ != Gravity::None));
    }
    const auto cowling = gravity_ == Gravity::Cowling;
    if (Fluid(layer)) {
      return gravity_ == Gravity::Full ? Of(FluidBlock4(c_, m))
                                       : Of(FluidBlock2(c_, m, cowling));
    }
    return gravity_ == Gravity::Full ? Of(SolidBlock6(c_, m))
                                     : Of(SolidBlock4(c_, m, cowling));
  }

  // Fixed steps of the classical fourth-order method, small near the centre.
  void Integrate(Int layer, double from, double to) {
    if (from == to) return;
    const auto h = std::min(2e-4, 0.02 * std::min(from, to));
    const auto steps = static_cast<int>(std::ceil(std::abs(to - from) / h));
    const auto dr = (to - from) / steps;
    for (auto n = 0; n < steps; ++n) {
      const auto r = from + n * dr;
      const auto a0 = At(layer, r);
      const auto a1 = At(layer, r + dr / 2);
      const auto a2 = At(layer, r + dr);
      for (auto& y : solutions_) {
        const auto k1 = Times(a0, y);
        const auto k2 = Times(a1, Plus(y, dr / 2, k1));
        const auto k3 = Times(a1, Plus(y, dr / 2, k2));
        const auto k4 = Times(a2, Plus(y, dr, k3));
        y = Plus(Plus(Plus(Plus(y, dr / 6, k1), dr / 3, k2), dr / 3, k3),
                 dr / 6, k4);
      }
    }
  }

  // The positions of the fluid variables among the solid ones, and those of
  // V and S.
  std::vector<std::size_t> Embedding() const {
    return gravity_ == Gravity::Full ? std::vector<std::size_t>{0, 2, 3, 5}
                                     : std::vector<std::size_t>{0, 2};
  }
  std::size_t SIndex() const { return gravity_ == Gravity::Full ? 4 : 3; }

  void ChangePhase(bool toFluid) {
    if (toFluid == inFluid_) return;
    const auto embedding = Embedding();
    if (toFluid) {
      // The combinations with no tangential traction, through the solution
      // with the most.
      const auto s = SIndex();
      auto pivot = std::size_t{0};
      for (auto i = std::size_t{1}; i < solutions_.size(); ++i) {
        if (std::abs(solutions_[i][s]) > std::abs(solutions_[pivot][s])) {
          pivot = i;
        }
      }
      // new_j = sum_k mixing[j][k] old_k.
      auto mixing = Matrix{};
      for (auto i = std::size_t{0}; i < solutions_.size(); ++i) {
        if (i == pivot) continue;
        auto row = Column(solutions_.size());
        row[i] = solutions_[pivot][s];
        row[pivot] = -solutions_[i][s];
        mixing.push_back(row);
      }
      auto mix = [&](const std::vector<Column>& old) {
        auto mixed = std::vector<Column>{};
        for (const auto& row : mixing) {
          auto y = Column(old.front().size());
          for (auto k = std::size_t{0}; k < old.size(); ++k) {
            y = Plus(y, row[k], old[k]);
          }
          mixed.push_back(y);
        }
        return mixed;
      };
      auto fluid = std::vector<Column>{};
      for (const auto& y : mix(solutions_)) {
        auto f = Column{};
        for (auto k : embedding) f.push_back(y[k]);
        fluid.push_back(f);
      }
      solutions_ = fluid;
      for (auto& [slot, kept] : kept_) kept.second = mix(kept.second);
    } else {
      // The fluid solutions with V = S = 0, and the solution that slips.
      const auto size = Size(false);
      auto solid = std::vector<Column>{};
      for (const auto& f : solutions_) {
        auto y = Column(size);
        for (auto k = std::size_t{0}; k < f.size(); ++k) y[embedding[k]] = f[k];
        solid.push_back(y);
      }
      auto slip = Column(size);
      slip[1] = 1;
      solid.push_back(slip);
      solutions_ = solid;
      for (auto& [slot, kept] : kept_) {
        kept.second.push_back(Column(kept.second.front().size()));
      }
    }
    inFluid_ = toFluid;
  }

  void Walk(Point to) {
    const auto up = to.layer >= here_.layer;
    for (auto layer = here_.layer;; layer += up ? 1 : -1) {
      const auto& L = model_.GetLayer(layer);
      const auto from =
          layer == here_.layer ? here_.r : (up ? L.Lower() : L.Upper());
      const auto end = layer == to.layer ? to.r : (up ? L.Upper() : L.Lower());
      ChangePhase(Fluid(layer));
      Integrate(layer, from, end);
      if (layer == to.layer) break;
    }
    here_ = to;
  }
};

// A response in physical variables, as the solver gives it.
struct Physical {
  Complex U, V, Phi, DPhi;
};

Physical ToPhysical(const Context& c, const Model& model, Point point,
                    Gravity gravity, bool fluid, const Column& y) {
  const auto r = point.r;
  const auto m = c(point.layer, r);
  const auto zeta = c.Zeta();
  const auto l = static_cast<double>(c.Degree());
  auto result = Physical{.U = y[0] / r, .V = {}, .Phi = {}, .DPhi = {}};
  auto rPhi = Complex{};
  auto psi = Complex{};
  if (!fluid) {
    result.V = y[1] / (zeta * r);
    if (gravity == Gravity::Full) {
      rPhi = y[2];
      psi = y[5];
    }
  } else {
    // In a fluid the tangential traction vanishes, which gives V.
    auto P = y[1];
    if (gravity == Gravity::Full) {
      rPhi = y[1];
      P = y[2];
      psi = y[3];
    }
    const auto g = gravity == Gravity::None ? 0.0 : m.g;
    result.V = (m.rho * (g * y[0] + rPhi) - P) / (c.Omega2() * m.rho * r * r);
  }
  if (gravity == Gravity::Full) {
    result.Phi = rPhi / r;
    const auto rho = model(point.layer, r).rho;
    result.DPhi = -4 * PiG * rho * result.U - (l + 1) * result.Phi / r +
                  4 * PiG * psi / r;
  }
  return result;
}

const auto kDepths = std::vector<double>{
    0, 1.5, 3, 10, 100, 1000, 2891, 2900, 3500, 5000, 5150, 5160, 5800, 6300};

class SolverAcrossBoundaries
    : public ::testing::TestWithParam<std::tuple<Gravity, int, double>> {};

TEST_P(SolverAcrossBoundaries, AgreeWithShooting) {
  const auto [gravity, l, frequency] = GetParam();
  const auto model = Model::ReadDeck(kPrem);
  const auto length = model.GetScales().Length();
  const auto w = 2 * std::numbers::pi * frequency * model.GetScales().Time();
  const auto omega = Complex(w, -0.02 * w);
  const auto rs = 1 - 30e3 / length;
  auto radii = std::vector<double>{};
  for (auto depth : kDepths) radii.push_back(1 - depth * 1e3 / length);
  const auto solver = RadialSolver(
      model, {rs}, radii, gravity, true,
      SolverOptions{.tolerance = 1e-10, .startDecay = 30, .stepFactor = 0.25});
  const auto responses = solver.Spheroidal(l, omega);
  const auto c = Context(model, l, omega, true);

  // Shooting from below to the source, and from above.
  const auto source = solver.Source(0);
  auto below = std::vector<std::pair<std::size_t, Point>>{};
  auto above = std::vector<std::pair<std::size_t, Point>>{};
  for (auto k = std::size_t{0}; k < radii.size(); ++k) {
    const auto point = solver.Receiver(k);
    (point.r < source.r ? below : above).push_back({k, point});
  }
  std::ranges::sort(below, {}, [](const auto& p) { return p.second.r; });
  std::ranges::sort(above, {}, [](const auto& p) { return -p.second.r; });
  auto fromBelow = Shooter(model, c, gravity, Kind::Spheroidal);
  fromBelow.StartBelow(1e-3);
  fromBelow.WalkTo(source, below);
  auto fromAbove = Shooter(model, c, gravity, Kind::Spheroidal);
  fromAbove.StartAbove();
  fromAbove.WalkTo(source, above);

  // The jumps y(above) - y(below) in (U, V, P, S).
  const auto& a = fromAbove.Solutions();
  const auto& b = fromBelow.Solutions();
  const auto n = a.front().size();
  const auto half = a.size();
  auto system = Matrix(n, Column(n));
  for (auto i = std::size_t{0}; i < n; ++i) {
    for (auto j = std::size_t{0}; j < half; ++j) {
      system[i][j] = a[j][i];
      system[i][half + j] = -b[j][i];
    }
  }
  const auto positions = gravity == Gravity::Full
                             ? std::array<std::size_t, 4>{0, 1, 3, 4}
                             : std::array<std::size_t, 4>{0, 1, 2, 3};
  for (auto jump : {0, 1, 2, 3}) {
    auto rhs = Column(n);
    rhs[positions[jump]] = 1;
    const auto x = Solve(system, rhs);
    const auto alpha = Column(x.begin(), x.begin() + half);
    const auto beta = Column(x.begin() + half, x.end());

    auto expected = std::vector<Physical>(radii.size());
    auto scale = Physical{};
    for (auto k = std::size_t{0}; k < radii.size(); ++k) {
      const auto point = solver.Receiver(k);
      const auto [fluid, y] = point.r < source.r ? fromBelow.Response(k, beta)
                                                 : fromAbove.Response(k, alpha);
      expected[k] = ToPhysical(c, model, point, gravity, fluid, y);
      scale.U = std::max(std::abs(scale.U), std::abs(expected[k].U));
      scale.V = std::max(std::abs(scale.V), std::abs(expected[k].V));
      scale.Phi = std::max(std::abs(scale.Phi), std::abs(expected[k].Phi));
      scale.DPhi = std::max(std::abs(scale.DPhi), std::abs(expected[k].DPhi));
    }
    for (auto k = std::size_t{0}; k < radii.size(); ++k) {
      const auto& got = responses(0, k);
      const auto& want = expected[k];
      const auto where = ::testing::Message() << "jump " << jump << ", depth "
                                              << kDepths[k] << " km";
      constexpr auto tolerance = 1e-7;
      EXPECT_NEAR(std::abs(got.U[jump] - want.U), 0,
                  tolerance * std::abs(scale.U))
          << where;
      EXPECT_NEAR(std::abs(got.V[jump] - want.V), 0,
                  tolerance * std::abs(scale.V))
          << where;
      if (gravity == Gravity::Full) {
        EXPECT_NEAR(std::abs(got.Phi[jump] - want.Phi), 0,
                    tolerance * std::abs(scale.Phi))
            << where;
        EXPECT_NEAR(std::abs(got.DPhi[jump] - want.DPhi), 0,
                    tolerance * std::abs(scale.DPhi))
            << where;
      }
    }
  }
}

INSTANTIATE_TEST_SUITE_P(Cases, SolverAcrossBoundaries,
                         ::testing::Combine(::testing::Values(Gravity::None,
                                                              Gravity::Cowling,
                                                              Gravity::Full),
                                            ::testing::Values(2, 5),
                                            ::testing::Values(2e-3, 4e-3)));

TEST(SolverAcrossBoundaries, GivesTheRadialResponseEverywhere) {
  const auto model = Model::ReadDeck(kPrem);
  const auto length = model.GetScales().Length();
  for (auto gravity : {Gravity::None, Gravity::Full}) {
    const auto w = 2 * std::numbers::pi * 3e-3 * model.GetScales().Time();
    const auto omega = Complex(w, -0.02 * w);
    const auto rs = 1 - 30e3 / length;
    auto radii = std::vector<double>{};
    for (auto depth : kDepths) radii.push_back(1 - depth * 1e3 / length);
    const auto solver = RadialSolver(
        model, {rs}, radii, gravity, true,
        SolverOptions{
            .tolerance = 1e-10, .startDecay = 30, .stepFactor = 0.25});
    const auto responses = solver.Radial(omega);
    const auto c = Context(model, 0, omega, true);
    const auto source = solver.Source(0);
    auto below = std::vector<std::pair<std::size_t, Point>>{};
    auto above = std::vector<std::pair<std::size_t, Point>>{};
    for (auto k = std::size_t{0}; k < radii.size(); ++k) {
      const auto point = solver.Receiver(k);
      (point.r < source.r ? below : above).push_back({k, point});
    }
    std::ranges::sort(below, {}, [](const auto& p) { return p.second.r; });
    std::ranges::sort(above, {}, [](const auto& p) { return -p.second.r; });
    auto fromBelow = Shooter(model, c, gravity, Kind::Radial);
    fromBelow.StartBelow(1e-3);
    fromBelow.WalkTo(source, below);
    auto fromAbove = Shooter(model, c, gravity, Kind::Radial);
    fromAbove.StartAbove();
    fromAbove.WalkTo(source, above);
    const auto& a = fromAbove.Solutions().front();
    const auto& b = fromBelow.Solutions().front();
    const auto system = Matrix{{a[0], -b[0]}, {a[1], -b[1]}};
    for (auto jump : {0, 1}) {
      auto rhs = Column(2);
      rhs[jump] = 1;
      const auto x = Solve(system, rhs);
      auto expected = Column(radii.size());
      auto scale = 0.0;
      for (auto k = std::size_t{0}; k < radii.size(); ++k) {
        const auto point = solver.Receiver(k);
        const auto [fluid, y] = point.r < source.r
                                    ? fromBelow.Response(k, {x[1]})
                                    : fromAbove.Response(k, {x[0]});
        expected[k] = y[0] / point.r;
        scale = std::max(scale, std::abs(expected[k]));
      }
      for (auto k = std::size_t{0}; k < radii.size(); ++k) {
        EXPECT_NEAR(std::abs(responses(0, k).U[jump] - expected[k]), 0,
                    1e-7 * scale)
            << "jump " << jump << ", depth " << kDepths[k] << " km";
      }
    }
  }
}

TEST(SolverAcrossBoundaries, GivesNoToroidalResponseOutsideTheSourceRegion) {
  const auto model = Model::ReadDeck(kPrem);
  const auto length = model.GetScales().Length();
  auto radii = std::vector<double>{};
  for (auto depth : kDepths) radii.push_back(1 - depth * 1e3 / length);
  // Integrations shared with other receivers stop at more radii: tightly
  // enough, that does not matter.
  const auto options =
      SolverOptions{.tolerance = 1e-10, .startDecay = 11, .stepFactor = 0.25};
  const auto solver = RadialSolver(model, 1 - 30e3 / length, 1 - 100e3 / length,
                                   Gravity::Full, true, options);
  const auto all = RadialSolver(model, {1 - 30e3 / length}, radii,
                                Gravity::Full, true, options);
  const auto omega = Complex(10, -0.2);
  const auto responses = all.Toroidal(3, omega);
  for (auto k = std::size_t{0}; k < radii.size(); ++k) {
    const auto point = all.Receiver(k);
    const auto inMantle =
        model.RegionIndex(point.layer) == model.RegionIndex(all.Source().layer);
    for (auto jump : {0, 1}) {
      if (inMantle) {
        EXPECT_NE(responses(0, k).W[jump], Complex{}) << kDepths[k];
      } else {
        EXPECT_EQ(responses(0, k).W[jump], Complex{}) << kDepths[k];
      }
    }
  }
  // The receivers outside the region leave those inside as they were.
  const auto alone = solver.Toroidal(3, omega)(0, 0);
  const auto index = static_cast<std::size_t>(
      std::ranges::find(kDepths, 100.0) - kDepths.begin());
  for (auto jump : {0, 1}) {
    EXPECT_NEAR(std::abs(responses(0, index).W[jump] - alone.W[jump]), 0,
                1e-6 * std::abs(alone.W[jump]));
  }
}

TEST(SolverAcrossBoundaries, ReachesCloseToTheCentre) {
  // Displacement at the centre comes from degree one alone: the responses
  // of higher degrees vanish there, as r^(l - 1).
  const auto model = Model::ReadDeck(kPrem);
  const auto length = model.GetScales().Length();
  const auto radii = std::vector<double>{1e-4, 1e-3, 1e-2};
  const auto solver =
      RadialSolver(model, {1 - 30e3 / length}, radii, Gravity::Full);
  const auto omega = Complex(12, -0.2);
  for (auto l : {1, 2, 3}) {
    const auto responses = solver.Spheroidal(l, omega);
    for (auto jump : {0, 1, 2, 3}) {
      const auto far = std::abs(responses(0, 2).U[jump]);
      const auto near = std::abs(responses(0, 0).U[jump]);
      ASSERT_TRUE(std::isfinite(near) && std::isfinite(far));
      if (l == 1) {
        EXPECT_NEAR(near, far, 1e-2 * far) << "jump " << jump;
      } else {
        EXPECT_LT(near, 2 * std::pow(1e-2, l - 1) * far)
            << "degree " << l << ", jump " << jump;
      }
    }
  }
  EXPECT_THROW((void)RadialSolver(model, 0.5, 0.0, Gravity::Full),
               std::invalid_argument);
}

TEST(SolverAcrossBoundaries, GivesTheRadialDerivatives) {
  // Against centred differences of the responses, in the mantle and in the
  // outer core (where only U has a derivative given).
  const auto model = Model::ReadDeck(kPrem);
  const auto length = model.GetScales().Length();
  const auto h = 1e-5;
  const auto options =
      SolverOptions{.tolerance = 1e-11, .startDecay = 11, .stepFactor = 0.25};
  for (auto gravity : {Gravity::None, Gravity::Full}) {
    for (auto depth : {800e3, 4000e3}) {
      const auto r = 1 - depth / length;
      const auto solver =
          RadialSolver(model, {1 - 30e3 / length}, {r - h, r, r + h}, gravity,
                       true, options);
      const auto fluid = model.GetLayer(solver.Receiver(1).layer).IsFluid();
      const auto omega = Complex(20, -0.3);
      auto near = [&](Complex derivative, Complex below, Complex above) {
        const auto difference = (above - below) / (2 * h);
        EXPECT_NEAR(std::abs(derivative - difference), 0,
                    1e-5 * std::abs(difference))
            << "depth " << depth / 1e3 << " km";
      };
      const auto spheroidal = solver.Spheroidal(4, omega);
      const auto toroidal = solver.Toroidal(4, omega);
      const auto radial = solver.Radial(omega);
      for (auto j : {0, 1, 2, 3}) {
        near(spheroidal(0, 1).DU[j], spheroidal(0, 0).U[j],
             spheroidal(0, 2).U[j]);
        if (!fluid) {
          near(spheroidal(0, 1).DV[j], spheroidal(0, 0).V[j],
               spheroidal(0, 2).V[j]);
        }
      }
      for (auto j : {0, 1}) {
        near(radial(0, 1).DU[j], radial(0, 0).U[j], radial(0, 2).U[j]);
        if (!fluid) {
          near(toroidal(0, 1).DW[j], toroidal(0, 0).W[j], toroidal(0, 2).W[j]);
        }
      }
    }
  }
}

TEST(SolverAcrossBoundaries, GivesEachSideOfADiscontinuity) {
  // Receivers at a discontinuity, in the layer on each side, against
  // receivers just inside each: at the core-mantle boundary, between solid
  // and fluid, and at 670 km, between two solids. Listed in both orders, the
  // integrations must pass them in the order of travel.
  const auto model = Model::ReadDeck(kPrem);
  const auto length = model.GetScales().Length();
  const auto h = 1e-9;
  for (auto radius : {3480e3, 5701e3}) {
    const auto r = radius / length;
    const auto below = model.LayerIndex(r, Side::Below);
    const auto above = model.LayerIndex(r, Side::Above);
    ASSERT_NE(below, above);
    for (auto reversed : {false, true}) {
      auto points = std::vector<Point>{
          {below, r - h}, {below, r}, {above, r}, {above, r + h}};
      if (reversed) std::ranges::reverse(points);
      const auto solver =
          RadialSolver(model, {1 - 30e3 / length}, points, Gravity::Full);
      auto index = [&](Int layer, double x) {
        for (auto i = std::size_t{0}; i < points.size(); ++i) {
          if (points[i].layer == layer && points[i].r == x) return i;
        }
        return points.size();
      };
      const auto omega = Complex(15, -0.2);
      const auto spheroidal = solver.Spheroidal(3, omega);
      const auto toroidal = solver.Toroidal(3, omega);
      for (auto [layer, x, near] :
           {std::tuple{below, r, r - h}, std::tuple{above, r, r + h}}) {
        const auto& at = spheroidal(0, index(layer, x));
        const auto& by = spheroidal(0, index(layer, near));
        const auto& atT = toroidal(0, index(layer, x));
        const auto& byT = toroidal(0, index(layer, near));
        for (auto j : {0, 1, 2, 3}) {
          for (auto [a, b] :
               {std::pair{at.U[j], by.U[j]}, std::pair{at.V[j], by.V[j]},
                std::pair{at.DU[j], by.DU[j]},
                std::pair{at.Phi[j], by.Phi[j]}}) {
            EXPECT_NEAR(std::abs(a - b), 0, 1e-6 * std::abs(b) + 1e-300)
                << radius / 1e3 << " km, layer " << layer << ", jump " << j;
          }
        }
        for (auto j : {0, 1}) {
          EXPECT_NEAR(std::abs(atT.W[j] - byT.W[j]), 0,
                      1e-6 * std::abs(byT.W[j]) + 1e-300);
        }
      }
    }
  }
}

TEST(SolverAcrossBoundaries, IsFiniteOnEverySideOfEveryDiscontinuity) {
  // With receivers on both sides of every discontinuity, the integration
  // from below may start at one of them, in the layer above: those below it
  // must be left out, not reached by integrating back across the boundary.
  const auto model = Model::ReadDeck(kPrem);
  const auto length = model.GetScales().Length();
  auto points = std::vector<Point>{};
  for (auto layer = Int{0}; layer < model.NumberOfLayers(); ++layer) {
    const auto& here = model.GetLayer(layer);
    if (here.Lower() > 0) points.push_back({layer, here.Lower()});
    points.push_back({layer, here.Upper()});
  }
  const auto solver =
      RadialSolver(model, {1 - 30e3 / length}, points, Gravity::Full);
  for (auto l : {1, 5, 11, 20, 40}) {
    for (auto f : {1.0, 4.0, 10.0, 25.0}) {
      const auto omega = Complex(f, -0.05);
      const auto spheroidal = solver.Spheroidal(l, omega);
      const auto toroidal = solver.Toroidal(l, omega);
      for (auto i = std::size_t{0}; i < points.size(); ++i) {
        for (auto j : {0, 1, 2, 3}) {
          EXPECT_TRUE(std::isfinite(std::abs(spheroidal(0, i).U[j])) &&
                      std::isfinite(std::abs(spheroidal(0, i).DU[j])))
              << "degree " << l << ", frequency " << f << ", point " << i;
        }
        EXPECT_TRUE(std::isfinite(std::abs(toroidal(0, i).W[0])));
      }
    }
  }
}

}  // namespace
