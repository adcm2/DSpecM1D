#pragma once

// The pieces from which the radial kernels are built: the mesh in radius, the
// products of the forward and adjoint fields at one degree and radius, and
// the accumulation of the products at one frequency into the kernels. They
// are kept apart from the computation so that each can be tested on fields
// given by hand.

#include <GaussQuad/All>
#include <YSpec/Kernels.hpp>
#include <YSpec/Model.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <ranges>
#include <span>
#include <stdexcept>
#include <vector>

#include "Coefficients.hpp"
#include "KernelParts.hpp"
#include "OrderArray.hpp"
#include "SourceWeights.hpp"

namespace YSpec {

// The mesh in radius, non-dimensional: elements between the discontinuities
// of the model and the given breaks, divided to be no larger than asked, each
// with the nodes of a Gauss-Legendre rule. Elements never cross a
// discontinuity, a source or a receiver, where the fields or their
// derivatives jump.
class Mesh {
 public:
  Mesh(const Model& model, const RadialMesh& options,
       std::span<const double> breaks)
      : order_{options.order} {
    if (order_ < 1 || !(options.maxElement > 0)) {
      throw std::invalid_argument(
          "the radial mesh needs at least one node in an element, and "
          "elements of positive size");
    }
    const auto rule = GaussQuad::GaussLegendreQuadrature1D<double>(order_);
    const auto& x = rule.Points();
    const auto& w = rule.Weights();
    const auto maxElement = options.maxElement / model.GetScales().Length();
    // Edges closer than this are one.
    const auto close = 1e-9 * model.SurfaceRadius();
    for (auto layer : std::views::iota(Int{0}, model.NumberOfLayers())) {
      const auto& here = model.GetLayer(layer);
      auto edges = std::vector{here.Lower(), here.Upper()};
      for (auto r : breaks) {
        if (r > edges.front() && r < edges.back()) edges.push_back(r);
      }
      std::ranges::sort(edges);
      auto lower = edges.front();
      for (auto upper : edges | std::views::drop(1)) {
        if (upper - lower < close) continue;
        const auto pieces = std::max(
            1,
            static_cast<int>(std::ceil((upper - lower) / maxElement - 1e-9)));
        const auto half = (upper - lower) / (2 * pieces);
        for (auto piece : std::views::iota(0, pieces)) {
          const auto centre = lower + (2 * piece + 1) * half;
          halves_.push_back(half);
          for (auto j : std::views::iota(0, order_)) {
            r_.push_back(centre + half * x[j]);
            weights_.push_back(half * w[j]);
            layers_.push_back(layer);
          }
        }
        lower = upper;
      }
    }
    // The integrals from each node to the top of its element of the
    // Lagrange polynomials through the nodes, on [-1, 1], by the same rule
    // mapped onto [x_i, 1], which is exact for them.
    auto basis = [&](int j, double t) {
      auto value = 1.0;
      for (auto k : std::views::iota(0, order_)) {
        if (k != j) value *= (t - x[k]) / (x[j] - x[k]);
      }
      return value;
    };
    above_.assign(static_cast<std::size_t>(order_ * order_), 0.0);
    for (auto i : std::views::iota(0, order_)) {
      const auto half = (1 - x[i]) / 2;
      for (auto j : std::views::iota(0, order_)) {
        auto sum = 0.0;
        for (auto q : std::views::iota(0, order_)) {
          sum += w[q] * basis(j, x[i] + half * (x[q] + 1));
        }
        above_[static_cast<std::size_t>(i * order_ + j)] = half * sum;
      }
    }
  }

  std::size_t Size() const { return r_.size(); }
  const std::vector<double>& Radii() const { return r_; }
  double Radius(std::size_t i) const { return r_[i]; }
  // The weight of a node in the integral over radius.
  double Weight(std::size_t i) const { return weights_[i]; }
  Int Layer(std::size_t i) const { return layers_[i]; }

  // The integral of f from each node to the top of the model.
  std::vector<double> IntegralAbove(std::span<const double> f) const {
    const auto order = static_cast<std::size_t>(order_);
    auto result = std::vector<double>(Size());
    auto tail = 0.0;
    for (auto element = halves_.size(); element-- > 0;) {
      const auto first = element * order;
      for (auto i : std::views::iota(std::size_t{0}, order)) {
        auto partial = 0.0;
        for (auto j : std::views::iota(std::size_t{0}, order)) {
          partial += above_[i * order + j] * f[first + j];
        }
        result[first + i] = tail + halves_[element] * partial;
      }
      for (auto j : std::views::iota(std::size_t{0}, order)) {
        tail += weights_[first + j] * f[first + j];
      }
    }
    return result;
  }

 private:
  int order_;
  std::vector<double> r_;
  std::vector<double> weights_;
  std::vector<Int> layers_;
  std::vector<double> halves_;  // half the width of each element
  std::vector<double> above_;
};

// The sum over orders of (-1)^m a[m] b[-m]. For fields of spin weights N and
// -N, expanded in X^N_lm exp(i m phi) and X^-N_lm exp(i m phi), this is the
// integral of their product over the unit sphere but for a factor (-1)^N,
// which the products below fold in; for the displacement, with u expanded as
// U Y r + V grad1 Y - W r x grad1 Y, the integral of u . u° is
// Pair(U, U°) + zeta^2 (Pair(V, V°) + Pair(W, W°)).
inline Complex Pair(const OrderArray<Complex>& a, const OrderArray<Complex>& b,
                    int l) {
  auto sum = Complex{};
  for (auto m : OrderArray<Complex>::Orders(l)) {
    const auto term = a[m] * b[-m];
    sum += m % 2 == 0 ? term : -term;
  }
  return sum;
}

// The coefficients of a function of a field's coefficients, order by order.
template <class Function>
OrderArray<Complex> Each(int l, Function f) {
  auto result = OrderArray<Complex>{};
  for (auto m : OrderArray<Complex>::Orders(l)) result[m] = f(m);
  return result;
}

// The horizontal trace of the strain.
inline OrderArray<Complex> Trace(const Coefficients& c, int l, double r) {
  const auto zeta2 = l * (l + 1.0);
  return Each(l, [&](int m) { return (2.0 * c.u[m] - zeta2 * c.v[m]) / r; });
}

// The consoidal and toroidal parts of 2 e_rh, the shear strain on a
// horizontal plane.
inline OrderArray<Complex> Shear(const Coefficients& c, int l, double r) {
  return Each(l, [&](int m) { return c.dv[m] + (c.u[m] - c.v[m]) / r; });
}
inline OrderArray<Complex> Twist(const Coefficients& c, int l, double r) {
  return Each(l, [&](int m) { return c.dw[m] - c.w[m] / r; });
}

// The integrals over a sphere of the products of the forward and adjoint
// fields that the kernels need, summed over degree. With h the horizontal
// trace of the strain,
//
//   A: h h°,   C: e_rr e°_rr,   F: e_rr h° + h e°_rr,
//   L: 4 e_rh . e°_rh,   N: 2 e_hh : e°_hh - 2 h h°,
//
// the coefficients of the moduli in the elastic energy density, and
//
//   kinetic: s . s°,   radial: s_r s°_r,
//   coupling: r (s_h . grad_h s°_r + s°_h . grad_h s_r),
//   potential: s° . grad phi + s . grad phi°.
//
// The radial products are also kept apart at the degrees with a potential,
// where the gravitational terms include the density through Poisson's
// equation.
struct Products {
  Complex A, C, F, L, N;
  Complex kinetic;
  Complex radial;
  Complex radialWithPotential;
  Complex coupling;
  Complex potential;

  // They are linear in the adjoint field.
  Products& operator+=(const Products& p) {
    const auto mine = Members();
    const auto theirs = p.Members();
    for (auto i : std::views::iota(std::size_t{0}, mine.size())) {
      *mine[i] += *theirs[i];
    }
    return *this;
  }
  friend Products operator*(Complex factor, Products p) {
    for (auto* member : p.Members()) *member *= factor;
    return p;
  }

 private:
  std::array<Complex*, 10> Members() {
    return {&A,        &C,        &F,      &L,
            &N,        &kinetic,  &radial, &radialWithPotential,
            &coupling, &potential};
  }
  std::array<const Complex*, 10> Members() const {
    return {&A,        &C,        &F,      &L,
            &N,        &kinetic,  &radial, &radialWithPotential,
            &coupling, &potential};
  }
};

// Adds the products at one degree, at radius r. In a fluid the shear terms
// are left zero: the radial derivative of the tangential displacement is not
// given there, and the fluid has no shear modulus to perturb.
inline void Add(Products& p, const Coefficients& a, const Coefficients& b,
                int l, double r, bool solid, bool potential) {
  const auto zeta2 = l * (l + 1.0);
  auto pair = [l](const OrderArray<Complex>& x, const OrderArray<Complex>& y) {
    return Pair(x, y, l);
  };
  const auto ha = Trace(a, l, r);
  const auto hb = Trace(b, l, r);
  const auto hh = pair(ha, hb);
  const auto uu = pair(a.u, b.u);
  const auto tangential = pair(a.v, b.v) + pair(a.w, b.w);
  p.A += hh;
  p.C += pair(a.du, b.du);
  p.F += pair(a.du, hb) + pair(ha, b.du);
  if (solid) {
    p.L += zeta2 * (pair(Shear(a, l, r), Shear(b, l, r)) +
                    pair(Twist(a, l, r), Twist(b, l, r)));
    p.N += zeta2 * (zeta2 - 2) / (r * r) * tangential - hh;
  }
  p.kinetic += uu + zeta2 * tangential;
  p.radial += uu;
  p.coupling += zeta2 * (pair(a.v, b.u) + pair(a.u, b.v));
  if (potential) {
    p.radialWithPotential += uu;
    p.potential += pair(b.u, a.dphi) + pair(a.u, b.dphi) +
                   zeta2 * (pair(b.v, a.phi) + pair(a.v, b.phi)) / r;
  }
}

// The jumps in the forward field's displacement at the source, at radius r,
// from the weights of the unit jumps in rU, zeta r V and zeta r W.
inline Coefficients Jumps(const SourceWeights& weights, int l, double r) {
  auto jumps = Coefficients{};
  if (l == 0) {
    jumps.u[0] = weights.radial[0] / r;
    return jumps;
  }
  const auto zeta = std::sqrt(l * (l + 1.0));
  for (auto m : OrderArray<Complex>::Orders(l)) {
    jumps.u[m] = weights.spheroidal[m][0] / r;
    jumps.v[m] = weights.spheroidal[m][1] / (zeta * r);
    jumps.w[m] = weights.toroidal[m][0] / (zeta * r);
  }
  return jumps;
}

// Adds the products at one degree of the part of the forward field's strain
// concentrated on the source's radius r, the jumps in U, V and W times a delta
// function in r, with the adjoint field there.
inline void AddAtSource(Products& p, const Coefficients& jumps,
                        const Coefficients& b, int l, double r) {
  const auto zeta2 = l * (l + 1.0);
  p.C += Pair(jumps.u, b.du, l);
  p.F += Pair(jumps.u, Trace(b, l, r), l);
  p.L += zeta2 *
         (Pair(jumps.v, Shear(b, l, r), l) + Pair(jumps.w, Twist(b, l, r), l));
}

// What is accumulated at each node: the kernels of the parameters, and after
// them the two angular integrals from which the adjoint potential follows,
// M_rr and the integrand of its radial integral.
inline constexpr auto KernelParameters =
    static_cast<std::size_t>(RadialParameter::Count);
inline constexpr auto RadialMomentSlot = KernelParameters;
inline constexpr auto FluxSlot = KernelParameters + 1;

// Adds the products at one frequency, with its weight in the sum over
// frequency, into the kernels: the sum is of the real parts, and each modulus
// is dispersed by its quality factor. The gravitational terms are those of the
// weight, through g, and where there is a potential those through its
// gradient and, in the gradient of gravity, the density.
inline void Accumulate(RadialAccumulator& sum, const Products& p,
                       const Material& m, double r, double weight,
                       Complex dispersion, Complex omega2, bool gravity) {
  auto add = [&](std::size_t slot, Complex value) {
    sum[slot] += weight * value.real();
  };
  auto at = [](RadialParameter q) { return static_cast<std::size_t>(q); };
  add(at(RadialParameter::A), -(1.0 + m.qA * dispersion) * p.A);
  add(at(RadialParameter::C), -(1.0 + m.qA * dispersion) * p.C);
  add(at(RadialParameter::F), -(1.0 + m.qF * dispersion) * p.F);
  add(at(RadialParameter::L), -(1.0 + m.qMu * dispersion) * p.L);
  add(at(RadialParameter::N), -(1.0 + m.qMu * dispersion) * p.N);
  add(at(RadialParameter::QA), -dispersion * (m.A * p.A + m.C * p.C));
  add(at(RadialParameter::QF), -dispersion * m.F * p.F);
  add(at(RadialParameter::QMu), -dispersion * (m.L * p.L + m.N * p.N));
  auto weightTerms = Complex{};
  if (gravity) {
    const auto flux = (p.coupling - 4.0 * p.radial) / r;
    weightTerms = m.g * flux + 4 * PiG * m.rho * p.radialWithPotential;
    add(RadialMomentSlot, m.rho * p.radialWithPotential);
    add(FluxSlot, m.rho * flux);
  }
  add(at(RadialParameter::Rho), omega2 * p.kinetic - weightTerms - p.potential);
}

}  // namespace YSpec
