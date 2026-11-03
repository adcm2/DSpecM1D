#pragma once

/**
 * @file SourceWeights.hpp
 * @brief The jumps a point source imposes on the radial solutions.
 *
 * @details In a frame with the source at the north pole only orders
 * @f$|m| \le 2@f$ are excited. For each degree and order the source imposes a
 * jump on the displacement and traction at the source radius; the radial
 * solvers give the response to unit jumps, and the response to the source is
 * their combination with the weights computed here. The time dependence is a
 * step, hence the factor @f$1/(i\omega)@f$. The formulas for a moment tensor
 * are those of the original YSpec.
 *
 * A moment tensor imposes a jump in displacement, @f$M\cdot\hat{r}@f$ over
 * the modulus and @f$r^2@f$, and in traction; a force @f$\mathbf{f}@f$ only
 * one in traction, @f$-\mathbf{f}/r^2@f$. The tractions are scaled as the
 * displacements are, so the force's weights are those of the first part of
 * the moment tensor's, with @f$M_{rr}/C@f$, @f$M_{r\theta}/L@f$ and
 * @f$M_{r\phi}/L@f$ replaced by @f$-f_r@f$, @f$-f_\theta@f$ and
 * @f$-f_\phi@f$, and placed in the jumps in traction.
 */

#include <YSpec/Source.hpp>
#include <array>
#include <cmath>
#include <complex>
#include <numbers>

#include "Equations.hpp"
#include "OrderArray.hpp"

namespace YSpec {

/** @brief For each order, the weights of the two unit toroidal jumps. */
using ToroidalWeights = OrderArray<std::array<Complex, 2>>;

/**
 * @brief For each order, the weights of the four unit spheroidal jumps, in
 * @f$(U, V, P, S)@f$.
 */
using SpheroidalWeights = OrderArray<std::array<Complex, 4>>;

/** @brief The weights of the two unit radial jumps. */
using RadialWeights = std::array<Complex, 2>;

/** @brief The weights of a source at one degree, of whichever jumps it has. */
struct SourceWeights {
  ToroidalWeights toroidal;      ///< For l >= 1.
  SpheroidalWeights spheroidal;  ///< For l >= 1.
  RadialWeights radial;          ///< For l = 0.
};

namespace SourceDetails {

// sqrt((2l+1)/4pi) / (i omega r^2), the factor common to every jump.
inline Complex Factor(const Context& c, double r) {
  const auto nu = std::sqrt((2.0 * c.Degree() + 1) / (4 * std::numbers::pi));
  return nu / (Complex(0, 1) * c.Omega() * r * r);
}

}  // namespace SourceDetails

/**
 * @brief The toroidal weights for degree @f$l \ge 1@f$, in the jumps
 * (W, T).
 * @param c The degree and frequency.
 * @param m The material at the source.
 * @param moment The moment tensor, non-dimensional.
 * @param force The force, non-dimensional.
 */
inline ToroidalWeights Toroidal(const Context& c, const ComplexMaterial& m,
                                const MomentTensor& moment,
                                const Force& force) {
  using namespace std::complex_literals;
  const auto rs = m.r;
  // sqrt((l-1)(l+2)), the factor of the terms of order +-2.
  const auto rootZeta2m2 = std::sqrt(std::max(c.Zeta2() - 2, 0.0));
  const auto factor = SourceDetails::Factor(c, rs);
  auto w = ToroidalWeights{};
  w[-1][0] = 0.5 * rs * (-moment.rp + 1i * moment.rt) / m.L * factor;
  w[1][0] = 0.5 * rs * (moment.rp + 1i * moment.rt) / m.L * factor;
  w[-2][1] = 0.25 * rootZeta2m2 *
             (-1i * moment.tt + 1i * moment.pp + 2.0 * moment.tp) * factor;
  w[2][1] = 0.25 * rootZeta2m2 *
            (1i * moment.tt - 1i * moment.pp + 2.0 * moment.tp) * factor;
  w[-1][1] += 0.5 * rs * (force.p - 1i * force.t) * factor;
  w[1][1] += 0.5 * rs * (-force.p - 1i * force.t) * factor;
  return w;
}

/**
 * @brief The spheroidal weights for degree @f$l \ge 1@f$.
 */
inline SpheroidalWeights Spheroidal(const Context& c, const ComplexMaterial& m,
                                    const MomentTensor& moment,
                                    const Force& force) {
  using namespace std::complex_literals;
  const auto rs = m.r;
  const auto zeta = c.Zeta();
  // sqrt((l-1)(l+2)), the factor of the terms of order +-2.
  const auto rootZeta2m2 = std::sqrt(std::max(c.Zeta2() - 2, 0.0));
  const auto factor = SourceDetails::Factor(c, rs);
  const auto fOverC = m.F / m.C;
  auto w = SpheroidalWeights{};
  // The jumps are in (U, V, P, S).
  w[-2][3] =
      0.25 * rootZeta2m2 * (-moment.tt + moment.pp - 2.0i * moment.tp) * factor;
  w[-1][1] = 0.5 * rs * (moment.rt + 1i * moment.rp) / m.L * factor;
  w[0][0] = rs * moment.rr / m.C * factor;
  w[0][2] = (-moment.tt - moment.pp + 2.0 * moment.rr * fOverC) * factor;
  w[0][3] = (0.5 * zeta * (moment.tt + moment.pp) - zeta * moment.rr * fOverC) *
            factor;
  w[1][1] = 0.5 * rs * (-moment.rt + 1i * moment.rp) / m.L * factor;
  w[2][3] =
      0.25 * rootZeta2m2 * (-moment.tt + moment.pp + 2.0i * moment.tp) * factor;
  w[0][2] += -rs * force.r * factor;
  w[-1][3] += 0.5 * rs * (-force.t - 1i * force.p) * factor;
  w[1][3] += 0.5 * rs * (force.t - 1i * force.p) * factor;
  return w;
}

/**
 * @brief The radial (l = 0) weights.
 */
inline RadialWeights Radial(const Context& c, const ComplexMaterial& m,
                            const MomentTensor& moment, const Force& force) {
  const auto rs = m.r;
  const auto factor = SourceDetails::Factor(c, rs);
  return {
      rs * moment.rr / m.C * factor,
      (-moment.tt - moment.pp + 2.0 * moment.rr * m.F / m.C - rs * force.r) *
          factor};
}

/**
 * @brief The weights of a source at the pole, at the degree of @p c.
 */
inline SourceWeights PoleWeights(const Context& c, const ComplexMaterial& m,
                                 const MomentTensor& moment,
                                 const Force& force) {
  auto weights = SourceWeights{};
  if (c.Degree() == 0) {
    weights.radial = Radial(c, m, moment, force);
  } else {
    weights.toroidal = Toroidal(c, m, moment, force);
    weights.spheroidal = Spheroidal(c, m, moment, force);
  }
  return weights;
}

/**
 * @brief The weights, in orders @f$|m| \le 2@f$ only, of an impulsive force
 * away from the pole.
 *
 * @details A force @f$\mathbf{f}\,\delta(\mathbf{x}-\mathbf{x}_0)@f$
 * excites each vector spherical harmonic in proportion to its complex
 * conjugate at @f$\mathbf{x}_0@f$, and at the pole that gives the weights
 * above. Elsewhere every order is excited, but a field with orders
 * @f$|m| \le 2@f$ alone, as from a source at the pole, meets only those in a
 * product integrated over the sphere. With
 * @f$f_\pm = f_\theta \mp i f_\phi@f$ and @f$X^N = X^N_{lm}(\theta_0)@f$
 * the weights are
 *
 * @f[ P: -r_0 f_r X^0, \qquad
 *     S: -\tfrac{r_0}{2}\left(f_- X^{-1} - f_+ X^{+1}\right), \qquad
 *     T: -\tfrac{i r_0}{2}\left(f_+ X^{+1} + f_- X^{-1}\right), @f]
 *
 * each times @f$e^{-im\phi_0}/r_0^2@f$. The force is impulsive: there is no
 * factor @f$1/(i\omega)@f$, and its spectrum multiplies the weights.
 * @param c The degree and frequency.
 * @param r0 The radius of the force.
 * @param force Its components in (r, theta, phi) at the point.
 * @param x The functions @f$X^N_{lm}(\theta_0)@f$, N = -2, ..., 2, at each
 * order @f$|m| \le \min(l, 2)@f$.
 * @param phi0 The longitude of the point.
 */
inline SourceWeights ForceWeights(const Context& c, double r0,
                                  const Force& force,
                                  const OrderArray<std::array<double, 5>>& x,
                                  double phi0) {
  using namespace std::complex_literals;
  const auto plus = Complex(force.t, -force.p);
  const auto minus = Complex(force.t, force.p);
  auto result = SourceWeights{};
  const auto l = c.Degree();
  for (auto m : OrderArray<Complex>::Orders(l)) {
    const auto e = std::exp(-1i * (m * phi0)) / (r0 * r0);
    const auto& X = x[m];
    if (l == 0) {
      result.radial[1] = -r0 * force.r * X[2] * e;
      continue;
    }
    result.spheroidal[m][2] = -r0 * force.r * X[2] * e;
    result.spheroidal[m][3] = -0.5 * r0 * (minus * X[1] - plus * X[3]) * e;
    result.toroidal[m][1] = -0.5i * r0 * (plus * X[3] + minus * X[1]) * e;
  }
  return result;
}

}  // namespace YSpec
