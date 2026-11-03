#pragma once

/**
 * @file Bessel.hpp
 * @brief The ratio of spherical Bessel functions that the regular solutions
 * in a homogeneous sphere are built from.
 */

#include <NumericConcepts/Numeric.hpp>
#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>

namespace YSpec {

/**
 * @brief The ratio @f$q_l(z) = x\, j_{l+1}(x) / j_l(x)@f$ with @f$z = x^2@f$.
 *
 * @details @f$j_l@f$ is the spherical Bessel function of the first kind. The
 * ratio is an even function of @f$x@f$, so it depends only on @f$z@f$, and it
 * is analytic there: no square root, and so no choice of branch, is needed for
 * a complex wavenumber. The recurrence for @f$j_l@f$ gives
 *
 * @f[ q_l(z) = \frac{z}{2l + 3 - q_{l+1}(z)}, @f]
 *
 * a continued fraction evaluated by the modified Lentz method (Thompson &
 * Barnett 1986, J. Comput. Phys. 64, 490). It converges for every @f$z@f$ off
 * the poles, in about @f$\sqrt{|z|}@f$ terms when @f$|z| \gg l^2@f$ and in a
 * few otherwise, which is where the starting solutions use it.
 *
 * @param l The degree, non-negative.
 * @param z The square of the argument.
 * @throws std::runtime_error if the fraction fails to converge.
 */
template <NumericConcepts::RealOrComplex Scalar>
Scalar BesselRatio(int l, Scalar z) {
  using Real = NumericConcepts::RemoveComplex<Scalar>;
  constexpr auto maxTerms = 100000;
  constexpr auto tolerance = 4 * std::numeric_limits<Real>::epsilon();
  constexpr auto tiny = std::numeric_limits<Real>::min() * 1e10;

  // The denominator 2l+3 - z/(2l+5 - z/(2l+7 - ...)), with f, c and d as in
  // Lentz's method: the approximant and its two running ratios.
  auto f = Scalar(2 * l + 3);
  auto c = f;
  auto d = Scalar(0);
  for (auto k = 1; k <= maxTerms; ++k) {
    const auto b = Scalar(2 * (l + k) + 3);
    d = b - z * d;
    if (std::abs(d) < tiny) d = tiny;
    c = b - z / c;
    if (std::abs(c) < tiny) c = tiny;
    d = Scalar(1) / d;
    const auto delta = c * d;
    f *= delta;
    if (std::abs(delta - Scalar(1)) < tolerance) return z / f;
  }
  throw std::runtime_error("BesselRatio: continued fraction did not converge");
}

}  // namespace YSpec
