#pragma once

#include <GSHTrans/GSHTrans.hpp>
#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <limits>
#include <numbers>

// Pins the value convention of the Wigner class.
//
// `Wigner`'s stored value at upper index N, degree l and order m is
//
//     sqrt((2l+1)/(4 pi)) * d^l_{Nm}(theta),
//
// where d^l_{Nm} = P^N_{lm}(cos theta) is the generalised Legendre function of
// Dahlen & Tromp (1998) eq. (C.115) -- the upper index is the FIRST subscript
// of d.
//
// This matters because d^l_{mN} = (-1)^{N-m} d^l_{Nm}: the two candidate
// conventions differ by a sign on exactly those entries with N - m odd, and
// nothing in the class's indexing distinguishes them. The convention becomes
// observable once the reality reduction sums a phase (-1)^N over an orbit of
// stored components, so a silent transposition here would surface as wrong
// numbers a long way from its cause.
//
// The l = 1 row of D&T (C.115) is written out in full below. It discriminates:
// the four entries with N - m odd change sign under transposition. The N = 0
// row is checked independently by CheckLegendre against std::sph_legendre,
// which pins the orthonormal scaling and the Condon-Shortley phase; this file
// pins the index order.

template <std::floating_point Real>
int CheckWignerConvention() {
  using namespace GSHTrans;
  using Int = std::ptrdiff_t;

  // d^1_{Nm}(theta) = P^N_{1m}(cos theta), Dahlen & Tromp (1998) eq. (C.115).
  auto DahlenTromp = [](Int N, Int m, Real theta) -> Real {
    const auto c = std::cos(theta);
    const auto s = std::sin(theta);
    const auto half = static_cast<Real>(1) / static_cast<Real>(2);
    const auto rootHalf = std::numbers::sqrt2_v<Real> * half;
    if (N == -1) {
      if (m == -1) return half * (1 + c);
      if (m == 0) return -rootHalf * s;
      return half * (1 - c);
    }
    if (N == 0) {
      if (m == -1) return rootHalf * s;
      if (m == 0) return c;
      return -rootHalf * s;
    }
    if (m == -1) return half * (1 - c);
    if (m == 0) return rootHalf * s;
    return half * (1 + c);
  };

  // Include the two boundaries and a value either side of pi/2. The two are
  // not alike: at zero the half-angle sine is exactly zero and
  // WignerDetails::Arguments flags the pole, while at pi the half-angle
  // cosine is 6e-17 and not zero, so nothing is flagged and the general path
  // has to be right there unaided.
  const auto angles =
      std::array<Real, 6>{static_cast<Real>(0),   static_cast<Real>(0.3),
                          static_cast<Real>(0.7), std::numbers::pi_v<Real> / 2,
                          static_cast<Real>(2.5), std::numbers::pi_v<Real>};

  constexpr auto eps = 100 * std::numeric_limits<Real>::epsilon();
  constexpr Int lMax = 1;

  // Counts the entries at which the two conventions actually differ, i.e.
  // where this test has discriminating power. Guards the oracle itself: at
  // theta = 0 and theta = pi every entry with N - m odd vanishes in both
  // conventions, so an angle list containing only those would pass while
  // testing nothing about the index order.
  auto discriminating = 0;

  for (auto theta : angles) {
    auto d = Wigner<Real, All, All, Single>(lMax, lMax, lMax, theta);
    const auto scale = std::sqrt(static_cast<Real>(3)) *
                       std::numbers::inv_sqrtpi_v<Real> / static_cast<Real>(2);

    for (Int N = -1; N <= 1; N++) {
      auto view = d[N, 0];
      for (Int m = -1; m <= 1; m++) {
        const auto expected = scale * DahlenTromp(N, m, theta);
        if (std::abs(view[lMax][m] - expected) > eps) return 1;

        // d^l_{mN} = (-1)^{N-m} d^l_{Nm}, so the transposed convention differs
        // exactly where N - m is odd and the value is nonzero.
        if (std::abs(expected - scale * DahlenTromp(m, N, theta)) > eps) {
          discriminating++;
        }
      }
    }
  }

  // Four entries per angle away from the poles: (N, m) with |N - m| = 1.
  if (discriminating == 0) return 2;

  return 0;
}
