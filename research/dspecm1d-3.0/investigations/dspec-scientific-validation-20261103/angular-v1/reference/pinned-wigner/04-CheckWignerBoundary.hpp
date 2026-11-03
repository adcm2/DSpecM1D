#pragma once

#include <GSHTrans/GSHTrans.hpp>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <vector>

// The values Wigner.hpp generates by recursion, checked against closed forms.
//
// Two paths. The seed row at l = |n| comes from an exact binomial row, and
// the boundary orders m = +-l at every higher degree come from a one-term
// recursion in l. The closed forms WignerMinOrder and WignerMaxOrder in
// WignerDetails are the definition the recursions are answerable to, and
// this is what they are kept for.
//
// Returns the largest relative discrepancy over the table. The comparison is
// relative above a floor, because the boundary values fall off like
// sin(theta/2)^{2l} and most of a large table is far below anything a relative
// comparison means something for; there all that is asked is that both paths
// be negligible.
constexpr std::ptrdiff_t CheckWignerBoundaryDegree = 128;

// What the comparison should be held to.
//
// The discrepancy grows linearly in the degree, and it is the closed form that
// is drifting rather than the recursion: the closed form exponentiates a
// logarithm of size O(l log 4), so it loses bits in proportion to the degree,
// while each recursion step is a multiplication by a factor of order one. Two
// pieces of evidence for reading it that way. The measured worst case is very
// nearly the same multiple of epsilon in double as in long double -- about
// 1000 -- so it tracks the precision rather than sitting at a fixed absolute
// size, and it scales with lMax rather than with the number of entries
// compared.
template <GSHTrans::RealFloatingPoint Real>
constexpr auto CheckWignerBoundaryTolerance() {
  return 20 * static_cast<Real>(CheckWignerBoundaryDegree) *
         std::numeric_limits<Real>::epsilon();
}

template <GSHTrans::RealFloatingPoint Real>
auto CheckWignerBoundary() {
  using namespace GSHTrans;
  using Int = std::ptrdiff_t;

  constexpr Int lMax = CheckWignerBoundaryDegree;
  constexpr Int nMax = 3;

  // Angles either side of pi/2, and two close to the poles, where the
  // recursion's factor is smallest and its seed is nearest the exact 0 or 1
  // that the closed form special-cases.
  const auto angles =
      std::vector<Real>{static_cast<Real>(0.02), static_cast<Real>(0.4),
                        static_cast<Real>(1.0),  std::numbers::pi_v<Real> / 2,
                        static_cast<Real>(2.3),  static_cast<Real>(3.12)};

  auto wigner = Wigner<Real, All, All, Multiple>(lMax, lMax, nMax, angles);

  const auto floor = static_cast<Real>(1e-100);
  auto worst = static_cast<Real>(0);
  auto Compare = [&worst, floor](Real actual, Real expected) {
    const auto scale = std::max(std::abs(expected), floor);
    worst = std::max(worst, std::abs(actual - expected) / scale);
  };

  for (auto n : wigner.UpperIndices()) {
    for (auto iTheta : wigner.AngleIndices()) {
      auto arg = WignerDetails::Arguments<Real>(angles[iTheta]);
      auto d = wigner[n, iTheta];

      for (auto l : d.Degrees()) {
        // The stored values carry the orthonormalisation, formed here exactly
        // as Compute forms it so that the comparison is of the recursion and
        // not of two spellings of the same constant.
        const auto factor = std::numbers::inv_sqrtpi_v<Real> /
                            static_cast<Real>(2) *
                            std::sqrt(static_cast<Real>(2 * l + 1));

        if (l == std::abs(n)) {
          // The seed row, over every order stored at that degree.
          for (auto m : d[l].Orders()) {
            const auto expected =
                n >= 0 ? WignerDetails::WignerMaxUpperIndex(l, m, arg)
                       : WignerDetails::WignerMinUpperIndex(l, m, arg);
            Compare(d[l][m], factor * expected);
          }
        } else {
          Compare(d[l][-l], factor * WignerDetails::WignerMinOrder(l, n, arg));
          Compare(d[l][l], factor * WignerDetails::WignerMaxOrder(l, n, arg));
        }
      }
    }
  }

  return worst;
}
