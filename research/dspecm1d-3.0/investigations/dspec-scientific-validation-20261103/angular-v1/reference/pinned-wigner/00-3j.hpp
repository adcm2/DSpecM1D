#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "NumericConcepts/Numeric.hpp"
#include "NumericConcepts/Ranges.hpp"

/**
 * @file 3j.hpp
 * @brief Wigner 3-j symbols by stable recursion over the @f$(m_1, m_3)@f$
 * plane.
 *
 * @details For fixed degrees (l1, l2, l3) the full table of symbols
 *
 *           ( l1  l2  l3 )
 *           ( m1  m2  m3 ),   with m2 = -(m1 + m3),
 *
 * is generated over the axes m1 in [-l1, l1] and m3 in [-l3, l3], one row at
 * a time, by the Schulten-Gordon algorithm.
 *
 * The order recurrence has a growing and a decaying solution, with a
 * classically allowed region between two forbidden ones. Recursing inward
 * from a forbidden end follows the growing solution and is stable; outward
 * follows the decaying one and loses the answer exponentially. So each row is
 * built inward from *both* ends, the two halves are matched where they
 * overlap, and the one remaining constant is fixed from the unitary property
 *
 *     (2 l1 + 1) sum_{m2} g(m2)^2 = 1.
 *
 * No closed-form seed is needed and no factorial is ever formed, so the
 * method holds at large degree: measured against an independent route
 * (column-permutation invariance, which runs the recursion along different
 * lines) it agrees to 1e-16 at (200,200,200) and 1e-15 at (1000,1000,1999).
 *
 * Recursing in one direction only -- from a closed-form seed at a corner --
 * is the obvious alternative and is not safe: it fails exponentially near
 * *stretched* triangles, where one degree approaches the sum of the other
 * two, which is the top of every coupling sum rather than an exotic corner.
 *
 * A second layout is offered beside the plain table -- see CouplingElement()
 * and FillCouplingMatrix -- in which the first order is negated and carries an
 * alternating phase. That is the arrangement in which the symbols appear in
 * normal-mode coupling matrices, and it is what Woodhouse's wig2.f tabulated,
 * so codes written against that routine want it. It is a phase and a
 * relabelling of the table below, not a separate algorithm.
 *
 * All degrees and orders are integers; there is no half-integer support.
 *
 * @see Schulten, K. and Gordon, R. G., J. Math. Phys. 16 (1975) 1961.
 */

namespace GSHTrans {

//--------------------------------------------------------------------------//
//                            Concepts and helpers                           //
//--------------------------------------------------------------------------//

/**
 * @brief Concept for a contiguous, sized range of real numbers that can be
 * read and written. Used for the in-place kernels that fill user-provided
 * storage.
 * @tparam R The range type to check.
 */
template <typename R>
concept RealContiguousWritableRange =
    NumericConcepts::RealRange<R> and NumericConcepts::RealWritableRange<R> and
    std::ranges::contiguous_range<R> and std::ranges::sized_range<R>;

/**
 * @brief Checks the triangle inequality |l1 - l3| <= l2 <= l1 + l3.
 * @param l1 The first degree.
 * @param l2 The second degree.
 * @param l3 The third degree.
 * @return true if (l1, l2, l3) can form a triangle.
 */
constexpr bool SatisfiesTriangle(int l1, int l2, int l3) {
  const auto ldel = l1 > l3 ? l1 - l3 : l3 - l1;
  return l2 >= ldel and l2 <= l1 + l3;
}

/**
 * @brief A closed integer interval [Min, Max] used to label an index axis
 * of a table of Wigner symbols (e.g. an order axis m in [-l, l]).
 * @details An Axis is itself a range over its values, so it can be used
 * directly in range-based for loops:
 * @code
 *   for (auto m1 : table.M1Axis()) { ... }
 * @endcode
 */
class Axis {
 public:
  /** @brief Constructs an empty axis. */
  constexpr Axis() = default;

  /**
   * @brief Constructs the axis [min, max].
   * @param min The smallest index on the axis.
   * @param max The largest index on the axis.
   */
  constexpr Axis(int min, int max) : min_{min}, max_{max} {}

  /** @brief Returns the smallest index on the axis. */
  constexpr auto Min() const { return min_; }

  /** @brief Returns the largest index on the axis. */
  constexpr auto Max() const { return max_; }

  /** @brief Returns the number of indices on the axis. */
  constexpr auto Size() const {
    return max_ >= min_ ? static_cast<std::size_t>(max_ - min_ + 1)
                        : std::size_t{0};
  }

  /** @brief Returns true if the index i lies on the axis. */
  constexpr auto Contains(int i) const { return i >= min_ and i <= max_; }

  /** @brief Returns the zero-based offset of the index i along the axis. */
  constexpr auto Offset(int i) const {
    assert(Contains(i));
    return static_cast<std::size_t>(i - min_);
  }

  /** @brief Returns a view of the indices on the axis. */
  constexpr auto Values() const {
    return std::views::iota(min_, std::max(min_, max_ + 1));
  }

  /** @brief Iterator to the first index on the axis. */
  constexpr auto begin() const { return Values().begin(); }

  /** @brief Iterator past the last index on the axis. */
  constexpr auto end() const { return Values().end(); }

  /** @brief Axes compare equal when their bounds agree. */
  constexpr bool operator==(const Axis&) const = default;

 private:
  int min_{0};
  int max_{-1};
};

/*------------------------------------------------------------------------*/
/*                        The computational kernels                        */
/*------------------------------------------------------------------------*/

namespace ThreeJDetails {

/**
 * @brief The two coefficients of the order recurrence.
 * @details With g(m2) = (l1 l2 l3; m1 m2 -m1-m2), the symbols satisfy
 *
 *     A(m2) g(m2) + B(m2) g(m2-1) + A(m2-1) g(m2-2) = 0
 *
 * with A and B as below. A(m2min) is identically zero -- that is the lower
 * endpoint of the run -- which is why the backward pass takes its coefficient
 * one step up from the index it is filling and never forms it.
 */
template <NumericConcepts::Real T>
T RecurrenceA([[maybe_unused]] int l1, int l2, int l3, int m1, int m2) {
  const auto m3 = -m1 - m2;
  return std::sqrt(static_cast<T>(l2 - m2 + 1) * static_cast<T>(l2 + m2) *
                   static_cast<T>(l3 + m3 + 1) * static_cast<T>(l3 - m3));
}

template <NumericConcepts::Real T>
T RecurrenceB(int l1, int l2, int l3, int m1, int m2) {
  const auto m3 = -m1 - m2;
  return static_cast<T>(l1 + l2 + l3 + 1) * static_cast<T>(l2 + l3 - l1) -
         static_cast<T>(l2 - m2 + 1) * static_cast<T>(l3 + m3 + 1) -
         static_cast<T>(l2 + m2 - 1) * static_cast<T>(l3 - m3 - 1);
}

/** @brief The orders m2 for which g(m2) can be non-zero. */
inline std::pair<int, int> OrderRange(int l2, int l3, int m1) {
  return {std::max(-l2, -l3 - m1), std::min(l2, l3 - m1)};
}

/**
 * @brief One row of the (m1, m3) plane by the Schulten-Gordon algorithm.
 * @details Fills g[0 .. n-1] with g(m2) for m2 from m2min to m2max at fixed
 * m1, where n = m2max - m2min + 1.
 *
 * The recurrence has a growing and a decaying solution, and a classically
 * allowed region between two forbidden ones. Recursing *inward* from a
 * forbidden end follows the growing solution, so errors decay; outward
 * follows the decaying one, so they grow exponentially -- which is why a
 * one-directional scheme fails, and which extra precision delays rather than
 * cures.
 *
 * So the run comes inward from both ends and the two halves are matched where
 * they overlap. Each is correct up to its own scale; matching leaves one
 * unknown constant, and the unitary property
 *
 *     (2 l1 + 1) sum_{m2} g(m2)^2 = 1
 *
 * fixes it. No closed-form seed is needed and no factorial is ever formed.
 *
 * The forward run stops where the values first fall, and the backward run
 * comes up to meet it. Recursing towards larger values is what is stable, so
 * the quantity to watch is the values themselves -- and the first fall is the
 * edge of the allowed region where there is one, and the top of the hump in a
 * near-stretched row, which has none. Either way the three points the halves
 * are matched on sit where the row is largest, which is where a least-squares
 * match is best conditioned.
 *
 * SLATEC watches the recurrence coefficient |c1| instead and stops when it
 * first rises. That is a proxy, good in a row with an allowed region because
 * |c1| is smallest near the middle of one, and wrong in a near-stretched row,
 * where its minimum is nowhere near the hump: one half would then run
 * downhill, in the direction the unwanted solution grows. That costs digits
 * at best, and where the halves end more than 1e154 apart the join overflows
 * and the table cannot be built at all -- which is the case for more than
 * half of all (l, 2l, l) tables up to l = 1000.
 *
 * Schulten, K. and Gordon, R. G., J. Math. Phys. 16 (1975) 1961, and the
 * companion at 1971 for the semiclassical picture behind the turning points.
 * The control flow was checked against SLATEC's DRC3JM.
 */
template <NumericConcepts::Real T>
int SchultenGordonRow(int l1, int l2, int l3, int m1, std::span<T> g) {
  const auto [m2Min, m2Max] = OrderRange(l2, l3, m1);
  const auto n = m2Max - m2Min + 1;
  if (n < 1) return 0;

  std::fill_n(g.begin(), n, T{0});

  // One order: the recurrence has nowhere to run, and the value is the
  // convention's phase over the square root of the perimeter.
  if (n == 1) {
    g[0] = (std::abs(l2 - l3 - m1) % 2 == 0 ? T{1} : T{-1}) /
           std::sqrt(static_cast<T>(l1 + l2 + l3 + 1));
    return n;
  }

  const auto big = std::sqrt(std::numeric_limits<T>::max() / 20);
  const auto rootBig = std::sqrt(big);
  const auto small = 1 / big;
  const auto rootSmall = 1 / rootBig;

  const auto A = [&](int m2) { return RecurrenceA<T>(l1, l2, l3, m1, m2); };
  const auto B = [&](int m2) { return RecurrenceB<T>(l1, l2, l3, m1, m2); };

  // Forward from m2Min, while the values grow.
  g[0] = rootSmall;
  auto sumAll = small;
  auto sumForward = small;
  auto factor = T{0};
  auto c1 = T{0};
  auto value = T{0};
  auto last = 0;

  for (auto k = 1; k < n; ++k) {
    const auto m2 = m2Min + k;
    const auto older = factor;
    factor = A(m2);
    c1 = -B(m2) / factor;

    if (k == 1) {
      value = rootSmall * c1;  // the third term vanishes at the first step
    } else {
      value = c1 * g[k - 1] - (older / factor) * g[k - 2];
    }
    g[k] = value;
    sumForward = sumAll;
    sumAll += value * value;
    last = k;
    if (k == n - 1) break;

    if (std::abs(value) >= rootBig) {
      for (auto i = 0; i <= k; ++i) {
        if (std::abs(g[i]) < rootSmall) g[i] = 0;
        g[i] /= rootBig;
      }
      sumAll /= big;
      sumForward /= big;
      value /= rootBig;
    }
    // The first fall. Not before the third value: the match needs three
    // points, and a hump at the very start is two harmless steps downhill.
    if (k >= 2 && std::abs(g[k]) < std::abs(g[k - 1])) break;
  }

  auto rootTotal = std::sqrt(sumAll);

  // The sign the tail carries, tracked rather than read back from storage.
  //
  // The convention fixes the sign of g(m2Max), and SLATEC recovers it from
  // the computed array. That is not safe here: a near-stretched row at high
  // degree spans a dynamic range of 1e201, so the rescaling above flushes the
  // tail to zero and the sign with it, and whole rows come out negated with
  // every magnitude correct to rounding. Neither the completeness relation
  // nor the recurrence residual can see that.
  auto tailSign = (g[n - 1] >= 0) ? T{1} : T{-1};

  if (last < n - 1) {
    // Backward from m2Max, overlapping the forward run at three points.
    const auto x1 = g[last];
    const auto x2 = g[last - 1];
    const auto x3 = g[last - 2];

    g[n - 1] = rootSmall;
    auto sumBack = small;
    auto sumBackward = small;
    factor = 0;
    auto y1 = T{0}, y2 = T{0}, y3 = T{0};

    for (auto j = n - 2; j >= last - 2; --j) {
      const auto older = factor;
      factor = A(m2Min + j + 1);
      c1 = -B(m2Min + j + 2) / factor;

      const auto y = (j == n - 2) ? rootSmall * c1  // the third term vanishes
                                  : c1 * g[j + 1] - (older / factor) * g[j + 2];

      if (j == last - 2) {  // the match point: compare, do not store
        y3 = y;
        y2 = g[j + 1];
        y1 = g[j + 2];
        break;
      }

      g[j] = y;
      sumBackward = sumBack;
      sumBack += y * y;

      if (std::abs(y) >= rootBig) {
        for (auto i = j; i < n; ++i) {
          if (std::abs(g[i]) < rootSmall) g[i] = 0;
          g[i] /= rootBig;
        }
        sumBack /= big;
        sumBackward /= big;
      }
    }

    // Least squares over the three overlapping points, which is steadier than
    // matching on one of them.
    auto ratio = (x1 * y1 + x2 * y2 + x3 * y3) / (x1 * x1 + x2 * x2 + x3 * x3);

    // The root of the total and not the total. Each half was seeded at
    // rootSmall and rescaled as it grew, so the two can differ in scale by
    // about the whole range of T, and the ratio squared then overflows where
    // nothing else in the row does. hypot forms the root without the square.
    if (std::abs(ratio) >= 1) {
      for (auto i = 0; i <= last - 2; ++i) g[i] *= ratio;
      rootTotal =
          std::hypot(ratio * std::sqrt(sumForward), std::sqrt(sumBackward));
      tailSign = 1;  // the backward seed was positive and is untouched
    } else {
      ratio = 1 / ratio;
      for (auto i = last - 1; i < n; ++i) g[i] *= ratio;
      rootTotal =
          std::hypot(std::sqrt(sumForward), ratio * std::sqrt(sumBackward));
      tailSign = (ratio >= 0) ? T{1} : T{-1};
    }
  }

  auto norm = 1 / (std::sqrt(static_cast<T>(2 * l1 + 1)) * rootTotal);
  const auto wanted = (std::abs(l2 - l3 - m1) % 2 == 0) ? T{1} : T{-1};
  if (tailSign * wanted < 0) norm = -norm;
  for (auto i = 0; i < n; ++i) g[i] *= norm;
  return n;
}

/**
 * @brief The tolerance the recurrence residual is judged against.
 * @details Relative to the row's own largest value and to the size of the
 * recurrence's coefficients, and proportional to the number of steps, since
 * the rounding of a recursion accumulates along it.
 *
 * The constant is measured and not chosen. Over the stretched sweeps to
 * l = 1000, fat and flat tables and four hundred random triples, the largest
 * residual a row leaves is 1.9 of these units in double, 2.6 in long double
 * and 2.8 in single precision, the worst always on a flat row such as
 * (0, l, l). Fifty is that with room for another compiler's arithmetic.
 *
 * A much looser constant would not do: a thousand, for a long row in single
 * precision, comes to a tenth of the row's scale, and lets through tables
 * wrong by a third of theirs. A broken match shows up at order one in these
 * units times the row length, so there is a great deal of room between the
 * two and no reason to sit at the far end of it.
 */
template <NumericConcepts::Real T>
T ResidualTolerance(int steps) {
  return static_cast<T>(50 * (steps + 1)) * std::numeric_limits<T>::epsilon();
}

/**
 * @brief Checks a completed row against the recurrence that defines it.
 * @details This is the runtime self-check. The completeness relation, the
 * natural check for a one-directional recursion seeded from a closed form,
 * is nearly worthless against Schulten-Gordon, which normalises every row by
 * that very identity, so the sum is one by construction whatever the row
 * holds -- and in particular a mis-scaled join, which is this algorithm's
 * characteristic failure, passes it.
 *
 * The defining recurrence does not pass it. A bad match violates the relation
 * at the join, a bad rescale violates it locally, and a NaN propagates.
 *
 * What it cannot see, since the recurrence is homogeneous: an overall factor
 * or sign on the row. Scale is fixed by the normalisation and sign by the
 * convention, and both are properties of the code rather than of the input,
 * so both are pinned by tests instead.
 */
template <NumericConcepts::Real T>
bool RowSatisfiesRecurrence(int l1, int l2, int l3, int m1,
                            std::span<const T> g, int n) {
  if (n < 3) return true;
  const auto [m2Min, m2Max] = OrderRange(l2, l3, m1);
  (void)m2Max;

  auto scale = T{0};
  for (auto k = 0; k < n; ++k) scale = std::max(scale, std::abs(g[k]));
  if (not(scale > 0)) return false;

  const auto tolerance = ResidualTolerance<T>(n) * scale;
  for (auto k = 2; k < n; ++k) {
    const auto m2 = m2Min + k;
    const auto residual = RecurrenceA<T>(l1, l2, l3, m1, m2) * g[k] +
                          RecurrenceB<T>(l1, l2, l3, m1, m2) * g[k - 1] +
                          RecurrenceA<T>(l1, l2, l3, m1, m2 - 1) * g[k - 2];
    // The coefficients are of order l^2, so the residual is measured against
    // the scale they multiply rather than against the values alone.
    const auto local =
        std::max({std::abs(RecurrenceA<T>(l1, l2, l3, m1, m2)),
                  std::abs(RecurrenceB<T>(l1, l2, l3, m1, m2)),
                  std::abs(RecurrenceA<T>(l1, l2, l3, m1, m2 - 1))});
    if (not(std::abs(residual) <= tolerance * local)) return false;
  }
  return true;
}

/**
 * @brief Converts between the plain (m1, m3) table and the coupling layout.
 * @details The coupling layout is
 *
 *     c(m, mp) = (-1)^m ( l1     l2      l3 )
 *                       ( -m   m - mp    mp ),
 *
 * so plain(m1, m3) = (-1)^m1 c(-m1, m3): the rows are mirrored in m1 and
 * carry an alternating phase. It is the array Woodhouse's wig2.f returned.
 *
 * **It is its own inverse**, since mirroring twice is the identity and the
 * phase squares to one, so one routine serves both directions and there is no
 * second convention to keep in step.
 */
template <NumericConcepts::Real T>
void SwapCouplingConvention(const int l1, const int l3, std::span<T> a) {
  const auto columns = static_cast<std::size_t>(2 * l3 + 1);
  for (auto r = 0; r < l1; ++r) {
    const auto phase = ((r + l1) % 2 == 0) ? T{1} : T{-1};
    auto* lower = a.data() + static_cast<std::size_t>(r) * columns;
    auto* upper = a.data() + static_cast<std::size_t>(2 * l1 - r) * columns;
    for (std::size_t j = 0; j < columns; ++j) {
      const auto tmp = lower[j];
      lower[j] = phase * upper[j];
      upper[j] = phase * tmp;
    }
  }
}

/**
 * @brief Fills the whole (m1, m3) plane, row by row.
 * @details One Schulten-Gordon row recursion per m1, scattered into the
 * row-major (m1, m3) table. Each row is normalised independently, which is
 * stronger than one global normalisation would be: an error in one row cannot
 * leak into another.
 */
template <NumericConcepts::Real T>
void Wigner3jPlane(int l1, int l2, int l3, std::span<T> table) {
  // Every public route to a table comes through here, so this is where the
  // arguments are checked, and they are thrown on rather than asserted: the
  // rows are scattered into the table by index, so a table that is too small
  // is a write past the end of the caller's storage.
  if (l1 < 0 || l2 < 0 || l3 < 0) {
    throw std::invalid_argument(
        "The degrees of a Wigner 3-j symbol cannot be negative, and these "
        "are (" +
        std::to_string(l1) + ", " + std::to_string(l2) + ", " +
        std::to_string(l3) + ")");
  }
  const auto needed = static_cast<std::size_t>(2 * l1 + 1) *
                      static_cast<std::size_t>(2 * l3 + 1);
  if (table.size() < needed) {
    throw std::invalid_argument(
        "The table of 3-j symbols for l1 = " + std::to_string(l1) +
        " and l3 = " + std::to_string(l3) + " has (2 l1 + 1)(2 l3 + 1) = " +
        std::to_string(needed) + " entries, and the storage given holds " +
        std::to_string(table.size()));
  }
  std::fill(table.begin(), table.end(), T{0});
  if (not SatisfiesTriangle(l1, l2, l3)) return;

  // The precision the rows are *computed* in, which is at least double
  // whatever they are stored in. A three-term recurrence accumulates about
  // n^2 epsilon along a row of n values, and in single precision that is two
  // or three digits left at l = 450 -- a flat row such as l1 = 0 shows it
  // plainly. Single precision is worth having for what is stored and moved,
  // and buys nothing in a recursion that runs once, so a row is computed in
  // double and narrowed as it is scattered. The row is a scratch buffer in
  // any case, so this costs the difference in its arithmetic and no memory.
  using Work = std::conditional_t<(std::numeric_limits<T>::digits <
                                   std::numeric_limits<double>::digits),
                                  double, T>;

  const auto columns = 2 * l3 + 1;
  auto row = std::vector<Work>(static_cast<std::size_t>(2 * l2 + 2));

  for (auto m1 = -l1; m1 <= l1; ++m1) {
    const auto [m2Min, m2Max] = OrderRange(l2, l3, m1);
    if (m2Max < m2Min) continue;
    const auto n =
        SchultenGordonRow<Work>(l1, l2, l3, m1, std::span<Work>(row));

    if (not RowSatisfiesRecurrence<Work>(
            l1, l2, l3, m1, std::span<const Work>(row.data(), n), n)) {
      throw std::runtime_error(
          "Wigner 3-j symbols: the recurrence is not satisfied for degrees (" +
          std::to_string(l1) + ", " + std::to_string(l2) + ", " +
          std::to_string(l3) + ") at m1 = " + std::to_string(m1) +
          ". The values are not to be trusted.");
    }

    for (auto k = 0; k < n; ++k) {
      const auto m3 = -m1 - (m2Min + k);
      if (std::abs(m3) > l3) continue;
      table[static_cast<std::size_t>(m1 + l1) * columns + (m3 + l3)] =
          static_cast<T>(row[k]);
    }
  }
}

}  // namespace ThreeJDetails

/*------------------------------------------------------------------------*/
/*                         Free kernel functions                           */
/*------------------------------------------------------------------------*/

/**
 * @brief Fills user-provided storage with the table of Wigner 3j symbols
 * (l1 l2 l3; m1, -(m1+m3), m3) in row-major order, m1 being the row axis
 * and m3 the column axis. Allocates nothing but one row of scratch.
 * @tparam Range A contiguous, sized, readable and writable range of reals.
 * @param l1 The first degree.
 * @param l2 The second degree.
 * @param l3 The third degree.
 * @param table Storage of at least (2 l1 + 1) x (2 l3 + 1) values.
 * @throws std::invalid_argument if a degree is negative or the storage is too
 * small; std::runtime_error if a row fails its self-check.
 */
template <RealContiguousWritableRange Range>
void FillWigner3jMatrix(int l1, int l2, int l3, Range&& table) {
  using T = std::ranges::range_value_t<Range>;
  auto span = std::span<T>(std::ranges::data(table), std::ranges::size(table));
  ThreeJDetails::Wigner3jPlane<T>(l1, l2, l3, span);
}

/**
 * @brief Fills user-provided storage with the coupling layout
 * (-1)^m (l1 l2 l3; -m, m-mp, mp) in row-major order, m being the row axis
 * and mp the column axis: the first order negated, with an alternating phase.
 * This is the arrangement the symbols appear in within normal-mode coupling
 * matrices, and element for element it is the array a(m+l1+1, mp+l3+1)
 * returned by the Fortran routine wig2(l1, l2, l3, a, id1) -- the memory
 * order differs, Fortran being column-major, so a buffer exchanged with
 * Fortran directly is the transpose. Allocates nothing but one row of
 * scratch.
 * @tparam Range A contiguous, sized, readable and writable range of reals.
 * @param l1 The first degree.
 * @param l2 The second degree.
 * @param l3 The third degree.
 * @param table Storage of at least (2 l1 + 1) x (2 l3 + 1) values.
 * @throws std::invalid_argument if a degree is negative or the storage is too
 * small; std::runtime_error if a row fails its self-check.
 */
template <RealContiguousWritableRange Range>
void FillCouplingMatrix(int l1, int l2, int l3, Range&& table) {
  using T = std::ranges::range_value_t<Range>;
  auto span = std::span<T>(std::ranges::data(table), std::ranges::size(table));
  ThreeJDetails::Wigner3jPlane<T>(l1, l2, l3, span);
  ThreeJDetails::SwapCouplingConvention<T>(l1, l3, span);
}

/*------------------------------------------------------------------------*/
/*                            Wigner3jMatrix                               */
/*------------------------------------------------------------------------*/

/**
 * @brief Table of Wigner 3j symbols for fixed degrees (l1, l2, l3) over
 * the order axes m1 in [-l1, l1] and m3 in [-l3, l3], with the middle
 * order fixed by the selection rule m2 = -(m1 + m3).
 * @details Storage is a single contiguous std::vector in row-major order
 * (m1 slowest, m3 fastest). Element access is by quantum numbers rather
 * than raw offsets; queries outside the axes return zero, which is the
 * mathematically consistent value of the symbol. The class models
 * std::ranges::random_access_range over its values, and therefore also
 * satisfies NumericConcepts::RealRange.
 * @tparam T The floating-point type used for storage. Rows are computed in
 * at least double precision and narrowed to T.
 */
template <NumericConcepts::Real T>
class Wigner3jMatrix {
 public:
  using value_type = T;  ///< The stored value type.

  /** @brief Default constructor: the trivial table (0 0 0; 0 0 0) = 1. */
  Wigner3jMatrix() : Wigner3jMatrix(0, 0, 0) {}

  /**
   * @brief Computes the table for degrees (l1, l2, l3).
   * @param l1 The first degree.
   * @param l2 The second degree.
   * @param l3 The third degree.
   */
  Wigner3jMatrix(int l1, int l2, int l3)
      : l1_{l1}, l2_{l2}, l3_{l3}, data_(TableSize(l1, l3)) {
    ThreeJDetails::Wigner3jPlane<T>(l1_, l2_, l3_, std::span<T>(data_));
  }

  /** @brief Returns the first degree. */
  auto L1() const { return l1_; }

  /** @brief Returns the second degree. */
  auto L2() const { return l2_; }

  /** @brief Returns the third degree. */
  auto L3() const { return l3_; }

  /** @brief Returns the degrees as an array {l1, l2, l3}. */
  auto Degrees() const { return std::array{l1_, l2_, l3_}; }

  /** @brief Returns the axis for the first order, m1 in [-l1, l1]. */
  auto M1Axis() const { return Axis(-l1_, l1_); }

  /** @brief Returns the axis for the third order, m3 in [-l3, l3]. */
  auto M3Axis() const { return Axis(-l3_, l3_); }

  /**
   * @brief Returns the symbol (l1 l2 l3; m1, -(m1+m3), m3).
   * @param m1 The first order.
   * @param m3 The third order.
   * @return The value of the symbol; zero if an order is out of range.
   */
  auto operator()(int m1, int m3) const {
    if (not(M1Axis().Contains(m1) and M3Axis().Contains(m3))) return T{0};
    return data_[M1Axis().Offset(m1) * M3Axis().Size() + M3Axis().Offset(m3)];
  }

  /**
   * @brief Returns the symbol (l1 l2 l3; m1 m2 m3) with all three orders
   * given explicitly.
   * @param m1 The first order.
   * @param m2 The second order.
   * @param m3 The third order.
   * @return The value of the symbol; zero if the selection rules
   * m1 + m2 + m3 = 0 and |mi| <= li are not met.
   */
  auto operator()(int m1, int m2, int m3) const {
    if (m1 + m2 + m3 != 0 or std::abs(m2) > l2_) return T{0};
    return (*this)(m1, m3);
  }

  /**
   * @brief Returns the symbol in the coupling layout,
   * (-1)^m (l1 l2 l3; -m, m-mp, mp): the first order negated, with an
   * alternating phase.
   * @details This is how the symbols enter a normal-mode coupling matrix, and
   * it is the quantity the Fortran routine wig2.f tabulated, so it is offered
   * for codes written against that convention. It is a relabelling of the
   * plain table rather than a different calculation.
   * @param m The row order, m in [-l1, l1].
   * @param mp The column order, mp in [-l3, l3].
   */
  auto CouplingElement(int m, int mp) const {
    const auto phase = (m % 2 == 0) ? T{1} : T{-1};
    return phase * (*this)(-m, mp);
  }

  /** @brief Returns a view of the row of symbols with fixed m1. */
  auto Row(int m1) const {
    return std::span<const T>(
        data_.data() + M1Axis().Offset(m1) * M3Axis().Size(), M3Axis().Size());
  }

  /** @brief Returns a flat view of the table, row-major in (m1, m3). */
  auto Data() const { return std::span<const T>(data_); }

  /** @brief Iterator to the start of the flattened table. */
  auto begin() const { return data_.cbegin(); }

  /** @brief Iterator past the end of the flattened table. */
  auto end() const { return data_.cend(); }

  /** @brief Returns the total number of stored values. */
  auto size() const { return data_.size(); }

 private:
  int l1_{0};
  int l2_{0};
  int l3_{0};
  std::vector<T> data_;

  // The storage is sized before the table is filled, so a negative degree has
  // to be met here: cast to a size it is an allocation of nearly everything,
  // and the complaint would be about memory.
  static std::size_t TableSize(int l1, int l3) {
    if (l1 < 0 || l3 < 0) {
      throw std::invalid_argument(
          "The degrees of a Wigner 3-j symbol cannot be negative");
    }
    return static_cast<std::size_t>(2 * l1 + 1) *
           static_cast<std::size_t>(2 * l3 + 1);
  }
};

/*------------------------------------------------------------------------*/
/*                             Wigner3jStack                               */
/*------------------------------------------------------------------------*/

/**
 * @brief Table of Wigner 3j symbols for fixed outer degrees (l1, l3) over
 * three index axes: the middle degree l2, and the orders m1 and m3, with
 * m2 = -(m1 + m3) fixed by the selection rule.
 * @details Stored as one Wigner3jMatrix per value of l2. By default l2
 * runs over the full triangle range [|l1 - l3|, l1 + l3], but any
 * non-negative range may be requested; matrices outside the triangle
 * range are identically zero. The middle degree is the degree axis because
 * that is the common pattern in applications (e.g. coupling
 * through a structural degree l2 in normal-mode seismology); other
 * arrangements follow from the column-permutation symmetries of the
 * symbols. The stack is itself a range over its matrices.
 * @tparam T The floating-point type used for storage, as for Wigner3jMatrix.
 */
template <NumericConcepts::Real T>
class Wigner3jStack {
 public:
  using value_type = T;  ///< The stored value type.

  /**
   * @brief Computes tables for all l2 in the triangle range
   * [|l1 - l3|, l1 + l3].
   * @param l1 The first degree.
   * @param l3 The third degree.
   */
  Wigner3jStack(int l1, int l3)
      : Wigner3jStack(l1, l3, std::abs(l1 - l3), l1 + l3) {}

  /**
   * @brief Computes tables for all l2 in [l2Min, l2Max].
   * @param l1 The first degree.
   * @param l3 The third degree.
   * @param l2Min The smallest middle degree.
   * @param l2Max The largest middle degree.
   */
  Wigner3jStack(int l1, int l3, int l2Min, int l2Max)
      : l1_{l1}, l3_{l3}, l2Axis_{l2Min, l2Max} {
    if (l1 < 0 || l3 < 0 || l2Min < 0 || l2Max < l2Min) {
      throw std::invalid_argument(
          "A stack of 3-j tables needs non-negative degrees and a middle "
          "degree running upwards, and was given l1 = " +
          std::to_string(l1) + ", l3 = " + std::to_string(l3) + ", l2 from " +
          std::to_string(l2Min) + " to " + std::to_string(l2Max));
    }
    matrices_.reserve(l2Axis_.Size());
    for (auto l2 : l2Axis_) matrices_.emplace_back(l1_, l2, l3_);
  }

  /** @brief Returns the first degree. */
  auto L1() const { return l1_; }

  /** @brief Returns the third degree. */
  auto L3() const { return l3_; }

  /** @brief Returns the axis for the middle degree l2. */
  auto L2Axis() const { return l2Axis_; }

  /** @brief Returns the axis for the first order, m1 in [-l1, l1]. */
  auto M1Axis() const { return Axis(-l1_, l1_); }

  /** @brief Returns the axis for the third order, m3 in [-l3, l3]. */
  auto M3Axis() const { return Axis(-l3_, l3_); }

  /**
   * @brief Returns the table of symbols for a given middle degree.
   * @param l2 The middle degree; must lie on L2Axis().
   */
  const Wigner3jMatrix<T>& Matrix(int l2) const {
    assert(l2Axis_.Contains(l2));
    return matrices_[l2Axis_.Offset(l2)];
  }

  /** @brief Equivalent to Matrix(l2), allowing stack[l2](m1, m3). */
  const Wigner3jMatrix<T>& operator[](int l2) const { return Matrix(l2); }

  /**
   * @brief Returns the symbol (l1 l2 l3; m1, -(m1+m3), m3).
   * @param l2 The middle degree.
   * @param m1 The first order.
   * @param m3 The third order.
   * @return The value of the symbol; zero if l2 is off the stored axis or
   * an order is out of range.
   */
  auto operator()(int l2, int m1, int m3) const {
    if (not l2Axis_.Contains(l2)) return T{0};
    return Matrix(l2)(m1, m3);
  }

  /**
   * @brief Returns the symbol (l1 l2 l3; m1 m2 m3) with all three orders
   * given explicitly.
   * @param l2 The middle degree.
   * @param m1 The first order.
   * @param m2 The second order.
   * @param m3 The third order.
   */
  auto operator()(int l2, int m1, int m2, int m3) const {
    if (not l2Axis_.Contains(l2)) return T{0};
    return Matrix(l2)(m1, m2, m3);
  }

  /**
   * @brief Returns the symbol in the coupling layout
   * (-1)^m (l1 l2 l3; -m, m-mp, mp) for a given middle degree.
   * @param l2 The middle degree.
   * @param m The row order.
   * @param mp The column order.
   */
  auto CouplingElement(int l2, int m, int mp) const {
    if (not l2Axis_.Contains(l2)) return T{0};
    return Matrix(l2).CouplingElement(m, mp);
  }

  /** @brief Iterator to the first stored matrix (smallest l2). */
  auto begin() const { return matrices_.cbegin(); }

  /** @brief Iterator past the last stored matrix. */
  auto end() const { return matrices_.cend(); }

  /** @brief Returns the number of stored matrices. */
  auto size() const { return matrices_.size(); }

 private:
  int l1_{0};
  int l3_{0};
  Axis l2Axis_;
  std::vector<Wigner3jMatrix<T>> matrices_;
};

/*------------------------------------------------------------------------*/
/*                          Single-symbol helper                           */
/*------------------------------------------------------------------------*/

/**
 * @brief Convenience function returning a single symbol
 * (l1 l2 l3; m1 m2 m3).
 * @details Builds the full (m1, m3) table internally, at a cost of
 * O((2 l1 + 1)(2 l3 + 1)); when many symbols with shared degrees are
 * needed, construct a Wigner3jMatrix or Wigner3jStack once and query it
 * instead.
 * @tparam T The floating-point type; defaults to double.
 */
template <NumericConcepts::Real T = double>
T Wigner3jSymbol(int l1, int l2, int l3, int m1, int m2, int m3) {
  if (m1 + m2 + m3 != 0) return T{0};
  if (std::abs(m1) > l1 or std::abs(m2) > l2 or std::abs(m3) > l3) return T{0};
  if (not SatisfiesTriangle(l1, l2, l3)) return T{0};
  return Wigner3jMatrix<T>(l1, l2, l3)(m1, m3);
}

/*------------------------------------------------------------------------*/

static_assert(NumericConcepts::RealRange<Wigner3jMatrix<double>>);
static_assert(NumericConcepts::RealRange<Wigner3jMatrix<float>>);

}  // namespace GSHTrans
