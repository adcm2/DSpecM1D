#pragma once

/**
 * @file OrderArray.hpp
 * @brief A value for each of the orders a point source at the pole excites.
 */

#include <algorithm>
#include <array>
#include <ranges>

namespace YSpec {

/**
 * @brief A value for each order @f$m = -2, \dots, 2@f$, indexed by @f$m@f$.
 */
template <class T>
class OrderArray {
 public:
  /// The largest order excited by a point source at the pole.
  static constexpr int MaxOrder = 2;

  /** @brief The value at order @p m. */
  T& operator[](int m) { return values_[m + MaxOrder]; }
  /** @brief The value at order @p m. */
  const T& operator[](int m) const { return values_[m + MaxOrder]; }

  /** @brief The orders that exist at degree @p l: @f$|m| \le \min(l, 2)@f$. */
  static auto Orders(int l) {
    const auto mMax = std::min(l, MaxOrder);
    return std::views::iota(-mMax, mMax + 1);
  }

 private:
  std::array<T, 2 * MaxOrder + 1> values_{};
};

}  // namespace YSpec
