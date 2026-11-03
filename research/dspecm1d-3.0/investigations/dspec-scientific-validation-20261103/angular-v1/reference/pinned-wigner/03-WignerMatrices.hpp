#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "Concepts.hpp"
#include "Indexing.hpp"
#include "OpenMP.hpp"
#include "Utility.hpp"
#include "Views.hpp"
#include "Wigner.hpp"

// The pragmas below are for a compiler that was asked for OpenMP. One that was
// not ignores them, which is right, and warns that it has, which in a
// header-only library is a warning in the caller's build. So the warning is
// turned off for this file alone, and only where OpenMP is off: see OpenMP.hpp.
#ifndef _OPENMP
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif

namespace GSHTrans {

/// The same d-function values as Wigner, laid out transform-major:
/// [n][m][l][theta] rather than [n][theta][(l, m)].
///
/// The matrix kernel writes the Legendre stage as one matrix product per
/// order,
///
///     f^n_{lm} = sum_i D^(n,m)_{li} b^(m)_i,
///     D^(n,m)_{li} = X^n_{lm}(theta_i),
///
/// and that needs D^(n,m) contiguous in (l, theta) at fixed (n, m). The table
/// Wigner builds is contiguous in (l, m) at fixed (n, theta), which is the
/// transpose of what is wanted. This class holds the other one.
///
/// In the full layout the total size is the same as Wigner's, which is worth
/// stating because it is not obvious: summing (lMax - max(|n|,|m|) + 1) over m
/// gives exactly the same count as summing (2 min(l, mMax) + 1) over l. The
/// same triangle, read the other way. At lMax = 256 and |n| = 2 that is 66,045
/// values per colatitude either way. The reflected layout (see Reflected())
/// stores the non-negative orders only, so roughly half.
///
/// -- Why a separate class, and not a layout flag on Wigner.
///
/// The two have unrelated interfaces. Wigner hands out a (l, m) block for one
/// (n, theta), through ConstGSHView; this hands out an (l, theta) matrix for
/// one (n, m), as a plain span a BLAS call can take. A single class serving
/// both would carry two interfaces and would have to branch in the accessor
/// the inner loop calls, which is the one place that cannot afford it. The
/// grid holds its Wigner table in an optional -- empty on a generating grid --
/// and holds this in a second optional beside it.
///
/// -- Why the values are generated rather than transposed.
///
/// Building a Wigner table and transposing it would need both live at once:
/// 1.3 GB at lMax = 256, nMax = 2, to end with 648 MB. So the recursion is run
/// directly into per-thread scratch, one (n, theta) block at a time, and
/// scattered into place. The recursion, the seeds, the evaluation order and
/// the orthonormalisation are Wigner's own -- WignerDetails::ComputeBlock is
/// called with the same arguments -- so the values are bit-identical to the
/// stored path by construction rather than by tolerance. The test checks
/// exactly that.
///
/// -- The write pattern, named because it is the one cost here.
///
/// A block computed for one (n, theta) scatters across every matrix, at stride
/// NumberOfAngles() in the degree. That is one cache line per value written in
/// the worst case. It is paid once, at construction and in parallel, where
/// the loop kernel's table is written in the order it is computed. Blocking
/// over colatitudes would remove it, and is not done.
template <RealFloatingPoint Real_, OrderIndexRange MRange_ = All,
          IndexRange NRange_ = All>
class WignerMatrices {
 public:
  using Int = std::ptrdiff_t;  ///< Signed index type used throughout.
  using Real = Real_;          ///< The precision.
  /// Whether all orders are stored, or only the non-negative ones.
  using MRange = MRange_;
  using NRange = NRange_;  ///< Which upper indices are covered.

  WignerMatrices() = default;

  /// Every order MRange admits.
  template <std::ranges::range Range>
  requires RealFloatingPoint<std::ranges::range_value_t<Range>>
  static auto Full(Int lMax, Int mMax, Int nMax, Range &&theta) {
    return WignerMatrices(lMax, mMax, nMax, theta, false);
  }

  /// Non-negative orders only, the negative ones recovered from the reflection
  /// (see Sign()). Requires colatitudes symmetric about pi/2, and checks it.
  template <std::ranges::range Range>
  requires RealFloatingPoint<std::ranges::range_value_t<Range>>
  static auto Reflected(Int lMax, Int mMax, Int nMax, Range &&theta) {
    return WignerMatrices(lMax, mMax, nMax, theta, true);
  }

  /**
   * @brief The whole table, computed at construction.
   * @param lMax The largest degree.
   * @param mMax The largest order.
   * @param nMax The largest upper index.
   * @param theta The colatitudes.
   * @param reflected Whether only the non-negative orders are stored, the
   * rest being reached through Sign(); see Reflected().
   * @throws std::invalid_argument if the degree is negative or above
   * MaxSafeDegree, if the order or upper index exceeds the degree, if a
   * colatitude lies outside @f$[0, \pi]@f$, or if @p reflected is asked of
   * colatitudes not symmetric about @f$\pi/2@f$.
   */
  template <std::ranges::random_access_range Range>
  requires std::ranges::sized_range<Range> &&
               RealFloatingPoint<std::ranges::range_value_t<Range>>
  WignerMatrices(Int lMax, Int mMax, Int nMax, Range &&theta,
                 bool reflected = false)
      : lMax_{lMax},
        mMax_{mMax},
        nMax_{nMax},
        nTheta_(static_cast<Int>(std::ranges::size(theta))),
        reflected_{reflected} {
    if (lMax < 0)
      throw std::invalid_argument("Maximum degree must not be negative");
    WignerDetails::CheckSafeDegree<Real>(lMax);
    if (mMax < 0 || mMax > lMax) {
      throw std::invalid_argument(
          "Maximum order must lie between zero and the maximum degree");
    }
    if (std::abs(nMax) > lMax) {
      throw std::invalid_argument(
          "Maximum upper index must not exceed the maximum degree");
    }

    for (auto angle : theta) {
      if (!(angle >= 0 && angle <= std::numbers::pi_v<Real>)) {
        throw std::invalid_argument("A colatitude must lie in [0, pi], and " +
                                    std::to_string(angle) + " does not");
      }
    }

    // The reflection is a statement about theta and pi - theta both being
    // sample points, so it is only available on a grid that has them. A
    // Gauss-Legendre grid does -- its nodes are symmetric to one ulp of pi and
    // its weights exactly -- but this class takes any angles at all, so it
    // asks rather than assumes. Getting this wrong would produce a table that
    // is quietly the wrong values for half the orders.
    if (reflected_) {
      const auto pi = std::numbers::pi_v<Real>;
      for (auto i = Int{0}; i < nTheta_; i++) {
        const auto mirror = theta[nTheta_ - 1 - i];
        if (std::abs(theta[i] + mirror - pi) >
            static_cast<Real>(64) * std::numeric_limits<Real>::epsilon() * pi) {
          throw std::invalid_argument(
              "The reflected layout stores non-negative orders only and "
              "recovers the rest from d^l_{nm}(pi - theta), so it needs "
              "colatitudes symmetric about pi / 2");
        }
      }
    }

    // One offset per (n, m), and the matrices laid end to end in that order.
    offset_.reserve(
        static_cast<std::size_t>(NumberOfUpperIndices() * NumberOfOrders()));
    auto size = std::size_t{0};
    for (auto n : UpperIndices()) {
      for (auto m : Orders()) {
        offset_.push_back(size);
        size += static_cast<std::size_t>(NumberOfDegrees(n, m) * nTheta_);
      }
    }
    data_ = std::vector<Real>(size);

    ComputeAll(theta);
  }

  // Degrees, orders, upper indices and angles. The upper-index accessors
  // match Wigner's exactly, since a grid hands both the same NRange.
  /** @brief The largest degree stored. */
  auto MaxDegree() const { return lMax_; }
  /** @brief The largest order stored. */
  auto MaxOrder() const { return mMax_; }

  /// Zero when reflected, whatever MRange: the negative orders are not
  /// stored and are reached through Sign() instead.
  auto MinOrder() const {
    if (reflected_) return Int{0};
    if constexpr (std::same_as<MRange, All>) {
      return -mMax_;
    } else {
      return Int{0};
    }
  }

  /** @brief Whether only the non-negative orders are stored. */
  auto IsReflected() const { return reflected_; }

  /// The reflection itself, as one function so that no caller writes the sign
  /// out by hand:
  ///
  ///     d^l_{nm}(pi - theta) = (-1)^{l+n} d^l_{n,-m}(theta)
  ///
  /// verified against this library's own values to 3.8e-15 on values of order
  /// one. Equivalently, the matrix at -m is the matrix at +m with its columns
  /// reversed and this sign applied to row l.
  static constexpr Real Sign(Int l, Int n) {
    return ((l + n) % 2 == 0) ? Real{1} : Real{-1};
  }

  /// Where the mirror of colatitude i lives.
  auto MirrorAngle(Int iTheta) const { return nTheta_ - 1 - iTheta; }

  /** @brief Every order stored. */
  auto Orders() const {
    return std::ranges::views::iota(MinOrder(), MaxOrder() + 1);
  }

  /** @brief How many orders are stored. */
  auto NumberOfOrders() const { return MaxOrder() - MinOrder() + 1; }

  /** @brief The smallest upper index covered. */
  auto MinUpperIndex() const {
    if constexpr (std::same_as<NRange, All>) {
      return -nMax_;
    } else if constexpr (std::same_as<NRange, NonNegative>) {
      return Int{0};
    } else {
      return nMax_;
    }
  }

  /** @brief The largest upper index covered. */
  auto MaxUpperIndex() const { return nMax_; }

  /** @brief Every upper index covered. */
  auto UpperIndices() const {
    return std::ranges::views::iota(MinUpperIndex(), MaxUpperIndex() + 1);
  }

  /** @brief How many upper indices are covered. */
  auto NumberOfUpperIndices() const {
    return MaxUpperIndex() - MinUpperIndex() + 1;
  }

  /** @brief How many colatitudes the table holds. */
  auto NumberOfAngles() const { return nTheta_; }

  /// The lowest degree present at (n, m), and how many there are.
  ///
  /// A d-function vanishes identically unless l >= |n| and l >= |m|, so the
  /// matrix at (n, m) starts at max(|n|, |m|) and its height falls linearly in
  /// |m|. That is the source of the load imbalance a parallel matrix kernel
  /// has to divide work for rather than count orders.
  auto MinDegree(Int n, Int m) const {
    assert(std::abs(n) <= nMax_);
    assert(m >= MinOrder() && m <= MaxOrder());
    return std::max(std::abs(n), std::abs(m));
  }

  /** @brief How many degrees are present at @p n and @p m. */
  auto NumberOfDegrees(Int n, Int m) const {
    return lMax_ - MinDegree(n, m) + 1;
  }

  /** @brief Every degree present at @p n and @p m. */
  auto Degrees(Int n, Int m) const {
    return std::ranges::views::iota(MinDegree(n, m), lMax_ + 1);
  }

  /// The matrix D^(n,m): NumberOfDegrees(n, m) rows by NumberOfAngles()
  /// columns, row-major, so the value at (l, iTheta) is at
  ///
  ///     (l - MinDegree(n, m)) * NumberOfAngles() + iTheta
  ///
  /// and the leading dimension a BLAS call wants is NumberOfAngles().
  ///
  /// The index is left to the caller rather than wrapped in an accessor,
  /// deliberately: the layout is the whole content of this class, and a test
  /// that computes the index itself is testing the layout rather than trusting
  /// a member that could be wrong in the same way twice.
  auto operator[](Int n, Int m) const {
    return std::span<const Real>(
        data_.data() + offset_[OffsetIndex(n, m)],
        static_cast<std::size_t>(NumberOfDegrees(n, m) * nTheta_));
  }

 private:
  Int lMax_{};
  Int mMax_{};
  Int nMax_{};
  Int nTheta_{};
  bool reflected_{false};

  std::vector<Real> data_;
  std::vector<std::size_t> offset_;

  auto OffsetIndex(Int n, Int m) const {
    assert(n >= MinUpperIndex() && n <= MaxUpperIndex());
    assert(m >= MinOrder() && m <= MaxOrder());
    return static_cast<std::size_t>((n - MinUpperIndex()) * NumberOfOrders() +
                                    (m - MinOrder()));
  }

  template <std::ranges::range Range>
  requires RealFloatingPoint<std::ranges::range_value_t<Range>>
  void ComputeAll(Range &&thetaRange) {
    // Recursed in at least double and narrowed at the scatter below; see
    // WignerRecursionReal. The block is computed into scratch in any case, so
    // for single precision the only cost is a wider scratch type.
    using Work = WignerRecursionReal<Real>;
    const auto [sqrtInt, sqrtIntInv] =
        WignerDetails::PreComputeRecursionTables<Real>(lMax_, mMax_, nMax_);
    const auto sqrtIntView = std::span<const Work>(sqrtInt);
    const auto sqrtIntInvView = std::span<const Work>(sqrtIntInv);

    // The largest (n, theta) block is the one at the smallest |n|, since a
    // block starts at degree |n|. One scratch buffer of that size per thread.
    const auto scratchSize = static_cast<std::size_t>(
        GSHIndices<MRange>(lMax_, mMax_, SmallestUpperIndexModulus()).Size());

    // Flattened to an integer loop and decoded inside, and the region
    // suppressed when one is already open, for the reasons Wigner::ComputeAll
    // gives: OpenMP's canonical loop form, and exactly one level threads.
    const auto count = NumberOfUpperIndices() * nTheta_;
    const auto minUpperIndex = MinUpperIndex();
    const auto nTheta = nTheta_;
    [[maybe_unused]] const bool nested = OpenMP::InParallel();

    // The per-thread scratch is allocated in here, which can throw, and
    // nothing may leave a region: see ExceptionCapture, including for why a
    // thread whose allocation failed still goes through the loop.
    auto capture = Details::ExceptionCapture{};
#pragma omp parallel if (!nested)
    {
      auto scratch = std::vector<Work>{};
      capture.Run([&] { scratch.resize(scratchSize); });

#pragma omp for schedule(static)
      for (Int index = 0; index < count; index++) {
        if (capture.Failed()) continue;
        capture.Run([&] {
          const auto n = minUpperIndex + index / nTheta;
          const auto iTheta = index % nTheta;

          // Exactly the call Wigner makes, into scratch instead of into place.
          // Orthonormalisation happens inside it, so what is scattered below is
          // already the stored value sqrt((2l+1)/(4 pi)) d^l_{nm}.
          auto indices = GSHIndices<MRange>(lMax_, mMax_, n);
          WignerDetails::ComputeBlock(
              GSHView<Work, MRange>(lMax_, mMax_, n, scratch.data()), n,
              static_cast<Work>(thetaRange[iTheta]), sqrtIntView,
              sqrtIntInvView);

          // Scatter: (l, m) in the block goes to column iTheta of matrix (n,
          // m).
          const auto lowest = MinOrder();
          for (auto l : indices.Degrees()) {
            auto [blockOffset, sub] = indices.Index(l);
            for (auto m : sub.Orders()) {
              if (m < lowest) continue;
              data_[offset_[OffsetIndex(n, m)] +
                    static_cast<std::size_t>((l - MinDegree(n, m)) * nTheta +
                                             iTheta)] =
                  static_cast<Real>(scratch[blockOffset + sub.Index(m)]);
            }
          }
        });
      }
    }
    capture.Rethrow();
  }

  // The upper index of smallest modulus this table carries, which is the one
  // whose (n, theta) block is largest.
  auto SmallestUpperIndexModulus() const {
    if (MinUpperIndex() <= 0 && MaxUpperIndex() >= 0) return Int{0};
    return MinUpperIndex() > 0 ? MinUpperIndex() : MaxUpperIndex();
  }
};

}  // namespace GSHTrans

#ifndef _OPENMP
#pragma GCC diagnostic pop
#endif
