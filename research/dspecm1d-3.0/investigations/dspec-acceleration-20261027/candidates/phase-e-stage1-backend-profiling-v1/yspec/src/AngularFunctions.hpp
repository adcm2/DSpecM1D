#pragma once

/**
 * @file AngularFunctions.hpp
 * @brief The angular functions that sum the response over degree and order.
 *
 * @details The angular functions are formed from the generalised Legendre
 * functions of GSHTrans: with @f$X^N_{lm}@f$ their orthonormalised values and
 * @f$\zeta = \sqrt{l(l+1)}@f$,
 *
 * @f[ \partial_\theta X_{lm} = \tfrac{\zeta}{2}(X^{-1}_{lm} - X^{+1}_{lm}),
 *     \qquad \frac{m}{\sin\theta} X_{lm} =
 *     -\tfrac{\zeta}{2}(X^{+1}_{lm} + X^{-1}_{lm}), @f]
 *
 * which are regular at both poles.
 */

#include <YSpec/Types.hpp>
#include <array>
#include <cstddef>
#include <vector>

namespace YSpec {

/**
 * @brief The angular functions at one point, for one degree and order.
 */
struct Angular {
  double value;     ///< @f$X_{lm}(\theta)@f$.
  double dTheta;    ///< @f$\partial_\theta X_{lm}(\theta)@f$.
  double mOverSin;  ///< @f$m X_{lm}(\theta) / \sin\theta@f$.
};

/**
 * @brief The angular functions at a set of colatitudes, for degrees 0 to
 * lMax and orders @f$|m| \le \min(l, 2)@f$ (zero for the others).
 */
class AngularFunctions {
 public:
  /**
   * @param lMax The largest degree.
   * @param theta The colatitudes, in @f$[0, \pi]@f$.
   */
  AngularFunctions(int lMax, std::vector<double> theta);

  /** @brief The largest degree. */
  int MaxDegree() const { return lMax_; }
  /** @brief The number of points. */
  Int Size() const { return count_; }

  /** @brief The functions at point @p i, degree @p l, order @p m. */
  const Angular& operator()(Int i, int l, int m) const {
    return values_[Index(i, l, m)];
  }

 private:
  int lMax_;
  Int count_;
  std::vector<Angular> values_;

  template <class Real>
  void Fill(const std::vector<double>& theta);

  // Degree outermost, so that a sum over the receivers at one degree reads
  // contiguous memory.
  std::size_t Index(Int i, int l, int m) const {
    return (static_cast<std::size_t>(l) * count_ + i) * 5 + (m + 2);
  }
  Angular& Value(Int i, int l, int m) { return values_[Index(i, l, m)]; }
};

/**
 * @brief The generalised Legendre functions @f$X^N_{lm}@f$ themselves, for
 * @f$|N| \le 2@f$, which give the strain, at a set of colatitudes, for
 * degrees 0 to lMax and orders @f$|m| \le \min(l, 2)@f$ (zero for the
 * others, and for @f$|N| > l@f$).
 *
 * @details With @f$a_N = \sqrt{(l+N)(l-N+1)}@f$ and
 * @f$b_N = \sqrt{(l-N)(l+N+1)}@f$ they satisfy
 *
 * @f[ \partial_\theta X^N = \tfrac12(a_N X^{N-1} - b_N X^{N+1}), \qquad
 *     \frac{m - N\cos\theta}{\sin\theta} X^N =
 *     -\tfrac12(a_N X^{N-1} + b_N X^{N+1}). @f]
 */
class SpinFunctions {
 public:
  /** @brief The values for N = -2, ..., 2. */
  using Values = std::array<double, 5>;

  /**
   * @param lMax The largest degree.
   * @param theta The colatitudes, in @f$[0, \pi]@f$.
   */
  SpinFunctions(int lMax, std::vector<double> theta);

  /** @brief The largest degree. */
  int MaxDegree() const { return lMax_; }
  /** @brief The number of points. */
  Int Size() const { return count_; }

  /** @brief The functions at point @p i, degree @p l, order @p m. */
  const Values& operator()(Int i, int l, int m) const {
    return values_[Index(i, l, m)];
  }

 private:
  int lMax_;
  Int count_;
  std::vector<Values> values_;

  template <class Real>
  void Fill(const std::vector<double>& theta);

  std::size_t Index(Int i, int l, int m) const {
    return (static_cast<std::size_t>(l) * count_ + i) * 5 + (m + 2);
  }
};

}  // namespace YSpec
