#pragma once

/**
 * @file RungeKutta.hpp
 * @brief Explicit Runge-Kutta integration of small complex linear systems,
 * with step-size control and a running logarithmic scale.
 *
 * @details The radial equations are linear and homogeneous, and their
 * solutions grow or decay exponentially through evanescent regions by far
 * more than a double can hold. A solution is therefore carried as a vector of
 * unit norm and the logarithm of the factor it has been divided by
 * (ScaledVector). Quantities formed from several solutions combine their
 * logarithms, and only a physical result is ever exponentiated.
 *
 * The systems must be linear in the solution: normalising the solution after a
 * step rescales the stage that a first-same-as-last method carries over.
 *
 * A method is a Butcher tableau with an embedded lower-order solution, given
 * as a type (DormandPrince54, CashKarp54). With step-size control on, the
 * local error is measured relative to the norm of the solution, which is the
 * natural scale for a linear system whose overall size is arbitrary. With it
 * off, an interval is divided into equal steps no longer than a given length,
 * which is how the original YSpec integrates.
 */

#include <YSpec/Types.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <concepts>
#include <cstddef>
#include <limits>
#include <numeric>
#include <ranges>
#include <stdexcept>

namespace YSpec {

/** @brief The state of a system of N complex equations. */
template <std::size_t N>
using Vector = std::array<Complex, N>;

/** @brief The Euclidean norm of a complex vector. */
template <std::size_t N>
double Norm(const Vector<N>& y) {
  auto sum = 0.0;
  for (const auto& v : y) sum += std::norm(v);
  return std::sqrt(sum);
}

/**
 * @brief @f$v e^{s}@f$, formed so that it overflows or underflows only if the
 * result does.
 */
inline Complex Rescale(Complex v, double s) {
  const auto magnitude = std::abs(v);
  if (magnitude == 0 || !std::isfinite(magnitude)) return v * std::exp(s);
  return v / magnitude * std::exp(std::log(magnitude) + s);
}

/**
 * @brief A complex number held as a value times @f$e^{s}@f$.
 */
struct ScaledComplex {
  Complex value;        ///< The mantissa.
  double logScale = 0;  ///< The logarithm of the scale.

  /** @brief The number itself; finite only if it is representable. */
  Complex Get() const { return Rescale(value, logScale); }
};

/** @brief The product of two scaled numbers. */
inline ScaledComplex operator*(const ScaledComplex& a, const ScaledComplex& b) {
  return {a.value * b.value, a.logScale + b.logScale};
}
/** @brief The quotient of two scaled numbers. */
inline ScaledComplex operator/(const ScaledComplex& a, const ScaledComplex& b) {
  return {a.value / b.value, a.logScale - b.logScale};
}

/**
 * @brief A complex vector held as a value times @f$e^{s}@f$.
 */
template <std::size_t N>
struct ScaledVector {
  Vector<N> value{};    ///< The mantissa, of unit norm after normalisation.
  double logScale = 0;  ///< The logarithm of the scale.

  /** @brief Divides the mantissa by its norm, moving it into the scale. */
  void Normalise() {
    const auto norm = Norm(value);
    if (norm == 0 || !std::isfinite(norm)) return;
    for (auto& v : value) v /= norm;
    logScale += std::log(norm);
  }

  /** @brief A component, with the vector's scale. */
  ScaledComplex operator[](std::size_t i) const { return {value[i], logScale}; }
};

/**
 * @brief A Butcher tableau with an embedded error estimate.
 */
template <class T>
concept ButcherTableau = requires {
  { T::Stages } -> std::convertible_to<std::size_t>;
  { T::ErrorOrder } -> std::convertible_to<int>;
  { T::FirstSameAsLast } -> std::convertible_to<bool>;
  T::c;
  T::a;
  T::b;
  T::e;
};

/**
 * @brief The fifth-order method of Dormand & Prince (1980) with its embedded
 * fourth-order error estimate. The last stage is the first of the next step.
 */
struct DormandPrince54 {
  /// The number of stages.
  static constexpr std::size_t Stages = 7;
  /// The order of the embedded, lower-order solution.
  static constexpr int ErrorOrder = 4;
  /// Whether the last stage is the first of the next step.
  static constexpr bool FirstSameAsLast = true;
  /// The nodes.
  static constexpr std::array<double, Stages> c = {
      0.0, 1.0 / 5, 3.0 / 10, 4.0 / 5, 8.0 / 9, 1.0, 1.0};
  /// The Runge-Kutta matrix, by rows.
  static constexpr std::array<std::array<double, Stages>, Stages> a = {{
      {},
      {1.0 / 5},
      {3.0 / 40, 9.0 / 40},
      {44.0 / 45, -56.0 / 15, 32.0 / 9},
      {19372.0 / 6561, -25360.0 / 2187, 64448.0 / 6561, -212.0 / 729},
      {9017.0 / 3168, -355.0 / 33, 46732.0 / 5247, 49.0 / 176, -5103.0 / 18656},
      {35.0 / 384, 0.0, 500.0 / 1113, 125.0 / 192, -2187.0 / 6784, 11.0 / 84},
  }};
  /// The weights of the higher-order solution.
  static constexpr std::array<double, Stages> b = {
      35.0 / 384,     0.0,       500.0 / 1113, 125.0 / 192,
      -2187.0 / 6784, 11.0 / 84, 0.0};
  /// The fifth-order weights less the fourth-order ones.
  static constexpr std::array<double, Stages> e = {
      71.0 / 57600,      0.0,        -71.0 / 16695, 71.0 / 1920,
      -17253.0 / 339200, 22.0 / 525, -1.0 / 40};
};

/**
 * @brief The fifth-order method of Cash & Karp (1990), with its embedded
 * fourth-order error estimate. The original YSpec uses it with fixed steps.
 */
struct CashKarp54 {
  /// The number of stages.
  static constexpr std::size_t Stages = 6;
  /// The order of the embedded, lower-order solution.
  static constexpr int ErrorOrder = 4;
  /// Whether the last stage is the first of the next step.
  static constexpr bool FirstSameAsLast = false;
  /// The nodes.
  static constexpr std::array<double, Stages> c = {0.0,     1.0 / 5, 3.0 / 10,
                                                   3.0 / 5, 1.0,     7.0 / 8};
  /// The Runge-Kutta matrix, by rows.
  static constexpr std::array<std::array<double, Stages>, Stages> a = {{
      {},
      {1.0 / 5},
      {3.0 / 40, 9.0 / 40},
      {3.0 / 10, -9.0 / 10, 6.0 / 5},
      {-11.0 / 54, 5.0 / 2, -70.0 / 27, 35.0 / 27},
      {1631.0 / 55296, 175.0 / 512, 575.0 / 13824, 44275.0 / 110592,
       253.0 / 4096},
  }};
  /// The weights of the higher-order solution.
  static constexpr std::array<double, Stages> b = {
      37.0 / 378, 0.0, 250.0 / 621, 125.0 / 594, 0.0, 512.0 / 1771};
  /// The fifth-order weights less the fourth-order ones.
  static constexpr std::array<double, Stages> e = {
      37.0 / 378 - 2825.0 / 27648,
      0.0,
      250.0 / 621 - 18575.0 / 48384,
      125.0 / 594 - 13525.0 / 55296,
      -277.0 / 14336,
      512.0 / 1771 - 1.0 / 4};
};

/**
 * @brief Settings of an integration.
 */
struct IntegrationOptions {
  /// Whether to control the step size; if not, steps are equal and no longer
  /// than the suggested step.
  bool adaptive = true;
  /// The local error allowed per step, relative to the norm of the solution.
  double tolerance = 1e-8;
  /// The largest number of steps over one interval.
  std::size_t maxSteps = 1000000;
};

/**
 * @brief A right-hand side @f$f(r, y)@f$ for a system of N equations.
 */
template <class F, std::size_t N>
concept RightHandSide =
    std::regular_invocable<F, double, const Vector<N>&> &&
    std::same_as<std::invoke_result_t<F, double, const Vector<N>&>, Vector<N>>;

/** @brief K solutions of one system of N equations. */
template <std::size_t N, std::size_t K>
using Vectors = std::array<Vector<N>, K>;

/**
 * @brief A right-hand side applied to K solutions of one system at once, so
 * that what they share -- the coefficients at a radius -- is formed once.
 */
template <class F, std::size_t N, std::size_t K>
concept BlockRightHandSide =
    std::regular_invocable<F, double, const Vectors<N, K>&> &&
    std::same_as<std::invoke_result_t<F, double, const Vectors<N, K>&>,
                 Vectors<N, K>>;

namespace RungeKuttaDetails {

// y += a x, for one solution or several.
template <std::size_t N>
void AddScaled(Vector<N>& y, double a, const Vector<N>& x) {
  for (auto i = std::size_t{0}; i < N; ++i) y[i] += a * x[i];
}
template <std::size_t N, std::size_t K>
void AddScaled(Vectors<N, K>& y, double a, const Vectors<N, K>& x) {
  for (auto k = std::size_t{0}; k < K; ++k) AddScaled(y[k], a, x[k]);
}

// The error of a step relative to the solution; for several solutions, the
// largest, each relative to its own size. Zero solutions stay zero and do not
// count.
template <std::size_t N>
double ErrorRatio(const Vector<N>& e, const Vector<N>& y0,
                  const Vector<N>& y1) {
  const auto norm = std::max(Norm(y0), Norm(y1));
  return norm > 0 ? Norm(e) / norm : 0.0;
}
template <std::size_t N, std::size_t K>
double ErrorRatio(const Vectors<N, K>& e, const Vectors<N, K>& y0,
                  const Vectors<N, K>& y1) {
  auto ratio = 0.0;
  for (auto k = std::size_t{0}; k < K; ++k) {
    ratio = std::max(ratio, ErrorRatio(e[k], y0[k], y1[k]));
  }
  return ratio;
}

// Divides the last stage of a step by the norms its solutions were divided
// by: the systems are linear.
template <std::size_t N>
void Divide(Vector<N>& k, double norm) {
  for (auto& v : k) v /= norm;
}
template <std::size_t N, std::size_t K>
void Divide(Vectors<N, K>& k, const std::array<double, K>& norms) {
  for (auto i = std::size_t{0}; i < K; ++i) Divide(k[i], norms[i]);
}

}  // namespace RungeKuttaDetails

/**
 * @brief An explicit Runge-Kutta integrator for one tableau.
 */
template <ButcherTableau Tableau = DormandPrince54>
class RungeKutta {
 public:
  /** @brief An integrator with the given settings. */
  explicit RungeKutta(IntegrationOptions options = {}) : options_{options} {}

  /** @brief The settings. */
  const IntegrationOptions& Options() const { return options_; }

  /**
   * @brief Integrates @f$y' = f(r, y)@f$ from @p r0 to @p r1, in either
   * direction.
   * @param f The right-hand side.
   * @param r0 The initial radius.
   * @param r1 The final radius.
   * @param y On entry the solution at @p r0, on exit at @p r1, normalised
   * after every step.
   * @param step On entry the suggested length of the first step (its sign is
   * ignored); on exit the length suggested for a following interval.
   * @throws std::runtime_error if the step count is exceeded or the step
   * size underflows.
   */
  template <std::size_t N, RightHandSide<N> F>
  void Integrate(F&& f, double r0, double r1, ScaledVector<N>& y,
                 double& step) const {
    // The systems are linear, so a zero solution stays zero.
    if (Norm(y.value) == 0) return;
    auto normalise = [&y](Vector<N>& value) {
      y.value = value;
      const auto norm = Norm(y.value);
      y.Normalise();
      value = y.value;
      return norm > 0 && std::isfinite(norm) ? norm : 1.0;
    };
    auto value = y.value;
    Run(f, r0, r1, value, step, normalise);
  }

  /**
   * @brief Integrates K solutions of one system together, with one sequence
   * of steps; the error of each is controlled relative to its own size, and
   * each keeps its own scale.
   * @param f The right-hand side for all K solutions at once.
   * @param r0 The initial radius.
   * @param r1 The final radius.
   * @param y On entry the solutions at @p r0, on exit at @p r1.
   * @param step As for Integrate.
   */
  template <std::size_t N, std::size_t K, BlockRightHandSide<N, K> F>
  void IntegrateTogether(F&& f, double r0, double r1,
                         std::array<ScaledVector<N>, K>& y,
                         double& step) const {
    auto value = Vectors<N, K>{};
    auto zero = true;
    for (auto k = std::size_t{0}; k < K; ++k) {
      value[k] = y[k].value;
      zero = zero && Norm(value[k]) == 0;
    }
    if (zero) return;
    auto normalise = [&y](Vectors<N, K>& v) {
      auto norms = std::array<double, K>{};
      for (auto k = std::size_t{0}; k < K; ++k) {
        y[k].value = v[k];
        norms[k] = Norm(v[k]);
        y[k].Normalise();
        v[k] = y[k].value;
        // A zero solution is left as it is.
        if (!(norms[k] > 0 && std::isfinite(norms[k]))) norms[k] = 1;
      }
      return norms;
    };
    Run(f, r0, r1, value, step, normalise);
  }

 private:
  IntegrationOptions options_;

  // The integration itself, for one solution (State a Vector) or several (a
  // Vectors). normalise(y) normalises the solution, records the factors in
  // the caller's scales and returns them.
  template <class State, class F, class Normalise>
  void Run(F& f, double r0, double r1, State& y, double& step,
           Normalise& normalise) const {
    using namespace RungeKuttaDetails;
    const auto length = r1 - r0;
    if (length == 0) return;
    const auto direction = length > 0 ? 1.0 : -1.0;
    step = std::abs(step);
    if (!(step > 0)) step = std::abs(length);

    if (!options_.adaptive) {
      // Equal steps, as many as needed to keep each below the suggestion.
      const auto steps = static_cast<std::size_t>(std::abs(length) / step) + 1;
      const auto h = length / static_cast<double>(steps);
      auto k0 = f(r0, y);
      for (auto i : std::views::iota(std::size_t{0}, steps)) {
        const auto r = r0 + static_cast<double>(i) * h;
        auto error = State{};
        Step(f, r, h, y, k0, error);
        normalise(y);
        k0 = f(r + h, y);
      }
      return;
    }

    auto r = r0;
    auto k0 = f(r, y);
    auto h = direction * std::min(step, std::abs(length));
    constexpr auto exponent = 1.0 / (Tableau::ErrorOrder + 1);
    for (auto count = std::size_t{0}; count < options_.maxSteps; ++count) {
      // The last step ends exactly at r1. The step it replaces is remembered,
      // since a step shortened to fit says little about the next interval.
      const auto remaining = r1 - r;
      const auto last = std::abs(h) >= std::abs(remaining);
      const auto proposed = h;
      if (last) h = remaining;

      auto trial = y;
      auto error = State{};
      auto kLast = Step(f, r, h, trial, k0, error);
      const auto ratio = ErrorRatio(error, y, trial) / options_.tolerance;

      if (ratio <= 1) {
        r = last ? r1 : r + h;
        y = trial;
        const auto norms = normalise(y);
        if constexpr (Tableau::FirstSameAsLast) {
          // The last stage was evaluated at the unnormalised solution.
          Divide(kLast, norms);
          k0 = kLast;
        } else {
          k0 = f(r, y);
        }
        const auto factor =
            ratio > 0 ? std::min(5.0, 0.9 * std::pow(ratio, -exponent)) : 5.0;
        if (last) {
          step = std::abs(h) * factor;
          if (std::abs(h) < std::abs(proposed)) {
            step = std::max(step, std::abs(proposed));
          }
          return;
        }
        h *= factor;
      } else {
        const auto factor =
            std::isfinite(ratio)
                ? std::max(0.2, 0.9 * std::pow(ratio, -exponent))
                : 0.2;
        h *= factor;
        if (std::abs(h) < 1e-14 * std::max(1.0, std::abs(r))) {
          throw std::runtime_error("RungeKutta: step size underflow");
        }
      }
    }
    throw std::runtime_error("RungeKutta: too many steps");
  }

  // Takes one step of length h from r, replacing y by the higher-order
  // solution and setting error to the error estimate. k0 is f(r, y). Returns
  // the last stage.
  template <class State, class F>
  static State Step(F& f, double r, double h, State& y, const State& k0,
                    State& error) {
    using RungeKuttaDetails::AddScaled;
    constexpr auto S = Tableau::Stages;
    auto k = std::array<State, S>{};
    k[0] = k0;
    for (auto s = std::size_t{1}; s < S; ++s) {
      auto ys = y;
      for (auto j = std::size_t{0}; j < s; ++j) {
        const auto a = Tableau::a[s][j];
        if (a != 0) AddScaled(ys, h * a, k[j]);
      }
      k[s] = f(r + Tableau::c[s] * h, ys);
    }
    error = State{};
    for (auto s = std::size_t{0}; s < S; ++s) {
      if (Tableau::b[s] != 0) AddScaled(y, h * Tableau::b[s], k[s]);
      if (Tableau::e[s] != 0) AddScaled(error, h * Tableau::e[s], k[s]);
    }
    return k[S - 1];
  }
};

}  // namespace YSpec
