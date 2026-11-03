#pragma once

/**
 * @file TimeSeries.hpp
 * @brief The frequency grid, and the transformation of spectra to time
 * series.
 *
 * @details The spectra are computed at complex frequencies
 * @f$\omega - i\epsilon@f$, which damps the response by @f$e^{-\epsilon t}@f$
 * so that the time series it is transformed to does not wrap around. The
 * damping is chosen to reduce the response by @f$e^{-5}@f$ over the length of
 * output asked for, the spacing of the frequencies follows from it, and the
 * number of samples is the next power of two. After the transformation the
 * damping is undone by multiplying by @f$e^{\epsilon t}@f$. This is the scheme
 * of the original YSpec.
 */

#include <YSpec/Calculation.hpp>
#include <YSpec/Types.hpp>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <numbers>
#include <optional>
#include <ranges>
#include <stdexcept>
#include <vector>

namespace YSpec {

/**
 * @brief The frequencies at which the spectra are computed, and the time
 * series they transform to. All quantities are non-dimensional.
 */
class FrequencyGrid {
 public:
  /**
   * @param fMin The lowest frequency wanted (cycles per unit time).
   * @param fMax The highest frequency wanted.
   * @param dt The time step.
   * @param length The length of the time series wanted.
   * @param decay The damping over that length, in e-folds.
   * @throws std::invalid_argument if the Nyquist frequency is below @p fMax
   * or the range holds fewer than two frequencies.
   */
  FrequencyGrid(double fMin, double fMax, double dt, double length,
                double decay = 5)
      : dt_{dt}, epsilon_{decay / length} {
    if (0.5 / dt < fMax) {
      throw std::invalid_argument(
          "the Nyquist frequency is below the highest frequency: reduce the "
          "time step");
    }
    // The spacing that resolves the damping, rounded down so that the number
    // of samples is a power of two.
    const auto spacing = epsilon_ / (2 * std::numbers::pi);
    const auto minimum = static_cast<std::size_t>(1 / (spacing * dt));
    samples_ = 1;
    while (samples_ <= minimum) samples_ *= 2;
    df_ = 1 / (static_cast<double>(samples_) * dt);
    first_ = std::max<std::size_t>(
                 static_cast<std::size_t>(std::floor(fMin / df_)), 2) -
             1;
    last_ = static_cast<std::size_t>(std::floor(fMax / df_)) + 1;
    if (last_ <= first_) {
      throw std::invalid_argument("too few frequencies between the limits");
    }
  }

  /** @brief The frequency spacing (cycles per unit time). */
  double Spacing() const { return df_; }
  /** @brief The imaginary part of the angular frequency, @f$-\epsilon@f$. */
  double Damping() const { return epsilon_; }
  /** @brief The time step. */
  double TimeStep() const { return dt_; }
  /** @brief The number of samples in the time series. */
  std::size_t Samples() const { return samples_; }
  /** @brief The index of the first frequency computed (zero frequency is 0). */
  std::size_t First() const { return first_; }
  /** @brief The index of the last frequency computed. */
  std::size_t Last() const { return last_; }
  /** @brief The number of frequencies computed. */
  std::size_t Count() const { return last_ - first_ + 1; }

  /** @brief The complex angular frequency of the j-th computed frequency. */
  Complex Omega(std::size_t j) const {
    return {2 * std::numbers::pi * static_cast<double>(first_ + j) * df_,
            -epsilon_};
  }

  /** @brief The frequency (cycles) of the j-th computed frequency. */
  double Frequency(std::size_t j) const {
    return static_cast<double>(first_ + j) * df_;
  }

 private:
  double dt_;
  double epsilon_;
  std::size_t samples_ = 0;
  double df_ = 0;
  std::size_t first_ = 0;
  std::size_t last_ = 0;
};

/** @brief The value of a taper at a frequency, in the taper's units. */
inline double TaperValue(const Taper& t, double f) {
  constexpr auto pi = std::numbers::pi;
  if (f < t.f1) return 0;
  if (f < t.f2) return 0.5 * (1 - std::cos(pi * (f - t.f1) / (t.f2 - t.f1)));
  if (f < t.f3) return 1;
  if (f < t.f4) return 0.5 * (1 - std::cos(pi * (t.f4 - f) / (t.f4 - t.f3)));
  return 0;
}

/**
 * @brief Transforms spectra to time series, with the damping undone.
 * @param grid The frequency grid.
 * @param taper A taper applied to the spectra, in the grid's units, if any.
 * @param spectra The spectra, @p count of them one after another, each
 * holding the grid's computed frequencies in order.
 * @param count The number of spectra.
 * @return The time series, one after another, each of grid.Samples()
 * samples.
 */
std::vector<double> InverseTransform(const FrequencyGrid& grid,
                                     const std::optional<Taper>& taper,
                                     const std::vector<Complex>& spectra,
                                     std::size_t count);

}  // namespace YSpec
