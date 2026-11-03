#pragma once

/**
 * @file Units.hpp
 * @brief The scales with which every quantity is made non-dimensional.
 *
 * @details A length, a density and the gravitational constant fix a time
 * scale, and with it the scale of every other quantity. The time scale is
 * chosen as @f$T = 1/\sqrt{\pi G \rho_0}@f$, so that @f$\pi G = 1@f$ in the
 * non-dimensional equations. This is the normalisation of the original YSpec
 * and of Minos, which keeps intermediate values directly comparable with
 * theirs.
 */

#include <cmath>
#include <numbers>

namespace YSpec {

/// The gravitational constant in SI units: the value used by Minos and by the
/// original YSpec, kept so that results can be compared with theirs.
inline constexpr auto GravitationalConstant = 6.6723e-11;

/// The value of @f$\pi G@f$ in the non-dimensional equations.
inline constexpr auto PiG = 1.0;

/**
 * @brief Scales for length, density and the quantities they imply.
 *
 * @details A dimensional quantity is divided by the matching scale to make it
 * non-dimensional, and multiplied by it to restore it.
 */
class Scales {
 public:
  /**
   * @brief Scales for a given length and density.
   * @param length The length scale in metres; by default the Earth's radius.
   * @param density The density scale in kg/m^3; by default the Earth's mean
   * density.
   */
  constexpr explicit Scales(double length = 6.371e6, double density = 5.515e3)
      : length_{length}, density_{density} {}

  /** @brief The length scale (m). */
  constexpr auto Length() const { return length_; }
  /** @brief The density scale (kg/m^3). */
  constexpr auto Density() const { return density_; }
  /** @brief The time scale @f$1/\sqrt{\pi G \rho_0}@f$ (s). */
  auto Time() const {
    return 1 / std::sqrt(std::numbers::pi * GravitationalConstant * density_);
  }
  /** @brief The frequency scale (rad/s), the inverse of the time scale. */
  auto Frequency() const { return 1 / Time(); }
  /** @brief The velocity scale (m/s). */
  auto Velocity() const { return length_ / Time(); }
  /** @brief The acceleration scale (m/s^2). */
  auto Acceleration() const { return Velocity() / Time(); }
  /** @brief The scale of elastic moduli (Pa). */
  auto Modulus() const { return density_ * Velocity() * Velocity(); }
  /** @brief The scale of a force (N). */
  auto Force() const { return Modulus() * length_ * length_; }
  /** @brief The scale of a moment tensor (N m). */
  auto Moment() const { return Force() * length_; }

 private:
  double length_;
  double density_;
};

}  // namespace YSpec
