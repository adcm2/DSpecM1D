#pragma once

/**
 * @file Calculation.hpp
 * @brief What a calculation computes, apart from its sources and receivers.
 *
 * @details Everything is in SI units: frequencies in Hz, times in seconds.
 */

#include <cstddef>
#include <functional>
#include <optional>

namespace YSpec {

/** @brief How gravity enters the equations. */
enum class Gravity {
  None,     ///< No gravity at all.
  Cowling,  ///< The equilibrium gravity, but no perturbation of the potential.
  Full,     ///< Full self-gravitation.
};

/** @brief The quantity the seismograms record. */
enum class Quantity { Displacement, Velocity, Acceleration };

/** @brief The spherical harmonic degrees summed. */
struct Degrees {
  int min = 0;  ///< The lowest degree.
  int max = 0;  ///< The highest degree.
};

/** @brief The band of frequencies computed (Hz). */
struct FrequencyBand {
  double min = 0;  ///< The lowest frequency.
  double max = 0;  ///< The highest frequency.
};

/** @brief The time series computed. */
struct Sampling {
  double duration = 0;  ///< The length of the seismograms (s).
  double timeStep = 0;  ///< The interval between samples (s).
};

/**
 * @brief A taper applied to the spectra: a raised cosine from zero at @c f1
 * to one at @c f2, one to @c f3, and a raised cosine back to zero at @c f4
 * (Hz).
 */
struct Taper {
  double f1 = 0;  ///< Where the taper starts to rise.
  double f2 = 0;  ///< Where it reaches one.
  double f3 = 0;  ///< Where it starts to fall.
  double f4 = 0;  ///< Where it reaches zero.
};

/** @brief Settings of the radial integrations. */
struct SolverOptions {
  /// The local error allowed per step, relative to the size of the solution.
  /// The default makes seismograms accurate to a few parts in 10^4 of their
  /// peak, about ten times better than the original YSpec; 1e-4 matches it at
  /// less cost.
  double tolerance = 1e-5;
  /// The decay, in e-folds, below a turning point at which to start.
  double startDecay = 11.0;
  /// The first step in a layer is this many times @f$\alpha/\omega@f$, for the
  /// largest P-wave speed in the layer.
  double stepFactor = 0.25;
};

/**
 * @brief What a calculation computes.
 *
 * @details For acceleration, @c corrections adds the change in gravity seen
 * by an instrument that moves with the ground and the tilt of its axis
 * ("potential and tilt corrections").
 */
struct Calculation {
  Gravity gravity = Gravity::Full;             ///< How gravity enters.
  bool attenuation = true;                     ///< Whether to attenuate.
  Quantity quantity = Quantity::Displacement;  ///< What is recorded.
  bool corrections = false;    ///< Potential and tilt corrections.
  Degrees degrees = {};        ///< The degrees summed.
  FrequencyBand band = {};     ///< The frequencies computed.
  Sampling sampling = {};      ///< The time series.
  std::optional<Taper> taper;  ///< A taper applied to the spectra, if any.
  SolverOptions solver = {};   ///< Settings of the radial integrations.
};

/**
 * @brief Called as frequencies are completed, with the number done and the
 * number in all; called from one thread at a time.
 */
using Progress = std::function<void(std::size_t, std::size_t)>;

}  // namespace YSpec
