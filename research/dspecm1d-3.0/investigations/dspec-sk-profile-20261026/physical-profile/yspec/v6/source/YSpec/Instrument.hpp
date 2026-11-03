#pragma once

/**
 * @file Instrument.hpp
 * @brief The response of an instrument, as poles and zeros.
 */

#include <filesystem>
#include <vector>

#include "Calculation.hpp"
#include "Types.hpp"

namespace YSpec {

/**
 * @brief The response of an instrument to ground motion, as poles and zeros.
 *
 * @details In terms of the Laplace variable @f$s = i\omega@f$, with
 * @f$\omega@f$ in radians per second, the response is
 * @f[ H(s) = c\,\frac{\prod_k (s - z_k)}{\prod_k (s - p_k)}, @f]
 * the ratio of the instrument's output to the ground motion it responds to,
 * which is in SI units: displacement in m, velocity in m/s or acceleration
 * in m/s^2. This is the convention of SAC's pole-zero files.
 */
struct Instrument {
  std::vector<Complex> zeros;               ///< The zeros (rad/s).
  std::vector<Complex> poles;               ///< The poles (rad/s).
  double constant = 1;                      ///< The constant @f$c@f$.
  Quantity input = Quantity::Displacement;  ///< The motion responded to.

  /** @brief The response at an angular frequency (rad/s), real or complex. */
  Complex Response(Complex omega) const {
    const auto s = Complex(0, 1) * omega;
    auto response = Complex(constant);
    for (const auto& zero : zeros) response *= s - zero;
    for (const auto& pole : poles) response /= s - pole;
    return response;
  }
};

/**
 * @brief Reads a SAC pole-zero file.
 *
 * @details The file holds the keywords ZEROS, POLES and CONSTANT. A count
 * follows ZEROS and POLES, and the lines after them hold the real and
 * imaginary parts of up to that many zeros or poles: those not listed are
 * zero. Lines starting with '*' are comments; one reading "INPUT UNIT : M"
 * (or M/S, or M/S**2), as in the files of the IRIS DMC, gives the motion the
 * instrument responds to, which is otherwise displacement.
 *
 * @throws std::runtime_error naming the file and the line of a fault,
 * including a file that holds more than one response.
 */
Instrument ReadSacPoleZero(const std::filesystem::path& path);

}  // namespace YSpec
