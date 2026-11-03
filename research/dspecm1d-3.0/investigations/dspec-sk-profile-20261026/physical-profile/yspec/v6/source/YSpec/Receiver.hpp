#pragma once

/**
 * @file Receiver.hpp
 * @brief A receiver.
 */

#include <optional>

#include "Instrument.hpp"
#include "Location.hpp"

namespace YSpec {

/**
 * @brief A receiver: where on the surface, how deep (m), and the instrument
 * that records there, if any.
 *
 * @details Without an instrument the seismograms are the ground motion, in SI
 * units; with one they are its output, applied to every component.
 */
struct Receiver {
  Location location;  ///< The point on the surface above.
  double depth = 0;   ///< The depth (m).
  std::optional<Instrument> instrument = {};  ///< The instrument, if any.
};

}  // namespace YSpec
