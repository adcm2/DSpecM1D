#pragma once

/**
 * @file Synthetics.hpp
 * @brief Synthetic seismograms for a set of sources and receivers.
 *
 * @details For each frequency and degree the radial equations are solved for
 * unit jumps at the source depth, and the responses are then combined with
 * each source's moment tensor and summed over degree and order at each
 * receiver. The solves depend only on the depths of the source and the
 * receiver, so sources at one depth, wherever they are and whatever their
 * moment tensors, share them, as do receivers at one depth.
 *
 * The seismograms record displacement, velocity or acceleration, in SI units,
 * in the vertical (up), north and east directions. For acceleration the
 * change in gravity seen by an instrument that moves with the ground, and the
 * tilt of its axis, can be included ("potential and tilt corrections").
 */

#include <cstddef>
#include <span>

#include "Calculation.hpp"
#include "Model.hpp"
#include "Receiver.hpp"
#include "Seismograms.hpp"
#include "Source.hpp"

namespace YSpec {

/**
 * @brief Computes synthetic seismograms.
 * @param model The earth model.
 * @param calculation What to compute.
 * @param sources The sources.
 * @param receivers The receivers.
 * @param progress Called as frequencies are completed, if given.
 * @throws std::invalid_argument if the calculation is inconsistent, or the
 * sources do not lie in one solid region of the model.
 */
Seismograms ComputeSeismograms(const Model& model,
                               const Calculation& calculation,
                               std::span<const Source> sources,
                               std::span<const Receiver> receivers,
                               const Progress& progress = {});

/**
 * @brief Computes synthetic seismograms of strain.
 *
 * @details As ComputeSeismograms, but recording the six components of the
 * strain tensor at each receiver, in spherical coordinates there (r up,
 * @f$\theta@f$ south, @f$\phi@f$ east), or their rate or second derivative.
 * The potential and tilt corrections do not apply. In a fluid, where the
 * radial derivative of the tangential displacement is not given, rt and rp
 * are zero.
 * @param model The earth model.
 * @param calculation What to compute.
 * @param sources The sources.
 * @param receivers The receivers.
 * @param progress Called as frequencies are completed, if given.
 * @throws std::invalid_argument as ComputeSeismograms does.
 */
Strains ComputeStrains(const Model& model, const Calculation& calculation,
                       std::span<const Source> sources,
                       std::span<const Receiver> receivers,
                       const Progress& progress = {});

}  // namespace YSpec
