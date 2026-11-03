#pragma once

/**
 * @file Seismograms.hpp
 * @brief Seismograms, of displacement or of strain, for every pair of a
 * source and a receiver.
 */

#include <cstddef>
#include <span>
#include <vector>

namespace YSpec {

/** @brief The components of a seismogram of displacement. */
enum Component { Vertical = 0, North = 1, East = 2 };

/**
 * @brief The components of a seismogram of strain, in spherical coordinates
 * at the receiver: r up, @f$\theta@f$ south and @f$\phi@f$ east, as for a
 * moment tensor.
 */
enum class StrainComponent { rr = 0, tt = 1, pp = 2, rt = 3, rp = 4, tp = 5 };

/**
 * @brief Traces for a set of sources and receivers, each with a number of
 * components, in SI units of the quantity recorded or, at a receiver with an
 * instrument, in the units of its output.
 *
 * @details Sample @f$i@f$ is at time @f$i\,\Delta t@f$ after the origin time.
 * Each trace is contiguous.
 *
 * @tparam ComponentType The enumeration of the components.
 * @tparam ComponentCount Their number.
 */
template <class ComponentType, std::size_t ComponentCount>
class Traces {
 public:
  /** @brief Zero traces of a given size. */
  Traces(std::size_t sources, std::size_t receivers, std::size_t samples,
         double timeStep)
      : sources_{sources},
        receivers_{receivers},
        samples_{samples},
        timeStep_{timeStep},
        data_(sources * receivers * ComponentCount * samples) {}

  /** @brief The number of sources. */
  std::size_t Sources() const { return sources_; }
  /** @brief The number of receivers. */
  std::size_t Receivers() const { return receivers_; }
  /** @brief The number of samples in each trace. */
  std::size_t Samples() const { return samples_; }
  /** @brief The interval between samples (s). */
  double TimeStep() const { return timeStep_; }

  /** @brief One component at one receiver from one source. */
  std::span<double> Trace(std::size_t source, std::size_t receiver,
                          ComponentType component) {
    return {data_.data() + Offset(source, receiver, component), samples_};
  }
  /** @brief One component at one receiver from one source. */
  std::span<const double> Trace(std::size_t source, std::size_t receiver,
                                ComponentType component) const {
    return {data_.data() + Offset(source, receiver, component), samples_};
  }

  /** @brief A sample. */
  double& operator()(std::size_t source, std::size_t receiver,
                     ComponentType component, std::size_t i) {
    return data_[Offset(source, receiver, component) + i];
  }
  /** @brief A sample. */
  double operator()(std::size_t source, std::size_t receiver,
                    ComponentType component, std::size_t i) const {
    return data_[Offset(source, receiver, component) + i];
  }

 private:
  std::size_t sources_;
  std::size_t receivers_;
  std::size_t samples_;
  double timeStep_;
  std::vector<double> data_;

  std::size_t Offset(std::size_t source, std::size_t receiver,
                     ComponentType component) const {
    return ((source * receivers_ + receiver) * ComponentCount +
            static_cast<std::size_t>(component)) *
           samples_;
  }
};

/** @brief Seismograms of displacement: vertical, north and east. */
using Seismograms = Traces<Component, 3>;

/** @brief Seismograms of strain: its six independent components. */
using Strains = Traces<StrainComponent, 6>;

}  // namespace YSpec
