#pragma once

// The spectra of seismograms, before they become time series, and what turns
// them into time series. The kernels need both: the adjoint source of a
// measurement is formed from the spectrum of the synthetic, and reaches the
// spectra through the same transform.

#include <YSpec/Calculation.hpp>
#include <YSpec/Instrument.hpp>
#include <YSpec/Model.hpp>
#include <YSpec/Receiver.hpp>
#include <YSpec/Source.hpp>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <ranges>
#include <span>
#include <vector>

#include "Coefficients.hpp"
#include "TimeSeries.hpp"

namespace YSpec {

// The spectra for every source, receiver and component, at the grid's
// frequencies.
class Spectra {
 public:
  Spectra(std::size_t sources, std::size_t receivers, std::size_t components,
          std::size_t frequencies)
      : receivers_{receivers},
        components_{components},
        frequencies_{frequencies},
        values_(sources * receivers * components * frequencies) {}

  Complex& operator()(std::size_t source, std::size_t receiver,
                      std::size_t component, std::size_t frequency) {
    return values_[Index(source, receiver, component, frequency)];
  }
  Complex operator()(std::size_t source, std::size_t receiver,
                     std::size_t component, std::size_t frequency) const {
    return values_[Index(source, receiver, component, frequency)];
  }
  const std::vector<Complex>& Values() const { return values_; }

 private:
  std::size_t receivers_;
  std::size_t components_;
  std::size_t frequencies_;
  std::vector<Complex> values_;

  std::size_t Index(std::size_t source, std::size_t receiver,
                    std::size_t component, std::size_t frequency) const {
    return ((source * receivers_ + receiver) * components_ + component) *
               frequencies_ +
           frequency;
  }
};

// The spectra of seismograms of displacement (or the quantity the
// calculation records) at the grid's computed frequencies, in (up, north,
// east) and through the receivers' instruments, without the taper,
// non-dimensional.
Spectra SeismogramSpectra(const Model& model, const Calculation& calculation,
                          std::span<const Source> sources,
                          std::span<const Receiver> receivers,
                          const Progress& progress = {});

// The response of an instrument to the quantity computed, at the grid's
// frequencies (non-dimensional, with the unit of time given): its response to
// the motion it responds to, converted by a power of i omega.
inline std::vector<Complex> InstrumentResponse(const Instrument& instrument,
                                               const FrequencyGrid& grid,
                                               double time, Quantity computed) {
  const auto power = Derivatives(instrument.input) - Derivatives(computed);
  auto response = std::vector<Complex>{};
  for (auto j : std::views::iota(std::size_t{0}, grid.Count())) {
    const auto omega = grid.Omega(j) / time;
    response.push_back(instrument.Response(omega) *
                       std::pow(Complex(0, 1) * omega, power));
  }
  return response;
}

// The number of samples of the seismograms: those of the length asked for,
// or of the transform if that is shorter.
inline std::size_t SampleCount(const FrequencyGrid& grid,
                               const Sampling& sampling) {
  return std::min<std::size_t>(
      grid.Samples(),
      static_cast<std::size_t>(sampling.duration / sampling.timeStep + 1e-9) +
          1);
}

}  // namespace YSpec
