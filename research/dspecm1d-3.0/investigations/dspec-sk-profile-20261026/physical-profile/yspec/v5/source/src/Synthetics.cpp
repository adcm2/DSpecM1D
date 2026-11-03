#include <YSpec/Synthetics.hpp>
#include <algorithm>
#include <atomic>
#include <cmath>
#include <exception>
#include <ranges>
#include <stdexcept>
#include <vector>

#include "AngularFunctions.hpp"
#include "Coefficients.hpp"
#include "Geometry.hpp"
#include "Kinds.hpp"
#include "OrderArray.hpp"
#include "Report.hpp"
#include "Spectra.hpp"
#include "TimeSeries.hpp"

namespace YSpec {

namespace {

// The spectra of seismograms of a kind, Displacement or Strain, turned to
// the receivers' frames and through their instruments.
template <class Kind>
Spectra ComputeSpectra(const Model& model, const Calculation& calculation,
                       std::span<const Source> sources,
                       std::span<const Receiver> receivers,
                       const Progress& progress) {
  constexpr auto components = Kind::Components;
  const auto& scales = model.GetScales();
  const auto& degrees = calculation.degrees;
  if (degrees.min < 0 || degrees.max < degrees.min) {
    throw std::invalid_argument("the degrees must satisfy 0 <= min <= max");
  }
  if (sources.empty() || receivers.empty()) {
    throw std::invalid_argument("there must be a source and a receiver");
  }

  const auto grid = MakeGrid(calculation, scales);

  // One solver for every source depth and every receiver depth.
  auto depths = std::vector<double>{};
  for (const auto& receiver : receivers) depths.push_back(receiver.depth);
  const auto solves = DepthSolves(model, calculation, sources, depths);

  const auto views =
      Views<typename Kind::Table>(sources, receivers, degrees.max);
  const auto receiverCount = receivers.size();
  const auto frequencyCount = grid.Count();
  auto spectra =
      Spectra(sources.size(), receiverCount, components, frequencyCount);
  auto done = std::atomic<std::size_t>{0};
  auto error = std::exception_ptr{};
  // Once a frequency fails the rest are skipped, and the error is thrown
  // when the loop ends.
  auto failed = std::atomic<bool>{false};

#pragma omp parallel for schedule(dynamic)
  for (long index = 0; index < static_cast<long>(frequencyCount); ++index) {
    if (failed) continue;
    try {
      const auto frequency = static_cast<std::size_t>(index);
      const auto omega = grid.Omega(frequency);
      // The sums at this frequency, gathered locally so that threads at
      // neighbouring frequencies do not write to the same cache lines.
      auto sums =
          std::vector<typename Kind::Sums>(sources.size() * receiverCount);
      for (auto l = degrees.min; l <= degrees.max; ++l) {
        ForEachSource(
            solves, model, calculation, sources, l, omega,
            [&](std::size_t source, std::size_t receiverDepth,
                const Coefficients& coefficients) {
              const auto point = solves.solver.Receiver(receiverDepth);
              const auto fluid = model.GetLayer(point.layer).IsFluid();
              for (auto receiver :
                   solves.receiverDepths[receiverDepth].members) {
                Kind::Add(sums[source * receiverCount + receiver], coefficients,
                          l, point.r, fluid, views, source, receiver);
              }
            });
      }
      for (auto source : std::views::iota(std::size_t{0}, sources.size())) {
        for (auto receiver : std::views::iota(std::size_t{0}, receiverCount)) {
          const auto& sum = sums[source * receiverCount + receiver];
          for (auto i : std::views::iota(std::size_t{0}, components)) {
            spectra(source, receiver, i, frequency) = sum[i];
          }
        }
      }
      const auto finished = ++done;
      Report(progress, finished, frequencyCount);
    } catch (...) {
      failed = true;
#pragma omp critical(YSpecError)
      if (!error) error = std::current_exception();
    }
  }
  if (error) std::rethrow_exception(error);

  // Turned to the receiver's frame, and through the instruments, if any:
  // both are linear and the same at every time, so they act on the spectra.
  // The instruments are evaluated at the same complex frequencies, which is
  // exact: the damping passes through a convolution unchanged.
  for (auto receiver : std::views::iota(std::size_t{0}, receiverCount)) {
    const auto& instrument = receivers[receiver].instrument;
    const auto response =
        instrument ? InstrumentResponse(*instrument, grid, scales.Time(),
                                        calculation.quantity)
                   : std::vector<Complex>(frequencyCount, 1.0);
    for (auto source : std::views::iota(std::size_t{0}, sources.size())) {
      const auto rotation = views.Of(source, receiver).geometry.rotation;
      for (auto frequency : std::views::iota(std::size_t{0}, frequencyCount)) {
        auto sum = typename Kind::Sums{};
        for (auto i : std::views::iota(std::size_t{0}, components)) {
          sum[i] = spectra(source, receiver, i, frequency);
        }
        const auto turned = Kind::ToReceiver(sum, rotation);
        for (auto i : std::views::iota(std::size_t{0}, components)) {
          spectra(source, receiver, i, frequency) =
              turned[i] * response[frequency];
        }
      }
    }
  }

  return spectra;
}

// The seismograms of a kind, Displacement or Strain.
template <class Kind>
typename Kind::Output Compute(const Model& model,
                              const Calculation& calculation,
                              std::span<const Source> sources,
                              std::span<const Receiver> receivers,
                              const Progress& progress) {
  constexpr auto components = Kind::Components;
  const auto spectra =
      ComputeSpectra<Kind>(model, calculation, sources, receivers, progress);
  const auto& scales = model.GetScales();
  const auto grid = MakeGrid(calculation, scales);
  const auto receiverCount = receivers.size();

  // To time series, in SI units or those of the instruments.
  const auto series =
      InverseTransform(grid, MakeTaper(calculation, scales), spectra.Values(),
                       sources.size() * receiverCount * components);
  const auto samples = SampleCount(grid, calculation.sampling);
  const auto unit = Kind::Unit(scales, calculation.quantity);
  auto output = typename Kind::Output(sources.size(), receiverCount, samples,
                                      calculation.sampling.timeStep);
  for (auto source : std::views::iota(std::size_t{0}, sources.size())) {
    for (auto receiver : std::views::iota(std::size_t{0}, receiverCount)) {
      for (auto component : std::views::iota(std::size_t{0}, components)) {
        const auto first =
            ((source * receiverCount + receiver) * components + component) *
            grid.Samples();
        const auto trace = std::span{series}.subspan(first, samples);
        auto out =
            output.Trace(source, receiver,
                         static_cast<typename Kind::ComponentType>(component));
        for (auto i : std::views::iota(std::size_t{0}, samples)) {
          out[i] = trace[i] * unit;
        }
      }
    }
  }
  return output;
}

}  // namespace

Spectra SeismogramSpectra(const Model& model, const Calculation& calculation,
                          std::span<const Source> sources,
                          std::span<const Receiver> receivers,
                          const Progress& progress) {
  return ComputeSpectra<Displacement>(model, calculation, sources, receivers,
                                      progress);
}

Seismograms ComputeSeismograms(const Model& model,
                               const Calculation& calculation,
                               std::span<const Source> sources,
                               std::span<const Receiver> receivers,
                               const Progress& progress) {
  return Compute<Displacement>(model, calculation, sources, receivers,
                               progress);
}

Strains ComputeStrains(const Model& model, const Calculation& calculation,
                       std::span<const Source> sources,
                       std::span<const Receiver> receivers,
                       const Progress& progress) {
  // The corrections are for instruments that measure motion.
  auto strain = calculation;
  strain.corrections = false;
  return Compute<Strain>(model, strain, sources, receivers, progress);
}

}  // namespace YSpec
