#include <YSpec/Snapshots.hpp>
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <exception>
#include <numbers>
#include <ranges>
#include <stdexcept>
#include <vector>

#include "AngularFunctions.hpp"
#include "Coefficients.hpp"
#include "OrderArray.hpp"
#include "Report.hpp"

namespace YSpec {

namespace {

using namespace std::complex_literals;

// Complex values for every source, receiver depth, distance, order and
// component, laid out in that order, with any number of values for each.
class Table {
 public:
  Table(std::size_t sources, std::size_t depths, std::size_t distances,
        std::size_t each)
      : depths_{depths},
        distances_{distances},
        each_{each},
        values_(sources * depths * distances * 5 * 3 * each) {}

  // The first of the values for a source, depth, distance, order and
  // component.
  std::size_t Index(std::size_t source, std::size_t depth, std::size_t distance,
                    int m, int component) const {
    const auto order = static_cast<std::size_t>(m + 2);
    return ((((source * depths_ + depth) * distances_ + distance) * 5 + order) *
                3 +
            static_cast<std::size_t>(component)) *
           each_;
  }
  Complex& operator[](std::size_t i) { return values_[i]; }
  Complex operator[](std::size_t i) const { return values_[i]; }
  std::vector<Complex>& Values() { return values_; }
  const std::vector<Complex>& Values() const { return values_; }

 private:
  std::size_t depths_;
  std::size_t distances_;
  std::size_t each_;
  std::vector<Complex> values_;
};

void Check(const Calculation& calculation, std::span<const Source> sources,
           const SnapshotGrid& grid, std::span<const double> times) {
  const auto& degrees = calculation.degrees;
  if (degrees.min < 0 || degrees.max < degrees.min) {
    throw std::invalid_argument("the degrees must satisfy 0 <= min <= max");
  }
  if (sources.empty()) throw std::invalid_argument("there must be a source");
  if (grid.depths.empty() || grid.distances.empty() || grid.azimuths.empty() ||
      times.empty()) {
    throw std::invalid_argument(
        "a snapshot needs a depth, a distance, an azimuth and a time");
  }
  for (auto distance : grid.distances) {
    if (!(distance >= 0 && distance <= 180)) {
      throw std::invalid_argument(
          "the distances of a snapshot must lie between 0 and 180 degrees");
    }
  }
  for (auto time : times) {
    if (!(time >= 0 && time <= calculation.sampling.duration)) {
      throw std::invalid_argument(
          "the times of a snapshot must lie between zero and the duration");
    }
  }
}

}  // namespace

Snapshots ComputeSnapshots(const Model& model, const Calculation& calculation,
                           std::span<const Source> sources,
                           const SnapshotGrid& grid,
                           std::span<const double> times,
                           const Progress& progress) {
  Check(calculation, sources, grid, times);
  const auto& scales = model.GetScales();
  const auto frequencies = MakeGrid(calculation, scales);
  const auto taper = MakeTaper(calculation, scales);
  const auto solves = DepthSolves(model, calculation, sources, grid.depths);

  // The angular functions at the distances.
  auto radians = std::vector<double>{};
  for (auto distance : grid.distances) {
    radians.push_back(distance * std::numbers::pi / 180);
  }
  const auto maxDegree = calculation.degrees.max;
  const auto angular = AngularFunctions(maxDegree, radians);

  // The times, non-dimensional.
  auto instants = std::vector<double>{};
  for (auto time : times) instants.push_back(time / scales.Time());

  const auto sourceCount = sources.size();
  const auto depthCount = solves.receiverDepths.size();
  const auto distanceCount = grid.distances.size();
  const auto timeCount = instants.size();
  const auto frequencyCount = frequencies.Count();

  // For each source, depth, distance, order and component, the sum over
  // frequency at each time of the spectrum times exp(i omega t).
  auto sums = Table(sourceCount, depthCount, distanceCount, timeCount);
  auto done = std::atomic<std::size_t>{0};
  auto error = std::exception_ptr{};
  // Once a frequency fails the rest are skipped, and the error is thrown
  // when the loop ends.
  auto failed = std::atomic<bool>{false};

  // The coefficients at one frequency, for each pair of a source and a depth,
  // each degree and each order: u, v and w.
  const auto minDegree = calculation.degrees.min;
  const auto degreeCount = static_cast<std::size_t>(maxDegree - minDegree + 1);
  const auto pairCount = sourceCount * depthCount;
  auto coefficientIndex = [&](std::size_t pair, int l) {
    return (pair * degreeCount + static_cast<std::size_t>(l - minDegree)) * 15;
  };

#pragma omp parallel
  {
    auto local = Table(sourceCount, depthCount, distanceCount, timeCount);
    auto spectrum = Table(sourceCount, depthCount, distanceCount, 1);
    // Orders that a degree does not have stay zero.
    auto byDegree = std::vector<Complex>(pairCount * degreeCount * 15);
#pragma omp for schedule(dynamic)
    for (long index = 0; index < static_cast<long>(frequencyCount); ++index) {
      if (failed) continue;
      try {
        const auto frequency = static_cast<std::size_t>(index);
        const auto omega = frequencies.Omega(frequency);
        for (auto l = minDegree; l <= maxDegree; ++l) {
          ForEachSource(
              solves, model, calculation, sources, l, omega,
              [&](std::size_t source, std::size_t depth,
                  const Coefficients& coefficients) {
                auto* at =
                    &byDegree[coefficientIndex(source * depthCount + depth, l)];
                for (auto m : OrderArray<Complex>::Orders(l)) {
                  const auto order = 3 * (m + 2);
                  at[order] = coefficients.u[m];
                  at[order + 1] = coefficients.v[m];
                  at[order + 2] = coefficients.w[m];
                }
              });
        }

        // The sums over degree, a block of distances at a time, so that the
        // block's angular functions stay in the cache for every pair, and
        // the sums for the block in registers or close by.
        constexpr auto block = std::size_t{8};
        for (auto first = std::size_t{0}; first < distanceCount;
             first += block) {
          const auto count = std::min(block, distanceCount - first);
          for (auto pair : std::views::iota(std::size_t{0}, pairCount)) {
            auto sum = std::array<Complex, block * 15>{};
            for (auto l = minDegree; l <= maxDegree; ++l) {
              const auto* c = &byDegree[coefficientIndex(pair, l)];
              for (auto b : std::views::iota(std::size_t{0}, count)) {
                // The five orders are contiguous.
                const auto* f = &angular(static_cast<Int>(first + b), l, -2);
                auto* at = &sum[b * 15];
                for (auto order = 0; order < 5; ++order) {
                  const auto u = c[3 * order];
                  const auto v = c[3 * order + 1];
                  const auto w = c[3 * order + 2];
                  const auto& a = f[order];
                  at[3 * order + Up] += u * a.value;
                  at[3 * order + Away] += v * a.dTheta + TimesI(w) * a.mOverSin;
                  at[3 * order + Transverse] +=
                      TimesI(v) * a.mOverSin - w * a.dTheta;
                }
              }
            }
            auto* values = &spectrum[(pair * distanceCount + first) * 15];
            std::ranges::copy(std::span{sum}.first(count * 15), values);
          }
        }

        // The spectrum at this frequency, at each time.
        const auto f = frequencies.Frequency(frequency);
        const auto weight = taper ? TaperValue(*taper, f) : 1.0;
        auto phase = std::vector<Complex>{};
        for (auto t : instants) {
          phase.push_back(weight * std::exp(2i * std::numbers::pi * f * t));
        }
        const auto& values = spectrum.Values();
        auto& sum = local.Values();
        for (auto i : std::views::iota(std::size_t{0}, values.size())) {
          if (values[i] == Complex{}) continue;
          for (auto k : std::views::iota(std::size_t{0}, timeCount)) {
            sum[i * timeCount + k] += values[i] * phase[k];
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
#pragma omp critical(YSpecSnapshotSum)
    {
      auto& total = sums.Values();
      const auto& part = local.Values();
      for (auto i : std::views::iota(std::size_t{0}, total.size())) {
        total[i] += part[i];
      }
    }
  }
  if (error) std::rethrow_exception(error);

  // The sums over the orders at each azimuth. The spectra are of the
  // positive frequencies of a real series, hence twice the real part; the
  // spacing makes the sum an integral, and the damping is undone.
  auto snapshots = Snapshots(sourceCount, {times.begin(), times.end()}, grid);
  const auto unit = Unit(scales, calculation.quantity);
  auto depthGroup = std::vector<std::size_t>(grid.depths.size());
  for (auto group : std::views::iota(std::size_t{0}, depthCount)) {
    for (auto depth : solves.receiverDepths[group].members) {
      depthGroup[depth] = group;
    }
  }
  for (auto source : std::views::iota(std::size_t{0}, sourceCount)) {
    for (auto k : std::views::iota(std::size_t{0}, timeCount)) {
      const auto scale = 2 * frequencies.Spacing() *
                         std::exp(frequencies.Damping() * instants[k]) * unit;
      for (auto depth : std::views::iota(std::size_t{0}, grid.depths.size())) {
        const auto group = depthGroup[depth];
        for (auto a : std::views::iota(std::size_t{0}, grid.azimuths.size())) {
          // The longitude in the frame with the source at the pole.
          const auto phi = std::numbers::pi * (1 - grid.azimuths[a] / 180);
          for (auto i : std::views::iota(std::size_t{0}, distanceCount)) {
            for (auto component : {Up, Away, Transverse}) {
              auto value = Complex{};
              for (auto m : OrderArray<Complex>::Orders(2)) {
                value += sums[sums.Index(source, group, i, m, component) + k] *
                         std::exp(1i * (m * phi));
              }
              snapshots(source, k, depth, i, a, component) =
                  scale * value.real();
            }
          }
        }
      }
    }
  }
  return snapshots;
}

}  // namespace YSpec
