#include "TimeSeries.hpp"

#include <algorithm>
#include <cmath>
#include <ranges>
#include <stdexcept>

#ifdef YSPEC_WITH_POCKETFFT
#include <pocketfft_hdronly.h>
#else
#include <FFTWpp/Ranges>
#endif

namespace YSpec {

namespace {

// The inverse real transforms of count spectra, each of bins frequencies from
// zero, to count series of samples each, unnormalised, by whichever library
// the build has.
#ifdef YSPEC_WITH_POCKETFFT
using ComplexStorage = std::vector<Complex>;
using RealStorage = std::vector<double>;

void Inverse(ComplexStorage& spectrum, RealStorage& result, int samples,
             int bins, int count) {
  const auto rows = static_cast<std::size_t>(count);
  const auto columns = static_cast<std::size_t>(samples);
  const auto complexStride = static_cast<std::ptrdiff_t>(sizeof(Complex));
  const auto realStride = static_cast<std::ptrdiff_t>(sizeof(double));
  pocketfft::c2r<double>({rows, columns}, {bins * complexStride, complexStride},
                         {samples * realStride, realStride}, 1,
                         pocketfft::BACKWARD, spectrum.data(), result.data(),
                         1.0);
}
#else
using ComplexStorage = FFTWpp::vector<Complex>;
using RealStorage = FFTWpp::vector<double>;

void Inverse(ComplexStorage& spectrum, RealStorage& result, int samples,
             int bins, int count) {
  using namespace FFTWpp;
  const auto complexLayout =
      Ranges::Layout(1, std::vector{bins}, count, std::vector{bins}, 1, bins);
  const auto realLayout = Ranges::Layout(1, std::vector{samples}, count,
                                         std::vector{samples}, 1, samples);
  // Planned by estimate, which leaves the spectrum as it is.
  auto plan = Ranges::Plan(Ranges::View(spectrum, complexLayout),
                           Ranges::View(result, realLayout), Estimate);
  plan.Execute();
}
#endif

}  // namespace

std::vector<double> InverseTransform(const FrequencyGrid& grid,
                                     const std::optional<Taper>& taper,
                                     const std::vector<Complex>& spectra,
                                     std::size_t count) {
  const auto samples = static_cast<int>(grid.Samples());
  // The frequencies of a real series of that length, zero to Nyquist.
  const auto bins = samples / 2 + 1;
  const auto frequencies = grid.Count();
  if (spectra.size() != count * frequencies) {
    throw std::invalid_argument("InverseTransform: spectra of the wrong size");
  }
  auto spectrum = ComplexStorage(count * static_cast<std::size_t>(bins));
  auto result = RealStorage(count * static_cast<std::size_t>(samples));

  // Frequencies above the Nyquist frequency are dropped.
  const auto last = std::min<std::size_t>(grid.Last(), bins - 1);
  for (auto series : std::views::iota(std::size_t{0}, count)) {
    for (auto i = grid.First(); i <= last; ++i) {
      const auto f = static_cast<double>(i) * grid.Spacing();
      const auto weight = taper ? TaperValue(*taper, f) : 1.0;
      spectrum[series * bins + i] =
          weight * spectra[series * frequencies + (i - grid.First())];
    }
  }
  Inverse(spectrum, result, samples, bins, static_cast<int>(count));

  // The sum approximates the inverse Fourier integral with the frequency
  // spacing, and the damping is undone.
  const auto timeStep = grid.TimeStep();
  const auto spacing = grid.Spacing();
  for (auto series : std::views::iota(std::size_t{0}, count)) {
    for (auto i : std::views::iota(0, samples)) {
      result[series * samples + i] *=
          spacing * std::exp(grid.Damping() * i * timeStep);
    }
  }
  return {result.begin(), result.end()};
}

}  // namespace YSpec
