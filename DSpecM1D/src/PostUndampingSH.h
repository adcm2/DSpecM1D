#ifndef DSPECM1D_POST_UNDAMPING_SH_H
#define DSPECM1D_POST_UNDAMPING_SH_H

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>
#include <Eigen/Core>
#include "SRInfo.h"
#include "frequency_info/PostprocessFunctions.h"

namespace SPARSESPEC {

/// Experimental SH result before any additional bandpass.
struct PostUndampingSHResult {
  Eigen::MatrixXd timeSeries;
  double sampleIntervalSeconds = 0.0;
  std::size_t solvedPairs = 0;
  std::size_t fullGridPairs = 0;
  double radialSeconds = 0.0;
  double synthesisSeconds = 0.0;
};

/// Component and slowness taper for the opt-in mixed-mode phase window.
enum class PostUndampingWindowComponent { Transverse, Vertical };

struct PostUndampingWindowOptions {
  PostUndampingWindowComponent component =
      PostUndampingWindowComponent::Transverse;
  // Ordered [outer-low, inner-low, inner-high, outer-high] in s/deg.
  std::array<double, 4> slownessCorners{};
  // Omitted retains the full requested solve grid. A supplied guard enables
  // pruning in the expanded outer interval.
  std::optional<double> preSolveGuardSPerDegree;
};

/// Experimental single-component trace synthesized from selected degree bins.
struct PostUndampingWindowResult {
  Eigen::MatrixXd timeSeries;
  double sampleIntervalSeconds = 0.0;
  std::size_t solvedPairs = 0;
  std::size_t fullGridPairs = 0;
  double radialSeconds = 0.0;
  double synthesisSeconds = 0.0;
};

namespace detail {

struct SHDegreeSpectrum {
  int degree = 0;
  int firstBin = 0;
  // Source-centred horizontal components only, over the retained bin suffix.
  Eigen::MatrixXcd horizontal;
};

struct SHDegreeSpectra {
  std::vector<SHDegreeSpectrum> degrees;
  std::size_t solvedPairs = 0;
  std::size_t fullGridPairs = 0;
};

struct WindowDegreeSpectrum {
  int degree = 0;
  int firstBin = 0;
  int lastBin = 0;
  // Receiver-major selected component over the compact retained span.
  // Original FFT coordinates are recovered from firstBin.
  Eigen::MatrixXcd components;
};

struct WindowDegreeSpectra {
  std::vector<WindowDegreeSpectrum> degrees;
  PostUndampingWindowOptions options;
  std::size_t solvedPairs = 0;
  std::size_t fullGridPairs = 0;
};

inline constexpr double shMaxSlowness = 26.0;

inline std::pair<double, double>
validatePostUndampingWindowOptions(const PostUndampingWindowOptions &options) {
  for (std::size_t i = 0; i < options.slownessCorners.size(); ++i) {
    if (!std::isfinite(options.slownessCorners[i]) ||
        (i > 0 && options.slownessCorners[i] <= options.slownessCorners[i - 1]))
      throw std::invalid_argument(
          "Slowness corners must be finite and strictly ordered.");
  }
  if (!options.preSolveGuardSPerDegree)
    return {-std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity()};
  const double guard = *options.preSolveGuardSPerDegree;
  if (!std::isfinite(guard) || guard < 0.0)
    throw std::invalid_argument(
        "Pre-solve guard must be finite and nonnegative.");
  const double low = std::max(0.0, options.slownessCorners[0] - guard);
  const double high = options.slownessCorners[3] + guard;
  if (!std::isfinite(high))
    throw std::invalid_argument("Expanded slowness window overflows.");
  return {low, high};
}

inline double
degreeSlownessSPerDegree(int degree, double frequencyHz) {
  return (static_cast<double>(degree) + 0.5) / (360.0 * frequencyHz);
}

inline double
shFirstSolveFrequencyHz(int degree) {
  return (static_cast<double>(degree) + 0.5) / (360.0 * shMaxSlowness);
}

inline Eigen::MatrixXd
synthesizePostUndampingSH(const SpectraSolver::FreqFull &freq, SRInfo &srInfo,
                          int numReceivers, double normFactor,
                          const SHDegreeSpectra &spectra) {
  if (freq.nt() < 2 || freq.nt() > std::numeric_limits<int>::max() / 2 ||
      !std::isfinite(freq.tout()) || freq.tout() <= 0.0) {
    throw std::invalid_argument("Invalid FFT grid for post-undamping SH.");
  }

  // Halving the sampling interval preserves the original solve-bin spacing.
  const int fftSize = 2 * freq.nt();
  const int numBins = fftSize / 2 + 1;
  const double dt = freq.dt() * 0.5;
  const double dfHz = freq.df() / freq.timeNorm();
  const int samples = static_cast<int>(std::min(
      static_cast<double>(fftSize), std::floor(freq.tout() / dt + 0.5) + 1.0));
  const int horizontalRows = 2 * numReceivers;
  Eigen::MatrixXcd accumulated = Eigen::MatrixXcd::Zero(horizontalRows, numBins);
  processfunctions::RealVector time(fftSize);
  processfunctions::ComplexVector bins(numBins);
  auto inverse = FFTWpp::Ranges::Plan(FFTWpp::Ranges::View(bins),
                                      FFTWpp::Ranges::View(time),
                                      FFTWpp::Estimate);
  auto forward = FFTWpp::Ranges::Plan(FFTWpp::Ranges::View(time),
                                      FFTWpp::Ranges::View(bins),
                                      FFTWpp::Estimate);

  constexpr double rampFraction = 0.025;
  const double pi = std::acos(-1.0);
  std::vector<double> undamp(samples);
  for (int i = 0; i < samples; ++i) {
    const double t = i * dt;
    const double ramp =
        std::clamp((t / freq.tout() - (1.0 - rampFraction)) / rampFraction,
                   0.0, 1.0);
    undamp[i] = std::exp(freq.ep() * t) * (0.5 + 0.5 * std::cos(pi * ramp));
  }

  std::vector<double> weights(numBins, 0.0);
  for (const auto &degree : spectra.degrees) {
    if (degree.horizontal.cols() == 0)
      continue;
    for (int i = 1; i < numBins; ++i) {
      const double p = (static_cast<double>(degree.degree) + 0.5) /
                       (360.0 * i * dfHz);
      weights[i] = filters::hannref(p, 0.0, 0.01, 23.5, shMaxSlowness);
    }
    for (int row = 0; row < horizontalRows; ++row) {
      std::fill(bins.begin(), bins.end(), std::complex<double>{});
      for (Eigen::Index i = 0; i < degree.horizontal.cols(); ++i)
        bins[degree.firstBin + i] = degree.horizontal(row, i);
      inverse.Execute();
      for (int i = 0; i < samples; ++i)
        time[i] *= undamp[i];
      std::fill(time.begin() + samples, time.end(), 0.0);
      forward.Execute();
      for (int i = 1; i < numBins; ++i)
        accumulated(row, i) += bins[i] * weights[i];
    }
  }

  // FFTW's transforms are unnormalized. The synthesis uses normalized df
  // like freq2time(), and applies output units once, after summing degrees.
  const double scale = (freq.df() * normFactor) / fftSize;
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(3 * numReceivers, samples);
  for (int receiver = 0; receiver < numReceivers; ++receiver) {
    const double baz = srInfo.backazimuths(receiver);
    const double sinBaz = std::sin(baz);
    const double cosBaz = std::cos(baz);
    for (int i = 0; i < numBins; ++i)
      bins[i] = -cosBaz * accumulated(2 * receiver, i) +
                sinBaz * accumulated(2 * receiver + 1, i);
    inverse.Execute();
    for (int i = 0; i < samples; ++i)
      result(3 * receiver + 1, i) = scale * time[i];
    for (int i = 0; i < numBins; ++i)
      bins[i] = sinBaz * accumulated(2 * receiver, i) +
                cosBaz * accumulated(2 * receiver + 1, i);
    inverse.Execute();
    for (int i = 0; i < samples; ++i)
      result(3 * receiver + 2, i) = scale * time[i];
  }
  return result;
}

inline Eigen::MatrixXd
synthesizePostUndampingWindow(const SpectraSolver::FreqFull &freq,
                              int numReceivers, double normFactor,
                              const WindowDegreeSpectra &spectra) {
  if (freq.nt() < 2 || freq.nt() > std::numeric_limits<int>::max() / 2 ||
      !std::isfinite(freq.tout()) || freq.tout() <= 0.0) {
    throw std::invalid_argument("Invalid FFT grid for post-undamping window.");
  }
  const int fftSize = 2 * freq.nt();
  const int numBins = fftSize / 2 + 1;
  const double dt = freq.dt() * 0.5;
  const double dfHz = freq.df() / freq.timeNorm();
  const int samples = static_cast<int>(std::min(
      static_cast<double>(fftSize), std::floor(freq.tout() / dt + 0.5) + 1.0));
  Eigen::MatrixXcd accumulated = Eigen::MatrixXcd::Zero(numReceivers, numBins);
  processfunctions::RealVector time(fftSize);
  processfunctions::ComplexVector bins(numBins);
  auto inverse = FFTWpp::Ranges::Plan(FFTWpp::Ranges::View(bins),
                                      FFTWpp::Ranges::View(time),
                                      FFTWpp::Estimate);
  auto forward = FFTWpp::Ranges::Plan(FFTWpp::Ranges::View(time),
                                      FFTWpp::Ranges::View(bins),
                                      FFTWpp::Estimate);
  constexpr double rampFraction = 0.025;
  const double pi = std::acos(-1.0);
  std::vector<double> undamp(samples);
  for (int i = 0; i < samples; ++i) {
    const double t = i * dt;
    const double ramp = std::clamp(
        (t / freq.tout() - (1.0 - rampFraction)) / rampFraction, 0.0, 1.0);
    undamp[i] = std::exp(freq.ep() * t) * (0.5 + 0.5 * std::cos(pi * ramp));
  }

  // The capture is receiver-major by degree in `components`; do each receiver
  // independently so FFT plans and damping operations are reused.
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(numReceivers, samples);
  const double scale = (freq.df() * normFactor) / fftSize;
  for (int receiver = 0; receiver < numReceivers; ++receiver) {
    for (const auto &degree : spectra.degrees) {
      if (degree.components.rows() != numReceivers ||
          degree.components.cols() != degree.lastBin - degree.firstBin)
        throw std::invalid_argument("Window degree spectrum has the wrong grid.");
      std::fill(bins.begin(), bins.end(), std::complex<double>{});
      for (int i = degree.firstBin; i < degree.lastBin; ++i) {
        if (i <= 0 || i >= numBins)
          continue;
        bins[i] = degree.components(receiver, i - degree.firstBin);
      }
      inverse.Execute();
      for (int i = 0; i < samples; ++i)
        time[i] *= undamp[i];
      std::fill(time.begin() + samples, time.end(), 0.0);
      forward.Execute();
      for (int i = 1; i < numBins; ++i) {
        const double p = degreeSlownessSPerDegree(degree.degree, i * dfHz);
        const auto &c = spectra.options.slownessCorners;
        accumulated(receiver, i) +=
            bins[i] * filters::hannref(p, c[0], c[1], c[2], c[3]);
      }
    }
    for (int i = 0; i < numBins; ++i)
      bins[i] = accumulated(receiver, i);
    inverse.Execute();
    for (int i = 0; i < samples; ++i)
      result(receiver, i) = scale * time[i];
  }
  return result;
}

} // namespace detail
} // namespace SPARSESPEC

#endif // DSPECM1D_POST_UNDAMPING_SH_H
