#ifndef DSPECM1D_OUTPUT_WRITERS_H
#define DSPECM1D_OUTPUT_WRITERS_H

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <DSpecM1D/FrequencyTools>
#include "InputParametersNew.h"

namespace DSpecM {

/**
 * @brief Writes a two-way frequency-domain comparison file.
 *
 * The output is a semicolon-delimited table containing the converted
 * frequency axis and the real, imaginary, and magnitude values for each
 * of the three component rows in the primary and reference spectra.
 */
inline void
writeFrequencyComparison(const std::string &path,
                         const std::vector<double> &angularFrequency,
                         const SpectraSolver::FreqFull &freq, double nval,
                         const Eigen::MatrixXcd &primary,
                         const Eigen::MatrixXcd &reference,
                         int extraFrequencies = 100, int precision = 16) {
  std::ofstream file(path);
  if (!file) {
    throw std::runtime_error("Error: unable to open output file: " + path);
  }

  constexpr double twoPi = 2.0 * 3.14159265358979323846;
  const std::size_t requested =
      static_cast<std::size_t>(std::max(0, freq.i2() + extraFrequencies));
  const std::size_t maxCount =
      std::min({requested, angularFrequency.size(),
                static_cast<std::size_t>(primary.cols()),
                static_cast<std::size_t>(reference.cols())});

  file << std::fixed << std::setprecision(precision);
  for (std::size_t idx = 0; idx < maxCount; ++idx) {
    file << (angularFrequency[idx] * nval * 1000.0 / twoPi) << ';'
         << primary(0, static_cast<Eigen::Index>(idx)).real() << ';'
         << primary(0, static_cast<Eigen::Index>(idx)).imag() << ';'
         << std::abs(primary(0, static_cast<Eigen::Index>(idx))) << ';'
         << primary(1, static_cast<Eigen::Index>(idx)).real() << ';'
         << primary(1, static_cast<Eigen::Index>(idx)).imag() << ';'
         << std::abs(primary(1, static_cast<Eigen::Index>(idx))) << ';'
         << primary(2, static_cast<Eigen::Index>(idx)).real() << ';'
         << primary(2, static_cast<Eigen::Index>(idx)).imag() << ';'
         << std::abs(primary(2, static_cast<Eigen::Index>(idx))) << ';'
         << reference(0, static_cast<Eigen::Index>(idx)).real() << ';'
         << reference(0, static_cast<Eigen::Index>(idx)).imag() << ';'
         << std::abs(reference(0, static_cast<Eigen::Index>(idx))) << ';'
         << reference(1, static_cast<Eigen::Index>(idx)).real() << ';'
         << reference(1, static_cast<Eigen::Index>(idx)).imag() << ';'
         << std::abs(reference(1, static_cast<Eigen::Index>(idx))) << ';'
         << reference(2, static_cast<Eigen::Index>(idx)).real() << ';'
         << reference(2, static_cast<Eigen::Index>(idx)).imag() << ';'
         << std::abs(reference(2, static_cast<Eigen::Index>(idx))) << '\n';
  }
}

/**
 * @brief Writes a three-way frequency-domain comparison file.
 */
inline void
writeFrequencyComparison(const std::string &path,
                         const std::vector<double> &angularFrequency,
                         const SpectraSolver::FreqFull &freq, double nval,
                         const Eigen::MatrixXcd &primary,
                         const Eigen::MatrixXcd &reference,
                         const Eigen::MatrixXcd &secondaryReference,
                         int extraFrequencies = 100, int precision = 16) {
  std::ofstream file(path);
  if (!file) {
    throw std::runtime_error("Error: unable to open output file: " + path);
  }

  constexpr double twoPi = 2.0 * 3.14159265358979323846;
  const std::size_t requested =
      static_cast<std::size_t>(std::max(0, freq.i2() + extraFrequencies));
  const std::size_t maxCount =
      std::min({requested, angularFrequency.size(),
                static_cast<std::size_t>(primary.cols()),
                static_cast<std::size_t>(reference.cols()),
                static_cast<std::size_t>(secondaryReference.cols())});

  file << std::fixed << std::setprecision(precision);
  for (std::size_t idx = 0; idx < maxCount; ++idx) {
    file << (angularFrequency[idx] * nval * 1000.0 / twoPi) << ';'
         << primary(0, static_cast<Eigen::Index>(idx)).real() << ';'
         << primary(0, static_cast<Eigen::Index>(idx)).imag() << ';'
         << std::abs(primary(0, static_cast<Eigen::Index>(idx))) << ';'
         << primary(1, static_cast<Eigen::Index>(idx)).real() << ';'
         << primary(1, static_cast<Eigen::Index>(idx)).imag() << ';'
         << std::abs(primary(1, static_cast<Eigen::Index>(idx))) << ';'
         << primary(2, static_cast<Eigen::Index>(idx)).real() << ';'
         << primary(2, static_cast<Eigen::Index>(idx)).imag() << ';'
         << std::abs(primary(2, static_cast<Eigen::Index>(idx))) << ';'
         << reference(0, static_cast<Eigen::Index>(idx)).real() << ';'
         << reference(0, static_cast<Eigen::Index>(idx)).imag() << ';'
         << std::abs(reference(0, static_cast<Eigen::Index>(idx))) << ';'
         << reference(1, static_cast<Eigen::Index>(idx)).real() << ';'
         << reference(1, static_cast<Eigen::Index>(idx)).imag() << ';'
         << std::abs(reference(1, static_cast<Eigen::Index>(idx))) << ';'
         << reference(2, static_cast<Eigen::Index>(idx)).real() << ';'
         << reference(2, static_cast<Eigen::Index>(idx)).imag() << ';'
         << std::abs(reference(2, static_cast<Eigen::Index>(idx))) << ';'
         << secondaryReference(0, static_cast<Eigen::Index>(idx)).real() << ';'
         << secondaryReference(0, static_cast<Eigen::Index>(idx)).imag() << ';'
         << std::abs(secondaryReference(0, static_cast<Eigen::Index>(idx)))
         << ';' << secondaryReference(1, static_cast<Eigen::Index>(idx)).real()
         << ';' << secondaryReference(1, static_cast<Eigen::Index>(idx)).imag()
         << ';'
         << std::abs(secondaryReference(1, static_cast<Eigen::Index>(idx)))
         << ';' << secondaryReference(2, static_cast<Eigen::Index>(idx)).real()
         << ';' << secondaryReference(2, static_cast<Eigen::Index>(idx)).imag()
         << ';'
         << std::abs(secondaryReference(2, static_cast<Eigen::Index>(idx)))
         << '\n';
  }
}

/**
 * @brief Convenience overload that derives frequency metadata from
 * `InputParametersNew`.
 */
inline void
writeFrequencyComparison(const std::string &path,
                         const InputParametersNew &paramsNew,
                         const Eigen::MatrixXcd &primary,
                         const Eigen::MatrixXcd &reference,
                         int extraFrequencies = 100, int precision = 16) {
  const auto &freq = paramsNew.freqFull();
  const auto vecW = freq.w();
  const double nval = 1.0 / paramsNew.earthModel().TimeNorm();
  writeFrequencyComparison(path, vecW, freq, nval, primary, reference,
                           extraFrequencies, precision);
}

/**
 * @brief Convenience overload that derives frequency metadata from
 * `InputParametersNew` for three-way comparisons.
 */
inline void
writeFrequencyComparison(const std::string &path,
                         const InputParametersNew &paramsNew,
                         const Eigen::MatrixXcd &primary,
                         const Eigen::MatrixXcd &reference,
                         const Eigen::MatrixXcd &secondaryReference,
                         int extraFrequencies = 100, int precision = 16) {
  const auto &freq = paramsNew.freqFull();
  const auto vecW = freq.w();
  const double nval = 1.0 / paramsNew.earthModel().TimeNorm();
  writeFrequencyComparison(path, vecW, freq, nval, primary, reference,
                           secondaryReference, extraFrequencies, precision);
}

/**
 * @brief Writes a two-way time-domain comparison file.
 */
inline void
writeTimeComparison(const std::string &path,
                    const SpectraSolver::FreqFull &freq, double timeNorm,
                    double tOutMinutes, const Eigen::MatrixXd &primary,
                    const Eigen::MatrixXd &reference, int precision = 16) {
  std::ofstream file(path);
  if (!file) {
    throw std::runtime_error("Error: unable to open output file: " + path);
  }

  const std::size_t maxCount =
      std::min(static_cast<std::size_t>(primary.cols()),
               static_cast<std::size_t>(reference.cols()));

  file << std::fixed << std::setprecision(precision);
  for (std::size_t idx = 0; idx < maxCount; ++idx) {
    const double currentTime = idx * freq.dt() * timeNorm;
    file << currentTime << ';' << primary(0, static_cast<Eigen::Index>(idx))
         << ';' << primary(1, static_cast<Eigen::Index>(idx)) << ';'
         << primary(2, static_cast<Eigen::Index>(idx)) << ';'
         << reference(0, static_cast<Eigen::Index>(idx)) << ';'
         << reference(1, static_cast<Eigen::Index>(idx)) << ';'
         << reference(2, static_cast<Eigen::Index>(idx)) << '\n';

    if (currentTime > tOutMinutes * 60.0) {
      break;
    }
  }
}

/**
 * @brief Convenience overload that derives the time metadata from
 * `InputParametersNew`.
 */
inline void
writeTimeComparison(const std::string &path,
                    const InputParametersNew &paramsNew,
                    const Eigen::MatrixXd &primary,
                    const Eigen::MatrixXd &reference, int precision = 16) {
  writeTimeComparison(
      path, paramsNew.freqFull(), paramsNew.earthModel().TimeNorm(),
      paramsNew.inputParameters().t_out(), primary, reference, precision);
}

/**
 * @brief Writes a standalone time-domain seismogram table.
 *
 * The first column is time in seconds. Each remaining column corresponds to
 * one row from the input matrix, making the helper suitable for single- or
 * multi-receiver outputs without assuming a comparison layout.
 */
inline void
writeTimeSeriesOnGrid(const std::string &path, double timeStep,
                      double timeNorm, double tOutMinutes,
                      const Eigen::MatrixXd &timeSeries, int precision = 16) {
  std::ofstream file(path);
  if (!file) {
    throw std::runtime_error("Error: unable to open output file: " + path);
  }

  file << std::fixed << std::setprecision(precision);
  for (Eigen::Index idx = 0; idx < timeSeries.cols(); ++idx) {
    const double currentTime =
        static_cast<double>(idx) * timeStep * timeNorm;
    file << currentTime;
    for (Eigen::Index row = 0; row < timeSeries.rows(); ++row) {
      file << ';' << timeSeries(row, idx);
    }
    file << '\n';

    if (currentTime > tOutMinutes * 60.0) {
      break;
    }
  }
}

inline void
writeTimeSeries(const std::string &path, const SpectraSolver::FreqFull &freq,
                double timeNorm, double tOutMinutes,
                const Eigen::MatrixXd &timeSeries, int precision = 16) {
  writeTimeSeriesOnGrid(path, freq.dt(), timeNorm, tOutMinutes, timeSeries,
                        precision);
}

/**
 * @brief Convenience overload that derives the time metadata from
 * `InputParametersNew`.
 */
inline void
writeTimeSeries(const std::string &path, const InputParametersNew &paramsNew,
                const Eigen::MatrixXd &timeSeries, int precision = 16) {
  writeTimeSeries(path, paramsNew.freqFull(), paramsNew.earthModel().TimeNorm(),
                  paramsNew.inputParameters().t_out(), timeSeries, precision);
}

/**
 * @brief Writes the unfiltered raw time series with a strict duration cutoff.
 *
 * Each output row is a time sample. The first column is physical time in
 * seconds; remaining columns are grouped receiver-major, with Z, N, E in each
 * receiver block. Unlike the legacy writer, the duration is inclusive and a
 * sample beyond it is never written.
 */
inline void
writeRawTimeSeries(const std::string &path, double sampleIntervalSeconds,
                   double durationSeconds, const Eigen::MatrixXd &timeSeries) {
  if (!(sampleIntervalSeconds > 0.0) || !(durationSeconds >= 0.0))
    throw std::invalid_argument("Raw time output requires a positive sample interval and non-negative duration.");
  std::ofstream file(path);
  if (!file)
    throw std::runtime_error("Error: unable to open output file: " + path);

  file << std::scientific
       << std::setprecision(std::numeric_limits<double>::max_digits10);
  for (Eigen::Index sample = 0; sample < timeSeries.cols(); ++sample) {
    double time = static_cast<double>(sample) * sampleIntervalSeconds;
    if (time > durationSeconds) {
      const double oneUlpLimit = std::nextafter(
          durationSeconds, std::numeric_limits<double>::infinity());
      const double twoUlpLimit = std::nextafter(
          oneUlpLimit, std::numeric_limits<double>::infinity());
      if (time > twoUlpLimit)
        break;
      // Keep the exported axis strictly within the requested duration even
      // when multiplication rounds an exact endpoint one or two ULPs high.
      time = durationSeconds;
    }
    file << time;
    for (Eigen::Index row = 0; row < timeSeries.rows(); ++row)
      file << ';' << timeSeries(row, sample);
    file << '\n';
  }
}

/**
 * @brief Writes damped raw complex spectra on the physical one-sided grid.
 *
 * The first column is frequency in Hz, followed by real and imaginary values
 * for each receiver's Z, N, E components. The supplied amplitude scale
 * converts solver output to observable SI units times seconds.
 */
inline void
writeRawFrequencySeries(const std::string &path,
                        const SpectraSolver::FreqFull &freq,
                        double timeNorm, double amplitudeScale,
                        const Eigen::MatrixXcd &spectrum) {
  const Eigen::Index expectedBins = freq.nt() / 2 + 1;
  if (spectrum.cols() != expectedBins)
    throw std::invalid_argument(
        "Raw frequency output must contain the full one-sided frequency grid.");
  if (spectrum.rows() % 3 != 0)
    throw std::invalid_argument(
        "Raw frequency output rows must contain receiver-major Z/N/E triples.");

  std::ofstream file(path);
  if (!file)
    throw std::runtime_error("Error: unable to open output file: " + path);

  file << std::scientific
       << std::setprecision(std::numeric_limits<double>::max_digits10);
  for (Eigen::Index bin = 0; bin < expectedBins; ++bin) {
    file << (freq.f(static_cast<int>(bin)) / timeNorm);
    for (Eigen::Index row = 0; row < spectrum.rows(); ++row) {
      const auto value = spectrum(row, bin) * amplitudeScale;
      file << ';' << value.real() << ';' << value.imag();
    }
    file << '\n';
  }
}

}   // namespace DSpecM

#endif   // DSPECM1D_OUTPUT_WRITERS_H
