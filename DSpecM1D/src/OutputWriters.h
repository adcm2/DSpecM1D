#ifndef DSPECM1D_OUTPUT_WRITERS_H
#define DSPECM1D_OUTPUT_WRITERS_H

#include <algorithm>
#include <fstream>
#include <iomanip>
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
writeTimeSeries(const std::string &path, const SpectraSolver::FreqFull &freq,
                double timeNorm, double tOutMinutes,
                const Eigen::MatrixXd &timeSeries, int precision = 16) {
  std::ofstream file(path);
  if (!file) {
    throw std::runtime_error("Error: unable to open output file: " + path);
  }

  file << std::fixed << std::setprecision(precision);
  for (Eigen::Index idx = 0; idx < timeSeries.cols(); ++idx) {
    const double currentTime =
        static_cast<double>(idx) * freq.dt() * timeNorm;
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

}   // namespace DSpecM

#endif   // DSPECM1D_OUTPUT_WRITERS_H
