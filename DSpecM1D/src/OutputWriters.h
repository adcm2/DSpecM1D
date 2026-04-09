#ifndef DSPECM1D_OUTPUT_WRITERS_H
#define DSPECM1D_OUTPUT_WRITERS_H

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <SpectraSolver/FF>
#include "InputParametersNew.h"

namespace DSpecM {

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

inline void
writeTimeComparison(const std::string &path,
                    const InputParametersNew &paramsNew,
                    const Eigen::MatrixXd &primary,
                    const Eigen::MatrixXd &reference, int precision = 16) {
  writeTimeComparison(
      path, paramsNew.freqFull(), paramsNew.earthModel().TimeNorm(),
      paramsNew.inputParameters().t_out(), primary, reference, precision);
}

}   // namespace DSpecM

#endif   // DSPECM1D_OUTPUT_WRITERS_H
