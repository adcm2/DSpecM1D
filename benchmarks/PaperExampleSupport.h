#ifndef DSPECM1D_PAPER_EXAMPLE_SUPPORT_H
#define DSPECM1D_PAPER_EXAMPLE_SUPPORT_H

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <DSpecM1D/FrequencyTools>
#include <DSpecM1D/All>

namespace PaperExamples {

/**
 * @brief Returns the SI scaling implied by the legacy `output_type` flag.
 *
 * This mirrors `InputParametersNew::normFactor()` for examples that still use
 * the legacy `InputParameters` workflow directly.
 */
template <class ModelType>
double
legacyNormFactor(const InputParameters &params, const ModelType &model) {
  if (params.output_type() == 0) {
    return model.LengthNorm();
  }
  if (params.output_type() == 1) {
    return model.LengthNorm() / model.TimeNorm();
  }
  if (params.output_type() == 2) {
    return model.LengthNorm() / (model.TimeNorm() * model.TimeNorm());
  }
  return 1.0;
}

/**
 * @brief Builds the standard filtering options used by the paper examples.
 */
inline DSpecM::FilterOptions
makeFilterOptions(double finalTaper, double preTaper = 0.05, int passes = 1,
                  bool enforceRealSignal = false) {
  DSpecM::FilterOptions options;
  options.preTaper = preTaper;
  options.finalTaper = finalTaper;
  options.passes = passes;
  options.enforceRealSignal = enforceRealSignal;
  return options;
}

/**
 * @brief Resolves the in-repo YSpec reference path from an output prefix.
 */
inline std::string
yspecReferencePath(const std::string &buildDir, const std::string &outputPrefix) {
  return buildDir + "data/reference/yspec/" +
         std::filesystem::path(outputPrefix).filename().string() + ".1";
}

/**
 * @brief Loads paired YSpec and Mineos comparison traces sharing one Mineos
 * base path.
 */
inline DSpecM::ReferenceTimeSeries
loadReferenceTimeSeriesWithMineosBase(const std::string &yspecPath,
                                      const std::string &mineosBase, int ncols,
                                      double mineosScale = 1e-9) {
  return DSpecM::loadReferenceTimeSeries(yspecPath, mineosBase + "LHZ.ASC",
                                         mineosBase + "LHN.ASC",
                                         mineosBase + "LHE.ASC", ncols,
                                         mineosScale);
}

/**
 * @brief Writes a four-way frequency comparison used by the paper examples.
 *
 * The output contains the converted frequency axis followed by the real,
 * imaginary, and magnitude values for the three components of:
 * 1. the primary DSpecM1D result
 * 2. the first reference result
 * 3. the second reference result
 * 4. the third reference result
 */
inline void
writeFourWayFrequencyComparison(const std::string &path,
                                const std::vector<double> &angularFrequency,
                                const SpectraSolver::FreqFull &freq,
                                double nval, const Eigen::MatrixXcd &primary,
                                const Eigen::MatrixXcd &reference,
                                const Eigen::MatrixXcd &secondaryReference,
                                const Eigen::MatrixXcd &tertiaryReference,
                                int extraFrequencies = 100,
                                int precision = 16) {
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
                static_cast<std::size_t>(secondaryReference.cols()),
                static_cast<std::size_t>(tertiaryReference.cols())});

  file << std::fixed << std::setprecision(precision);
  for (std::size_t idx = 0; idx < maxCount; ++idx) {
    const auto eidx = static_cast<Eigen::Index>(idx);
    file << (angularFrequency[idx] * nval * 1000.0 / twoPi) << ';'
         << primary(0, eidx).real() << ';' << primary(0, eidx).imag() << ';'
         << std::abs(primary(0, eidx)) << ';' << primary(1, eidx).real() << ';'
         << primary(1, eidx).imag() << ';' << std::abs(primary(1, eidx)) << ';'
         << primary(2, eidx).real() << ';' << primary(2, eidx).imag() << ';'
         << std::abs(primary(2, eidx)) << ';' << reference(0, eidx).real()
         << ';' << reference(0, eidx).imag() << ';'
         << std::abs(reference(0, eidx)) << ';' << reference(1, eidx).real()
         << ';' << reference(1, eidx).imag() << ';'
         << std::abs(reference(1, eidx)) << ';' << reference(2, eidx).real()
         << ';' << reference(2, eidx).imag() << ';'
         << std::abs(reference(2, eidx)) << ';'
         << secondaryReference(0, eidx).real() << ';'
         << secondaryReference(0, eidx).imag() << ';'
         << std::abs(secondaryReference(0, eidx)) << ';'
         << secondaryReference(1, eidx).real() << ';'
         << secondaryReference(1, eidx).imag() << ';'
         << std::abs(secondaryReference(1, eidx)) << ';'
         << secondaryReference(2, eidx).real() << ';'
         << secondaryReference(2, eidx).imag() << ';'
         << std::abs(secondaryReference(2, eidx)) << ';'
         << tertiaryReference(0, eidx).real() << ';'
         << tertiaryReference(0, eidx).imag() << ';'
         << std::abs(tertiaryReference(0, eidx)) << ';'
         << tertiaryReference(1, eidx).real() << ';'
         << tertiaryReference(1, eidx).imag() << ';'
         << std::abs(tertiaryReference(1, eidx)) << ';'
         << tertiaryReference(2, eidx).real() << ';'
         << tertiaryReference(2, eidx).imag() << ';'
         << std::abs(tertiaryReference(2, eidx)) << '\n';
  }
}

/**
 * @brief Writes a four-way time-domain comparison used by the paper examples.
 *
 * @param timeOffsetSeconds Optional offset applied to the output time axis,
 * useful for source-time-function-shifted comparisons.
 */
inline void
writeFourWayTimeComparison(const std::string &path,
                           const SpectraSolver::FreqFull &freq,
                           double timeNorm, double tOutMinutes,
                           const Eigen::MatrixXd &primary,
                           const Eigen::MatrixXd &reference,
                           const Eigen::MatrixXd &secondaryReference,
                           const Eigen::MatrixXd &tertiaryReference,
                           double timeOffsetSeconds = 0.0,
                           int precision = 16) {
  std::ofstream file(path);
  if (!file) {
    throw std::runtime_error("Error: unable to open output file: " + path);
  }

  const std::size_t maxCount =
      std::min({static_cast<std::size_t>(primary.cols()),
                static_cast<std::size_t>(reference.cols()),
                static_cast<std::size_t>(secondaryReference.cols()),
                static_cast<std::size_t>(tertiaryReference.cols())});

  file << std::fixed << std::setprecision(precision);
  for (std::size_t idx = 0; idx < maxCount; ++idx) {
    const auto eidx = static_cast<Eigen::Index>(idx);
    const double currentTime = idx * freq.dt() * timeNorm;

    file << (currentTime + timeOffsetSeconds) << ';' << primary(0, eidx) << ';'
         << primary(1, eidx) << ';' << primary(2, eidx) << ';'
         << reference(0, eidx) << ';' << reference(1, eidx) << ';'
         << reference(2, eidx) << ';' << secondaryReference(0, eidx) << ';'
         << secondaryReference(1, eidx) << ';' << secondaryReference(2, eidx)
         << ';' << tertiaryReference(0, eidx) << ';'
         << tertiaryReference(1, eidx) << ';' << tertiaryReference(2, eidx)
         << '\n';

    if (currentTime > tOutMinutes * 60.0) {
      break;
    }
  }
}

}   // namespace PaperExamples

#endif   // DSPECM1D_PAPER_EXAMPLE_SUPPORT_H
