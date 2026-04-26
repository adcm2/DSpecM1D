#ifndef DSPECM1D_SIGNAL_FILTERING_H
#define DSPECM1D_SIGNAL_FILTERING_H

#include <algorithm>
#include <stdexcept>
#include <Eigen/Core>
#include <DSpecM1D/FrequencyTools>

namespace DSpecM {

/**
 * @brief Options controlling the standard DSpecM1D filtering pipeline.
 */
struct FilterOptions {
  /// Hann taper width used during the initial time-to-frequency pass.
  double preTaper = 0.05;
  /// Hann taper width used for the final frequency-domain output.
  double finalTaper = 0.5;
  /// Number of filter passes to apply.
  int passes = 1;
  /// Whether to enforce a purely real time-domain signal.
  bool enforceRealSignal = false;
};

/**
 * @brief Bundle containing both filtered time-domain and frequency-domain
 * results.
 */
struct FilterResult {
  Eigen::MatrixXd timeSeries;
  Eigen::MatrixXcd frequencySeries;
};

/**
 * @brief Applies the standard DSpecM1D filtering pipeline to raw frequency
 * spectra.
 *
 * @param rawFrequency Raw complex spectra arranged as rows = channels,
 * columns = frequencies.
 * @param freq Frequency helper used by the FFT/filter routines.
 * @param options Filtering options controlling tapers and pass count.
 * @return Filtered time-domain and frequency-domain products.
 */
inline FilterResult
applyFilter(const Eigen::MatrixXcd &rawFrequency, SpectraSolver::FreqFull &freq,
            const FilterOptions &options = {}) {
  if (options.passes < 1) {
    throw std::invalid_argument("FilterOptions.passes must be >= 1");
  }
  if (rawFrequency.cols() != freq.nt() / 2 + 1) {
    throw std::invalid_argument(
        "rawFrequency must have exactly freq.nt() / 2 + 1 columns");
  }

  FilterResult out;
  out.timeSeries = processfunctions::freq2time(rawFrequency, freq);

  Eigen::MatrixXcd workingFreq = rawFrequency;
  for (int pass = 0; pass < options.passes; ++pass) {
    workingFreq =
        processfunctions::fulltime2freq(out.timeSeries, freq, options.preTaper);
    out.timeSeries = processfunctions::filtfreq2time(workingFreq, freq,
                                                     options.enforceRealSignal);
    workingFreq = processfunctions::fulltime2freq(out.timeSeries, freq,
                                                  options.finalTaper);
  }

  out.frequencySeries = std::move(workingFreq);
  return out;
}

/**
 * @brief Applies the standard DSpecM1D filtering pipeline to raw time series.
 *
 * @param rawTime Raw time-domain input arranged as rows = channels,
 * columns = samples.
 * @param freq Frequency helper used by the FFT/filter routines.
 * @param options Filtering options controlling tapers and pass count.
 * @return Filtered time-domain and frequency-domain products.
 */
inline FilterResult
applyFilter(const Eigen::MatrixXd &rawTime, SpectraSolver::FreqFull &freq,
            const FilterOptions &options = {}) {
  if (options.passes < 1) {
    throw std::invalid_argument("FilterOptions.passes must be >= 1");
  }
  if (rawTime.cols() != freq.nt()) {
    throw std::invalid_argument(
        "rawTime must have exactly freq.nt() columns");
  }

  FilterResult out;
  out.timeSeries = rawTime;

  Eigen::MatrixXcd workingFreq;
  for (int pass = 0; pass < options.passes; ++pass) {
    workingFreq =
        processfunctions::fulltime2freq(out.timeSeries, freq, options.preTaper);
    out.timeSeries = processfunctions::filtfreq2time(workingFreq, freq,
                                                     options.enforceRealSignal);
    workingFreq = processfunctions::fulltime2freq(out.timeSeries, freq,
                                                  options.finalTaper);
  }

  out.frequencySeries = std::move(workingFreq);
  return out;
}

}   // namespace DSpecM

#endif   // DSPECM1D_SIGNAL_FILTERING_H
