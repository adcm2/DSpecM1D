#ifndef DSPECM1D_SIGNAL_FILTERING_H
#define DSPECM1D_SIGNAL_FILTERING_H

#include <algorithm>
#include <stdexcept>
#include <Eigen/Core>
#include <SpectraSolver/FF>

namespace DSpecM {

struct FilterOptions {
  double preTaper = 0.05;
  double finalTaper = 0.5;
  int passes = 1;
  bool enforceRealSignal = false;
};

struct FilterResult {
  Eigen::MatrixXd timeSeries;
  Eigen::MatrixXcd frequencySeries;
};

inline FilterResult
applyFilter(const Eigen::MatrixXcd &rawFrequency, SpectraSolver::FreqFull &freq,
            const FilterOptions &options = {}) {
  if (options.passes < 1) {
    throw std::invalid_argument("FilterOptions.passes must be >= 1");
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

inline FilterResult
applyFilter(const Eigen::MatrixXd &rawTime, SpectraSolver::FreqFull &freq,
            const FilterOptions &options = {}) {
  if (options.passes < 1) {
    throw std::invalid_argument("FilterOptions.passes must be >= 1");
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
