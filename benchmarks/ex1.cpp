#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

// Standard library includes
#include "config.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <complex>
#include <cmath>
#include <filesystem>

// Project-specific includes
#include <PlanetaryModel/All>
#include <DSpecM1D/Timer>
#include <DSpecM1D/All>
#include <SpectraSolver/FF>

int
main() {
  using Complex = std::complex<double>;
  using MATRIX = Eigen::MatrixXcd;

  Timer timer1;

  // --- 0. Initial Setup ---
  // get paths required for input parameters and Earth model
  std::string param_path =
      std::string(PROJECT_BUILD_DIR) + "data/params/ex1.txt";
  DSpecM::InputParametersNew paramsNew(param_path);

  // --- 4. Setup PREM and Frequency Class ---
  timer1.start();

  auto &myff = paramsNew.freqFull();

  // Added a stop here since you called start() again immediately after
  timer1.stop("Total time for PREM and Frequency setup");

  // --- 5. Source Information & Spectrum Generation ---
  DSpecM::SparseFSpec mytest;

  timer1.start();
  MATRIX vec_raw = mytest.spectra(paramsNew);
  timer1.stop("Total time for sparse frequency spectrum");

  // --- 6. Normalization ---
  double norm_factor = paramsNew.normFactor();
  vec_raw *= norm_factor;

  // --- 7. Signal Filtering ---
  double hann_w = 0.5;
  DSpecM::FilterOptions filterOptions;
  filterOptions.preTaper = 0.05;
  filterOptions.finalTaper = hann_w;
  filterOptions.passes = 1;
  filterOptions.enforceRealSignal = false;

  auto filtered = DSpecM::applyFilter(vec_raw, myff, filterOptions);
  auto &vec_filt_t = filtered.timeSeries;
  auto &a_filt = filtered.frequencySeries;

  // --- 8. Read and Process YSpec Data ---
  std::string yspec_path =
      std::string(PROJECT_BUILD_DIR) + "../../YSpec/output/yspec.lf.out.1";
  Eigen::MatrixXd yspec_t =
      DSpecM::loadYSpecTimeSeries(yspec_path, vec_filt_t.cols());

  auto filteredYSpec = DSpecM::applyFilter(yspec_t, myff, filterOptions);
  auto &vec_filt_t_yspec = filteredYSpec.timeSeries;
  auto &a_filt_yspec = filteredYSpec.frequencySeries;

  // --- 9. Output Frequency Spectrum ---
  std::string pathtofile =
      std::string(PROJECT_BUILD_DIR) + "../plotting/outputs/ex1_w.out";
  DSpecM::writeFrequencyComparison(pathtofile, paramsNew, a_filt, a_filt_yspec);

  // --- 10. Output MinEOS Time Series ---
  std::string pathtofile_mineos =
      std::string(PROJECT_BUILD_DIR) + "../plotting/outputs/ex1_t.out";
  DSpecM::writeTimeComparison(pathtofile_mineos, paramsNew, vec_filt_t,
                              vec_filt_t_yspec);

  return 0;
}