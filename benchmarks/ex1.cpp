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
  using MatrixC = Eigen::MatrixXcd;

  Timer timer1;

  // --- 0. Initial Setup ---
  // get paths required for input parameters and Earth model
  std::string paramPath =
      std::string(PROJECT_BUILD_DIR) + "data/params/ex1.txt";
  DSpecM::InputParametersNew paramsNew(paramPath);

  // --- 4. Setup PREM and Frequency Class ---
  timer1.start();

  auto &myff = paramsNew.freqFull();

  // Added a stop here since you called start() again immediately after
  timer1.stop("Total time for PREM and Frequency setup");

  // --- 5. Source Information & Spectrum Generation ---
  DSpecM::SparseFSpec specSolver;

  timer1.start();
  MatrixC vecRaw = specSolver.spectra(paramsNew);
  timer1.stop("Total time for sparse frequency spectrum");

  // --- 6. Normalization ---
  double normFactor = paramsNew.normFactor();
  vecRaw *= normFactor;

  // --- 7. Signal Filtering ---
  double hannW = 0.5;
  DSpecM::FilterOptions filterOptions;
  filterOptions.preTaper = 0.05;
  filterOptions.finalTaper = hannW;
  filterOptions.passes = 1;
  filterOptions.enforceRealSignal = false;

  auto filtered = DSpecM::applyFilter(vecRaw, myff, filterOptions);
  auto &vecFiltT = filtered.timeSeries;
  auto &aFilt = filtered.frequencySeries;

  // --- 8. Read and Process YSpec Data ---
  std::string yspec_path =
      std::string(PROJECT_BUILD_DIR) + "../../YSpec/output/yspec.lf.out.1";
  Eigen::MatrixXd yspec_t =
      DSpecM::loadYSpecTimeSeries(yspec_path, vecFiltT.cols());

  auto filteredYSpec = DSpecM::applyFilter(yspec_t, myff, filterOptions);
  auto &vecFiltTYSpec = filteredYSpec.timeSeries;
  auto &aFiltYSpec = filteredYSpec.frequencySeries;

  // --- 9. Read and Process SpecNM Data ---
  std::string specnm_path = std::string(PROJECT_BUILD_DIR) +
                            "../../specnm/outputs/"
                            "seismogram_lf_traces_semicolon.txt";
  Eigen::MatrixXd specnm_t =
      DSpecM::loadSpecnmTimeSeries(specnm_path, vecFiltT.cols());

  auto filteredSpecnm = DSpecM::applyFilter(specnm_t, myff, filterOptions);
  auto &vecFiltTSpecnm = filteredSpecnm.timeSeries;
  auto &aFiltSpecnm = filteredSpecnm.frequencySeries;

  // --- 10. Output Frequency Spectrum ---
  std::string pathToFile =
      std::string(PROJECT_BUILD_DIR) + "../plotting/outputs/ex1_w.out";
  DSpecM::writeFrequencyComparison(pathToFile, paramsNew, aFilt, aFiltYSpec,
                                   aFiltSpecnm);

  // --- 11. Output MinEOS Time Series ---
  std::string pathToFileMineos =
      std::string(PROJECT_BUILD_DIR) + "../plotting/outputs/ex1_t.out";
  DSpecM::writeTimeComparison(pathToFileMineos, paramsNew, vecFiltT,
                              vecFiltTSpecnm);

  return 0;
}
