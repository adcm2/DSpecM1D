#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#include "config.h"
#include "PaperExampleSupport.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <complex>
#include <cmath>
#include <filesystem>

#include <DSpecM1D/ModelInput>
#include <DSpecM1D/Timer>
#include <DSpecM1D/All>
#include <DSpecM1D/FrequencyTools>

int
main() {
  using Complex = std::complex<double>;
  using MatrixC = Eigen::MatrixXcd;

  Timer timer1;

  // Build the preferred high-level workflow context from the example input.
  std::string paramPath =
      std::string(PROJECT_BUILD_DIR) + "data/params/ex1.txt";
  InputParametersNew paramsNew(paramPath);

  // Reuse the pre-built frequency helper carried by InputParametersNew.
  timer1.start();
  auto &myff = paramsNew.freqFull();
  timer1.stop("Total time for PREM and Frequency setup");

  // Solve for the raw complex spectra using the preferred release-facing API.
  SPARSESPEC::SparseFSpec specSolver;
  timer1.start();
  MatrixC vecRaw = specSolver.spectra(paramsNew);
  timer1.stop("Total time for sparse frequency spectrum");

  // Convert the non-dimensional spectra into the requested SI output type.
  double normFactor = paramsNew.normFactor();
  vecRaw *= normFactor;

  // Apply the standard DSpecM1D filtering pipeline to the modelled spectra.
  double hannW = 0.5;
  auto filterOptions = PaperExamples::makeFilterOptions(hannW);
  auto filtered = DSpecM::applyFilter(vecRaw, myff, filterOptions);
  auto &vecFiltT = filtered.timeSeries;
  auto &aFilt = filtered.frequencySeries;

  // Load the paper reference traces and run them through the same filter chain.
  std::string yspecPath =
      std::string(PROJECT_BUILD_DIR) + "data/reference/yspec/yspec.lf.out.1";
  Eigen::MatrixXd yspecTime =
      DSpecM::loadYSpecTimeSeries(yspecPath, vecFiltT.cols());
  auto filteredYSpec = DSpecM::applyFilter(yspecTime, myff, filterOptions);
  auto &vecFiltTYSpec = filteredYSpec.timeSeries;
  auto &aFiltYSpec = filteredYSpec.frequencySeries;

  std::string specnmPath = std::string(PROJECT_BUILD_DIR) +
                           "data/reference/specnm/"
                           "seismogram_lf_traces_semicolon.txt";
  Eigen::MatrixXd specnmTime =
      DSpecM::loadSpecnmTimeSeries(specnmPath, vecFiltT.cols());
  auto filteredSpecnm = DSpecM::applyFilter(specnmTime, myff, filterOptions);
  auto &vecFiltTSpecnm = filteredSpecnm.timeSeries;
  auto &aFiltSpecnm = filteredSpecnm.frequencySeries;

  // Write the frequency- and time-domain comparison products used for plotting.
  std::string pathToFile =
      std::string(PROJECT_BUILD_DIR) + "../plotting/outputs/ex1_w.out";
  DSpecM::writeFrequencyComparison(pathToFile, paramsNew, aFilt, aFiltYSpec,
                                   aFiltSpecnm);

  std::string pathToMineosFile =
      std::string(PROJECT_BUILD_DIR) + "../plotting/outputs/ex1_t.out";
  DSpecM::writeTimeComparison(pathToMineosFile, paramsNew, vecFiltT,
                              vecFiltTSpecnm);

  return 0;
}
