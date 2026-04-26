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
#include <algorithm>

#include <DSpecM1D/ModelInput>
#include <DSpecM1D/Timer>
#include <DSpecM1D/All>
#include <DSpecM1D/FrequencyTools>

int
main() {
  using Complex = std::complex<double>;
  using MatrixC = Eigen::MatrixXcd;

  Timer timer1;

  // Legacy explicit setup: parse the example input, build the model, and
  // prepare the source/receiver geometry used in the paper comparison.
  std::string paramPath =
      std::string(PROJECT_BUILD_DIR) + "data/params/ex7.txt";
  InputParameters params(paramPath);
  std::string earthModelPath =
      std::string(PROJECT_BUILD_DIR) + "data/" + params.earth_model();

  SRInfo srInfo(params);
  auto cmt = SourceInfo::EarthquakeCMT(params);

  prem_norm<double> normClass;
  auto prem = EarthModels::ModelInput(earthModelPath, normClass);

  int nq = 6;
  double dt = params.time_step_sec();
  double tout = params.t_out() / 60.0;
  double df0 = 1.0;
  double wtb = 0.05;
  double t1 = 0.0;
  double t2 = tout;
  int qex = 1;

  // Build the frequency helper used by both the forward solve and the
  // comparison-data filtering path.
  timer1.start();
  SpectraSolver::FreqFull myff(params.f1(), params.f2(), params.f11(),
                               params.f12(), params.f21(), params.f22(), dt,
                               tout, df0, wtb, t1, t2, qex, prem.TimeNorm());
  auto vecW = myff.w();
  timer1.stop("Total time for reading PREM and setting up frequency class");

  // Compute the raw DSpecM1D spectra for this benchmark configuration.
  SPARSESPEC::SparseFSpec specSolver;
  timer1.start();
  MatrixC vecRaw = specSolver.spectra(myff, prem, cmt, params, nq, srInfo,
                                      params.relative_error());
  timer1.stop("Total time for sparse frequency spectrum");

  // Convert the spectra into the requested SI quantity before filtering.
  double normFactor = PaperExamples::legacyNormFactor(params, prem);
  vecRaw *= normFactor;

  // Use the standard filtering pipeline for both the model result and the
  // comparison traces so they are compared on the same footing.
  double hannW = 0.2;
  auto filterOptions = PaperExamples::makeFilterOptions(hannW);
  auto filtered = DSpecM::applyFilter(vecRaw, myff, filterOptions);
  auto &vecFiltT = filtered.timeSeries;
  auto &aFilt = filtered.frequencySeries;

  std::string yspecPath = PaperExamples::yspecReferencePath(
      std::string(PROJECT_BUILD_DIR), params.output_prefix());
  std::string mineosBase =
      std::string(PROJECT_BUILD_DIR) +
      "data/reference/mineos/noheader/Syndat.2000014:23:37:10.TLY.";
  auto referenceSeries = PaperExamples::loadReferenceTimeSeriesWithMineosBase(
      yspecPath, mineosBase, vecFiltT.cols());
  auto filteredYSpec =
      DSpecM::applyFilter(referenceSeries.yspecTime, myff, filterOptions);
  auto &vecFiltTYSpec = filteredYSpec.timeSeries;
  auto &aFiltYSpec = filteredYSpec.frequencySeries;

  auto filteredMineos =
      DSpecM::applyFilter(referenceSeries.mineosTime, myff, filterOptions);
  auto &vecFiltTMineos = filteredMineos.timeSeries;
  auto &aFiltMineos = filteredMineos.frequencySeries;

  std::string specnmPath = std::string(PROJECT_BUILD_DIR) +
                           "data/reference/specnm/"
                           "seismogram_ex7.txt";
  Eigen::MatrixXd specnmTime =
      DSpecM::loadSpecnmTimeSeries(specnmPath, vecFiltT.cols());

  auto filteredSpecnm = DSpecM::applyFilter(specnmTime, myff, filterOptions);
  auto &vecFiltTSpecnm = filteredSpecnm.timeSeries;
  auto &aFiltSpecnm = filteredSpecnm.frequencySeries;

  // Write the four-way paper comparison in both frequency and time domains.
  std::string pathToFreqFile =
      std::string(PROJECT_BUILD_DIR) + "../plotting/outputs/ex7_w.out";
  double nval = 1.0 / prem.TimeNorm();
  PaperExamples::writeFourWayFrequencyComparison(
      pathToFreqFile, vecW, myff, nval, aFilt, aFiltYSpec, aFiltMineos,
      aFiltSpecnm);

  std::string pathToTimeFile =
      std::string(PROJECT_BUILD_DIR) + "../plotting/outputs/ex7_t.out";
  PaperExamples::writeFourWayTimeComparison(
      pathToTimeFile, myff, prem.TimeNorm(), params.t_out(), vecFiltT,
      vecFiltTYSpec, vecFiltTMineos, vecFiltTSpecnm);

  return 0;
}
