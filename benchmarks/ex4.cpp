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

#include <DSpecM1D/ModelInput>
#include <DSpecM1D/Timer>
#include <DSpecM1D/All>
#include <DSpecM1D/FrequencyTools>

int
main() {
  using MatrixC = Eigen::MatrixXcd;

  Timer timer1;

  // ex4 writes a record section, so it keeps the setup explicit even though
  // the filtering stage now reuses the shared helper path.
  prem_norm<double> normClass;
  auto timeNorm = normClass.TimeNorm();
  std::string paramPath =
      std::string(PROJECT_BUILD_DIR) + "data/params/ex4.txt";
  InputParameters params(paramPath);
  std::string earthModelPath =
      std::string(PROJECT_BUILD_DIR) + "data/" + params.earth_model();

  SRInfo srInfo(params);

  double dt = params.time_step_sec();
  double tout = params.t_out() / 60.0;
  double df0 = 1.0;
  double wtb = 0.05;
  double t1 = 0.0;
  double t2 = tout + 1.0;
  int nq = 5;
  int qex = 1;

  // Build the explicit frequency helper for this record-section run.
  SpectraSolver::FreqFull myff(params.f1(), params.f2(), params.f11(),
                               params.f12(), params.f21(), params.f22(), dt,
                               tout, df0, wtb, t1, t2, qex, timeNorm);

  auto prem = EarthModels::ModelInput(earthModelPath, normClass);
  auto cmt = SourceInfo::EarthquakeCMT(params);

  // Solve for the raw spectra and then filter them back to the time domain.
  SPARSESPEC::SparseFSpec specSolver;
  timer1.start();
  MatrixC vecRaw = specSolver.spectra(myff, prem, cmt, params, nq, srInfo,
                                      params.relative_error());
  timer1.stop("Total time for sparse frequency spectrum");

  auto filterOptions = PaperExamples::makeFilterOptions(0.5, 0.01);
  auto filtered = DSpecM::applyFilter(vecRaw, myff, filterOptions);
  auto &vecFiltT = filtered.timeSeries;

  // ex4 is an acceleration record section, so scale the filtered traces
  // accordingly before writing the output grid.
  double accelNorm = prem.LengthNorm() / (prem.TimeNorm() * prem.TimeNorm());
  vecFiltT *= accelNorm;

  // Write one row per time sample, with the three components for each receiver.
  std::string pathToFile =
      std::string(PROJECT_BUILD_DIR) + "../plotting/outputs/ex4.out";
  std::ofstream file(pathToFile);
  if (!file) {
    std::cerr << "Error: unable to open output file: " << pathToFile << "\n";
    return 1;
  }

  file << std::fixed << std::setprecision(22);

  for (std::size_t idx = 0; idx < static_cast<std::size_t>(vecFiltT.cols());
       ++idx) {
    double tSec = idx * myff.dt() * prem.TimeNorm();
    file << tSec;

    for (int jidx = 0; jidx < params.num_receivers(); ++jidx) {
      file << ';' << vecFiltT(3 * jidx + 0, idx) << ';'
           << vecFiltT(3 * jidx + 1, idx) << ';' << vecFiltT(3 * jidx + 2, idx);
    }
    file << '\n';

    if (tSec > params.t_out() * 60.0) {
      break;
    }
  }
  file.close();

  return 0;
}
