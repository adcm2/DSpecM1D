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
#include <algorithm>

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

  // --- 1. Read Inputs & Earth Model ---
  // get paths required for input parameters and Earth model
  std::string paramPath =
      std::string(PROJECT_BUILD_DIR) + "data/params/ex7.txt";
  InputParameters params(paramPath);
  std::string earthModelPath =
      std::string(PROJECT_BUILD_DIR) + "data/" + params.earth_model();

  SRInfo srInfo(params);
  auto cmt = SourceInfo::EarthquakeCMT(params);

  prem_norm<double> norm_class;
  auto prem = EarthModels::ModelInput(earthModelPath, norm_class, "true");

  // --- 2. Frequency Solver Parameters ---
  int NQ = 6;
  double dt = params.time_step_sec();
  double tout = params.t_out() / 60.0;
  double df0 = 1.0;
  double wtb = 0.05;
  double t1 = 0.0;
  double t2 = tout;
  int qex = 1;

  // --- 3. Setup Frequency Class ---
  timer1.start();
  SpectraSolver::FreqFull myff(params.f1(), params.f2(), params.f11(),
                               params.f12(), params.f21(), params.f22(), dt,
                               tout, df0, wtb, t1, t2, qex, prem.TimeNorm());
  auto vecW = myff.w();
  timer1.stop("Total time for reading PREM and setting up frequency class");

  // --- 4. Compute Sparse Frequency Spectrum ---
  SPARSESPEC::SparseFSpec specSolver;

  timer1.start();
  MatrixC vecRaw = specSolver.spectra(myff, prem, cmt, params, NQ, srInfo,
                                      params.relative_error());
  timer1.stop("Total time for sparse frequency spectrum");

  // --- 5. Normalization ---
  double normFactor = 1.0;
  double accel_norm = prem.LengthNorm() / (prem.TimeNorm() * prem.TimeNorm());

  if (params.output_type() == 0) {
    normFactor = prem.LengthNorm();
  } else if (params.output_type() == 1) {
    normFactor = prem.LengthNorm() / prem.TimeNorm();
  } else if (params.output_type() == 2) {
    normFactor = accel_norm;
  }
  vecRaw *= normFactor;

  // --- 6. Base Responses (Filter) ---
  double hannW = 0.2;
  auto vecR2TB = processfunctions::freq2time(vecRaw, myff);
  auto aFilt0 = processfunctions::fulltime2freq(vecR2TB, myff, 0.05);
  auto vecFiltT = processfunctions::filtfreq2time(aFilt0, myff, false);
  auto aFilt = processfunctions::fulltime2freq(vecFiltT, myff, hannW);

  //////////////////////////////////////////////////////////////////////////////
  // --- 7. Read and Process YSpec Data ---
  std::string yspec_path = std::string(PROJECT_BUILD_DIR) + "../../YSpec/" +
                           params.output_prefix() + ".1";
  YSPECREADER::DataColumns yspec_data(yspec_path);

  std::size_t maxcoly = std::min(static_cast<std::size_t>(vecR2TB.cols()),
                                 yspec_data.getColumn1().size());

  Eigen::MatrixXd yspec_t = Eigen::MatrixXd::Zero(3, vecR2TB.cols());
  for (std::size_t idx = 0; idx < maxcoly; ++idx) {
    yspec_t(0, idx) = yspec_data.getColumn2()[idx];
    yspec_t(1, idx) = yspec_data.getColumn3()[idx];
    yspec_t(2, idx) = yspec_data.getColumn4()[idx];
  }

  auto a_filt_yspec0 = processfunctions::fulltime2freq(yspec_t, myff, 0.05);
  auto vecFiltTYSpec =
      processfunctions::filtfreq2time(a_filt_yspec0, myff, false);
  auto aFiltYSpec = processfunctions::fulltime2freq(vecFiltTYSpec, myff, hannW);

  // --- 8. Read and Process MinEOS Data ---
  std::string mineos_base =
      std::string(PROJECT_BUILD_DIR) +
      "../../mineos/DEMO/MYEX/Syndat_ASC_NOHEADER/Syndat.2000014:23:37:10.TLY.";
  MINEOSREADER::DataColumns mineos_data(mineos_base + "LHZ.ASC");
  MINEOSREADER::DataColumns mineos_data1(mineos_base + "LHN.ASC");
  MINEOSREADER::DataColumns mineos_data2(mineos_base + "LHE.ASC");

  std::size_t maxcol = std::min(static_cast<std::size_t>(vecR2TB.cols()),
                                mineos_data.getColumn1().size());

  Eigen::MatrixXd mineos_t = Eigen::MatrixXd::Zero(3, vecR2TB.cols());
  for (std::size_t idx = 0; idx < maxcol; ++idx) {
    mineos_t(0, idx) = mineos_data.getColumn2()[idx] * 1e-9;
    mineos_t(1, idx) = mineos_data1.getColumn2()[idx] * 1e-9;
    mineos_t(2, idx) = mineos_data2.getColumn2()[idx] * 1e-9;
  }

  auto a_mineos_0 = processfunctions::fulltime2freq(mineos_t, myff, 0.05);
  auto vec_filt_t_mineos =
      processfunctions::filtfreq2time(a_mineos_0, myff, false);
  auto a_filt_mineos =
      processfunctions::fulltime2freq(vec_filt_t_mineos, myff, hannW);

  // --- 9. Read and Process SpecNM Data ---
  DSpecM::FilterOptions filterOptions;
  filterOptions.preTaper = 0.05;
  filterOptions.finalTaper = hannW;
  filterOptions.passes = 1;
  filterOptions.enforceRealSignal = false;
  std::string specnm_path = std::string(PROJECT_BUILD_DIR) +
                            "../../specnm/outputs/"
                            "seismogram_ex7.txt";
  Eigen::MatrixXd specnm_t =
      DSpecM::loadSpecnmTimeSeries(specnm_path, vecFiltT.cols());

  auto filteredSpecnm = DSpecM::applyFilter(specnm_t, myff, filterOptions);
  auto &vecFiltTSpecnm = filteredSpecnm.timeSeries;
  auto &aFiltSpecnm = filteredSpecnm.frequencySeries;

  //////////////////////////////////////////////////////////////////////////////
  // --- 9. Output Frequency Spectrum ---
  std::string ptf_w =
      std::string(PROJECT_BUILD_DIR) + "../plotting/outputs/ex7_w.out";
  std::ofstream file_w(ptf_w);
  if (!file_w) {
    std::cerr << "Error: unable to open output file_w: " << ptf_w << "\n";
    return 1;
  }

  double nval = 1.0 / prem.TimeNorm();
  file_w << std::fixed << std::setprecision(16);

  for (std::size_t idx = 0; idx < myff.i2() + 100; ++idx) {
    file_w << (vecW[idx] * nval * 1000.0 / TWO_PI) << ';'
           << aFilt(0, idx).real() << ';' << aFilt(0, idx).imag() << ';'
           << std::abs(aFilt(0, idx)) << ';' << aFilt(1, idx).real() << ';'
           << aFilt(1, idx).imag() << ';' << std::abs(aFilt(1, idx)) << ';'
           << aFilt(2, idx).real() << ';' << aFilt(2, idx).imag() << ';'
           << std::abs(aFilt(2, idx)) << ';' << aFiltYSpec(0, idx).real() << ';'
           << aFiltYSpec(0, idx).imag() << ';' << std::abs(aFiltYSpec(0, idx))
           << ';' << aFiltYSpec(1, idx).real() << ';'
           << aFiltYSpec(1, idx).imag() << ';' << std::abs(aFiltYSpec(1, idx))
           << ';' << aFiltYSpec(2, idx).real() << ';'
           << aFiltYSpec(2, idx).imag() << ';' << std::abs(aFiltYSpec(2, idx))
           << ';' << a_filt_mineos(0, idx).real() << ';'
           << a_filt_mineos(0, idx).imag() << ';'
           << std::abs(a_filt_mineos(0, idx)) << ';'
           << a_filt_mineos(1, idx).real() << ';'
           << a_filt_mineos(1, idx).imag() << ';'
           << std::abs(a_filt_mineos(1, idx)) << ';'
           << a_filt_mineos(2, idx).real() << ';'
           << a_filt_mineos(2, idx).imag() << ';'
           << std::abs(a_filt_mineos(2, idx)) << ';'
           << aFiltSpecnm(0, idx).real() << ';' << aFiltSpecnm(0, idx).imag()
           << ';' << std::abs(aFiltSpecnm(0, idx)) << ';'
           << aFiltSpecnm(1, idx).real() << ';' << aFiltSpecnm(1, idx).imag()
           << ';' << std::abs(aFiltSpecnm(1, idx)) << ';'
           << aFiltSpecnm(2, idx).real() << ';' << aFiltSpecnm(2, idx).imag()
           << ';' << std::abs(aFiltSpecnm(2, idx)) << '\n';
  }
  file_w.close();

  // --- 10. Output Time Series ---
  std::string ptf_t =
      std::string(PROJECT_BUILD_DIR) + "../plotting/outputs/ex7_t.out";
  std::ofstream file_t(ptf_t);
  if (!file_t) {
    std::cerr << "Error: unable to open output file_t: " << ptf_t << "\n";
    return 1;
  }

  file_t << std::fixed << std::setprecision(16);

  for (std::size_t idx = 0;
       idx < static_cast<std::size_t>(vec_filt_t_mineos.cols()); ++idx) {
    auto tval = idx * myff.dt() * prem.TimeNorm();
    file_t << tval << ';' << vecFiltT(0, idx) << ';' << vecFiltT(1, idx) << ';'
           << vecFiltT(2, idx) << ';' << vecFiltTYSpec(0, idx) << ';'
           << vecFiltTYSpec(1, idx) << ';' << vecFiltTYSpec(2, idx) << ';'
           << vec_filt_t_mineos(0, idx) << ';' << vec_filt_t_mineos(1, idx)
           << ';' << vec_filt_t_mineos(2, idx) << ';' << vecFiltTSpecnm(0, idx)
           << ';' << vecFiltTSpecnm(1, idx) << ';' << vecFiltTSpecnm(2, idx)
           << '\n';

    if (tval > params.t_out() * 60.0) {
      break;
    }
  }
  file_t.close();

  return 0;
}