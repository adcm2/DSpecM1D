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
  InputParameters params(param_path);
  SRInfo sr_info(params);
  std::string earth_model_path =
      std::string(PROJECT_BUILD_DIR) + "data/" + params.earth_model();

  // --- 2. Spectral Element Method (SEM) Parameters ---
  int lval = params.lmax();
  int NQ = 5;
  double maxstep = 0.05;

  // --- 3. Frequency Solver Parameters ---
  double dt = params.time_step_sec();
  double tout = params.t_out() / 60.0;
  double df0 = 1.0;
  double wtb = 0.05;
  double t1 = 0.0;
  double t2 = tout;
  int qex = 1;

  // --- 4. Setup PREM and Frequency Class ---
  timer1.start();

  prem_norm<double> norm_class;
  auto prem = EarthModels::ModelInput(earth_model_path, norm_class, "true");

  SpectraSolver::FreqFull myff(params.f1(), params.f2(), params.f11(),
                               params.f12(), params.f21(), params.f22(), dt,
                               tout, df0, wtb, t1, t2, qex, prem.TimeNorm());

  auto vec_w = myff.w();

  // Added a stop here since you called start() again immediately after
  timer1.stop("Total time for PREM and Frequency setup");

  // --- 5. Setup SEM Class ---
  timer1.start();
  std::cout << "Setting up SEM class...\n";
  Full1D::specsem sem(prem, maxstep, NQ, lval);
  timer1.stop("Total time for setting up SEM class");

  // --- 6. Source Information & Spectrum Generation ---
  auto cmt = SourceInfo::EarthquakeCMT(params);
  SPARSESPEC::Sparse_F_Spec mytest;

  timer1.start();
  MATRIX vec_raw = mytest.Spectra(myff, prem, cmt, params, NQ, sr_info,
                                  params.relative_error());
  timer1.stop("Total time for sparse frequency spectrum");

  // --- 7. Normalization ---
  double norm_factor = 1.0;
  double accel_norm = prem.LengthNorm() / (prem.TimeNorm() * prem.TimeNorm());

  if (params.output_type() == 0) {
    norm_factor = prem.LengthNorm();
  } else if (params.output_type() == 1) {
    norm_factor = prem.LengthNorm() / prem.TimeNorm();
  } else if (params.output_type() == 2) {
    norm_factor = accel_norm;
  }
  vec_raw *= norm_factor;

  // --- 8. Signal Filtering ---
  double hann_w = 0.5;
  auto vec_r2t_b = processfunctions::freq2time(vec_raw, myff);
  auto a_filt0 = processfunctions::fulltime2freq(vec_r2t_b, myff, 0.05);
  auto vec_filt_t = processfunctions::filtfreq2time(a_filt0, myff, false);
  auto a_filt = processfunctions::fulltime2freq(vec_filt_t, myff, hann_w);

  // --- 9. Read and Process YSpec Data ---
  std::string yspec_path =
      std::string(PROJECT_BUILD_DIR) + "../../YSpec/output/yspec.lf.out.1";
  YSPECREADER::DataColumns yspec_data(yspec_path);

  // Safe boundary checking for columns
  std::size_t maxcoly = std::min(static_cast<std::size_t>(vec_r2t_b.cols()),
                                 yspec_data.getColumn1().size());

  Eigen::MatrixXd yspec_t = Eigen::MatrixXd::Zero(3, vec_r2t_b.cols());
  for (std::size_t idx = 0; idx < maxcoly; ++idx) {
    // Replaced += with = since the matrix is already zero-initialized
    yspec_t(0, idx) = yspec_data.getColumn2()[idx];
    yspec_t(1, idx) = yspec_data.getColumn3()[idx];
    yspec_t(2, idx) = yspec_data.getColumn4()[idx];
  }

  auto a_filt_yspec0 = processfunctions::fulltime2freq(yspec_t, myff, 0.05);
  auto vec_filt_t_yspec =
      processfunctions::filtfreq2time(a_filt_yspec0, myff, false);
  auto a_filt_yspec =
      processfunctions::fulltime2freq(vec_filt_t_yspec, myff, hann_w);

  // --- 10. Output Frequency Spectrum ---
  std::string pathtofile =
      std::string(PROJECT_BUILD_DIR) + "../plotting/outputs/ex1_w.out";
  std::ofstream file(pathtofile);
  if (!file) {
    std::cerr << "Error: unable to open output file: " << pathtofile << "\n";
    return 1;
  }

  double nval = 1.0 / prem.TimeNorm();
  file << std::fixed << std::setprecision(16);

  for (std::size_t idx = 0; idx < myff.i2() + 100; ++idx) {
    file << (vec_w[idx] * nval * 1000.0 / TWO_PI) << ';'
         << a_filt(0, idx).real() << ';' << a_filt(0, idx).imag() << ';'
         << std::abs(a_filt(0, idx)) << ';' << a_filt(1, idx).real() << ';'
         << a_filt(1, idx).imag() << ';' << std::abs(a_filt(1, idx)) << ';'
         << a_filt(2, idx).real() << ';' << a_filt(2, idx).imag() << ';'
         << std::abs(a_filt(2, idx)) << ';' << a_filt_yspec(0, idx).real()
         << ';' << a_filt_yspec(0, idx).imag() << ';'
         << std::abs(a_filt_yspec(0, idx)) << ';' << a_filt_yspec(1, idx).real()
         << ';' << a_filt_yspec(1, idx).imag() << ';'
         << std::abs(a_filt_yspec(1, idx)) << ';' << a_filt_yspec(2, idx).real()
         << ';' << a_filt_yspec(2, idx).imag() << ';'
         << std::abs(a_filt_yspec(2, idx)) << '\n';
  }
  file.close();

  // --- 11. Output MinEOS Time Series ---
  std::string pathtofile_mineos =
      std::string(PROJECT_BUILD_DIR) + "../plotting/outputs/ex1_t.out";
  std::ofstream file_mineos(pathtofile_mineos);
  if (!file_mineos) {
    std::cerr << "Error: unable to open output file: " << pathtofile_mineos
              << "\n";
    return 1;
  }

  file_mineos << std::fixed << std::setprecision(16);

  for (std::size_t idx = 0;
       idx < static_cast<std::size_t>(vec_filt_t_yspec.cols()); ++idx) {
    double current_time = idx * myff.dt() * prem.TimeNorm();

    file_mineos << current_time << ';' << vec_filt_t(0, idx) << ';'
                << vec_filt_t(1, idx) << ';' << vec_filt_t(2, idx) << ';'
                << vec_filt_t_yspec(0, idx) << ';' << vec_filt_t_yspec(1, idx)
                << ';' << vec_filt_t_yspec(2, idx) << '\n';

    if (current_time > params.t_out() * 60.0) {
      break;
    }
  }
  file_mineos.close();

  return 0;
}