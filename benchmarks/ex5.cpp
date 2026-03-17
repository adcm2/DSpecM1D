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
#include <vector>

// Project-specific includes
#include <PlanetaryModel/All>
#include <DSpecM1D/Timer>
#include <DSpecM1D/All>
#include <SpectraSolver/FF>

// constexpr double PI = 3.1415926535897932;
// constexpr double TWO_PI = 2.0 * PI;

// template <typename FLOAT> class prem_norm {
// public:
//   prem_norm() = default;

//   FLOAT LengthNorm() const { return _length_norm; }
//   FLOAT MassNorm() const { return _mass_norm; }
//   FLOAT TimeNorm() const { return _time_norm; }

// private:
//   FLOAT _length_norm = 1000.0;
//   FLOAT _mass_norm = 5515.0 * std::pow(_length_norm, 3.0);
//   FLOAT _time_norm = 1.0 / std::sqrt(PI * 6.67230e-11 * 5515.0);
// };

int
main() {
  using MATRIX = Eigen::MatrixXcd;
  Timer timer1;

  // --- 1. Read Inputs & Earth Model ---
  // get paths required for input parameters and Earth model
  std::string param_path =
      std::string(PROJECT_BUILD_DIR) + "data/params/ex5.txt";
  InputParameters params(param_path);
  std::string earth_model_path =
      std::string(PROJECT_BUILD_DIR) + "data/" + params.earth_model();

  prem_norm<double> norm_class;
  auto prem = EarthModels::ModelInput(earth_model_path, norm_class, "true");
  auto cmt = SourceInfo::EarthquakeCMT(params);

  // --- 2. Parameters ---
  int lval = params.lmax();
  int NQ = 6;
  int qex = 1;
  double maxstep = 0.05;

  double dt = params.time_step_sec();
  double tout = params.t_out() / 60.0;
  double df0 = 1.0;
  double wtb = 0.05;
  double t1 = 0.0;
  double t2 = tout;

  // --- 3. Setup Frequency Class ---
  SpectraSolver::FreqFull myff(params.f1(), params.f2(), params.f11(),
                               params.f12(), params.f21(), params.f22(), dt,
                               tout, df0, wtb, t1, t2, qex, prem.TimeNorm());
  auto vec_w = myff.w();

  // --- 4. Setup Convergence Steps ---
  int nsteps = 50;
  double step_0 = 2.0 * 0.63 / myff.f22();
  std::cout << "Initial step: " << step_0 << "\n";

  std::vector<double> vec_step;
  vec_step.reserve(nsteps);
  for (int idx = 0; idx < nsteps; ++idx) {
    vec_step.push_back(
        step_0 / std::pow(10.0, 2.0 * idx / static_cast<double>(nsteps - 1)));
  }

  // --- 5. Normalization ---
  double norm_factor = 1.0;
  double accel_norm = prem.LengthNorm() / (prem.TimeNorm() * prem.TimeNorm());

  if (params.output_type() == 0) {
    norm_factor = prem.LengthNorm();
  } else if (params.output_type() == 1) {
    norm_factor = prem.LengthNorm() / prem.TimeNorm();
  } else if (params.output_type() == 2) {
    norm_factor = accel_norm;
  }
  double hann_w = 0.2;

  // --- 6. Execute Iterative Convergence Test ---
  SPARSESPEC::SparseFSpec mytest;
  std::vector<Eigen::MatrixXcd> vec_final_w;
  std::vector<Eigen::MatrixXd> vec_final_t;
  vec_final_w.reserve(nsteps);
  vec_final_t.reserve(nsteps);

  for (int idx = 0; idx < nsteps; ++idx) {
    maxstep = vec_step[idx];
    int nskip = 3 * static_cast<int>(
                        std::floor(maxstep / ((vec_w[1] - vec_w[0]) * 0.003))) +
                1;

    Full1D::SEM sem(prem, maxstep, NQ, lval);
    std::cout << "\nDoing step: " << maxstep << ", nskip: " << nskip << "\n";

    MATRIX vec_raw = mytest.spectra(myff, sem, prem, cmt, params, nskip);
    vec_raw *= norm_factor;

    // Process responses
    auto vec_r2t_b = processfunctions::freq2time(vec_raw, myff);
    auto a_filt0 = processfunctions::fulltime2freq(vec_r2t_b, myff, 0.05);
    auto vec_filt_t = processfunctions::filtfreq2time(a_filt0, myff, false);
    auto a_filt = processfunctions::fulltime2freq(vec_filt_t, myff, hann_w);

    vec_final_w.push_back(a_filt);
    vec_final_t.push_back(vec_filt_t);
  }

  // --- 7. Error Calculation ---
  std::size_t num_err_steps = vec_step.size() - 1;
  std::vector<std::vector<double>> vec_err_t(num_err_steps,
                                             std::vector<double>(3, 0.0));
  std::vector<std::vector<double>> vec_err_w(num_err_steps,
                                             std::vector<double>(3, 0.0));
  std::vector<std::vector<double>> vec_l2_err_t(num_err_steps,
                                                std::vector<double>(3, 0.0));
  std::vector<std::vector<double>> vec_l2_err_w(num_err_steps,
                                                std::vector<double>(3, 0.0));

  int idxout = vec_final_t[0].cols() - 1;
  for (int idx = 0; idx < vec_final_t[0].cols(); ++idx) {
    auto tval = idx * myff.dt() * prem.TimeNorm();
    if (tval > t2 * 3600.0) {
      idxout = idx;
      break;
    }
  }

  for (int idx = 0; idx < nsteps; ++idx) {
    vec_final_t[idx]
        .block(0, idxout, 3, vec_final_t[idx].cols() - idxout)
        .setZero();
  }

  // "Exact" solution is the one with the smallest step
  auto idxback = nsteps - 1;
  auto vec_t_ex = vec_final_t[idxback];
  auto vec_w_ex = vec_final_w[idxback];

  // Reference Norms Setup
  auto tmp_t = 100.0 / idxout;
  auto mv_t_1 = tmp_t / vec_t_ex.row(0).lpNorm<Eigen::Infinity>();
  auto mv_t_2 = tmp_t / vec_t_ex.row(1).lpNorm<Eigen::Infinity>();
  auto mv_t_3 = tmp_t / vec_t_ex.row(2).lpNorm<Eigen::Infinity>();

  auto mv_t_1_l2 = 100.0 / vec_t_ex.row(0).norm();
  auto mv_t_2_l2 = 100.0 / vec_t_ex.row(1).norm();
  auto mv_t_3_l2 = 100.0 / vec_t_ex.row(2).norm();

  auto wcols = vec_final_w[idxback].cols();
  auto tmp_w = 100.0 / wcols;
  auto mv_w_1 = tmp_w / vec_w_ex.row(0).lpNorm<Eigen::Infinity>();
  auto mv_w_2 = tmp_w / vec_w_ex.row(1).lpNorm<Eigen::Infinity>();
  auto mv_w_3 = tmp_w / vec_w_ex.row(2).lpNorm<Eigen::Infinity>();

  auto mv_w_1_l2 = 100.0 / vec_w_ex.row(0).norm();
  auto mv_w_2_l2 = 100.0 / vec_w_ex.row(1).norm();
  auto mv_w_3_l2 = 100.0 / vec_w_ex.row(2).norm();

  for (std::size_t idx = 0; idx < num_err_steps; ++idx) {
    auto vec_t_err = vec_final_t[idx] - vec_t_ex;
    auto vec_w_err = vec_final_w[idx] - vec_w_ex;

    // L1 norms (scaled by exact solution's L-infinity norm)
    vec_err_t[idx][0] = vec_t_err.row(0).lpNorm<1>() * mv_t_1;
    vec_err_t[idx][1] = vec_t_err.row(1).lpNorm<1>() * mv_t_2;
    vec_err_t[idx][2] = vec_t_err.row(2).lpNorm<1>() * mv_t_3;
    vec_err_w[idx][0] = vec_w_err.row(0).lpNorm<1>() * mv_w_1;
    vec_err_w[idx][1] = vec_w_err.row(1).lpNorm<1>() * mv_w_2;
    vec_err_w[idx][2] = vec_w_err.row(2).lpNorm<1>() * mv_w_3;

    // L2 norms
    vec_l2_err_t[idx][0] = vec_t_err.row(0).norm() * mv_t_1_l2;
    vec_l2_err_t[idx][1] = vec_t_err.row(1).norm() * mv_t_2_l2;
    vec_l2_err_t[idx][2] = vec_t_err.row(2).norm() * mv_t_3_l2;
    vec_l2_err_w[idx][0] = vec_w_err.row(0).norm() * mv_w_1_l2;
    vec_l2_err_w[idx][1] = vec_w_err.row(1).norm() * mv_w_2_l2;
    vec_l2_err_w[idx][2] = vec_w_err.row(2).norm() * mv_w_3_l2;
  }

  // --- 8. Outputs ---
  double nval = 1.0 / prem.TimeNorm();

  // 8a. Output Frequency Series
  std::string ptf_w = std::string(PROJECT_BUILD_DIR) +
                      "../plotting/outputs/ex5_w_NQ" + std::to_string(NQ) +
                      "_step.out";
  std::ofstream file_w(ptf_w);
  if (!file_w) {
    std::cerr << "Error: unable to open output file_w: " << ptf_w << "\n";
    return 1;
  }

  file_w << std::fixed << std::setprecision(22);
  for (std::size_t idx = 0; idx < myff.i2() + 100; ++idx) {
    file_w << (vec_w[idx] * nval * 1000.0 / TWO_PI);
    for (std::size_t idx2 = 0; idx2 < vec_step.size(); ++idx2) {
      file_w << ";" << vec_final_w[idx2](0, idx).real() << ';'
             << vec_final_w[idx2](0, idx).imag() << ';'
             << std::abs(vec_final_w[idx2](0, idx)) << ';'
             << vec_final_w[idx2](1, idx).real() << ';'
             << vec_final_w[idx2](1, idx).imag() << ';'
             << std::abs(vec_final_w[idx2](1, idx)) << ';'
             << vec_final_w[idx2](2, idx).real() << ';'
             << vec_final_w[idx2](2, idx).imag() << ';'
             << std::abs(vec_final_w[idx2](2, idx));
    }
    file_w << '\n';
  }
  file_w.close();

  // 8b. Output Time Series
  std::string ptf_t = std::string(PROJECT_BUILD_DIR) +
                      "../plotting/outputs/ex5_t_NQ" + std::to_string(NQ) +
                      "_step.out";
  std::ofstream file_t(ptf_t);
  if (!file_t) {
    std::cerr << "Error: unable to open output file_t: " << ptf_t << "\n";
    return 1;
  }

  file_t << std::fixed << std::setprecision(22);
  for (std::size_t idx = 0; idx < static_cast<std::size_t>(idxout); ++idx) {
    file_t << (idx * myff.dt()) * prem.TimeNorm();
    for (std::size_t idx2 = 0; idx2 < static_cast<std::size_t>(nsteps);
         ++idx2) {
      file_t << ";" << vec_final_t[idx2](0, idx) << ";"
             << vec_final_t[idx2](1, idx) << ";" << vec_final_t[idx2](2, idx);
    }
    file_t << '\n';
  }
  file_t.close();

  // 8c. Output Error Matrix (L1-based)
  std::string ptf_err = std::string(PROJECT_BUILD_DIR) +
                        "../plotting/outputs/ex5_NQ" + std::to_string(NQ) +
                        "_step_error_" +
                        std::to_string(static_cast<int>(params.f22())) + ".out";
  std::ofstream file_err(ptf_err);
  if (!file_err) {
    std::cerr << "Error: unable to open output file_err: " << ptf_err << "\n";
    return 1;
  }

  file_err << std::fixed << std::setprecision(22);
  for (std::size_t idx = 0; idx < num_err_steps; ++idx) {
    file_err << vec_step[idx] << ";" << vec_err_t[idx][0] << ";"
             << vec_err_t[idx][1] << ";" << vec_err_t[idx][2] << ";"
             << vec_err_w[idx][0] << ";" << vec_err_w[idx][1] << ";"
             << vec_err_w[idx][2] << "\n";
  }
  file_err.close();

  // 8d. Output L2 Error Matrix
  std::string ptf_l2 = std::string(PROJECT_BUILD_DIR) +
                       "../plotting/outputs/ex5_NQ" + std::to_string(NQ) +
                       "_step_error_l2_" +
                       std::to_string(static_cast<int>(params.f22())) + ".out";
  std::ofstream file_l2(ptf_l2);
  if (!file_l2) {
    std::cerr << "Error: unable to open output file_l2: " << ptf_l2 << "\n";
    return 1;
  }

  file_l2 << std::fixed << std::setprecision(22);
  for (std::size_t idx = 0; idx < num_err_steps; ++idx) {
    file_l2 << vec_step[idx] << ";" << vec_l2_err_t[idx][0] << ";"
            << vec_l2_err_t[idx][1] << ";" << vec_l2_err_t[idx][2] << ";"
            << vec_l2_err_w[idx][0] << ";" << vec_l2_err_w[idx][1] << ";"
            << vec_l2_err_w[idx][2] << "\n";
  }
  file_l2.close();

  return 0;
}