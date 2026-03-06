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

// Global constants
// constexpr double PI = 3.1415926535897932;
// constexpr double TWO_PI = 2.0 * PI;
// constexpr std::complex<double> I(0.0, 1.0);   // Imaginary unit

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
  using Complex = std::complex<double>;
  using MATRIX = Eigen::MatrixXcd;

  Timer timer1;

  // --- 1. Read Inputs & Earth Model ---
  // get paths required for input parameters and Earth model
  std::string param_path =
      std::string(PROJECT_BUILD_DIR) + "data/params/ex2.txt";
  InputParameters params(param_path);
  std::string earth_model_path =
      std::string(PROJECT_BUILD_DIR) + "data/" + params.earth_model();

  SRInfo sr_info(params);
  auto cmt = SourceInfo::EarthquakeCMT(params);

  prem_norm<double> norm_class;
  auto prem = EarthModels::ModelInput(earth_model_path, norm_class, "true");

  // --- 2. Frequency Solver Parameters ---
  double dt = params.time_step_sec();
  double tout = params.t_out() / 60.0;
  double df0 = 1.0;
  double wtb = 0.05;
  double t1 = 0.0;
  double t2 = tout;
  int qex = 1;
  int NQ = 6;

  SpectraSolver::FreqFull myff(params.f1(), params.f2(), params.f11(),
                               params.f12(), params.f21(), params.f22(), dt,
                               tout, df0, wtb, t1, t2, qex, prem.TimeNorm());
  auto vec_w = myff.w();

  // --- 3. Compute Sparse Frequency Spectrum ---
  SPARSESPEC::Sparse_F_Spec spec;

  timer1.start();
  MATRIX vec_raw = spec.Spectra(myff, prem, cmt, params, NQ, sr_info,
                                params.relative_error());
  timer1.stop("Total time for sparse frequency spectrum");

  // --- 4. Normalization ---
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

  // --- 5. Base Responses (Pre-Filter) ---
  double hann_w = 0.2;
  auto vec_r2t_b = processfunctions::freq2time(vec_raw, myff);
  auto a_filt0 = processfunctions::fulltime2freq(vec_r2t_b, myff, 0.01);

  // --- 6. Source Time Function (STF) Convolution ---

  double hd = 1e-9;       // Half duration (set to near zero for testing)
  double dec_t = 1.628;   // Decay constant
  double kap_val = dec_t / (std::sqrt(PI) * hd);
  double st_time = 2.0 / kap_val;

  std::cout << "Source time function parameters:\n"
            << "Half duration (s): " << hd << "\n"
            << "Kappa: " << kap_val << "\n"
            << "df: " << myff.df() << "\n"
            << "Source time (s): " << st_time * prem.TimeNorm() << "\n";

  auto a_filt_stf0 = a_filt0;
  for (int idx = 0; idx < a_filt0.cols(); ++idx) {
    auto wval = myff.w(idx);
    Complex stf_factor =
        std::exp(-I * wval * st_time) *
        std::exp(-(1.0 / (4.0 * PI)) * std::pow(wval / kap_val, 2.0));
    a_filt_stf0.col(idx) *= stf_factor;
  }

  auto vec_filt_t = processfunctions::filtfreq2time(a_filt_stf0, myff, false);
  auto a_filt = processfunctions::fulltime2freq(vec_filt_t, myff, hann_w);

  //////////////////////////////////////////////////////////////////////////////

  // --- 7. Read and Process YSpec Data ---
  std::string yspec_path = std::string(PROJECT_BUILD_DIR) +
                           "../../YSpec/output/yspec.out.bolivia.mf.1";
  YSPECREADER::DataColumns yspec_data(yspec_path);

  // FIX: Added safe boundary checking here to prevent out-of-bounds crashes
  std::size_t maxcoly = std::min(static_cast<std::size_t>(vec_r2t_b.cols()),
                                 yspec_data.getColumn1().size());

  Eigen::MatrixXd yspec_t = Eigen::MatrixXd::Zero(3, vec_r2t_b.cols());
  for (std::size_t idx = 0; idx < maxcoly; ++idx) {
    yspec_t(0, idx) = yspec_data.getColumn2()[idx];
    yspec_t(1, idx) = yspec_data.getColumn3()[idx];
    yspec_t(2, idx) = yspec_data.getColumn4()[idx];
  }

  auto a_filt_yspec0 = processfunctions::fulltime2freq(yspec_t, myff, 0.01);
  auto a_yspec_stf0 = a_filt_yspec0;

  for (int idx = 1; idx < a_yspec_stf0.cols(); ++idx) {
    auto wval = myff.w(idx);
    Complex stf_factor =
        std::exp(-I * wval * st_time) *
        std::exp(-(1.0 / (4.0 * PI)) * std::pow(wval / kap_val, 2.0));
    a_yspec_stf0.col(idx) *= stf_factor;
  }

  auto vec_filt_t_yspec =
      processfunctions::filtfreq2time(a_yspec_stf0, myff, false);
  auto a_filt_yspec =
      processfunctions::fulltime2freq(vec_filt_t_yspec, myff, hann_w);

  // --- 8. Read and Process MinEOS Data ---
  std::string mineos_base =
      std::string(PROJECT_BUILD_DIR) +
      "../../mineos/DEMO/MYEX/Syndat_ASC_BOLIVIA/Syndat2.2000160: 0:33:16.TLY.";
  MINEOSREADER::DataColumns mineos_data(mineos_base + "LHZ.ASC");
  MINEOSREADER::DataColumns mineos_data1(mineos_base + "LHN.ASC");
  MINEOSREADER::DataColumns mineos_data2(mineos_base + "LHE.ASC");

  std::size_t maxcol = std::min(static_cast<std::size_t>(vec_r2t_b.cols()),
                                mineos_data.getColumn1().size());

  Eigen::MatrixXd mineos_t = Eigen::MatrixXd::Zero(3, vec_r2t_b.cols());
  for (std::size_t idx = 0; idx < maxcol; ++idx) {
    mineos_t(0, idx) = mineos_data.getColumn2()[idx] * 1e-9;
    mineos_t(1, idx) = mineos_data1.getColumn2()[idx] * 1e-9;
    mineos_t(2, idx) = mineos_data2.getColumn2()[idx] * 1e-9;
  }

  auto a_mineos_0 = processfunctions::fulltime2freq(mineos_t, myff, 0.01);
  auto a_mineos_stf0 = a_mineos_0;

  // Start at idx=1 to avoid division by zero (wval * wval)
  for (int idx = 1; idx < a_mineos_0.cols(); ++idx) {
    auto wval = myff.w(idx);
    Complex stf_factor =
        -1.0 / (wval * wval) * std::exp(-I * wval * st_time) *
        std::exp(-(1.0 / (4.0 * PI)) * std::pow(wval / kap_val, 2.0));
    a_mineos_stf0.col(idx) *= stf_factor;
  }

  a_mineos_stf0 *= prem.TimeNorm() * prem.TimeNorm();
  auto vec_filt_t_mineos =
      processfunctions::filtfreq2time(a_mineos_stf0, myff, false);
  auto a_filt_mineos =
      processfunctions::fulltime2freq(vec_filt_t_mineos, myff, hann_w);

  //////////////////////////////////////////////////////////////////////////////
  // --- 9. Output Frequency Spectrum ---
  std::string ptf_w =
      std::string(PROJECT_BUILD_DIR) + "../plotting/outputs/ex2_w.out";
  std::ofstream file_w(ptf_w);
  if (!file_w) {
    std::cerr << "Error: unable to open output file_w: " << ptf_w << "\n";
    return 1;
  }

  double nval = 1.0 / prem.TimeNorm();
  file_w << std::fixed << std::setprecision(16);

  for (std::size_t idx = 0; idx < myff.i2() + 100; ++idx) {
    file_w << (vec_w[idx] * nval * 1000.0 / TWO_PI) << ';'
           << a_filt(0, idx).real() << ';' << a_filt(0, idx).imag() << ';'
           << std::abs(a_filt(0, idx)) << ';' << a_filt(1, idx).real() << ';'
           << a_filt(1, idx).imag() << ';' << std::abs(a_filt(1, idx)) << ';'
           << a_filt(2, idx).real() << ';' << a_filt(2, idx).imag() << ';'
           << std::abs(a_filt(2, idx)) << ';' << a_filt_yspec(0, idx).real()
           << ';' << a_filt_yspec(0, idx).imag() << ';'
           << std::abs(a_filt_yspec(0, idx)) << ';'
           << a_filt_yspec(1, idx).real() << ';' << a_filt_yspec(1, idx).imag()
           << ';' << std::abs(a_filt_yspec(1, idx)) << ';'
           << a_filt_yspec(2, idx).real() << ';' << a_filt_yspec(2, idx).imag()
           << ';' << std::abs(a_filt_yspec(2, idx)) << ';'
           << a_filt_mineos(0, idx).real() << ';'
           << a_filt_mineos(0, idx).imag() << ';'
           << std::abs(a_filt_mineos(0, idx)) << ';'
           << a_filt_mineos(1, idx).real() << ';'
           << a_filt_mineos(1, idx).imag() << ';'
           << std::abs(a_filt_mineos(1, idx)) << ';'
           << a_filt_mineos(2, idx).real() << ';'
           << a_filt_mineos(2, idx).imag() << ';'
           << std::abs(a_filt_mineos(2, idx)) << '\n';
  }
  file_w.close();

  // --- 10. Output Time Series ---
  std::string ptf_t =
      std::string(PROJECT_BUILD_DIR) + "../plotting/outputs/ex2_t.out";
  std::ofstream file_t(ptf_t);
  if (!file_t) {
    std::cerr << "Error: unable to open output file_t: " << ptf_t << "\n";
    return 1;
  }

  file_t << std::fixed << std::setprecision(16);

  for (std::size_t idx = 0;
       idx < static_cast<std::size_t>(vec_filt_t_mineos.cols()); ++idx) {
    auto tval = idx * myff.dt() * prem.TimeNorm();
    file_t << (idx * myff.dt() - st_time) * prem.TimeNorm() << ';'
           << vec_filt_t(0, idx) << ';' << vec_filt_t(1, idx) << ';'
           << vec_filt_t(2, idx) << ';' << vec_filt_t_yspec(0, idx) << ';'
           << vec_filt_t_yspec(1, idx) << ';' << vec_filt_t_yspec(2, idx) << ';'
           << vec_filt_t_mineos(0, idx) << ';' << vec_filt_t_mineos(1, idx)
           << ';' << vec_filt_t_mineos(2, idx) << '\n';

    if (tval > t2 * 3600.0) {
      break;
    }
  }
  file_t.close();

  return 0;
}