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

// Project-specific includes
#include <PlanetaryModel/All>
#include <DSpecM1D/Timer>
#include <DSpecM1D/All>
#include <SpectraSolver/FF>

// constexpr double PI = 3.1415926535897932;

// Normalisation class for PREM - stores length, mass, and time scales
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

  // --- 1. Normalisation and Input Parameters ---
  prem_norm<double> norm_class;
  auto timenorm = norm_class.TimeNorm();

  // get paths required for input parameters and Earth model
  std::string param_path =
      std::string(PROJECT_BUILD_DIR) + "data/params/ex4.txt";
  InputParameters params(param_path);
  std::string earth_model_path =
      std::string(PROJECT_BUILD_DIR) + "data/" + params.earth_model();

  SRInfo srInfo(params);

  // --- 2. Frequency Solver Parameters ---
  double dt = params.time_step_sec();
  double tout = params.t_out() / 60.0;
  double df0 = 1.0;
  double wtb = 0.05;
  double t1 = 0.0;
  double t2 = tout + 1.0;
  int NQ = 5;
  int qex = 1;

  // --- 3. Setup Frequency Class ---
  SpectraSolver::FreqFull myff(params.f1(), params.f2(), params.f11(),
                               params.f12(), params.f21(), params.f22(), dt,
                               tout, df0, wtb, t1, t2, qex, timenorm);

  // --- 4. Earth Model and Source ---
  auto prem = EarthModels::ModelInput(earth_model_path, norm_class, "true");
  auto cmt = SourceInfo::EarthquakeCMT(params);

  // --- 5. Compute Raw Frequency-Domain Spectra ---
  SPARSESPEC::SparseFSpec mytest;

  timer1.start();
  MATRIX vec_raw = mytest.spectra(myff, prem, cmt, params, NQ, srInfo,
                                  params.relative_error());
  timer1.stop("Total time for sparse frequency spectrum");

  // --- 6. Post-processing (Freq -> Time -> Filter) ---
  auto vec_r2t_b = processfunctions::freq2time(vec_raw, myff);
  auto a_filt0 = processfunctions::fulltime2freq(vec_r2t_b, myff, 0.01);
  auto vec_filt_t = processfunctions::filtfreq2time(a_filt0, myff, false);

  // Dimensionalise to acceleration (m/s^2)
  double accel_norm = prem.LengthNorm() / (prem.TimeNorm() * prem.TimeNorm());
  vec_filt_t *= accel_norm;

  //////////////////////////////////////////////////////////////////////////////
  // --- 7. Write Record Section to File ---
  std::string pathtofile =
      std::string(PROJECT_BUILD_DIR) + "../plotting/outputs/ex4.out";
  std::ofstream file(pathtofile);

  // Safety check to ensure the file opened successfully
  if (!file) {
    std::cerr << "Error: unable to open output file: " << pathtofile << "\n";
    return 1;
  }

  file << std::fixed << std::setprecision(22);

  for (std::size_t idx = 0; idx < static_cast<std::size_t>(vec_filt_t.cols());
       ++idx) {
    double t_sec = idx * myff.dt() * prem.TimeNorm();
    file << t_sec;

    for (int jidx = 0; jidx < params.num_receivers(); ++jidx) {
      file << ';' << vec_filt_t(3 * jidx + 0, idx) << ';'
           << vec_filt_t(3 * jidx + 1, idx) << ';'
           << vec_filt_t(3 * jidx + 2, idx);
    }
    file << '\n';

    if (t_sec > params.t_out() * 60.0) {
      break;
    }
  }
  file.close();

  return 0;
}