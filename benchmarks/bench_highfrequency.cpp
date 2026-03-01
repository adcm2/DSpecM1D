#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE

#include <iostream>
#include <PlanetaryModel/All>
#include <new_coupling/Timer>

#include "sem_full.h"
#include "../SpectraSolver/SpectraSolver/FF"
#include "../SpectraSolver/SpectraSolver/src/ODE_Spectra/postprocessfunctions.h"
#include "read_station.h"
#include "input_parser.h"   // Use the new input parser
#include "read_yspec.h"
#include "full_spec.h"
#include "spectra_master.h"

template <typename FLOAT> class prem_norm {
public:
  FLOAT LengthNorm() const { return _length_norm; };
  FLOAT MassNorm() const { return _mass_norm; };
  FLOAT TimeNorm() const { return _time_norm; };

  prem_norm() {};

private:
  FLOAT _length_norm = 1000.0;
  FLOAT _mass_norm = 5515 * std::pow(_length_norm, 3.0);
  FLOAT _time_norm = 1 / std::sqrt(3.1415926535897932 * 6.67230e-11 * 5515);
};

int
main() {
  using Complex = std::complex<double>;
  using MATRIX = Eigen::MatrixXcd;
  using SMATRIX = Eigen::SparseMatrix<Complex>;
  using SLU = Eigen::SparseLU<SMATRIX, Eigen::COLAMDOrdering<int>>;
  Timer timer1;

  //////////////////////////////////////////////////////////////////////////////
  // Read all parameters from the input file
  InputParameters params("../YSpec/input_bench_hf.txt");

  // Earth model
  std::string earth_model_path = "../YSpec/" + params.earth_model();

  //////////////////////////////////////////////////////////////////////////////
  // Parameters of SEM
  int lval = params.lmax(), NQ = 5;
  double maxstep = 0.05;
  std::cout << "Enter max step size for SEM integration: \n";
  std::cin >> maxstep;
  const double twopi = 2.0 * 3.1415926535897932;

  // Frequency solver parameters
  double dt = params.time_step_sec(), tout = params.t_out() / 60.0, df0 = 1.0,
         wtb = 0.05, t1 = 0, t2 = tout;
  int qex = 1;

  //////////////////////////////////////////////////////////////////////////////
  // Initialize Earth model and frequency grid
  timer1.start();
  prem_norm<double> norm_class;
  auto prem = EarthModels::ModelInput(earth_model_path, norm_class, "true");

  SpectraSolver::FreqFull myff(params.f1(), params.f2(), params.f11(),
                               params.f12(), params.f21(), params.f22(), dt,
                               tout, df0, wtb, t1, t2, qex, prem.TimeNorm());
  auto vec_w = myff.w();
  std::cout << "Lowest frequencies: " << vec_w[0] << " " << vec_w[1] << " Hz\n";
  timer1.stop("Total time for reading PREM and setting up frequency class");

  std::cout << "Velocity norm: " << prem.VelocityNorm() << "\n";
  //////////////////////////////////////////////////////////////////////////////
  // Initialize SEM and compute frequency spectrum
  timer1.start();
  Full1D::specsem sem(prem, maxstep, NQ, lval);
  timer1.stop("Total time for setting up SEM class");

  auto cmt = SourceInfo::EarthquakeCMT(params);

  int nskip = 5 * std::floor(maxstep / ((vec_w[1] - vec_w[0]) * 0.003));
  timer1.start();
  SPARSESPEC::Sparse_F_Spec mytest;
  MATRIX vec_raw = mytest.FrequencySpectrum_TEST_SPECSEM(myff, sem, prem, cmt,
                                                         params, nskip);
  timer1.stop("Total time for sparse frequency spectrum");

  //////////////////////////////////////////////////////////////////////////////
  // Post-processing: Normalization and filtering
  auto accel_norm = prem.LengthNorm() / (prem.TimeNorm() * prem.TimeNorm());
  vec_raw *= accel_norm;

  double hann_w = 0.2;

  auto vec_r2t_b = processfunctions::freq2time(vec_raw, myff);
  auto a_filt0 = processfunctions::fulltime2freq(vec_r2t_b, myff, 0.05);
  auto vec_filt_t = processfunctions::filtfreq2time(a_filt0, myff, false);
  auto a_filt = processfunctions::fulltime2freq(vec_filt_t, myff, hann_w);

  //////////////////////////////////////////////////////////////////////////////
  // Output frequency-domain results
  std::string pathtofile = "./work/benchmarks/spectra_hf.out";
  std::ofstream file(pathtofile);
  if (!file) {
    std::cerr << "Error: unable to open output file: " << pathtofile << "\n";
    return 1;
  }

  double nval = 1.0 / prem.TimeNorm();
  const auto &vec_w_ref = vec_w;

  file.setf(std::ios::fixed);
  file << std::setprecision(22);

  for (std::size_t idx = 0; idx < myff.i2() + 100; ++idx) {
    file << (vec_w_ref[idx] * nval * 1000.0 / (twopi)) << ';'
         << a_filt(0, idx).real() << ';' << a_filt(0, idx).imag() << ';'
         << std::abs(a_filt(0, idx)) << ';' << a_filt(1, idx).real() << ';'
         << a_filt(1, idx).imag() << ';' << std::abs(a_filt(1, idx)) << ';'
         << a_filt(2, idx).real() << ';' << a_filt(2, idx).imag() << ';'
         << std::abs(a_filt(2, idx)) << '\n';
  }
  file.close();

  //////////////////////////////////////////////////////////////////////////////
  // Output time-domain comparison results
  std::string pathtofile_compare = "./work/benchmarks/time_hf.out";
  std::ofstream file_compare(pathtofile_compare);
  if (!file_compare) {
    std::cerr << "Error: unable to open output file: " << pathtofile_compare
              << "\n";
    return 1;
  }
  file_compare.setf(std::ios::fixed);
  file_compare << std::setprecision(22);

  for (std::size_t idx = 0; idx < vec_filt_t.cols(); ++idx) {
    auto tval = idx * myff.dt() * prem.TimeNorm();
    file_compare << tval << ';' << vec_filt_t(0, idx) << ";"
                 << vec_filt_t(1, idx) << ";" << vec_filt_t(2, idx) << '\n';
    if (tval > params.t_out() * 60.0) {
      break;
    }
  }
  file_compare.close();

  return 0;
}
#endif