#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#include <iostream>
#include <PlanetaryModel/All>
#include <new_coupling/Timer>
#include "sem_full.h"
// #include "sem_spheroidal_debug.h"
#include "../SpectraSolver/SpectraSolver/FF"
// #include "../SpectraSolver/SpectraSolver/src/ODE_Spectra/filter_base.h"
#include "../SpectraSolver/SpectraSolver/src/ODE_Spectra/postprocessfunctions.h"
#include "read_station.h"
#include "input_parser.h"   // Use the new input parser
#include "read_yspec.h"
#include "read_mineos.h"
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

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Read all parameters from the input file_w
  InputParameters params("../YSpec/input_bench_tor.txt");

  // earth model
  std::string cpath = params.earth_model();
  std::string earth_model_path = "../YSpec/" + params.earth_model();

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // parameters of sem
  int lval = params.lmax(), NQ = 5;
  bool toaug = false;
  double maxstep = 0.05;
  std::cout << "Max step: \n";
  std::cin >> maxstep;
  const double twopi = 2.0 * 3.1415926535897932;
  std::string pathpert;

  //////////////////////////////
  // frequency solver parameters
  double dt = params.time_step_sec(), tout = params.t_out() / 60.0, df0 = 1.0,
         wtb = 0.05, t1 = 0, t2 = tout;
  int qex = 1;
  Complex myi(0.0, 1.0);
  double droptol = 1e-4;
  auto tend = 0.6 * params.t_out();   // length of time series in seconds
  // std::cout << "tend: " << tend << "\n";
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // prem
  // auto prem = EarthModels::ModelInput(pathtoprem);
  timer1.start();
  prem_norm<double> norm_class;
  auto prem = EarthModels::ModelInput(earth_model_path, norm_class, "true");

  // frequency class
  SpectraSolver::FreqFull myff(params.f1(), params.f2(), params.f11(),
                               params.f12(), params.f21(), params.f22(), dt,
                               tout, df0, wtb, t1, t2, qex, prem.TimeNorm());
  auto vec_w = myff.w();

  timer1.stop("Total time for reading PREM and setting up frequency class");

  timer1.start();
  // sem class
  // Full1D::sem sem(prem, maxstep, NQ, lval);
  Full1D::specsem sem(prem, maxstep, NQ, lval);
  timer1.stop("Total time for setting up SEM class");

  // source information
  auto cmt = SourceInfo::EarthquakeCMT(params);

  //////////////////////////////////////////////////////////////////////////////

  // test new functionality
  SPARSESPEC::Sparse_F_Spec mytest;
  // int nskip = 20;
  // std::cout << "Approximate nskip: "
  //           << maxstep / ((vec_w[1] - vec_w[0]) * 0.003) << "\n";
  int nskip = 5 * std::floor(maxstep / ((vec_w[1] - vec_w[0]) * 0.003));
  omp_set_num_threads(10);
  std::cout << "Using " << omp_get_num_threads()
            << " threads for sparse frequency spectrum calculation\n";
  // set_omp_num_threads(10);
  timer1.start();
  MATRIX vec_raw = mytest.FrequencySpectrum_TEST_SPECSEM(myff, sem, prem, cmt,
                                                         params, nskip);
  timer1.stop("Total time for sparse frequency spectrum");

  // normalise
  double norm_factor;
  auto accel_norm = prem.LengthNorm() / (prem.TimeNorm() * prem.TimeNorm());
  if (params.output_type() == 0) {
    norm_factor = prem.LengthNorm();
  } else if (params.output_type() == 1) {
    norm_factor = prem.LengthNorm() / prem.TimeNorm();
  } else if (params.output_type() == 2) {
    norm_factor = accel_norm;
  }
  vec_raw *= norm_factor;

  double hann_w = 0.2;
  // process responses
  // don't do hann filter until you have time domain response etc
  auto vec_r2t_b = processfunctions::freq2time(vec_raw, myff);
  auto a_filt0 = processfunctions::fulltime2freq(vec_r2t_b, myff, 0.05);
  auto a_filt_stf0 = a_filt0;
  auto vec_filt_t = processfunctions::filtfreq2time(a_filt_stf0, myff, false);
  auto a_filt = processfunctions::fulltime2freq(vec_filt_t, myff, hann_w);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // testing read in of yspec
  // std::cout << "First output\n";
  std::string yspec_path = "../YSpec/output/yspec.out.tor.1";
  YSPECREADER::DataColumns yspec_data(yspec_path);

  Eigen::MatrixXd yspec_t = Eigen::MatrixXd::Zero(3, vec_r2t_b.cols());
  for (int idx = 0; idx < yspec_data.getColumn1().size(); ++idx) {
    yspec_t(0, idx) += yspec_data.getColumn2()[idx];
    yspec_t(1, idx) += yspec_data.getColumn3()[idx];
    yspec_t(2, idx) += yspec_data.getColumn4()[idx];
  }
  auto a_filt_yspec0 = processfunctions::fulltime2freq(yspec_t, myff, 0.05);
  auto a_yspec_stf0 = a_filt_yspec0;
  auto vec_filt_t_yspec =
      processfunctions::filtfreq2time(a_yspec_stf0, myff, false);
  auto a_filt_yspec =
      processfunctions::fulltime2freq(vec_filt_t_yspec, myff, hann_w);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // output
  // outputting result
  // std::cout << "Second output\n";
  std::string ptf_w = "./work/benchmarks/bench_w_tor.out";
  std::ofstream file_w(ptf_w);
  if (!file_w) {
    std::cerr << "Error: unable to open output file_w: " << ptf_w << "\n";
    return 1;
  }

  double nval = 1.0 / prem.TimeNorm();
  const auto &vec_w_ref = vec_w;             // alias
  const std::size_t nw = vec_w_ref.size();   // use actual vector size

  // Write header (optional)
  // file_w <<
  // "#freq_mHz;Z_re;Z_im;Z_abs;TH_re;TH_im;TH_abs;PH_re;PH_im;PH_abs\n";
  file_w.setf(std::ios::fixed);
  file_w << std::setprecision(22);

  for (std::size_t idx = 0; idx < myff.i2() + 100; ++idx) {
    // index bounds are safe since nw == vec_w.size() and a_filt has
    // compatible columns
    file_w << (vec_w_ref[idx] * nval * 1000.0 / (twopi)) << ';'
           << a_filt(0, idx).real() << ';' << a_filt(0, idx).imag() << ';'
           << std::abs(a_filt(0, idx)) << ';' << a_filt(1, idx).real() << ';'
           << a_filt(1, idx).imag() << ';' << std::abs(a_filt(1, idx)) << ';'
           << a_filt(2, idx).real() << ';' << a_filt(2, idx).imag() << ';'
           << std::abs(a_filt(2, idx)) << ";" << a_filt_yspec(0, idx).real()
           << ";" << a_filt_yspec(0, idx).imag() << ";"
           << std::abs(a_filt_yspec(0, idx)) << ";"
           << a_filt_yspec(1, idx).real() << ";" << a_filt_yspec(1, idx).imag()
           << ";" << std::abs(a_filt_yspec(1, idx)) << ";"
           << a_filt_yspec(2, idx).real() << ";" << a_filt_yspec(2, idx).imag()
           << ";" << std::abs(a_filt_yspec(2, idx)) << ';'
           << vec_raw(0, idx).real() << ';' << vec_raw(0, idx).imag() << ';'
           << std::abs(vec_raw(0, idx)) << ';' << vec_raw(1, idx).real() << ';'
           << vec_raw(1, idx).imag() << ';' << std::abs(vec_raw(1, idx)) << ';'
           << vec_raw(2, idx).real() << ';' << vec_raw(2, idx).imag() << ';'
           << std::abs(vec_raw(2, idx)) << '\n';
  }
  file_w.close();

  //////////////////////////////////////////////////////////////////////////////
  // output mineos time series
  std::string ptf_t = "./work/benchmarks/bench_t_tor.out";
  std::ofstream file_t(ptf_t);
  if (!file_t) {
    std::cerr << "Error: unable to open output file_t: " << ptf_t << "\n";
    return 1;
  }
  file_t.setf(std::ios::fixed);
  file_t << std::setprecision(22);
  for (std::size_t idx = 0; idx < vec_filt_t_yspec.cols(); ++idx) {
    // index bounds are safe since nw == vec_w.size() and a_filt has
    // compatible columns
    auto tval = idx * myff.dt() * prem.TimeNorm();
    file_t << (idx * myff.dt()) * prem.TimeNorm() << ";" << vec_filt_t(0, idx)
           << ";" << vec_filt_t(1, idx) << ";" << vec_filt_t(2, idx) << ";"
           << vec_filt_t_yspec(0, idx) << ";" << vec_filt_t_yspec(1, idx) << ";"
           << vec_filt_t_yspec(2, idx) << '\n';
    if (tval > t2 * 3600.0) {
      break;
    }
  }
  file_t.close();

  return 0;
}
#endif
