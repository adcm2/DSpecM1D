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
#include "SR_Info.h"

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
  // Read all parameters from the input file
  InputParameters params("../MyVersion/YSpec/input_bench_lf.txt");
  SRInfo sr_info(params);
  // earth model
  std::string cpath = params.earth_model();
  std::string earth_model_path = "../MyVersion/YSpec/" + params.earth_model();
  std::cout << "Path to earth model: " << earth_model_path << "\n";
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // parameters of sem
  int lval = params.lmax(), NQ = 5;
  bool toaug = false;
  double maxstep = 0.05;
  // std::cout << "Max step: \n";
  // std::cin >> maxstep;
  const double twopi = 2.0 * 3.1415926535897932;
  std::string pathpert;

  //////////////////////////////
  // frequency solver parameters
  double dt = params.time_step_sec(), tout = params.t_out() / 60.0, df0 = 1.0,
         wtb = 0.05, t1 = 0, t2 = tout;
  int qex = 4;
  Complex myi(0.0, 1.0);
  double droptol = 1e-4;

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
  std::cout << "Setting up SEM class...\n";

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
  nskip += 1;
  int num_chunks = 5;
  timer1.start();
  MATRIX vec_raw = mytest.FrequencySpectrum_RED(myff, prem, cmt, params, NQ,
                                                nskip, num_chunks, sr_info);

  // MATRIX vec_raw = mytest.FrequencySpectrum_TEST_SPECSEM(myff, sem, prem,
  // cmt,
  //                                                        params, nskip);
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

  double hann_w = 0.5;
  // std::cout << "Finished computing raw frequency spectrum.\n";
  // process responses
  // don't do hann filter until you have time domain response etc
  auto vec_r2t_b = processfunctions::freq2time(vec_raw, myff);
  // std::cout << "Finished converting raw frequency spectrum to time
  // domain.\n";
  auto a_filt0 = processfunctions::fulltime2freq(vec_r2t_b, myff, 0.05);
  // std::cout
  // << "Finished converting raw time domain response to frequency
  // domain.\n";
  auto vec_filt_t = processfunctions::filtfreq2time(a_filt0, myff, false);
  // std::cout << "Finished converting filtered frequency response back to time
  // "
  // "domain.\n";
  auto a_filt = processfunctions::fulltime2freq(vec_filt_t, myff, hann_w);

  // std::cout << "Finished processing response.\n";
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // testing read in of yspec
  // std::cout << "First output\n";
  std::string yspec_path = "../MyVersion/YSpec/output/yspec.lf.out.1";
  YSPECREADER::DataColumns yspec_data(yspec_path);
  std::size_t maxcoly;
  if (yspec_data.getColumn1().size() > vec_r2t_b.cols()) {
    maxcoly = vec_r2t_b.cols();
  } else {
    maxcoly = yspec_data.getColumn1().size();
  }
  Eigen::MatrixXd yspec_t = Eigen::MatrixXd::Zero(3, vec_r2t_b.cols());
  for (int idx = 0; idx < maxcoly; ++idx) {
    yspec_t(0, idx) += yspec_data.getColumn2()[idx];
    yspec_t(1, idx) += yspec_data.getColumn3()[idx];
    yspec_t(2, idx) += yspec_data.getColumn4()[idx];
  }
  // std::cout << "Finished reading YSpec output.\n";
  // std::cout << "The size of yspec_t: " << yspec_t.rows() << " x "
  //           << yspec_t.cols() << "\n";
  auto a_filt_yspec0 = processfunctions::fulltime2freq(yspec_t, myff, 0.05);
  auto vec_filt_t_yspec =
      processfunctions::filtfreq2time(a_filt_yspec0, myff, false);
  auto a_filt_yspec =
      processfunctions::fulltime2freq(vec_filt_t_yspec, myff, hann_w);
  // std::cout << "Finished processing YSpec output.\n";
  // read in mineos output
  // std::string mineos_path = "../mineos/DEMO/MYEX/Syndat_ASC_NOHEADER/"
  //                           "Syndat.2000014:23:37:10.TLY.LHZ.ASC";
  // MINEOSREADER::DataColumns mineos_data(mineos_path);
  // std::string mineos_path1 = "../mineos/DEMO/MYEX/Syndat_ASC_NOHEADER/"
  //                            "Syndat.2000014:23:37:10.TLY.LHN.ASC";
  // MINEOSREADER::DataColumns mineos_data1(mineos_path1);
  // std::string mineos_path2 = "../mineos/DEMO/MYEX/Syndat_ASC_NOHEADER/"
  //                            "Syndat.2000014:23:37:10.TLY.LHE.ASC";
  // MINEOSREADER::DataColumns mineos_data2(mineos_path2);

  // read in mineos output
  // std::string mineos_path = "../mineos/DEMO/MYEX/Syndat_ASC_BOLIVIA/"
  //                           "Syndat.2000160: 0:33:16.TLY.LHZ.ASC";
  // MINEOSREADER::DataColumns mineos_data(mineos_path);
  // std::string mineos_path1 = "../mineos/DEMO/MYEX/Syndat_ASC_BOLIVIA/"
  //                            "Syndat.2000160: 0:33:16.TLY.LHN.ASC";
  // MINEOSREADER::DataColumns mineos_data1(mineos_path1);
  // std::string mineos_path2 = "../mineos/DEMO/MYEX/Syndat_ASC_BOLIVIA/"
  //                            "Syndat.2000160: 0:33:16.TLY.LHE.ASC";
  // MINEOSREADER::DataColumns mineos_data2(mineos_path2);

  // std::size_t maxcol;
  // if (mineos_data.getColumn1().size() > vec_r2t_b.cols()) {
  //   maxcol = vec_r2t_b.cols();
  // } else {
  //   maxcol = mineos_data.getColumn1().size();
  // }
  Eigen::MatrixXd mineos_t = yspec_t;
  // for (int idx = 0; idx < maxcol; ++idx) {
  //   mineos_t(0, idx) += mineos_data.getColumn2()[idx] * 1e-9;
  //   mineos_t(1, idx) += mineos_data1.getColumn2()[idx] * 1e-9;
  //   mineos_t(2, idx) += mineos_data2.getColumn2()[idx] * 1e-9;
  // }
  // std::cout << "Finished reading Mineos output.\n";
  // std::cout << "The size of mineos_t: " << mineos_t.rows() << " x "
  //           << mineos_t.cols() << "\n";
  // auto a_mineos_00 = processfunctions::fulltime2freq(mineos_t, myff, 0.001);
  // for (int idx = 0; idx < a_mineos_00.cols(); ++idx) {
  //   auto fval = myff.df() * idx;
  //   if ((fval > myff.f22()) || (fval < myff.f11())) {
  //     a_mineos_00.col(idx) *= 0.0;
  //   }
  // }
  // auto vec_filt_t_mineos0 =
  //     processfunctions::filtfreq2time(a_mineos_00, myff, false);
  auto a_mineos_0 = processfunctions::fulltime2freq(mineos_t, myff, 0.05);
  auto vec_filt_t_mineos =
      processfunctions::filtfreq2time(a_mineos_0, myff, false);
  auto a_filt_mineos =
      processfunctions::fulltime2freq(vec_filt_t_mineos, myff, hann_w);
  // auto a_filt_mineos =
  // std::cout << "Finished processing Mineos output.\n";
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // output
  // outputting result
  // std::cout << "Second output\n";
  std::string pathtofile = "./work/benchmarks/bench_w_lf.out";
  std::ofstream file(pathtofile);
  if (!file) {
    std::cerr << "Error: unable to open output file: " << pathtofile << "\n";
    return 1;
  }

  double nval = 1.0 / prem.TimeNorm();
  const auto &vec_w_ref = vec_w;             // alias
  const std::size_t nw = vec_w_ref.size();   // use actual vector size

  // Write header (optional)
  // file <<
  // "#freq_mHz;Z_re;Z_im;Z_abs;TH_re;TH_im;TH_abs;PH_re;PH_im;PH_abs\n";
  file.setf(std::ios::fixed);
  file << std::setprecision(16);

  for (std::size_t idx = 0; idx < myff.i2() + 100; ++idx) {
    // index bounds are safe since nw == vec_w.size() and a_filt has
    // compatible columns
    file << (vec_w_ref[idx] * nval * 1000.0 / (twopi)) << ';'
         << a_filt(0, idx).real() << ';' << a_filt(0, idx).imag() << ';'
         << std::abs(a_filt(0, idx)) << ';' << a_filt(1, idx).real() << ';'
         << a_filt(1, idx).imag() << ';' << std::abs(a_filt(1, idx)) << ';'
         << a_filt(2, idx).real() << ';' << a_filt(2, idx).imag() << ';'
         << std::abs(a_filt(2, idx)) << ";" << a_filt_yspec(0, idx).real()
         << ";" << a_filt_yspec(0, idx).imag() << ";"
         << std::abs(a_filt_yspec(0, idx)) << ";" << a_filt_yspec(1, idx).real()
         << ";" << a_filt_yspec(1, idx).imag() << ";"
         << std::abs(a_filt_yspec(1, idx)) << ";" << a_filt_yspec(2, idx).real()
         << ";" << a_filt_yspec(2, idx).imag() << ";"
         << std::abs(a_filt_yspec(2, idx)) << ";"
         << a_filt_mineos(0, idx).real() << ";" << a_filt_mineos(0, idx).imag()
         << ";" << std::abs(a_filt_mineos(0, idx)) << ";"
         << a_filt_mineos(1, idx).real() << ";" << a_filt_mineos(1, idx).imag()
         << ";" << std::abs(a_filt_mineos(1, idx)) << ";"
         << a_filt_mineos(2, idx).real() << ";" << a_filt_mineos(2, idx).imag()
         << ";" << std::abs(a_filt_mineos(2, idx)) << '\n';
  }
  file.close();

  //////////////////////////////////////////////////////////////////////////////
  // Repeat output for raw data
  // std::cout << "Third output\n";
  std::string pathtofile2 = "./work/benchmarks/bench_raw_w_lf.out";
  std::ofstream file2(pathtofile2);
  if (!file2) {
    std::cerr << "Error: unable to open output file: " << pathtofile << "\n";
    return 1;
  }

  // Write header (optional)
  // file <<
  // "#freq_mHz;Z_re;Z_im;Z_abs;TH_re;TH_im;TH_abs;PH_re;PH_im;PH_abs\n";
  file2.setf(std::ios::fixed);
  file2 << std::setprecision(16);

  for (std::size_t idx = 0; idx < myff.i2() + 100; ++idx) {
    // index bounds are safe since nw == vec_w.size() and a_filt has
    // compatible columns
    file2 << (vec_w_ref[idx] * nval * 1000.0 / (twopi)) << ';'
          << vec_raw(0, idx).real() << ';' << vec_raw(0, idx).imag() << ';'
          << std::abs(vec_raw(0, idx)) << ';' << vec_raw(1, idx).real() << ';'
          << vec_raw(1, idx).imag() << ';' << std::abs(vec_raw(1, idx)) << ';'
          << vec_raw(2, idx).real() * (accel_norm * nval * myff.dt()) << ';'
          << vec_raw(2, idx).imag() * (accel_norm * nval * myff.dt()) << ';'
          << std::abs(vec_raw(2, idx)) * (accel_norm * nval * myff.dt())
          << '\n';
  }
  file2.close();

  //////////////////////////////////////////////////////////////////////////////
  // output mineos time series
  std::string pathtofile_mineos = "./work/benchmarks/bench_t_lf.out";
  std::ofstream file_mineos(pathtofile_mineos);
  if (!file_mineos) {
    std::cerr << "Error: unable to open output file: " << pathtofile_mineos
              << "\n";
    return 1;
  }
  file_mineos.setf(std::ios::fixed);
  file_mineos << std::setprecision(16);
  for (std::size_t idx = 0; idx < vec_filt_t_mineos.cols(); ++idx) {
    // index bounds are safe since nw == vec_w.size() and a_filt has
    // compatible columns
    auto tval = idx * myff.dt() * prem.TimeNorm();
    file_mineos << (idx * myff.dt() * prem.TimeNorm()) << ";"
                << vec_filt_t(0, idx) << ";" << vec_filt_t(1, idx) << ";"
                << vec_filt_t(2, idx) << ";" << vec_filt_t_yspec(0, idx) << ";"
                << vec_filt_t_yspec(1, idx) << ";" << vec_filt_t_yspec(2, idx)
                << ';' << vec_filt_t_mineos(0, idx) << ';'
                << vec_filt_t_mineos(1, idx) << ';' << vec_filt_t_mineos(2, idx)
                << '\n';
    if (tval > params.t_out() * 60.0) {
      break;
    }
  }
  file_mineos.close();

  return 0;
}
#endif