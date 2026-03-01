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
  InputParameters params("../YSpec/input_full.txt");

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
         wtb = 0.05, t1 = 0, t2 = tout + 1.0;
  int qex = 4;
  Complex myi(0.0, 1.0);
  double droptol = 1e-4;
  std::cout << "Drop tol: " << "\n";
  // std::cin >> droptol;
  // lval = params.lmax();

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // prem
  // auto prem = EarthModels::ModelInput(pathtoprem);

  prem_norm<double> norm_class;
  auto prem = EarthModels::ModelInput(earth_model_path, norm_class, "true");

  // frequency class
  SpectraSolver::FreqFull myff(params.f1(), params.f2(), params.f11(),
                               params.f12(), params.f21(), params.f22(), dt,
                               tout, df0, wtb, t1, t2, qex, prem.TimeNorm());
  auto vec_w = myff.w();
  // for (int idx = myff.i1(); idx < myff.i2(); ++idx) {
  //   if ((idx == myff.i1()) || (idx == myff.i2() - 1)) {
  //     std::cout << "Frequency index: " << idx << std::setprecision(10)
  //               << ", w: " << vec_w[idx] * 1000.0 / (prem.TimeNorm() * twopi)
  //               << "\n";
  //   }
  // }
  std::cout << "Frequency step: " << std::setprecision(16)
            << (vec_w[1] - vec_w[0]) * 1000.0 / (prem.TimeNorm() * twopi)
            << "\n";
  std::cout << std::setprecision(15) << "ep: " << myff.ep() << "\n";

  // sem class
  Full1D::sem sem(prem, maxstep, NQ, lval);

  // source information
  auto cmt = SourceInfo::EarthquakeCMT(params);

  //////////////////////////////////////////////////////////////////////////////

  // test new functionality
  std::cout << "Testing new functionality\n";
  SPARSESPEC::Sparse_F_Spec mytest;
  timer1.start();
  MATRIX vec_raw = mytest.FrequencySpectrum_TEST_CLEAN(myff, sem, prem, cmt,
                                                       params, droptol);
  timer1.stop("Total time for sparse frequency spectrum");

  // normalise
  auto accel_norm = prem.LengthNorm() / (prem.TimeNorm() * prem.TimeNorm());
  std::cout << "Acceleration norm: " << accel_norm << "\n";
  std::cout << "Time norm: " << prem.TimeNorm() << "\n";
  vec_raw *= accel_norm;

  // process responses
  auto vec_r2t_b = processfunctions::filtfreq2time(vec_raw, myff);

  //////////////////////////////////////////////////////////////////////////////
  // Repeat output for raw data
  // std::cout << "Third output\n";
  std::string pathtofile2 =
      "./work/spheroidal/groundresponse_raw_transverse.out";
  std::ofstream file2(pathtofile2);

  // Write header (optional)
  // file <<
  // "#freq_mHz;Z_re;Z_im;Z_abs;TH_re;TH_im;TH_abs;PH_re;PH_im;PH_abs\n";
  file2.setf(std::ios::fixed);
  file2 << std::setprecision(22);
  auto nval = 1.0 / prem.TimeNorm();
  for (std::size_t idx = 0; idx < myff.i2() + 100; ++idx) {
    // index bounds are safe since nw == vec_w.size() and a_filt has
    // compatible columns
    file2 << (vec_w[idx] * 1000.0 / (prem.TimeNorm() * twopi)) << ';'
          << vec_raw(0, idx).real() << ';' << vec_raw(0, idx).imag() << ';'
          << std::abs(vec_raw(0, idx)) << ';' << vec_raw(1, idx).real() << ';'
          << vec_raw(1, idx).imag() << ';' << std::abs(vec_raw(1, idx)) << ';'
          << vec_raw(2, idx).real() << ';' << vec_raw(2, idx).imag() << ';'
          << std::abs(vec_raw(2, idx)) << '\n';
  }
  file2.close();

  // output time series
  std::string pathtofile3 = "./work/spheroidal/check_time.out";
  std::ofstream file3(pathtofile3);
  file3.setf(std::ios::fixed);
  file3 << std::setprecision(22);
  for (std::size_t idx = 0; idx < vec_r2t_b.cols(); ++idx) {
    auto ct = idx * myff.dt() * prem.TimeNorm();
    file3 << ct << ';' << vec_r2t_b(0, idx) << ';' << vec_r2t_b(1, idx) << ';'
          << vec_r2t_b(2, idx) << '\n';
    if (ct > params.t_out() * 60.0) {
      break;
    }
  }

  file3.close();
  return 0;
}
