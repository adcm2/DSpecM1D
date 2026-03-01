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
  // norm class
  prem_norm<double> norm_class;
  auto timenorm = norm_class.TimeNorm();

  // Read all parameters from the input file for record section
  InputParameters params("../YSpec/input_record.txt");

  // earth model
  std::string earth_model_path = "../YSpec/" + params.earth_model();

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // frequency solver parameters
  double dt = params.time_step_sec(), tout = params.t_out() / 60.0, df0 = 1.0,
         wtb = 0.05, t1 = 0, t2 = tout + 1.0, f1 = params.f1(),
         f2 = params.f2(), f11 = params.f11(), f12 = params.f12(),
         f21 = params.f21(), f22 = params.f22();
  int lval = params.lmax(), NQ = 5;
  int qex = 1;
  double maxstep = 0.003;
  // std::cout << "Max step: \n";
  // std::cin >> maxstep;

  // frequency class
  SpectraSolver::FreqFull myff(f1, f2, f11, f12, f21, f22, dt, tout, df0, wtb,
                               t1, t2, qex, timenorm);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // prem
  auto prem = EarthModels::ModelInput(earth_model_path, norm_class, "true");

  // sem class
  Full1D::specsem sem(prem, maxstep, NQ, lval);

  // source information
  auto cmt = SourceInfo::EarthquakeCMT(params);

  //////////////////////////////////////////////////////////////////////////////
  // get raw spectra
  SPARSESPEC::Sparse_F_Spec mytest;
  int nskip = 3 * std::floor(maxstep / ((myff.w(1) - myff.w(0)) * 0.003));
  if (nskip == 0) {
    nskip = 1;
  }
  timer1.start();
  MATRIX vec_raw = mytest.FrequencySpectrum_TEST_SPECSEM(myff, sem, prem, cmt,
                                                         params, nskip);
  timer1.stop("Total time for sparse frequency spectrum");

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // process responses
  // procedure is freq -> time -> filter -> freq -> filter -> time
  // filter chosen is Hann window. In time domain taken with 0.01 at start and
  // end and in frequency domain with frequencies in the input file
  auto vec_r2t_b = processfunctions::freq2time(vec_raw, myff);
  auto a_filt0 = processfunctions::fulltime2freq(vec_r2t_b, myff, 0.01);
  auto vec_filt_t = processfunctions::filtfreq2time(a_filt0, myff, false);

  // Dimensionalise
  auto accel_norm = prem.LengthNorm() / (prem.TimeNorm() * prem.TimeNorm());
  vec_filt_t *= accel_norm;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // output record section
  std::string pathtofile = "./work/benchmarks/record_section.out";
  std::ofstream file(pathtofile);
  file.setf(std::ios::fixed);
  file << std::setprecision(22);

  for (std::size_t idx = 0; idx < vec_filt_t.cols(); ++idx) {
    file << idx * myff.dt() * prem.TimeNorm();
    for (auto jidx = 0; jidx < params.num_receivers(); ++jidx) {
      file << ';' << vec_filt_t(3 * jidx + 0, idx) << ';'
           << vec_filt_t(3 * jidx + 1, idx) << ';'
           << vec_filt_t(3 * jidx + 2, idx);
    }
    file << std::endl;
    if ((idx * myff.dt() * prem.TimeNorm()) > params.t_out() * 60.0) {
      break;
    }
  }
  file.close();

  return 0;
}
