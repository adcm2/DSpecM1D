#include <iostream>
#include <PlanetaryModel/All>
#include <new_coupling/Timer>
#include "spectra_master.h"
// #include "sem_spheroidal_debug.h"
#include "../SpectraSolver/SpectraSolver/FF"
// #include "../SpectraSolver/SpectraSolver/src/ODE_Spectra/filter_base.h"
// #include
// "../SpectraSolver/SpectraSolver/src/ODE_Spectra/postprocessfunctions.h"
// #include "read_station.h"
#include "input_parser.h"   // Use the new input parser
#include "start_element.h"

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
  std::string cpath = params.earth_model();
  std::string earth_model_path = "../YSpec/" + params.earth_model();

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // parameters of sem
  int nfreq = 5, nsolid = 1, lval = params.lmax(), NQ = 5;
  bool toaug = false;
  double maxstep = 0.001;
  const double twopi = 2.0 * 3.1415926535897932;
  std::string pathpert;

  //////////////////////////////
  // frequency solver parameters
  double f1 = params.f1(), f2 = params.f2(), dt = params.time_step_sec(),
         tout = params.t_out() / 60.0, df0 = 1.0, wtb = 0.05, t1 = 0,
         t2 = tout + 1.0;
  int qex = 4;
  Complex myi(0.0, 1.0);
  // lval = params.lmax();

  prem_norm<double> norm_class;
  auto prem = EarthModels::ModelInput(earth_model_path, norm_class, true);
  SpectraSolver::FreqFull myff(params.f1(), params.f2(), params.f11(),
                               params.f12(), params.f21(), params.f22(), dt,
                               tout, df0, wtb, t1, t2, qex, prem.TimeNorm());

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // prem
  // auto prem = EarthModels::ModelInput(pathtoprem);

  // sem class
  Full1D::specsem sem(prem, maxstep, NQ, lval);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  int l = 1;
  std::cout << "Enter l value to check start element: ";
  std::cin >> l;
  auto vec_w = myff.w();
  std::vector<std::vector<double>> vec_wout;
  std::vector<std::vector<double>> vec_rout;
  auto mfact = 6.0 * 2.0 * std::sqrt(2.0 / (l * (l + 1)));
  for (int idxl = 1; idxl < l + 1; ++idxl) {
    std::vector<double> wtmp, rtmp;
    for (int idx = 0; idx < myff.i2() + 1; ++idx) {
      wtmp.push_back(vec_w[idx]);
      auto tmp1 = SpectralTools::StartRadiusClean(sem, idxl, vec_w[idx], 0);
      auto tmp2 = SpectralTools::StartRadiusClean(sem, idxl, vec_w[idx], 1);
      auto tmp = std::min(tmp1, tmp2);
      rtmp.push_back(tmp);
    }
    vec_wout.push_back(wtmp);
    vec_rout.push_back(rtmp);
  }

  for (int idx = 1; idx < vec_wout.size(); ++idx) {
    // std::cout << vec_wout[idx] << ": " << vec_rout[idx] << ": "
    //           << vec_rout[idx - 1] - 0.5 * mfact *
    //                                      (vec_wout[idx] * vec_wout[idx] -
    //                                       vec_wout[idx - 1] * vec_wout[idx -
    //                                       1])
    //           << " "
    //           << vec_rout[idx - 1] -
    //                  16.0 * std::sqrt(2.0) *
    //                      std::log(vec_wout[idx] / vec_wout[idx - 1])
    //           << "\n";

    // std::cout << vec_wout[idx] << ": " << vec_rout[idx] << ": "
    //           << vec_rout[idx - 1] * vec_wout[idx - 1] / vec_wout[idx] <<
    //           "\n";
  }

  for (int idxl = 0; idxl < l; ++idxl) {
    std::string pathtofile = "./work/test/rad_start_";
    pathtofile = pathtofile + std::to_string(idxl + 2) + ".out";
    std::ofstream file(pathtofile);
    file << std::setprecision(16);

    for (std::size_t idx = 1; idx < vec_wout[idxl].size(); ++idx) {
      // index bounds are safe since nw == vec_w.size() and a_filt has
      // compatible columns
      file << vec_wout[idxl][idx] << ';'
           << vec_rout[idxl][idx] * prem.LengthNorm() / (6371000) << '\n';
    }
    file.close();
  }

  // for (int idxl = 0; idxl < l; ++idxl) {
  //   std::string pathtofile = "./work/test/rad_deriv_";
  //   pathtofile = pathtofile + std::to_string(idxl + 2) + ".out";
  //   std::ofstream file(pathtofile);
  //   file << std::setprecision(16);

  //   for (std::size_t idx = 1; idx < vec_wout[idxl].size(); ++idx) {
  //     // index bounds are safe since nw == vec_w.size() and a_filt has
  //     // compatible columns
  //     file << vec_wout[idxl][idx]  << ';'
  //          << vec_rout[idxl][idx] * prem.LengthNorm() / (6371000) << '\n';
  //   }
  //   file.close();
  // }

  // index calculation
  // check ltg_r
  //   for (int idxe = 0; idxe < sem.mesh().NE(); ++idxe) {
  //     for (int idxn = 0; idxn < NQ; ++idxn) {
  //       std::cout << "Element " << idxe << ", Node " << idxn << ", is solid?
  //       "
  //                 << sem.mesh().IsSolid(idxe)
  //                 << ", LtG_R(0) = " << sem.LtG_R(0, idxe, idxn) << "\n";
  //       std::cout << "Element " << idxe << ", Node " << idxn << ", is solid?
  //       "
  //                 << sem.mesh().IsSolid(idxe)
  //                 << ", LtG_R(1) = " << sem.LtG_R(1, idxe, idxn) << "\n";
  //     }
  //     std::cout << "-----------------------\n";
  //   }

  return 0;
}
