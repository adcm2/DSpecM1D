#include <iostream>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <Eigen/Core>
#include <Eigen/LU>
#include <PlanetaryModel/All>
#include "../SpectraSolver/SpectraSolver/FF"
#include "SourceInfo.h"
#include "MatrixIndices.h"
#include "toroidal_bench.h"
#include <new_coupling/Timer>
#include "sem_toroidal.h"
#include <Eigen/IterativeLinearSolvers>
#include <iomanip>
#include "../SpectraSolver/SpectraSolver/src/ODE_Spectra/filter_base.h"
#include "../SpectraSolver/SpectraSolver/src/ODE_Spectra/postprocessfunctions.h"

int
main() {
  using Complex = std::complex<double>;
  // parameters
  int nfreq = 5;
  int nsolid = 1;
  int pn = 6;
  double twopi = 2.0 * 3.1415926535897932;
  double maxstep;
  int lval;
  double sigshift;
  bool toaug;
  double maxfreq = 0.002;
  double theta_S = 0.4, phi_S = 3.1;
  std::string pathpert;
  std::cout << "Enter maxstep: \n";
  std::cin >> maxstep;
  std::cout << "Enter shift:\n";
  std::cin >> sigshift;
  std::cout << "Enter l:\n";
  std::cin >> lval;
  std::cout << "Enter # of frequencies:\n";
  std::cin >> nfreq;
  std::cout << "Enter solid layer:\n";
  std::cin >> nsolid;
  std::cout << "To aug or not?\n";
  std::cin >> toaug;
  std::cout << "Path to perturbed model:\n";
  std::cin >> pathpert;
  std::cout << "maxfreq:\n";
  std::cin >> maxfreq;

  //////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////
  // paths
  std::string pathtopremeig = "../mineos-1.0.2/OUTPUT/eprem_noocean_S_IC";
  std::string pathtoprem =
      "../mineos-1.0.2/DEMO/models/prem.200.no.noatten.txt";
  std::string pathtoprempert =
      "../mineos-1.0.2/DEMO/models/prem.200.no.noatten.perturb.txt";
  std::string pathtopremnew = "work/toroidal/prem_pert.out";
  //   std::string pathtocmt = "./examples/china_cmt_event";
  std::string pathtocmt = "./examples/bolivia_cmt_event";

  // models
  auto prem = EarthModels::ModelInput(pathtoprem);
  auto prempert = EarthModels::ModelInput(pathpert);
  Toroidal_Bench::spherical_model homosphere(1.0, 1.0, 1.0);

  // trying frequency solver
  double f1 = 0.2, f2 = 1000 * maxfreq, dt = 20, tout = 100, df0 = 1.0,
         wtb = 0.05, t1 = 0, t2 = 128;
  int qex = 4;

  // frequency solver
  SpectraSolver::FreqFull myff(f1, f2, dt, tout, df0, wtb, t1, t2, qex,
                               prem.TimeNorm());
  double nval = 1 / (prem.TimeNorm());
  std::cout << "f1: " << myff.f1() << ", f2: " << myff.f2()
            << ", df: " << myff.df() << ", nt: " << myff.nt() << "\n";
  // std::cout << myff.w()

  // test sem
  Timer timer1;
  timer1.start();
  Toroidal::sem sem(prem, maxstep, pn, lval, nsolid);
  //   Toroidal::sem sem2(prem, maxstep, pn, lval, nsolid);
  Toroidal::sem sem_p(prempert, maxstep, pn, lval, nsolid);
  timer1.stop("Toroidal::sem constructor");

  timer1.start();
  sem.CalculateEigenfrequenciesSeparate(nfreq, sigshift);
  //   sem2.CalculateEigenfrequenciesSeparate(nfreq, sigshift);
  sem_p.CalculateEigenfrequenciesSeparate(nfreq, sigshift);
  timer1.stop("Toroidal::sem CalculateEigenfrequenciesSeparate");

  timer1.start();
  sem.FindModesForCoupling(maxfreq);
  sem.augment_basis_calculate();
  //   sem2.augment_basis_calculate();
  sem_p.augment_basis_calculate();
  timer1.stop("Augmentation");

  // get source
  auto cmt = SourceInfo::EarthquakeCMT(pathtocmt);
  std::cout << "Got CMT\n";
  // frequency vector
  auto vec_w = myff.w();
  std::cout << "Got freq vector\n\n";

  //////////////////////////////////////////////////////////////////

  // get matrices and vectors
  Eigen::MatrixXcd km_0 = sem.NMC_KE(prem, lval, toaug);
  Eigen::MatrixXcd in_0 = sem.NMC_INERTIA(prem, lval, toaug);
  Eigen::MatrixXcd f_0 = sem.CalculateForceNMC(cmt, lval, toaug);
  Eigen::MatrixXcd theta_0 = sem.ReceiverVectorThetaSurface_NMCL(
      prem.OuterRadius(), theta_S, phi_S, lval, toaug);
  Eigen::MatrixXcd phi_0 = sem.ReceiverVectorPhiSurface_NMCL(
      prem.OuterRadius(), theta_S, phi_S, lval, toaug);

  // getting raw spectra
  std::cout << "Size of force: " << f_0.rows() << " " << f_0.cols()
            << std::endl;

  auto soltol = 1e-9;
  Eigen::MatrixXcd gm_0(2, vec_w.size());
  gm_0.setZero();
  for (int idxm2 = -lval; idxm2 < lval + 1; ++idxm2) {
    Eigen::MatrixXcd VR2(theta_0.rows(), 2);
    VR2.col(0) = theta_0.col(idxm2 + lval);
    VR2.col(1) = phi_0.col(idxm2 + lval);
    std::cout << "Working on m: " << idxm2 << "\n";
    gm_0 +=
        myff.Spectra_Raw_NoCoriolis_LU(km_0, in_0, VR2, f_0.col(idxm2 + lval));
  }
  std::cout << "Finished getting raw spectra 0\n";

  // processing
  auto r2t_0 = processfunctions::filtfreq2time(gm_0, myff);
  auto r2f_0 = processfunctions::fulltime2freq(r2t_0, myff);
  std::cout << "Finished processing 0\n";

  //////////////////////////////////////////////////////////////////

  // get matrices and vectors
  Eigen::MatrixXcd km_0p = sem.NMC_KE(prempert, lval, toaug);
  Eigen::MatrixXcd in_0p = sem.NMC_INERTIA(prempert, lval, toaug);
  Eigen::MatrixXcd f_0p = sem.CalculateForceNMC(cmt, lval, toaug);
  Eigen::MatrixXcd theta_0p = sem.ReceiverVectorThetaSurface_NMCL(
      prempert.OuterRadius(), theta_S, phi_S, lval, toaug);
  Eigen::MatrixXcd phi_0p = sem.ReceiverVectorPhiSurface_NMCL(
      prempert.OuterRadius(), theta_S, phi_S, lval, toaug);

  std::cout << "First few rows of km_0p:\n";
  std::cout << std::setprecision(15) << (km_0p - km_0).norm() << "\n\n";
  //   std::cout << "First few rows of in_0p:\n";
  //   std::cout << in_0p - in_0 << "\n\n";

  // look at matrices
  //   std::cout << "First few rows of km_0p:\n";
  //   std::cout << std::setprecision(3) << km_0p.block(0, 0, 5, 5) << "\n\n";
  //   std::cout << "First few rows of in_0p:\n";
  //   std::cout << in_0p.block(0, 0, 5, 5) << "\n\n";

  // getting raw spectra
  Eigen::MatrixXcd gm_0p(2, vec_w.size());
  gm_0p.setZero();
  for (int idxm2 = -lval; idxm2 < lval + 1; ++idxm2) {
    Eigen::MatrixXcd VR2(theta_0p.rows(), 2);
    VR2.col(0) = theta_0p.col(idxm2 + lval);
    VR2.col(1) = phi_0p.col(idxm2 + lval);
    std::cout << "Working on m: " << idxm2 << "\n";
    gm_0p += myff.Spectra_Raw_NoCoriolis_LU(km_0p, in_0p, VR2,
                                            f_0p.col(idxm2 + lval));
  }

  std::cout << "Finished getting raw spectra 0\n";
  std::cout << gm_0p.block(0, 0, 2, 5) << "\n";

  // processing
  auto r2t_0p = processfunctions::filtfreq2time(gm_0p, myff);
  auto r2f_0p = processfunctions::fulltime2freq(r2t_0p, myff);
  std::cout << "Finished processing 0\n";
  //////////////////////////////////////////////////////////////////

  // get matrices and vectors
  Eigen::MatrixXcd km_p = sem_p.NMC_KE(prempert, lval, toaug);
  Eigen::MatrixXcd in_p = sem_p.NMC_INERTIA(prempert, lval, toaug);
  Eigen::MatrixXcd f_p = sem_p.CalculateForceNMC(cmt, lval, toaug);
  Eigen::MatrixXcd theta_p = sem_p.ReceiverVectorThetaSurface_NMCL(
      prempert.OuterRadius(), theta_S, phi_S, lval, toaug);
  Eigen::MatrixXcd phi_p = sem_p.ReceiverVectorPhiSurface_NMCL(
      prempert.OuterRadius(), theta_S, phi_S, lval, toaug);
  std::cout << "Finished getting matrices p\n";

  // getting raw spectra
  //   auto soltol = 1e-9;
  Eigen::MatrixXcd gm_p(2, vec_w.size());
  gm_p.setZero();
  for (int idxm2 = -lval; idxm2 < lval + 1; ++idxm2) {
    Eigen::MatrixXcd VR2(theta_p.rows(), 2);
    VR2.col(0) = theta_p.col(idxm2 + lval);
    VR2.col(1) = phi_p.col(idxm2 + lval);
    gm_p +=
        myff.Spectra_Raw_NoCoriolis_LU(km_p, in_p, VR2, f_p.col(idxm2 + lval));
  }

  // processing
  auto r2t_p = processfunctions::filtfreq2time(gm_p, myff);
  auto r2f_p = processfunctions::fulltime2freq(r2t_p, myff);

  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  std::string pathtofile = "./work/toroidal/groundresponsel_comp.out";
  auto file = std::ofstream(pathtofile);
  for (int idx = 0; idx < vec_w.size(); ++idx) {
    // file << std::setprecision(16)
    //      << vec_w[idx] * nval * 1000 / (2.0 * 3.1415926535) << ";"
    //      << r2f_0(0, idx).real() << ";" << r2f_0(0, idx).imag() << ";"
    //      << std::abs(r2f_0(0, idx)) << ";" << r2f_0(1, idx).real() << ";"
    //      << r2f_0(1, idx).imag() << ";" << std::abs(r2f_0(1, idx)) <<
    //      std::endl;
    file << std::setprecision(16)
         << vec_w[idx] * nval * 1000 / (2.0 * 3.1415926535) << ";"
         << r2f_0(0, idx).real() << ";" << r2f_0(0, idx).imag() << ";"
         << std::abs(r2f_0(0, idx)) << ";" << r2f_0(1, idx).real() << ";"
         << r2f_0(1, idx).imag() << ";" << std::abs(r2f_0(1, idx)) << ";"
         << r2f_p(0, idx).real() << ";" << r2f_p(0, idx).imag() << ";"
         << std::abs(r2f_p(0, idx)) << ";" << r2f_p(1, idx).real() << ";"
         << r2f_p(1, idx).imag() << ";" << std::abs(r2f_p(1, idx)) << ";"
         << r2f_0p(0, idx).real() << ";" << r2f_0p(0, idx).imag() << ";"
         << std::abs(r2f_0p(0, idx)) << ";" << r2f_0p(1, idx).real() << ";"
         << r2f_0p(1, idx).imag() << ";" << std::abs(r2f_0p(1, idx))
         << std::endl;
  }
  file.close();

  std::string pathtofile2 = "./work/toroidal/groundresponsel_comp2.out";
  auto file2 = std::ofstream(pathtofile2);
  for (int idx = 0; idx < vec_w.size(); ++idx) {
    // file << std::setprecision(16)
    //      << vec_w[idx] * nval * 1000 / (2.0 * 3.1415926535) << ";"
    //      << r2f_0(0, idx).real() << ";" << r2f_0(0, idx).imag() << ";"
    //      << std::abs(r2f_0(0, idx)) << ";" << r2f_0(1, idx).real() << ";"
    //      << r2f_0(1, idx).imag() << ";" << std::abs(r2f_0(1, idx)) <<
    //      std::endl;
    file2 << std::setprecision(16)
          << vec_w[idx] * nval * 1000 / (2.0 * 3.1415926535) << ";"
          << gm_0(0, idx).real() << ";" << gm_0(0, idx).imag() << ";"
          << std::abs(gm_0(0, idx)) << ";" << gm_0(1, idx).real() << ";"
          << gm_0(1, idx).imag() << ";" << std::abs(gm_0(1, idx)) << ";"
          << gm_p(0, idx).real() << ";" << gm_p(0, idx).imag() << ";"
          << std::abs(gm_p(0, idx)) << ";" << gm_p(1, idx).real() << ";"
          << gm_p(1, idx).imag() << ";" << std::abs(gm_p(1, idx)) << ";"
          << gm_0p(0, idx).real() << ";" << gm_0p(0, idx).imag() << ";"
          << std::abs(gm_0p(0, idx)) << ";" << gm_0p(1, idx).real() << ";"
          << gm_0p(1, idx).imag() << ";" << std::abs(gm_0p(1, idx))
          << std::endl;
  }
  file2.close();

  return 0;
}
