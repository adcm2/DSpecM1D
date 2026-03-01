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

// #include <cmath>
// #include <Eigen/Core>
// #include <Spectra/MatOp/SparseGenMatProd.h>
// #include <Spectra/SymGEigsShiftSolver.h>
// #include <GSHTrans/Wigner>
// #include <Interpolation/Lagrange>
// #include "SourceInfo.h"
// #include <EarthMesh/All>
// #include "MatrixIndices.h"

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

  ///////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////
  // paths
  std::string pathtopremeig = "../mineos-1.0.2/OUTPUT/eprem_noocean_S_IC";
  std::string pathtoprem =
      "../mineos-1.0.2/DEMO/models/prem.200.no.noatten.txt";
  std::string pathtoprempert =
      "../mineos-1.0.2/DEMO/models/prem.200.no.noatten.perturb.txt";
  std::string pathtopremnew = "work/toroidal/prem_pert.out";
  // std::string pathtocmt = "./examples/china_cmt_event";
  std::string pathtocmt = "./examples/bolivia_cmt_event";
  // std::cout << pathpert << std::endl;

  // models
  auto prem = EarthModels::ModelInput(pathtoprem);
  auto prempert = EarthModels::ModelInput(pathtoprem);
  Toroidal_Bench::spherical_model homosphere(1.0, 1.0, 1.0);

  // trying frequency solver
  double f1 = 0.2, f2 = 1000 * maxfreq, dt = 20, tout = 100, df0 = 1.0,
         wtb = 0.05, t1 = 0, t2 = 128;
  int qex = 4;
  SpectraSolver::FreqFull myff(f1, f2, dt, tout, df0, wtb, t1, t2, qex,
                               prem.TimeNorm());
  std::cout << "f1: " << myff.f1() << ", f2: " << myff.f2()
            << ", df: " << myff.df() << ", nt: " << myff.nt() << "\n";
  // std::cout << myff.w()

  // test sem
  Timer timer1;
  timer1.start();
  Toroidal::sem sem(prem, maxstep, pn, lval, nsolid);
  //   Toroidal::sem sem2(homosphere, maxstep, pn, lval);
  timer1.stop("Toroidal::sem constructor");

  timer1.start();
  sem.CalculateEigenfrequencies(nfreq, sigshift);

  timer1.stop("Toroidal::sem CalculateEigenfrequencies");

  timer1.start();
  sem.CalculateEigenfrequenciesSeparate(nfreq, sigshift);
  timer1.stop("Toroidal::sem CalculateEigenfrequenciesSeparate");
  Eigen::VectorXd w2 = sem.efrequencies_gen();
  // sem.PrintModesUpToFreq(maxfreq);

  // auto eig_check = sem.efunctions_ref(2);

  // timer1.start();
  // sem2.CalculateEigenfrequencies(nfreq, sigshift);
  // timer1.stop("Toroidal::sem CalculateEigenfrequencies");
  // Eigen::VectorXd w1 = sem2.efrequencies_gen();

  auto cmt = SourceInfo::EarthquakeCMT(pathtocmt);

  timer1.start();
  Eigen::VectorXcd vec_force = sem.CalculateForce(cmt);
  timer1.stop("Toroidal::sem CalculateForce");

  // different frequencies
  timer1.start();
  Eigen::SparseMatrix<Complex> mat_ke =
      sem.GetStiffnessMatrix().cast<Complex>();

  Eigen::SparseMatrix<Complex> mat_inertia =
      sem.GetInertiaMatrix().cast<Complex>();
  mat_ke.makeCompressed();
  mat_inertia.makeCompressed();
  timer1.stop("Toroidal::sem GetStiffnessMatrix and InertiaMatrix");

  timer1.start();
  // get receiver vectors
  Eigen::VectorXcd vec_thetar =
      sem.ReceiverVectorThetaSurface(prem.OuterRadius(), theta_S, phi_S);
  Eigen::VectorXcd vec_phir =
      sem.ReceiverVectorPhiSurface(prem.OuterRadius(), theta_S, phi_S);
  timer1.stop("Toroidal::sem ReceiverVectorTheta and ReceiverVectorPhi");

  // vector of frequency values
  double nval = 1 / (prem.TimeNorm());
  // std::vector<double> vec_w(100, 0.0);
  // for (int idx = 0; idx < vec_w.size(); ++idx) {
  //   vec_w[idx] = idx * maxfreq / (nval * vec_w.size());
  // };
  // for (auto idx : vec_w) {
  //   std::cout << idx << std::endl;
  // };
  auto vec_w = myff.w();
  // for (auto &idx : vec_w) {
  //   idx *= 1.0 / nval;
  // };
  std::cout << vec_w[0] << " " << vec_w.back() << std::endl;
  std::cout << "Number of frequencies: " << vec_w.size() << std::endl;
  std::cout << "Low and high indices: " << myff.i1() << " " << myff.i2()
            << std::endl;
  // Complex ieps(0.0, -1.0e-2);
  Complex ieps = myff.ep() * Complex(0.0, -1.0);
  Complex myi = Complex(0.0, 1.0);
  Complex myvali = Complex(0.0, -0.2);

  // std::vector<Complex> vec_tm(vec_w.size(), 0.0), vec_pm(vec_w.size(), 0.0);
  std::vector<Complex> vec_tm2(vec_w.size(), 0.0), vec_pm2(vec_w.size(), 0.0);
  Eigen::MatrixXcd vec_tm2f(2, vec_w.size());
  vec_tm2f.setZero();
  std::vector<Complex> vec_tm3(vec_w.size(), 0.0), vec_pm3(vec_w.size(), 0.0);
  Eigen::SparseLU<Eigen::SparseMatrix<Complex>, Eigen::COLAMDOrdering<int>>
      solver;
  Eigen::VectorXcd vec_xstore;

  timer1.start();
  // for (int idxl = 1; idxl < lval + 1; ++idxl) {
  // get matrices for each l
  Eigen::SparseMatrix<Complex> mat_kel =
      sem.GetStiffnessMatrixL(lval).cast<Complex>();
  Eigen::SparseMatrix<Complex> mat_inertial =
      sem.GetInertiaMatrixL(lval).cast<Complex>();

  // get force vector for each l
  Eigen::MatrixXcd vec_force2 = sem.CalculateForce(cmt, lval);

  // get receiver vectors for each l
  Eigen::MatrixXcd vec_thetar2 =
      sem.ReceiverVectorThetaSurfaceL(prem.OuterRadius(), theta_S, phi_S, lval);
  Eigen::MatrixXcd vec_phir2 =
      sem.ReceiverVectorPhiSurfaceL(prem.OuterRadius(), theta_S, phi_S, lval);
  std::cout << "Number of rows and cols: \n"
            << vec_thetar2.rows() << " " << vec_thetar2.cols() << std::endl;

  ///////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////
  // checks for nmc
  auto eig_check = sem.efunctions_ref(lval);
  std::cout << "Size of eig_check: " << eig_check.rows() << " "
            << eig_check.cols() << std::endl;
  auto nrows = eig_check.rows();
  auto ncols = eig_check.cols();
  Eigen::MatrixXcd eig_check2 = eig_check(
      Eigen::seq(0, nrows - 1, 1), Eigen::seq(ncols - 1, 0, Eigen::fix<-1>));
  Eigen::MatrixXcd eig_checknorm = eig_check2;
  std::cout << "Check 1\n";
  for (int idx = 0; idx < eig_check2.cols(); ++idx) {
    eig_checknorm.col(idx) *=
        1.0 / (eig_check2.col(idx).norm() * eig_check2.col(idx).norm());
  }
  auto v22 = eig_check2.transpose() * vec_force2;
  auto mat_ke_check = eig_check2.transpose() * mat_kel * eig_check2;
  auto mat_inertia_check = eig_check2.transpose() * mat_inertial * eig_check2;
  std::cout << "Check 2\n";
  // std::cout << "\n\n" << mat_inertia_check << "\n\n";
  auto vec_theta_check = eig_check2.transpose() * vec_thetar2;
  auto vec_phi_check = eig_check2.transpose() * vec_phir2;
  std::cout << "Check 3\n";
  // forming vec_force3:
  Eigen::MatrixXcd vec_force3(vec_force2.rows(), vec_force2.cols());
  vec_force3.setZero();
  for (int idx = 0; idx < vec_force2.cols(); ++idx) {
    // for (int idx1 = 0; idx1 < vec_force2.rows(); ++idx1) {
    //   vec_force3(idx1, idx) += eig_check2(idx1, 0) * v22(0, idx);
    // }
    vec_force3.col(idx) = eig_checknorm * v22.col(idx);
  }
  std::cout << "Check 4\n";
  // matrix G, ie the transformation matrix
  // Eigen::MatrixXcd mat_G(sem.NumModesL(lval), sem.NumModesL(lval));
  // mat_G.setZero();
  // for (int idx = 0; idx < sem.NumModesL(lval); ++idx) {
  //   for (int idx1 = 0; idx1 < sem.NumModesL(lval); ++idx1) {
  //     mat_G(idx, idx1) = eig_check2.col(idx).transpose() *
  //     eig_check2.col(idx1);
  //   }
  // }
  // Eigen::FullPivLU<Eigen::MatrixXcd> lu(mat_G);
  Eigen::MatrixXcd vec_coeff_force(eig_check2.cols(), vec_force2.cols());
  for (int idx = 0; idx < eig_check2.cols(); ++idx) {
    // std::complex<double> nval =
    //     eig_check2.col(idx).transpose() * mat_inertial * eig_check2.col(idx);
    for (int idx1 = 0; idx1 < vec_force2.cols(); ++idx1) {
      vec_coeff_force(idx, idx1) =
          eig_check2.col(idx).transpose() * vec_force2.col(idx1);
      vec_coeff_force(idx, idx1) /= mat_inertia_check(idx, idx);
    }
  }
  std::cout << "Check 5\n";

  // auto vec_coeff_force = lu.solve(v22);
  Eigen::MatrixXcd vec_force4(vec_force2.rows(), vec_force2.cols());
  vec_force4.setZero();
  for (int idx = 0; idx < vec_coeff_force.cols(); ++idx) {
    vec_force4.col(idx) = mat_inertial * eig_check2 * vec_coeff_force.col(idx);
  }
  std::cout << "Check 6\n";
  ///////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////

  //   std::cout << "\n\n" << vec_thetar2.col(0) << "\n\n";
  // std::cout << myff.i1() << " " << myff.i2() << std::endl;
  // go through frequencies
  for (int idx = myff.i1(); idx < myff.i2(); ++idx) {
    Complex w = static_cast<Complex>(vec_w[idx]) + ieps;
    // Complex w = static_cast<Complex>(vec_w[idx]) + myvali;

    // if (idx == myff.i1()) {
    //   std::cout << std::setprecision(15) << idx << " " << w << "\n";
    // }
    Eigen::SparseMatrix<Complex> mat_w = -w * w * mat_inertial + mat_kel;
    mat_w.makeCompressed();
    solver.compute(mat_w);
    Eigen::MatrixXcd vec_fw = vec_force2 / (myi * w);
    Eigen::MatrixXcd vec_x = solver.solve(vec_fw);
    //   vec_tm2[idx] += -w * w * (vec_thetar2.cwiseProduct(vec_x)).sum();
    for (int j = 0; j < vec_thetar2.cols(); ++j) {
      // int j = 1;
      for (int i = 0; i < vec_thetar2.rows(); ++i) {
        vec_tm2[idx] += -w * w * vec_thetar2(i, j) * vec_x(i, j);
        vec_tm2f(0, idx) += -w * w * vec_thetar2(i, j) * vec_x(i, j);
        vec_pm2[idx] += -w * w * vec_phir2(i, j) * vec_x(i, j);
        vec_tm2f(1, idx) += -w * w * vec_phir2(i, j) * vec_x(i, j);
      }
    }
    //   vec_tm2[idx] += -w * w * vec_thetar2.transpose()*vec_x;
    //   vec_pm2[idx] += -w * w * (vec_phir2.cwiseProduct(vec_x)).sum();
    // if (std::abs(idx - myff.i2()) < 5) {
    //   std::cout << std::setprecision(15) << idx << " " << w.real() << " "
    //             << w.imag() << "\n";
    // }
  };
  // }
  timer1.stop("Toroidal::sem Solve");

  // timer1.start();
  sem.FindModesForCoupling(maxfreq);
  // test augmentation
  timer1.start();
  sem.augment_basis_calculate();
  timer1.stop("Toroidal::sem augment_basis_calculate");
  // bool to_augment = true;
  bool to_augment = false;
  bool aug2 = toaug;
  // std::cout << "Augment?\n";
  // std::cin >> aug2;
  // Eigen::MatrixXcd mat_nmc_ke = sem.NMC_KE(prem, to_augment);
  // Eigen::MatrixXcd mat_nmc_inertia = sem.NMC_INERTIA(prem, to_augment);
  // Eigen::VectorXcd vec_force_nmc = sem.NMC_FORCE(cmt, to_augment);
  // Eigen::VectorXcd vec_theta_nmc = sem.ReceiverVectorThetaSurfaceCoupling(
  //     prem.OuterRadius(), theta_S, phi_S, to_augment);
  // Eigen::VectorXcd vec_phi_nmc = sem.ReceiverVectorPhiSurfaceCoupling(
  //     prem.OuterRadius(), theta_S, phi_S, to_augment);

  // just at l = 2
  int idxl2 = lval;
  Eigen::MatrixXcd mat_nmc_ke2 = sem.NMC_KE(prem, idxl2, aug2);
  // std::cout << "Check 1\n";
  Eigen::MatrixXcd mat_nmc_inertia2 = sem.NMC_INERTIA(prem, idxl2, aug2);
  // std::cout << "Check 2\n";
  Eigen::MatrixXcd vec_force_nmc2 = sem.CalculateForceNMC(cmt, lval, aug2);
  // std::cout << "Check 3\n";
  Eigen::MatrixXcd vec_theta_nmc2 = sem.ReceiverVectorThetaSurface_NMCL(
      prem.OuterRadius(), theta_S, phi_S, idxl2, aug2);
  // std::cout << "Check 4\n";
  Eigen::MatrixXcd vec_phi_nmc2 = sem.ReceiverVectorPhiSurface_NMCL(
      prem.OuterRadius(), theta_S, phi_S, idxl2, aug2);
  // std::cout << "Check 5\n";

  auto soltol = 1e-9;
  // std::cout << vec_theta_nmc.cols() << "\n\n";
  // Eigen::MatrixXcd VR(vec_theta_nmc.rows(), 2);
  // VR.col(0) = vec_theta_nmc;
  // VR.col(1) = vec_phi_nmc;

  // auto vec_tm4 = myff.Spectra_Raw_NoCoriolis(mat_nmc_ke, mat_nmc_inertia, VR,
  //                                            vec_force_nmc, soltol);
  // timer1.stop("FreqFull Solve");
  // auto vec_response = myff.Spectra_Raw_NoCoriolis_LU(
  //     mat_nmc_ke, mat_nmc_inertia, VR, vec_force_nmc);
  std::cout << "Size of theta rec: " << vec_theta_nmc2.rows() << " "
            << vec_theta_nmc2.cols() << std::endl;
  std::cout << "Size of phi rec: " << vec_phi_nmc2.rows() << " "
            << vec_phi_nmc2.cols() << std::endl;
  Eigen::MatrixXcd vec_response2(2, vec_w.size());
  vec_response2.setZero();
  for (int idxm2 = -idxl2; idxm2 < idxl2 + 1; ++idxm2) {
    Eigen::MatrixXcd VR2(vec_theta_nmc2.rows(), 2);
    VR2.col(0) = vec_theta_nmc2.col(idxm2 + idxl2);
    VR2.col(1) = vec_phi_nmc2.col(idxm2 + idxl2);
    // std::cout << "Check 6: " << idxm2 << "\n";
    vec_response2 += myff.Spectra_Raw_NoCoriolis_LU(
        mat_nmc_ke2, mat_nmc_inertia2, VR2, vec_force_nmc2.col(idxm2 + idxl2));
  }

  timer1.stop("FreqFull Low Memory Solve");

  // filter
  auto vec_r2t = processfunctions::filtfreq2time(vec_response2, myff);
  auto vec_r2f = processfunctions::fulltime2freq(vec_r2t, myff);

  // std::cout << "Size of r2t: " << vec_r2t.rows() << " " << vec_r2t.cols() <<
  // std::endl; std::cout << "Size of r2f: " << vec_r2f.rows() << " " <<
  // vec_r2f.cols() << std::endl; std::cout << "Size of tm2f: " <<
  // vec_tm2f.rows() << " " << vec_tm2f.cols() << std::endl; std::cout << "Size
  // of pm2f: " << vec_pm2f.rows() << " " << vec_pm2f.cols() << std::endl;

  auto vec_r2t_b = processfunctions::filtfreq2time(vec_tm2f, myff);
  auto vec_r2f_b = processfunctions::fulltime2freq(vec_r2t_b, myff);
  // std::cout << "Check filter\n" << vec_r2f.rows() << " " << vec_r2f.cols() <<
  // " " << vec_response2.rows() << " " << vec_response2.cols() << "\n";
  ///////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////

  // check the matrices and vectors
  // Eigen::FullPivLU<Eigen::MatrixXcd> lu_check(mat_nmc_inertia2);
  // Eigen::MatrixXcd mat_check = lu_check.solve(mat_nmc_ke2);
  // std::cout << std::setprecision(3) << "\n\n"
  //           << mat_nmc_inertia2 << "\n\n"
  //           << mat_nmc_ke2 << "\n\n";

  ///////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////
  // check matrices and vectors

  {
    // auto vec_store = sem.modes_coupled();
    std::size_t ovidx = 0;
    auto totissue = 0;

    // std::cout << "Size of eig_check: " << eig_check.rows() << " "
    //           << eig_check.cols() << std::endl;
    // std::cout << "\n\n" << v22.rows() << " " << v22.cols() << std::endl;
    // std::cout << vec_force_nmc2.rows() << " " << vec_force_nmc2.cols()
    //           << std::endl;
    // std::cout << std::setprecision(3) << v22 << "\n\n" << std::endl;
    // std::cout << std::setprecision(3) << vec_force_nmc2 << "\n\n" <<
    // std::endl;
    auto v1v2norm = (vec_force_nmc2 - v22).norm();
    std::cout << "Force difference norm: " << v1v2norm << "\n";
    auto matkenorm = (mat_nmc_ke2 - mat_ke_check).norm();
    std::cout << "KE matrix difference norm: " << matkenorm << "\n";
    auto matinorm = (mat_nmc_inertia2 - mat_inertia_check).norm();
    std::cout << "Inertia matrix difference norm: " << matinorm << "\n";
    auto vecthetanorm = (vec_theta_nmc2 - vec_theta_check).norm();
    std::cout << "Theta receiver vector difference norm: " << vecthetanorm
              << "\n";
    auto vecphinorm = (vec_phi_nmc2 - vec_phi_check).norm();
    std::cout << "Phi receiver vector difference norm: " << vecphinorm << "\n";
  }

  ///////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  // output
  // outputting result
  std::string pathtofile = "./work/toroidal/groundresponsel.out";
  auto file = std::ofstream(pathtofile);
  for (int idx = 0; idx < vec_w.size(); ++idx) {
    file << std::setprecision(16)
         << vec_w[idx] * nval * 1000 / (2.0 * 3.1415926535) << ";"
         << vec_tm2[idx].real() << ";" << vec_tm2[idx].imag() << ";"
         << std::abs(vec_tm2[idx]) << ";" << vec_pm2[idx].real() << ";"
         << vec_pm2[idx].imag() << ";" << std::abs(vec_pm2[idx]) << ";"
         << vec_response2(0, idx).real() << ";" << vec_response2(0, idx).imag()
         << ";" << std::abs(vec_response2(0, idx)) << ";"
         << vec_response2(1, idx).real() << ";" << vec_response2(1, idx).imag()
         << ";" << std::abs(vec_response2(1, idx)) << std::endl;
  }
  file.close();

  std::string pathtofile2 = "./work/toroidal/groundresponsel_filt.out";
  auto file2 = std::ofstream(pathtofile2);
  for (int idx = 0; idx < vec_w.size(); ++idx) {
    file2 << std::setprecision(16)
          << vec_w[idx] * nval * 1000 / (2.0 * 3.1415926535) << ";"
          << vec_r2f_b(0, idx).real() << ";" << vec_r2f_b(0, idx).imag() << ";"
          << std::abs(vec_r2f_b(idx)) << ";" << vec_r2f_b(1, idx).real() << ";"
          << vec_r2f_b(1, idx).imag() << ";" << std::abs(vec_r2f_b(idx)) << ";"
          << vec_r2f(0, idx).real() << ";" << vec_r2f(0, idx).imag() << ";"
          << std::abs(vec_r2f(0, idx)) << ";" << vec_r2f(1, idx).real() << ";"
          << vec_r2f(1, idx).imag() << ";" << std::abs(vec_r2f(1, idx))
          << std::endl;
  }
  file2.close();

  return 0;
}
