#include <iostream>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <Eigen/Core>
#include <Eigen/LU>
#include <PlanetaryModel/All>
// #include <SpectraSolver/FF>
#include "../SpectraSolver/SpectraSolver/FF"
#include "SourceInfo.h"
#include "MatrixIndices.h"
#include "toroidal_bench.h"
#include <new_coupling/Timer>
#include "sem_toroidal.h"
#include <Eigen/IterativeLinearSolvers>
#include <iomanip>

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
  int pn = 8;
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

  ///////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////
  // paths
  std::string pathtopremeig =
      "/space/adcm2/mineos-1.0.2/OUTPUT/eprem_noocean_S_IC";
  std::string pathtoprem =
      "/space/adcm2/mineos-1.0.2/DEMO/models/prem.200.no.noatten.txt";
  std::string pathtoprempert =
      "/space/adcm2/mineos-1.0.2/DEMO/models/prem.200.no.noatten.perturb.txt";
  std::string pathtopremnew = "work/toroidal/prem_pert.out";
  std::string pathtocmt = "./examples/china_cmt_event";
  // std::cout << pathpert << std::endl;

  // models
  auto prem = EarthModels::ModelInput(pathtoprem);
  auto prempert = EarthModels::ModelInput(pathpert);
  Toroidal_Bench::spherical_model homosphere(1.0, 1.0, 1.0);

  // truying frequency solver
  double f1 = 0.2, f2 = 1000 * maxfreq, dt = 20, tout = 32, df0 = 2.0,
         wtb = 0.05, t1 = 0, t2 = 60;
  int qex = 2;
  SpectraSolver::FreqFull myff(f1, f2, dt, tout, df0, wtb, t1, t2, qex,
                               prem.TimeNorm());
  std::cout << "f1: " << myff.f1() << ", f2: " << myff.f2()
            << ", df: " << myff.df() << ", nt: " << myff.nt() << "\n";
  // std::cout << myff.w()

  // test sem
  Timer timer1;
  timer1.start();
  Toroidal::sem sem(prem, maxstep, pn, lval, nsolid);
  Toroidal::sem sem2(homosphere, maxstep, pn, lval);
  timer1.stop("Toroidal::sem constructor");

  timer1.start();
  sem.CalculateEigenfrequencies(nfreq, sigshift);

  timer1.stop("Toroidal::sem CalculateEigenfrequencies");

  timer1.start();
  sem.CalculateEigenfrequenciesSeparate(10, sigshift);
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

  std::vector<Complex> vec_tm(vec_w.size(), 0.0), vec_pm(vec_w.size(), 0.0);
  std::vector<Complex> vec_tm2(vec_w.size(), 0.0), vec_pm2(vec_w.size(), 0.0);
  std::vector<Complex> vec_tm3(vec_w.size(), 0.0), vec_pm3(vec_w.size(), 0.0);
  Eigen::SparseLU<Eigen::SparseMatrix<Complex>, Eigen::COLAMDOrdering<int>>
      solver;
  Eigen::VectorXcd vec_xstore;
  timer1.start();
  // {
  //   Eigen::SparseMatrix<Complex> mat_test = mat_inertia + mat_ke;
  //   mat_test.makeCompressed();
  //   solver.analyzePattern(mat_test);
  // }
  for (int idx = myff.i1(); idx < myff.i2(); ++idx) {
    Complex w = static_cast<Complex>(vec_w[idx]) + ieps;
    Eigen::SparseMatrix<Complex> mat_w = -w * w * mat_inertia + mat_ke;
    mat_w.makeCompressed();
    solver.compute(mat_w);
    Eigen::VectorXcd vec_fw = vec_force / (myi * w);
    Eigen::VectorXcd vec_x = solver.solve(vec_fw);
    // if (std::abs(idx - myff.i1()) < 5) {
    //   std::cout << idx << " " << vec_thetar.transpose() * vec_x << " "
    //             << vec_phir.transpose() * vec_x << "\n";
    // }

    // if (std::abs(idx - myff.i2()) < 5) {
    //   std::cout << std::setprecision(15) << idx << " " << w.real() << "  "
    //             << w.imag() << "\n";
    // }
    vec_tm[idx] = -w * w * vec_thetar.transpose() * vec_x;
    vec_pm[idx] = -w * w * vec_phir.transpose() * vec_x;
  };
  timer1.stop("Toroidal::sem Solve");

  {
    timer1.start();
    for (int idxl = 1; idxl < lval + 1; ++idxl) {
      // get matrices for each l
      Eigen::SparseMatrix<Complex> mat_kel =
          sem.GetStiffnessMatrixL(idxl).cast<Complex>();
      Eigen::SparseMatrix<Complex> mat_inertial =
          sem.GetInertiaMatrixL(idxl).cast<Complex>();

      // analyze pattern for each l
      // {
      //   Eigen::SparseMatrix<Complex> mat_testl = mat_inertial + mat_kel;
      //   mat_testl.makeCompressed();
      //   solver.analyzePattern(mat_testl);
      // }

      // get force vector for each l
      Eigen::MatrixXcd vec_force2 = sem.CalculateForce(cmt, idxl);

      // get receiver vectors for each l
      Eigen::MatrixXcd vec_thetar2 = sem.ReceiverVectorThetaSurfaceL(
          prem.OuterRadius(), theta_S, phi_S, idxl);
      Eigen::MatrixXcd vec_phir2 = sem.ReceiverVectorPhiSurfaceL(
          prem.OuterRadius(), theta_S, phi_S, idxl);

      // std::cout << myff.i1() << " " << myff.i2() << std::endl;
      // go through frequencies
      for (int idx = myff.i1(); idx < myff.i2(); ++idx) {
        Complex w = static_cast<Complex>(vec_w[idx]) + ieps;
        // if (idx == myff.i1()) {
        //   std::cout << std::setprecision(15) << idx << " " << w << "\n";
        // }
        Eigen::SparseMatrix<Complex> mat_w = -w * w * mat_inertial + mat_kel;
        mat_w.makeCompressed();
        solver.compute(mat_w);
        Eigen::MatrixXcd vec_fw = vec_force2 / (myi * w);
        Eigen::MatrixXcd vec_x = solver.solve(vec_fw);
        vec_tm2[idx] += -w * w * (vec_thetar2.cwiseProduct(vec_x)).sum();
        vec_pm2[idx] += -w * w * (vec_phir2.cwiseProduct(vec_x)).sum();
        // if (std::abs(idx - myff.i2()) < 5) {
        //   std::cout << std::setprecision(15) << idx << " " << w.real() << " "
        //             << w.imag() << "\n";
        // }
      };
    }
    timer1.stop("Toroidal::sem Solve");
  }

  // timer1.start();
  sem.FindModesForCoupling(maxfreq);
  // test augmentation
  timer1.start();
  sem.augment_basis_calculate();
  timer1.stop("Toroidal::sem augment_basis_calculate");
  // bool to_augment = true;
  bool to_augment = false;
  bool aug2 = false;
  Eigen::MatrixXcd mat_nmc_ke = sem.NMC_KE(prem, to_augment);
  Eigen::MatrixXcd mat_nmc_inertia = sem.NMC_INERTIA(prem, to_augment);
  Eigen::VectorXcd vec_force_nmc = sem.NMC_FORCE(cmt, to_augment);
  Eigen::VectorXcd vec_theta_nmc = sem.ReceiverVectorThetaSurfaceCoupling(
      prem.OuterRadius(), theta_S, phi_S, to_augment);
  Eigen::VectorXcd vec_phi_nmc = sem.ReceiverVectorPhiSurfaceCoupling(
      prem.OuterRadius(), theta_S, phi_S, to_augment);

  // just at l = 2
  int idxl2 = 2;
  Eigen::MatrixXcd mat_nmc_ke2 = sem.NMC_KE(prem, idxl2, aug2);
  Eigen::MatrixXcd mat_nmc_inertia2 = sem.NMC_INERTIA(prem, idxl2, aug2);
  Eigen::MatrixXcd vec_force_nmc2 = sem.CalculateForceNMC(cmt, idxl2, aug2);
  Eigen::MatrixXcd vec_theta_nmc2 = sem.ReceiverVectorThetaSurface_NMCL(
      prem.OuterRadius(), theta_S, phi_S, idxl2, aug2);
  Eigen::MatrixXcd vec_phi_nmc2 = sem.ReceiverVectorPhiSurface_NMCL(
      prem.OuterRadius(), theta_S, phi_S, idxl2, aug2);

  auto soltol = 1e-9;
  // std::cout << vec_theta_nmc.cols() << "\n\n";
  Eigen::MatrixXcd VR(vec_theta_nmc.rows(), 2);
  VR.col(0) = vec_theta_nmc;
  VR.col(1) = vec_phi_nmc;

  // auto vec_tm4 = myff.Spectra_Raw_NoCoriolis(mat_nmc_ke, mat_nmc_inertia, VR,
  //                                            vec_force_nmc, soltol);
  // timer1.stop("FreqFull Solve");
  auto vec_response = myff.Spectra_Raw_NoCoriolis_LU(
      mat_nmc_ke, mat_nmc_inertia, VR, vec_force_nmc);
  Eigen::MatrixXcd vec_response2 = vec_response;
  vec_response2.setZero();
  for (int idxm2 = -idxl2; idxm2 < idxl2 + 1; ++idxm2) {
    Eigen::MatrixXcd VR2(vec_theta_nmc2.rows(), 2);
    VR2.col(0) = vec_theta_nmc2.col(idxm2 + idxl2);
    VR2.col(1) = vec_phi_nmc2.col(idxm2 + idxl2);
    vec_response2 += myff.Spectra_Raw_No_Coriolis_Low_Memory(
        mat_nmc_ke2, mat_nmc_inertia2, VR2, vec_force_nmc2.col(idxm2 + idxl2),
        soltol);
  }

  timer1.stop("FreqFull Low Memory Solve");

  // auto vec_pm4 = myff.Spectra_Raw_NoCoriolis(mat_nmc_ke, mat_nmc_inertia,
  //                                            vec_phi_nmc, vec_force_nmc,
  //                                            soltol);

  ///////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  // output
  // outputting result
  std::string pathtofile = "./work/toroidal/groundresponse.out";
  auto file = std::ofstream(pathtofile);
  for (int idx = 0; idx < vec_w.size(); ++idx) {
    file << std::setprecision(16)
         << vec_w[idx] * nval * 1000 / (2.0 * 3.1415926535) << ";"
         << vec_tm[idx].real() << ";" << vec_tm[idx].imag() << ";"
         << std::abs(vec_tm[idx]) << ";" << vec_pm[idx].real() << ";"
         << vec_pm[idx].imag() << ";" << std::abs(vec_pm[idx]) << ";"
         << std::abs(vec_tm2[idx]) << ";" << std::abs(vec_pm2[idx]) << ";"
         << std::abs(vec_response(0, idx)) << ";"
         << std::abs(vec_response(1, idx)) << ";" << vec_response(0, idx).real()
         << ";" << vec_response(1, idx).real() << ";"
         << vec_response(0, idx).imag() << ";" << vec_response(1, idx).imag()
         << ";" << std::abs(vec_response2(0, idx)) << ";"
         << std::abs(vec_response2(1, idx)) << ";"
         << vec_response2(0, idx).real() << ";" << vec_response2(1, idx).real()
         << ";" << vec_response2(0, idx).imag() << ";"
         << vec_response2(1, idx).imag() << std::endl;
  }
  file.close();

  return 0;
}

/*
{
auto vec_store = sem.modes_coupled();
std::size_t ovidx = 0;
auto totissue = 0;
for (int idx = 0; idx < vec_store.size(); ++idx) {
  auto idxl = vec_store[idx].GetL();
  Eigen::SparseMatrix<Complex> mat_ke2 =
      sem.GetStiffnessMatrixL(idxl).cast<Complex>();
  Eigen::SparseMatrix<Complex> mat_inertia2 =
      sem.GetInertiaMatrixL(idxl).cast<Complex>();
  auto eig_check = sem.efunctions_ref(idxl);
  Eigen::MatrixXcd vec_force2 = sem.CalculateForce(cmt, idxl);
  auto v1 = vec_force_nmc.segment(ovidx, 2 * idxl + 1);
  auto v2 =
      eig_check.col(9 - vec_store[idx].GetN()).transpose() * vec_force2;

  auto v1v2norm = (v1 - v2).norm();
  if (v1v2norm > 1e-10) {
    std::cout << "Large force check norm for l = " << idxl
              << ", overtone = " << vec_store[idx].GetN() << ":\n"
              << v1v2norm << "\n\n";
    ++totissue;
  }

  auto matcheck =
      mat_nmc_ke.block(ovidx, ovidx, 2 * idxl + 1, 2 * idxl + 1);
  Complex matcheck2 = eig_check.col(9 - vec_store[idx].GetN()).transpose() *
                      mat_ke2 * eig_check.col(9 - vec_store[idx].GetN());
  auto matcheck3 = matcheck2 * Eigen::MatrixXcd::Identity(matcheck.rows(),
                                                          matcheck.cols());
  // std::cout << std::setprecision(15) << "KE check for l = " << idxl
  //           << ", overtone = " << vec_store[idx].GetN() << ":\n"
  //           << (matcheck - matcheck3).norm() << "\n\n";
  auto mchecknorm = (matcheck - matcheck3).norm();
  if (mchecknorm > 1e-10) {
    std::cout << "Large KE check norm for l = " << idxl
              << ", overtone = " << vec_store[idx].GetN() << ":\n"
              << mchecknorm << "\n\n";
    ++totissue;
  }

  auto matchecki =
      mat_nmc_inertia.block(ovidx, ovidx, 2 * idxl + 1, 2 * idxl + 1);
  Complex matchecki2 =
      eig_check.col(9 - vec_store[idx].GetN()).transpose() * mat_inertia2 *
      eig_check.col(9 - vec_store[idx].GetN());
  auto matchecki3 = matchecki2 * Eigen::MatrixXcd::Identity(
                                     matchecki.rows(), matchecki.cols());
  // std::cout << std::setprecision(15) << "Inertia check for l = " << idxl
  //           << ", overtone = " << vec_store[idx].GetN() << ":\n"
  //           << (matchecki - matchecki3).norm() << "\n\n";

  auto mcheckinorm = (matchecki - matchecki3).norm();
  if (mcheckinorm > 1e-10) {
    std::cout << "Large Inertia check norm for l = " << idxl
              << ", overtone = " << vec_store[idx].GetN() << ":\n"
              << mcheckinorm << "\n\n";
    ++totissue;
  }

  // check receiver vectors
  Eigen::VectorXcd vec_thetar2 = sem.ReceiverVectorThetaSurfaceL(
      prem.OuterRadius(), theta_S, phi_S, idxl);
  Eigen::VectorXcd vec_phir2 = sem.ReceiverVectorPhiSurfaceL(
      prem.OuterRadius(), theta_S, phi_S, idxl);
  auto v3 = vec_theta_nmc.segment(ovidx, 2 * idxl + 1);
  auto v4 =
      vec_thetar2.transpose() * eig_check.col(9 - vec_store[idx].GetN());
  // std::cout << std::setprecision(15)
  //           << "Receiver theta check f/*
{
auto vec_store = sem.modes_coupled();
std::size_t ovidx = 0;
auto totissue = 0;
for (int idx = 0; idx < vec_store.size(); ++idx) {
  auto idxl = vec_store[idx].GetL();
  Eigen::SparseMatrix<Complex> mat_ke2 =
      sem.GetStiffnessMatrixL(idxl).cast<Complex>();
  Eigen::SparseMatrix<Complex> mat_inertia2 =
      sem.GetInertiaMatrixL(idxl).cast<Complex>();
  auto eig_check = sem.efunctions_ref(idxl);
  Eigen::MatrixXcd vec_force2 = sem.CalculateForce(cmt, idxl);
  auto v1 = vec_force_nmc.segment(ovidx, 2 * idxl + 1);
  auto v2 =
      eig_check.col(9 - vec_store[idx].GetN()).transpose() * vec_force2;

  auto v1v2norm = (v1 - v2).norm();
  if (v1v2norm > 1e-10) {
    std::cout << "Large force check norm for l = " << idxl
              << ", overtone = " << vec_store[idx].GetN() << ":\n"
              << v1v2norm << "\n\n";
    ++totissue;
  }

  auto matcheck =
      mat_nmc_ke.block(ovidx, ovidx, 2 * idxl + 1, 2 * idxl + 1);
  Complex matcheck2 = eig_check.col(9 - vec_store[idx].GetN()).transpose() *
                      mat_ke2 * eig_check.col(9 - vec_store[idx].GetN());
  auto matcheck3 = matcheck2 * Eigen::MatrixXcd::Identity(matcheck.rows(),
                                                          matcheck.cols());
  // std::cout << std::setprecision(15) << "KE check for l = " << idxl
  //           << ", overtone = " << vec_store[idx].GetN() << ":\n"
  //           << (matcheck - matcheck3).norm() << "\n\n";
  auto mchecknorm = (matcheck - matcheck3).norm();
  if (mchecknorm > 1e-10) {
    std::cout << "Large KE check norm for l = " << idxl
              << ", overtone = " << vec_store[idx].GetN() << ":\n"
              << mchecknorm << "\n\n";
    ++totissue;
  }

  auto matchecki =
      mat_nmc_inertia.block(ovidx, ovidx, 2 * idxl + 1, 2 * idxl + 1);
  Complex matchecki2 =
      eig_check.col(9 - vec_store[idx].GetN()).transpose() * mat_inertia2 *
      eig_check.col(9 - vec_store[idx].GetN());
  auto matchecki3 = matchecki2 * Eigen::MatrixXcd::Identity(
                                     matchecki.rows(), matchecki.cols());
  // std::cout << std::setprecision(15) << "Inertia check for l = " << idxl
  //           << ", overtone = " << vec_store[idx].GetN() << ":\n"
  //           << (matchecki - matchecki3).norm() << "\n\n";

  auto mcheckinorm = (matchecki - matchecki3).norm();
  if (mcheckinorm > 1e-10) {
    std::cout << "Large Inertia check norm for l = " << idxl
              << ", overtone = " << vec_store[idx].GetN() << ":\n"
              << mcheckinorm << "\n\n";
    ++totissue;
  }

  // check receiver vectors
  Eigen::VectorXcd vec_thetar2 = sem.ReceiverVectorThetaSurfaceL(
      prem.OuterRadius(), theta_S, phi_S, idxl);
  Eigen::VectorXcd vec_phir2 = sem.ReceiverVectorPhiSurfaceL(
      prem.OuterRadius(), theta_S, phi_S, idxl);
  auto v3 = vec_theta_nmc.segment(ovidx, 2 * idxl + 1);
  auto v4 =
      vec_thetar2.transpose() * eig_check.col(9 - vec_store[idx].GetN());
  // std::cout << std::setprecision(15)
  //           << "Receiver theta check for l = " << idxl
  //           << ", overtone = " << vec_store[idx].GetN() << ":\n"
  //           << (v3 - v4).norm() << "\n\n";

  auto v3v4norm = (v3 - v4).norm();
  if (v3v4norm > 1e-10) {
    std::cout << "Large receiver theta check norm for l = " << idxl
              << ", overtone = " << vec_store[idx].GetN() << ":\n"
              << v3v4norm << "\n\n";
    ++totissue;
  }

  // check phi receiver vectors
  auto v5 = vec_phi_nmc.segment(ovidx, 2 * idxl + 1);
  auto v6 =
      vec_phir2.transpose() * eig_check.col(9 - vec_store[idx].GetN());
  auto v5v6norm = (v5 - v6).norm();
  if (v5v6norm > 1e-10) {
    std::cout << "Large receiver phi chec/*
{
auto vec_store = sem.modes_coupled();
std::size_t ovidx = 0;
auto totissue = 0;
for (int idx = 0; idx < vec_store.size(); ++idx) {
  auto idxl = vec_store[idx].GetL();
  Eigen::SparseMatrix<Complex> mat_ke2 =
      sem.GetStiffnessMatrixL(idxl).cast<Complex>();
  Eigen::SparseMatrix<Complex> mat_inertia2 =
      sem.GetInertiaMatrixL(idxl).cast<Complex>();
  auto eig_check = sem.efunctions_ref(idxl);
  Eigen::MatrixXcd vec_force2 = sem.CalculateForce(cmt, idxl);
  auto v1 = vec_force_nmc.segment(ovidx, 2 * idxl + 1);
  auto v2 =
      eig_check.col(9 - vec_store[idx].GetN()).transpose() * vec_force2;

  auto v1v2norm = (v1 - v2).norm();
  if (v1v2norm > 1e-10) {
    std::cout << "Large force check norm for l = " << idxl
              << ", overtone = " << vec_store[idx].GetN() << ":\n"
              << v1v2norm << "\n\n";
    ++totissue;
  }

  auto matcheck =
      mat_nmc_ke.block(ovidx, ovidx, 2 * idxl + 1, 2 * idxl + 1);
  Complex matcheck2 = eig_check.col(9 - vec_store[idx].GetN()).transpose() *
                      mat_ke2 * eig_check.col(9 - vec_store[idx].GetN());
  auto matcheck3 = matcheck2 * Eigen::MatrixXcd::Identity(matcheck.rows(),
                                                          matcheck.cols());
  // std::cout << std::setprecision(15) << "KE check for l = " << idxl
  //           << ", overtone = " << vec_store[idx].GetN() << ":\n"
  //           << (matcheck - matcheck3).norm() << "\n\n";
  auto mchecknorm = (matcheck - matcheck3).norm();
  if (mchecknorm > 1e-10) {
    std::cout << "Large KE check norm for l = " << idxl
              << ", overtone = " << vec_store[idx].GetN() << ":\n"
              << mchecknorm << "\n\n";
    ++totissue;
  }

  auto matchecki =
      mat_nmc_inertia.block(ovidx, ovidx, 2 * idxl + 1, 2 * idxl + 1);
  Complex matchecki2 =
      eig_check.col(9 - vec_store[idx].GetN()).transpose() * mat_inertia2 *
      eig_check.col(9 - vec_store[idx].GetN());
  auto matchecki3 = matchecki2 * Eigen::MatrixXcd::Identity(
                                     matchecki.rows(), matchecki.cols());
  // std::cout << std::setprecision(15) << "Inertia check for l = " << idxl
  //           << ", overtone = " << vec_store[idx].GetN() << ":\n"
  //           << (matchecki - matchecki3).norm() << "\n\n";

  auto mcheckinorm = (matchecki - matchecki3).norm();
  if (mcheckinorm > 1e-10) {
    std::cout << "Large Inertia check norm for l = " << idxl
              << ", overtone = " << vec_store[idx].GetN() << ":\n"
              << mcheckinorm << "\n\n";
    ++totissue;
  }

  // check receiver vectors
  Eigen::VectorXcd vec_thetar2 = sem.ReceiverVectorThetaSurfaceL(
      prem.OuterRadius(), theta_S, phi_S, idxl);
  Eigen::VectorXcd vec_phir2 = sem.ReceiverVectorPhiSurfaceL(
      prem.OuterRadius(), theta_S, phi_S, idxl);
  auto v3 = vec_theta_nmc.segment(ovidx, 2 * idxl + 1);
  auto v4 =
      vec_thetar2.transpose() * eig_check.col(9 - vec_store[idx].GetN());
  // std::cout << std::setprecision(15)
  //           << "Receiver theta check for l = " << idxl
  //           << ", overtone = " << vec_store[idx].GetN() << ":\n"
  //           << (v3 - v4).norm() << "\n\n";

  auto v3v4norm = (v3 - v4).norm();
  if (v3v4norm > 1e-10) {
    std::cout << "Large receiver theta check norm for l = " << idxl
              << ", overtone = " << vec_store[idx].GetN() << ":\n"
              << v3v4norm << "\n\n";
    ++totissue;
  }

  // check phi receiver vectors
  auto v5 = vec_phi_nmc.segment(ovidx, 2 * idxl + 1);
  auto v6 =
      vec_phir2.transpose() * eig_check.col(9 - vec_store[idx].GetN());
  auto v5v6norm = (v5 - v6).norm();
  if (v5v6norm > 1e-10) {
    std::cout << "Large receiver phi check norm for l = " << idxl
              << ", overtone = " << vec_store[idx].GetN() << ":\n"
              << v5v6norm << "\n\n";
    ++totissue;
  }
  // increment
  ovidx += 2 * idxl + 1;
}
if (totissue == 0) {
  std::cout << "All checks passed!\n\n";
} else {
  std::cout << "Total number of issues found: " << totissue << "\n\n";
}

// std::cout << "Force 2: \n"
//           << eig_check.col(9).transpose() *
//                  vec_force2;   // get force vector for each l
}
k norm for l = " << idxl
              << ", overtone = " << vec_store[idx].GetN() << ":\n"
              << v5v6norm << "\n\n";
++totissue;
}
// increment
ovidx += 2 * idxl + 1;
}
if (totissue == 0) {
  std::cout << "All checks passed!\n\n";
} else {
  std::cout << "Total number of issues found: " << totissue << "\n\n";
}

// std::cout << "Force 2: \n"
//           << eig_check.col(9).transpose() *
//                  vec_force2;   // get force vector for each l
}
 or l = " << idxl
    //           << ", overtone = " << vec_store[idx].GetN() << ":\n"
    //           << (v3 - v4).norm() << "\n\n";

    auto v3v4norm = (v3 - v4).norm();
if (v3v4norm > 1e-10) {
  std::cout << "Large receiver theta check norm for l = " << idxl
            << ", overtone = " << vec_store[idx].GetN() << ":\n"
            << v3v4norm << "\n\n";
  ++totissue;
}

// check phi receiver vectors
auto v5 = vec_phi_nmc.segment(ovidx, 2 * idxl + 1);
auto v6 = vec_phir2.transpose() * eig_check.col(9 - vec_store[idx].GetN());
auto v5v6norm = (v5 - v6).norm();
if (v5v6norm > 1e-10) {
  std::cout << "Large receiver phi check norm for l = " << idxl
            << ", overtone = " << vec_store[idx].GetN() << ":\n"
            << v5v6norm << "\n\n";
  ++totissue;
}
// increment
ovidx += 2 * idxl + 1;
}
if (totissue == 0) {
  std::cout << "All checks passed!\n\n";
} else {
  std::cout << "Total number of issues found: " << totissue << "\n\n";
}

// std::cout << "Force 2: \n"
//           << eig_check.col(9).transpose() *
//                  vec_force2;   // get force vector for each l
}
*/