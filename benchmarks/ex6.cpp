#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#include "config.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <complex>
#include <cmath>
#include <vector>

#include <PlanetaryModel/All>
#include <DSpecM1D/Timer>
#include <DSpecM1D/All>
#include <SpectraSolver/FF>

int
main() {
  using Complex = std::complex<double>;
  using SparseMatrixC = Eigen::SparseMatrix<Complex>;

  Timer timer1;

  // ex6 is intentionally low-level: it works directly with the SEM matrices to
  // study tidal forcing, buoyancy structure, and kinetic-energy partitioning.
  std::string paramPath =
      std::string(PROJECT_BUILD_DIR) + "data/params/ex6.txt";
  InputParameters params(paramPath);
  std::string earthModelPath =
      std::string(PROJECT_BUILD_DIR) + "data/" + params.earth_model();

  auto cmt = SourceInfo::EarthquakeCMT(params);

  int lval = params.lmax();
  int nq = 5;
  double maxstep = 0.05;

  std::cout << "Enter max step size for SEM integration: \n";
  std::cin >> maxstep;

  // Build the Earth model and SEM mesh explicitly so the benchmark can inspect
  // mesh-level quantities as well as solve-level quantities.
  timer1.start();
  prem_norm<double> normClass;
  auto prem = EarthModels::ModelInput(earthModelPath, normClass, "true");

  Full1D::SEM sem(prem, maxstep, nq, lval);
  auto mesh = sem.mesh();
  auto meshmodel = sem.meshModel();

  std::cout << sem.meshModel().Gravity(sem.mesh().NE() - 1,
                                       sem.mesh().NN() - 1) *
                   prem.AccelerationNorm()
            << "\n";

  // Restrict the buoyancy-frequency calculation to the fluid shell.
  std::vector<std::vector<double>> vecN2;
  auto fsb = mesh.FS_Boundaries();
  int idxLow = fsb[0] + 1;
  int idxUp = fsb[1] + 1;

  double maxN2 = 0.0;
  double minN2 = 0.0;

  for (int idxe = idxLow; idxe < idxUp; ++idxe) {
    std::vector<double> tmp;
    tmp.reserve(mesh.NN());
    for (int idxn = 0; idxn < mesh.NN(); ++idxn) {
      auto crad = mesh.NodeRadius(idxe, idxn);
      auto rho = meshmodel.Density(idxe, idxn);
      auto drho_dr = prem.Density(1).Derivative(crad);
      auto gval = meshmodel.Gravity(idxe, idxn);

      auto tn2 =
          -drho_dr * gval / rho - rho * gval * gval / prem.Kappa(1)(crad);
      tmp.push_back(tn2);

      if (tn2 > maxN2)
        maxN2 = tn2;
      if (tn2 < minN2)
        minN2 = tn2;
    }
    vecN2.push_back(tmp);
  }

  // Sample a few long periods relevant to the tidal forcing experiment.
  auto maxImag = std::sqrt(-minN2) / (TWO_PI * prem.TimeNorm());
  auto nHours = maxImag * 3600.0;
  auto nHPeriod = 1.0 / nHours;

  std::vector<double> vecTPeriods{12.0, 48.0, 168.0, 1680.0};
  std::vector<Complex> vecOmega;
  vecOmega.reserve(vecTPeriods.size());

  for (double period : vecTPeriods) {
    auto freq = 1.0 / (period * 3600.0);
    Complex omegaC = TWO_PI * freq * prem.TimeNorm();
    vecOmega.push_back(omegaC);
  }

  // Assemble the tidal-forcing right-hand side directly in the SEM basis.
  auto idxl = 3;
  SparseMatrixC keS = sem.hS(idxl).cast<Complex>();
  SparseMatrixC inS = sem.pS(idxl).cast<Complex>();

  Eigen::VectorXcd rhs = Eigen::VectorXcd::Zero(keS.rows());

  for (int idx = 0; idx < mesh.NE(); ++idx) {
    auto jacval = 2.0 / (mesh.EUR(idx) - mesh.ELR(idx));
    for (int idxn = 0; idxn < mesh.NN(); ++idxn) {
      auto crad = mesh.NodeRadius(idx, idxn);
      auto idx_u = sem.ltgS(0, idx, idxn);
      auto idx_v = sem.ltgS(1, idx, idxn);
      auto tmp = mesh.GLL().W(idxn) * jacval * std::pow(crad, idxl + 1);

      rhs(idx_u) += idxl * tmp;
      rhs(idx_v) += idxl * (idxl + 1.0) * tmp * std::sqrt(2.0);
    }
  }

  // Solve one complex linear system per forcing frequency.
  timer1.start();
  Eigen::SparseLU<SparseMatrixC, Eigen::COLAMDOrdering<int>> solver;
  std::vector<Eigen::VectorXcd> vecSolutions;
  vecSolutions.reserve(vecOmega.size());

  for (std::size_t idx = 0; idx < vecOmega.size(); ++idx) {
    Complex omegaC = vecOmega[idx];
    SparseMatrixC A = -omegaC * omegaC * inS + keS;
    A.makeCompressed();
    solver.compute(A);

    if (solver.info() != Eigen::Success) {
      std::cerr << "Decomposition failed for frequency index: " << idx << "\n";
      return 1;
    }
    vecSolutions.push_back(solver.solve(rhs));
  }
  timer1.stop("Total time for solving linear system");

  // Integrate the total and fluid-only kinetic energies from the SEM solution.
  std::vector<double> vecTotalKE, vecFluidKE;
  vecTotalKE.reserve(vecOmega.size());
  vecFluidKE.reserve(vecOmega.size());

  auto energynorm = prem.MassNorm() * prem.LengthNorm() * prem.LengthNorm() /
                    (prem.TimeNorm() * prem.TimeNorm());

  for (std::size_t idxf = 0; idxf < vecOmega.size(); ++idxf) {
    auto sol = vecSolutions[idxf];
    double total_ke = 0.0;
    double fluid_ke = 0.0;

    for (int idxe = 0; idxe < mesh.NE(); ++idxe) {
      for (int idxn = 0; idxn < mesh.NN(); ++idxn) {
        auto crad = mesh.NodeRadius(idxe, idxn);
        auto rho = meshmodel.Density(idxe, idxn);
        auto jacval = 2.0 / (mesh.EUR(idxe) - mesh.ELR(idxe));
        auto wval = mesh.GLL().W(idxn);

        auto uidx = sem.ltgS(0, idxe, idxn);
        auto vidx = sem.ltgS(1, idxe, idxn);
        auto uval = sol(uidx);
        auto vval = sol(vidx);

        auto ke_local = jacval * wval * rho *
                        (std::norm(uval) + std::norm(vval)) *
                        std::pow(crad, 2.0);
        total_ke += ke_local;

        if (idxe >= idxLow && idxe < idxUp) {
          fluid_ke += ke_local;
        }
      }
    }
    vecTotalKE.push_back(total_ke * energynorm);
    vecFluidKE.push_back(fluid_ke * energynorm);
  }

  // Write the derived buoyancy-frequency profile, the frequency summary, and
  // the radial SEM response sampled across the mesh.
  std::string n2OutPath = std::string(PROJECT_BUILD_DIR) +
                          "../plotting/outputs/ex6_N2_" +
                          std::to_string(maxstep) + ".out";
  std::ofstream n2File(n2OutPath);
  if (!n2File) {
    std::cerr << "Error: unable to open output file: " << n2OutPath << "\n";
    return 1;
  }

  n2File << std::fixed << std::setprecision(22);
  auto fnorm2 = 1.0 / (prem.TimeNorm() * prem.TimeNorm());

  for (int idxe = idxLow; idxe < idxUp; ++idxe) {
    for (int idxn = 0; idxn < mesh.NN(); ++idxn) {
      auto crad = mesh.NodeRadius(idxe, idxn) * prem.LengthNorm();
      n2File << crad << ";" << vecN2[idxe - idxLow][idxn] * fnorm2 << "\n";
    }
  }
  n2File.close();

  std::string freqOutPath = std::string(PROJECT_BUILD_DIR) +
                            "../plotting/outputs/ex6_w_" +
                            std::to_string(maxstep) + ".out";
  std::ofstream freqFile(freqOutPath);
  if (!freqFile) {
    std::cerr << "Error: unable to open output file: " << freqOutPath << "\n";
    return 1;
  }

  freqFile << std::fixed << std::setprecision(22);
  for (std::size_t idx = 0; idx < vecOmega.size(); ++idx) {
    freqFile << vecTPeriods[idx] / nHPeriod << ";"
             << vecFluidKE[idx] / vecTotalKE[idx] << "\n";
  }
  freqFile.close();

  std::string pathToFile = std::string(PROJECT_BUILD_DIR) +
                           "../plotting/outputs/ex6_radial_response_" +
                           std::to_string(maxstep) + ".out";
  std::ofstream file(pathToFile);
  if (!file) {
    std::cerr << "Error: unable to open output file: " << pathToFile << "\n";
    return 1;
  }

  file << std::fixed << std::setprecision(22);
  for (int idxe = 0; idxe < mesh.NE(); ++idxe) {
    for (int idxq = 0; idxq < nq; ++idxq) {
      auto crad = mesh.NodeRadius(idxe, idxq) * prem.LengthNorm();
      auto uidx = sem.ltgS(0, idxe, idxq);
      auto vidx = sem.ltgS(1, idxe, idxq);

      file << crad;
      for (std::size_t idxf = 0; idxf < vecOmega.size(); ++idxf) {
        auto sol = vecSolutions[idxf];
        file << ";" << sol(uidx).real() << ";" << sol(uidx).imag() << ";"
             << std::abs(sol(uidx)) << ";" << sol(vidx).real() << ";"
             << sol(vidx).imag() << ";" << std::abs(sol(vidx));
      }
      file << '\n';
    }
  }
  file.close();

  return 0;
}
