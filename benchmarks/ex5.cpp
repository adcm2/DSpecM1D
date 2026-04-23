#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#include "config.h"
#include "PaperExampleSupport.h"
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
  using MatrixC = Eigen::MatrixXcd;
  Timer timer1;

  // ex5 is a convergence study, so it keeps the solver setup and postprocessing
  // explicit while reusing the shared normalization/filter helpers.
  std::string paramPath =
      std::string(PROJECT_BUILD_DIR) + "data/params/ex5.txt";
  InputParameters params(paramPath);
  std::string earthModelPath =
      std::string(PROJECT_BUILD_DIR) + "data/" + params.earth_model();

  prem_norm<double> normClass;
  auto prem = EarthModels::ModelInput(earthModelPath, normClass, "true");
  auto cmt = SourceInfo::EarthquakeCMT(params);

  int lval = params.lmax();
  int nq = 6;
  int qex = 1;
  double maxstep = 0.05;

  double dt = params.time_step_sec();
  double tout = params.t_out() / 60.0;
  double df0 = 1.0;
  double wtb = 0.05;
  double t1 = 0.0;
  double t2 = tout;

  // Build the frequency helper once and reuse it for every convergence step.
  SpectraSolver::FreqFull myff(params.f1(), params.f2(), params.f11(),
                               params.f12(), params.f21(), params.f22(), dt,
                               tout, df0, wtb, t1, t2, qex, prem.TimeNorm());
  auto vecW = myff.w();

  // Sweep a logarithmically decreasing SEM step size to study convergence.
  int nSteps = 50;
  double step0 = 2.0 * 0.63 / myff.f22();
  std::cout << "Initial step: " << step0 << "\n";

  std::vector<double> vecStep;
  vecStep.reserve(nSteps);
  for (int idx = 0; idx < nSteps; ++idx) {
    vecStep.push_back(
        step0 / std::pow(10.0, 2.0 * idx / static_cast<double>(nSteps - 1)));
  }

  double normFactor = PaperExamples::legacyNormFactor(params, prem);
  double hannW = 0.2;
  auto filterOptions = PaperExamples::makeFilterOptions(hannW);

  SPARSESPEC::SparseFSpec specSolver;
  std::vector<Eigen::MatrixXcd> vecFinalW;
  std::vector<Eigen::MatrixXd> vecFinalT;
  vecFinalW.reserve(nSteps);
  vecFinalT.reserve(nSteps);

  // For each SEM step size, rebuild the mesh, solve, and store the filtered
  // frequency/time responses for later error analysis.
  for (int idx = 0; idx < nSteps; ++idx) {
    maxstep = vecStep[idx];
    int nskip = 3 * static_cast<int>(
                        std::floor(maxstep / ((vecW[1] - vecW[0]) * 0.003))) +
                1;

    Full1D::SEM sem(prem, maxstep, nq, lval);
    std::cout << "\nDoing step: " << maxstep << ", nskip: " << nskip << "\n";

    MatrixC vecRaw = specSolver.spectra(myff, sem, prem, cmt, params, nskip);
    vecRaw *= normFactor;
    auto filtered = DSpecM::applyFilter(vecRaw, myff, filterOptions);

    vecFinalW.push_back(filtered.frequencySeries);
    vecFinalT.push_back(filtered.timeSeries);
  }

  // Compare every run against the finest-step result treated as the reference.
  std::size_t numErrSteps = vecStep.size() - 1;
  std::vector<std::vector<double>> vecErrT(numErrSteps,
                                           std::vector<double>(3, 0.0));
  std::vector<std::vector<double>> vecErrW(numErrSteps,
                                           std::vector<double>(3, 0.0));
  std::vector<std::vector<double>> vecL2ErrT(numErrSteps,
                                             std::vector<double>(3, 0.0));
  std::vector<std::vector<double>> vecL2ErrW(numErrSteps,
                                             std::vector<double>(3, 0.0));

  int idxOut = vecFinalT[0].cols() - 1;
  for (int idx = 0; idx < vecFinalT[0].cols(); ++idx) {
    auto tval = idx * myff.dt() * prem.TimeNorm();
    if (tval > t2 * 3600.0) {
      idxOut = idx;
      break;
    }
  }

  for (int idx = 0; idx < nSteps; ++idx) {
    vecFinalT[idx]
        .block(0, idxOut, 3, vecFinalT[idx].cols() - idxOut)
        .setZero();
  }

  auto idxBack = nSteps - 1;
  auto vecTEx = vecFinalT[idxBack];
  auto vecWEx = vecFinalW[idxBack];

  auto tmpT = 100.0 / idxOut;
  auto mvT1 = tmpT / vecTEx.row(0).lpNorm<Eigen::Infinity>();
  auto mvT2 = tmpT / vecTEx.row(1).lpNorm<Eigen::Infinity>();
  auto mvT3 = tmpT / vecTEx.row(2).lpNorm<Eigen::Infinity>();

  auto mvT1L2 = 100.0 / vecTEx.row(0).norm();
  auto mvT2L2 = 100.0 / vecTEx.row(1).norm();
  auto mvT3L2 = 100.0 / vecTEx.row(2).norm();

  auto wCols = vecFinalW[idxBack].cols();
  auto tmpW = 100.0 / wCols;
  auto mvW1 = tmpW / vecWEx.row(0).lpNorm<Eigen::Infinity>();
  auto mvW2 = tmpW / vecWEx.row(1).lpNorm<Eigen::Infinity>();
  auto mvW3 = tmpW / vecWEx.row(2).lpNorm<Eigen::Infinity>();

  auto mvW1L2 = 100.0 / vecWEx.row(0).norm();
  auto mvW2L2 = 100.0 / vecWEx.row(1).norm();
  auto mvW3L2 = 100.0 / vecWEx.row(2).norm();

  for (std::size_t idx = 0; idx < numErrSteps; ++idx) {
    auto vecTErr = vecFinalT[idx] - vecTEx;
    auto vecWErr = vecFinalW[idx] - vecWEx;

    vecErrT[idx][0] = vecTErr.row(0).lpNorm<1>() * mvT1;
    vecErrT[idx][1] = vecTErr.row(1).lpNorm<1>() * mvT2;
    vecErrT[idx][2] = vecTErr.row(2).lpNorm<1>() * mvT3;
    vecErrW[idx][0] = vecWErr.row(0).lpNorm<1>() * mvW1;
    vecErrW[idx][1] = vecWErr.row(1).lpNorm<1>() * mvW2;
    vecErrW[idx][2] = vecWErr.row(2).lpNorm<1>() * mvW3;

    vecL2ErrT[idx][0] = vecTErr.row(0).norm() * mvT1L2;
    vecL2ErrT[idx][1] = vecTErr.row(1).norm() * mvT2L2;
    vecL2ErrT[idx][2] = vecTErr.row(2).norm() * mvT3L2;
    vecL2ErrW[idx][0] = vecWErr.row(0).norm() * mvW1L2;
    vecL2ErrW[idx][1] = vecWErr.row(1).norm() * mvW2L2;
    vecL2ErrW[idx][2] = vecWErr.row(2).norm() * mvW3L2;
  }

  // Write the stored convergence responses plus the derived L1/L2 error tables.
  double nval = 1.0 / prem.TimeNorm();

  std::string pathToFreqFile = std::string(PROJECT_BUILD_DIR) +
                               "../plotting/outputs/ex5_w_NQ" +
                               std::to_string(nq) + "_step.out";
  std::ofstream freqFile(pathToFreqFile);
  if (!freqFile) {
    std::cerr << "Error: unable to open output file_w: " << pathToFreqFile
              << "\n";
    return 1;
  }

  freqFile << std::fixed << std::setprecision(22);
  for (std::size_t idx = 0; idx < myff.i2() + 100; ++idx) {
    freqFile << (vecW[idx] * nval * 1000.0 / TWO_PI);
    for (std::size_t idx2 = 0; idx2 < vecStep.size(); ++idx2) {
      freqFile << ";" << vecFinalW[idx2](0, idx).real() << ';'
               << vecFinalW[idx2](0, idx).imag() << ';'
               << std::abs(vecFinalW[idx2](0, idx)) << ';'
               << vecFinalW[idx2](1, idx).real() << ';'
               << vecFinalW[idx2](1, idx).imag() << ';'
               << std::abs(vecFinalW[idx2](1, idx)) << ';'
               << vecFinalW[idx2](2, idx).real() << ';'
               << vecFinalW[idx2](2, idx).imag() << ';'
               << std::abs(vecFinalW[idx2](2, idx));
    }
    freqFile << '\n';
  }
  freqFile.close();

  std::string pathToTimeFile = std::string(PROJECT_BUILD_DIR) +
                               "../plotting/outputs/ex5_t_NQ" +
                               std::to_string(nq) + "_step.out";
  std::ofstream timeFile(pathToTimeFile);
  if (!timeFile) {
    std::cerr << "Error: unable to open output file_t: " << pathToTimeFile
              << "\n";
    return 1;
  }

  timeFile << std::fixed << std::setprecision(22);
  for (std::size_t idx = 0; idx < static_cast<std::size_t>(idxOut); ++idx) {
    timeFile << (idx * myff.dt()) * prem.TimeNorm();
    for (std::size_t idx2 = 0; idx2 < static_cast<std::size_t>(nSteps);
         ++idx2) {
      timeFile << ";" << vecFinalT[idx2](0, idx) << ";"
               << vecFinalT[idx2](1, idx) << ";" << vecFinalT[idx2](2, idx);
    }
    timeFile << '\n';
  }
  timeFile.close();

  std::string pathToErrFile =
      std::string(PROJECT_BUILD_DIR) + "../plotting/outputs/ex5_NQ" +
      std::to_string(nq) + "_step_error_" +
      std::to_string(static_cast<int>(params.f22())) + ".out";
  std::ofstream errFile(pathToErrFile);
  if (!errFile) {
    std::cerr << "Error: unable to open output file_err: " << pathToErrFile
              << "\n";
    return 1;
  }

  errFile << std::fixed << std::setprecision(22);
  for (std::size_t idx = 0; idx < numErrSteps; ++idx) {
    errFile << vecStep[idx] << ";" << vecErrT[idx][0] << ";" << vecErrT[idx][1]
            << ";" << vecErrT[idx][2] << ";" << vecErrW[idx][0] << ";"
            << vecErrW[idx][1] << ";" << vecErrW[idx][2] << "\n";
  }
  errFile.close();

  std::string pathToL2File =
      std::string(PROJECT_BUILD_DIR) + "../plotting/outputs/ex5_NQ" +
      std::to_string(nq) + "_step_error_l2_" +
      std::to_string(static_cast<int>(params.f22())) + ".out";
  std::ofstream l2File(pathToL2File);
  if (!l2File) {
    std::cerr << "Error: unable to open output file_l2: " << pathToL2File
              << "\n";
    return 1;
  }

  l2File << std::fixed << std::setprecision(22);
  for (std::size_t idx = 0; idx < numErrSteps; ++idx) {
    l2File << vecStep[idx] << ";" << vecL2ErrT[idx][0] << ";"
           << vecL2ErrT[idx][1] << ";" << vecL2ErrT[idx][2] << ";"
           << vecL2ErrW[idx][0] << ";" << vecL2ErrW[idx][1] << ";"
           << vecL2ErrW[idx][2] << "\n";
  }
  l2File.close();

  return 0;
}
