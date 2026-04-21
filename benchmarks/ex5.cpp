#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

// Standard library includes
#include "config.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <complex>
#include <cmath>
#include <vector>

// Project-specific includes
#include <PlanetaryModel/All>
#include <DSpecM1D/Timer>
#include <DSpecM1D/All>
#include <SpectraSolver/FF>

// constexpr double PI = 3.1415926535897932;
// constexpr double TWO_PI = 2.0 * PI;

// template <typename FLOAT> class prem_norm {
// public:
//   prem_norm() = default;

//   FLOAT LengthNorm() const { return _length_norm; }
//   FLOAT MassNorm() const { return _mass_norm; }
//   FLOAT TimeNorm() const { return _time_norm; }

// private:
//   FLOAT _length_norm = 1000.0;
//   FLOAT _mass_norm = 5515.0 * std::pow(_length_norm, 3.0);
//   FLOAT _time_norm = 1.0 / std::sqrt(PI * 6.67230e-11 * 5515.0);
// };

int
main() {
  using MatrixC = Eigen::MatrixXcd;
  Timer timer1;

  // --- 1. Read Inputs & Earth Model ---
  // get paths required for input parameters and Earth model
  std::string paramPath =
      std::string(PROJECT_BUILD_DIR) + "data/params/ex5.txt";
  InputParameters params(paramPath);
  std::string earthModelPath =
      std::string(PROJECT_BUILD_DIR) + "data/" + params.earth_model();

  prem_norm<double> normClass;
  auto prem = EarthModels::ModelInput(earthModelPath, normClass, "true");
  auto cmt = SourceInfo::EarthquakeCMT(params);

  // --- 2. Parameters ---
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

  // --- 3. Setup Frequency Class ---
  SpectraSolver::FreqFull myff(params.f1(), params.f2(), params.f11(),
                               params.f12(), params.f21(), params.f22(), dt,
                               tout, df0, wtb, t1, t2, qex, prem.TimeNorm());
  auto vecW = myff.w();

  // --- 4. Setup Convergence Steps ---
  int nSteps = 50;
  double step0 = 2.0 * 0.63 / myff.f22();
  std::cout << "Initial step: " << step0 << "\n";

  std::vector<double> vecStep;
  vecStep.reserve(nSteps);
  for (int idx = 0; idx < nSteps; ++idx) {
    vecStep.push_back(
        step0 / std::pow(10.0, 2.0 * idx / static_cast<double>(nSteps - 1)));
  }

  // --- 5. Normalization ---
  double normFactor = 1.0;
  double accelNorm = prem.LengthNorm() / (prem.TimeNorm() * prem.TimeNorm());

  if (params.output_type() == 0) {
    normFactor = prem.LengthNorm();
  } else if (params.output_type() == 1) {
    normFactor = prem.LengthNorm() / prem.TimeNorm();
  } else if (params.output_type() == 2) {
    normFactor = accelNorm;
  }
  double hannW = 0.2;

  // --- 6. Execute Iterative Convergence Test ---
  SPARSESPEC::SparseFSpec specSolver;
  std::vector<Eigen::MatrixXcd> vecFinalW;
  std::vector<Eigen::MatrixXd> vecFinalT;
  vecFinalW.reserve(nSteps);
  vecFinalT.reserve(nSteps);

  for (int idx = 0; idx < nSteps; ++idx) {
    maxstep = vecStep[idx];
    int nskip = 3 * static_cast<int>(
                        std::floor(maxstep / ((vecW[1] - vecW[0]) * 0.003))) +
                1;

    Full1D::SEM sem(prem, maxstep, nq, lval);
    std::cout << "\nDoing step: " << maxstep << ", nskip: " << nskip << "\n";

    MatrixC vecRaw = specSolver.spectra(myff, sem, prem, cmt, params, nskip);
    vecRaw *= normFactor;

    // Process responses
    auto vecR2TB = processfunctions::freq2time(vecRaw, myff);
    auto aFilt0 = processfunctions::fulltime2freq(vecR2TB, myff, 0.05);
    auto vecFiltT = processfunctions::filtfreq2time(aFilt0, myff, false);
    auto aFilt = processfunctions::fulltime2freq(vecFiltT, myff, hannW);

    vecFinalW.push_back(aFilt);
    vecFinalT.push_back(vecFiltT);
  }

  // --- 7. Error Calculation ---
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

  // "Exact" solution is the one with the smallest step
  auto idxBack = nSteps - 1;
  auto vecTEx = vecFinalT[idxBack];
  auto vecWEx = vecFinalW[idxBack];

  // Reference Norms Setup
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

    // L1 norms (scaled by exact solution's L-infinity norm)
    vecErrT[idx][0] = vecTErr.row(0).lpNorm<1>() * mvT1;
    vecErrT[idx][1] = vecTErr.row(1).lpNorm<1>() * mvT2;
    vecErrT[idx][2] = vecTErr.row(2).lpNorm<1>() * mvT3;
    vecErrW[idx][0] = vecWErr.row(0).lpNorm<1>() * mvW1;
    vecErrW[idx][1] = vecWErr.row(1).lpNorm<1>() * mvW2;
    vecErrW[idx][2] = vecWErr.row(2).lpNorm<1>() * mvW3;

    // L2 norms
    vecL2ErrT[idx][0] = vecTErr.row(0).norm() * mvT1L2;
    vecL2ErrT[idx][1] = vecTErr.row(1).norm() * mvT2L2;
    vecL2ErrT[idx][2] = vecTErr.row(2).norm() * mvT3L2;
    vecL2ErrW[idx][0] = vecWErr.row(0).norm() * mvW1L2;
    vecL2ErrW[idx][1] = vecWErr.row(1).norm() * mvW2L2;
    vecL2ErrW[idx][2] = vecWErr.row(2).norm() * mvW3L2;
  }

  // --- 8. Outputs ---
  double nval = 1.0 / prem.TimeNorm();

  // 8a. Output Frequency Series
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

  // 8b. Output Time Series
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

  // 8c. Output Error Matrix (L1-based)
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

  // 8d. Output L2 Error Matrix
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