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
#include <algorithm>

// Project-specific includes
#include <PlanetaryModel/All>
#include <DSpecM1D/Timer>
#include <DSpecM1D/All>
#include <SpectraSolver/FF>

// Global constants
// constexpr double PI = 3.1415926535897932;
// constexpr double TWO_PI = 2.0 * PI;
// constexpr std::complex<double> I(0.0, 1.0);   // Imaginary unit

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
  using Complex = std::complex<double>;
  using MatrixC = Eigen::MatrixXcd;

  Timer timer1;

  // --- 1. Read Inputs & Earth Model ---
  // get paths required for input parameters and Earth model
  std::string paramPath =
      std::string(PROJECT_BUILD_DIR) + "data/params/ex2.txt";
  InputParameters params(paramPath);
  std::string earthModelPath =
      std::string(PROJECT_BUILD_DIR) + "data/" + params.earth_model();

  SRInfo srInfo(params);
  auto cmt = SourceInfo::EarthquakeCMT(params);

  prem_norm<double> normClass;
  auto prem = EarthModels::ModelInput(earthModelPath, normClass, "true");

  // --- 2. Frequency Solver Parameters ---
  double dt = params.time_step_sec();
  double tout = params.t_out() / 60.0;
  double df0 = 1.0;
  double wtb = 0.05;
  double t1 = 0.0;
  double t2 = tout;
  int qex = 1;
  int nq = 6;

  SpectraSolver::FreqFull myff(params.f1(), params.f2(), params.f11(),
                               params.f12(), params.f21(), params.f22(), dt,
                               tout, df0, wtb, t1, t2, qex, prem.TimeNorm());
  auto vecW = myff.w();

  // --- 3. Compute Sparse Frequency Spectrum ---
  SPARSESPEC::SparseFSpec spec;

  timer1.start();
  MatrixC vecRaw = spec.spectra(myff, prem, cmt, params, nq, srInfo,
                                params.relative_error());
  timer1.stop("Total time for sparse frequency spectrum");

  // --- 4. Normalization ---
  double normFactor = 1.0;
  double accelNorm = prem.LengthNorm() / (prem.TimeNorm() * prem.TimeNorm());

  if (params.output_type() == 0) {
    normFactor = prem.LengthNorm();
  } else if (params.output_type() == 1) {
    normFactor = prem.LengthNorm() / prem.TimeNorm();
  } else if (params.output_type() == 2) {
    normFactor = accelNorm;
  }
  vecRaw *= normFactor;

  // --- 5. Base Responses (Pre-Filter) ---
  double hannW = 0.2;
  auto vecR2TB = processfunctions::freq2time(vecRaw, myff);
  auto aFilt0 = processfunctions::fulltime2freq(vecR2TB, myff, 0.01);

  // --- 6. Source Time Function (STF) Convolution ---

  double hd = 1e-9;        // Half duration (set to near zero for testing)
  double decayT = 1.628;   // Decay constant
  double kapVal = decayT / (std::sqrt(PI) * hd);
  double stTime = 2.0 / kapVal;

  std::cout << "Source time function parameters:\n"
            << "Half duration (s): " << hd << "\n"
            << "Kappa: " << kapVal << "\n"
            << "df: " << myff.df() << "\n"
            << "Source time (s): " << stTime * prem.TimeNorm() << "\n";

  auto aFiltStf0 = aFilt0;
  for (int idx = 0; idx < aFilt0.cols(); ++idx) {
    auto wval = myff.w(idx);
    Complex stfFactor =
        std::exp(-I * wval * stTime) *
        std::exp(-(1.0 / (4.0 * PI)) * std::pow(wval / kapVal, 2.0));
    aFiltStf0.col(idx) *= stfFactor;
  }

  auto vecFiltT = processfunctions::filtfreq2time(aFiltStf0, myff, false);
  auto aFilt = processfunctions::fulltime2freq(vecFiltT, myff, hannW);

  //////////////////////////////////////////////////////////////////////////////
  // --- 7. Signal Filtering ---
  DSpecM::FilterOptions filterOptions;
  filterOptions.preTaper = 0.05;
  filterOptions.finalTaper = hannW;
  filterOptions.passes = 1;
  filterOptions.enforceRealSignal = false;

  // --- 7. Read and Process YSpec Data ---
  std::string yspecPath = std::string(PROJECT_BUILD_DIR) +
                          "data/reference/yspec/yspec.out.bolivia.mf.1";
  YSPECREADER::DataColumns yspecData(yspecPath);

  // FIX: Added safe boundary checking here to prevent out-of-bounds crashes
  std::size_t maxColY = std::min(static_cast<std::size_t>(vecR2TB.cols()),
                                 yspecData.getColumn1().size());

  Eigen::MatrixXd yspecTime = Eigen::MatrixXd::Zero(3, vecR2TB.cols());
  for (std::size_t idx = 0; idx < maxColY; ++idx) {
    yspecTime(0, idx) = yspecData.getColumn2()[idx];
    yspecTime(1, idx) = yspecData.getColumn3()[idx];
    yspecTime(2, idx) = yspecData.getColumn4()[idx];
  }

  auto aFiltYSpec0 = processfunctions::fulltime2freq(yspecTime, myff, 0.01);
  auto aYSpecStf0 = aFiltYSpec0;

  for (int idx = 1; idx < aYSpecStf0.cols(); ++idx) {
    auto wval = myff.w(idx);
    Complex stfFactor =
        std::exp(-I * wval * stTime) *
        std::exp(-(1.0 / (4.0 * PI)) * std::pow(wval / kapVal, 2.0));
    aYSpecStf0.col(idx) *= stfFactor;
  }

  auto vecFiltTYSpec = processfunctions::filtfreq2time(aYSpecStf0, myff, false);
  auto aFiltYSpec = processfunctions::fulltime2freq(vecFiltTYSpec, myff, hannW);

  // --- 8. Read and Process MinEOS Data ---
  std::string mineosBase = std::string(PROJECT_BUILD_DIR) +
                           "data/reference/mineos/bolivia/"
                           "Syndat2.2000160: 0:33:16.TLY.";
  MINEOSREADER::DataColumns mineosDataZ(mineosBase + "LHZ.ASC");
  MINEOSREADER::DataColumns mineosDataN(mineosBase + "LHN.ASC");
  MINEOSREADER::DataColumns mineosDataE(mineosBase + "LHE.ASC");

  std::size_t maxCol = std::min(static_cast<std::size_t>(vecR2TB.cols()),
                                mineosDataZ.getColumn1().size());

  Eigen::MatrixXd mineosTime = Eigen::MatrixXd::Zero(3, vecR2TB.cols());
  for (std::size_t idx = 0; idx < maxCol; ++idx) {
    mineosTime(0, idx) = mineosDataZ.getColumn2()[idx] * 1e-9;
    mineosTime(1, idx) = mineosDataN.getColumn2()[idx] * 1e-9;
    mineosTime(2, idx) = mineosDataE.getColumn2()[idx] * 1e-9;
  }

  auto aMineos0 = processfunctions::fulltime2freq(mineosTime, myff, 0.01);
  auto aMineosStf0 = aMineos0;

  // Start at idx=1 to avoid division by zero (wval * wval)
  for (int idx = 1; idx < aMineos0.cols(); ++idx) {
    auto wval = myff.w(idx);
    Complex stfFactor =
        -1.0 / (wval * wval) * std::exp(-I * wval * stTime) *
        std::exp(-(1.0 / (4.0 * PI)) * std::pow(wval / kapVal, 2.0));
    aMineosStf0.col(idx) *= stfFactor;
  }

  aMineosStf0 *= prem.TimeNorm() * prem.TimeNorm();
  auto vecFiltTMineos =
      processfunctions::filtfreq2time(aMineosStf0, myff, false);
  auto aFiltMineos =
      processfunctions::fulltime2freq(vecFiltTMineos, myff, hannW);

  // --- 9. Read and Process SpecNM Data ---
  std::string specnmPath = std::string(PROJECT_BUILD_DIR) +
                           "data/reference/specnm/"
                           "seismogram_mf_traces_semicolon.txt";
  Eigen::MatrixXd specnmTime =
      DSpecM::loadSpecnmTimeSeries(specnmPath, vecFiltT.cols());

  auto filteredSpecnm = DSpecM::applyFilter(specnmTime, myff, filterOptions);
  auto &vecFiltTSpecnm = filteredSpecnm.timeSeries;
  auto &aFiltSpecnm = filteredSpecnm.frequencySeries;

  //////////////////////////////////////////////////////////////////////////////
  // --- 9. Output Frequency Spectrum ---
  std::string pathToFreqFile =
      std::string(PROJECT_BUILD_DIR) + "../plotting/outputs/ex2_w.out";
  std::ofstream freqFile(pathToFreqFile);
  if (!freqFile) {
    std::cerr << "Error: unable to open output file_w: " << pathToFreqFile
              << "\n";
    return 1;
  }

  double nval = 1.0 / prem.TimeNorm();
  freqFile << std::fixed << std::setprecision(16);

  for (std::size_t idx = 0; idx < myff.i2() + 100; ++idx) {
    freqFile << (vecW[idx] * nval * 1000.0 / TWO_PI) << ';'
             << aFilt(0, idx).real() << ';' << aFilt(0, idx).imag() << ';'
             << std::abs(aFilt(0, idx)) << ';' << aFilt(1, idx).real() << ';'
             << aFilt(1, idx).imag() << ';' << std::abs(aFilt(1, idx)) << ';'
             << aFilt(2, idx).real() << ';' << aFilt(2, idx).imag() << ';'
             << std::abs(aFilt(2, idx)) << ';' << aFiltYSpec(0, idx).real()
             << ';' << aFiltYSpec(0, idx).imag() << ';'
             << std::abs(aFiltYSpec(0, idx)) << ';' << aFiltYSpec(1, idx).real()
             << ';' << aFiltYSpec(1, idx).imag() << ';'
             << std::abs(aFiltYSpec(1, idx)) << ';' << aFiltYSpec(2, idx).real()
             << ';' << aFiltYSpec(2, idx).imag() << ';'
             << std::abs(aFiltYSpec(2, idx)) << ';'
             << aFiltMineos(0, idx).real() << ';' << aFiltMineos(0, idx).imag()
             << ';' << std::abs(aFiltMineos(0, idx)) << ';'
             << aFiltMineos(1, idx).real() << ';' << aFiltMineos(1, idx).imag()
             << ';' << std::abs(aFiltMineos(1, idx)) << ';'
             << aFiltMineos(2, idx).real() << ';' << aFiltMineos(2, idx).imag()
             << ';' << std::abs(aFiltMineos(2, idx)) << ";"
             << aFiltSpecnm(0, idx).real() << ";" << aFiltSpecnm(0, idx).imag()
             << ";" << std::abs(aFiltSpecnm(0, idx)) << ";"
             << aFiltSpecnm(0, idx).real() << ";" << aFiltSpecnm(1, idx).imag()
             << ";" << std::abs(aFiltSpecnm(1, idx)) << ";"
             << aFiltSpecnm(2, idx).real() << ";" << aFiltSpecnm(2, idx).imag()
             << ";" << std::abs(aFiltSpecnm(2, idx)) << '\n';
  }
  freqFile.close();

  // --- 10. Output Time Series ---
  std::string pathToTimeFile =
      std::string(PROJECT_BUILD_DIR) + "../plotting/outputs/ex2_t.out";
  std::ofstream timeFile(pathToTimeFile);
  if (!timeFile) {
    std::cerr << "Error: unable to open output file_t: " << pathToTimeFile
              << "\n";
    return 1;
  }

  timeFile << std::fixed << std::setprecision(16);

  for (std::size_t idx = 0;
       idx < static_cast<std::size_t>(vecFiltTMineos.cols()); ++idx) {
    auto tval = idx * myff.dt() * prem.TimeNorm();
    timeFile << (idx * myff.dt() - stTime) * prem.TimeNorm() << ';'
             << vecFiltT(0, idx) << ';' << vecFiltT(1, idx) << ';'
             << vecFiltT(2, idx) << ';' << vecFiltTYSpec(0, idx) << ';'
             << vecFiltTYSpec(1, idx) << ';' << vecFiltTYSpec(2, idx) << ';'
             << vecFiltTMineos(0, idx) << ';' << vecFiltTMineos(1, idx) << ';'
             << vecFiltTMineos(2, idx) << ";" << vecFiltTSpecnm(0, idx) << ";"
             << vecFiltTSpecnm(1, idx) << ";" << vecFiltTSpecnm(2, idx) << '\n';

    if (tval > t2 * 3600.0) {
      break;
    }
  }
  timeFile.close();

  return 0;
}
