#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#include "config.h"
#include "PaperExampleSupport.h"
#include <iostream>
#include <string>
#include <complex>
#include <cmath>

#include <PlanetaryModel/All>
#include <DSpecM1D/Timer>
#include <DSpecM1D/All>
#include <SpectraSolver/FF>

int
main() {
  using Complex = std::complex<double>;
  using MatrixC = Eigen::MatrixXcd;

  Timer timer1;

  // This example keeps the legacy explicit setup because it inserts a
  // source-time-function workflow between the raw solve and the final outputs.
  std::string paramPath =
      std::string(PROJECT_BUILD_DIR) + "data/params/ex2.txt";
  InputParameters params(paramPath);
  std::string earthModelPath =
      std::string(PROJECT_BUILD_DIR) + "data/" + params.earth_model();

  SRInfo srInfo(params);
  auto cmt = SourceInfo::EarthquakeCMT(params);

  prem_norm<double> normClass;
  auto prem = EarthModels::ModelInput(earthModelPath, normClass, "true");

  double dt = params.time_step_sec();
  double tout = params.t_out() / 60.0;
  double df0 = 1.0;
  double wtb = 0.05;
  double t1 = 0.0;
  double t2 = tout;
  int qex = 1;
  int nq = 6;

  // Build the frequency grid explicitly so the later STF processing can refer
  // to the same samples used by the forward solve.
  SpectraSolver::FreqFull myff(params.f1(), params.f2(), params.f11(),
                               params.f12(), params.f21(), params.f22(), dt,
                               tout, df0, wtb, t1, t2, qex, prem.TimeNorm());
  auto vecW = myff.w();

  // Solve for the raw DSpecM1D spectra.
  SPARSESPEC::SparseFSpec spec;
  timer1.start();
  MatrixC vecRaw = spec.spectra(myff, prem, cmt, params, nq, srInfo,
                                params.relative_error());
  timer1.stop("Total time for sparse frequency spectrum");

  // Convert the non-dimensional solver output into the requested SI quantity.
  double normFactor = PaperExamples::legacyNormFactor(params, prem);
  vecRaw *= normFactor;

  // Start from the unshifted model response before the STF is applied.
  double hannW = 0.2;
  auto vecR2TB = processfunctions::freq2time(vecRaw, myff);
  auto aFilt0 = processfunctions::fulltime2freq(vecR2TB, myff, 0.01);

  // Apply the STF in the frequency domain before the final filtered outputs.
  double hd = 1e-9;
  double decayT = 1.628;
  double kapVal = decayT / (std::sqrt(PI) * hd);
  double stTime = 2.0 / kapVal;
  const Complex imagUnit(0.0, 1.0);

  auto applyGaussianStf = [&](const Eigen::MatrixXcd &spectrum,
                              bool divideByOmegaSquared = false) {
    auto shifted = spectrum;
    const int startCol = divideByOmegaSquared ? 1 : 0;

    for (int idx = startCol; idx < shifted.cols(); ++idx) {
      auto wval = myff.w(idx);
      Complex stfFactor =
          std::exp(-imagUnit * wval * stTime) *
          std::exp(-(1.0 / (4.0 * PI)) * std::pow(wval / kapVal, 2.0));

      if (divideByOmegaSquared) {
        stfFactor *= -1.0 / (wval * wval);
      }

      shifted.col(idx) *= stfFactor;
    }

    return shifted;
  };

  std::cout << "Source time function parameters:\n"
            << "Half duration (s): " << hd << "\n"
            << "Kappa: " << kapVal << "\n"
            << "df: " << myff.df() << "\n"
            << "Source time (s): " << stTime * prem.TimeNorm() << "\n";

  // Apply the STF to the DSpecM1D response and then return to the filtered
  // time/frequency products used for the paper comparison.
  auto aFiltStf0 = applyGaussianStf(aFilt0);
  auto vecFiltT = processfunctions::filtfreq2time(aFiltStf0, myff, false);
  auto aFilt = processfunctions::fulltime2freq(vecFiltT, myff, hannW);

  // Load the external comparison traces that accompany the paper.
  auto filterOptions = PaperExamples::makeFilterOptions(hannW);
  std::string yspecPath = std::string(PROJECT_BUILD_DIR) +
                          "data/reference/yspec/yspec.out.bolivia.mf.1";
  std::string mineosBase = std::string(PROJECT_BUILD_DIR) +
                           "data/reference/mineos/bolivia/"
                           "Syndat2.2000160: 0:33:16.TLY.";
  auto referenceSeries = PaperExamples::loadReferenceTimeSeriesWithMineosBase(
      yspecPath, mineosBase, vecR2TB.cols());
  const auto &yspecTime = referenceSeries.yspecTime;
  const auto &mineosTime = referenceSeries.mineosTime;

  // YSpec uses the same STF shift as the DSpecM1D result.
  auto aYSpecStf0 =
      applyGaussianStf(processfunctions::fulltime2freq(yspecTime, myff, 0.01));
  auto vecFiltTYSpec = processfunctions::filtfreq2time(aYSpecStf0, myff, false);
  auto aFiltYSpec = processfunctions::fulltime2freq(vecFiltTYSpec, myff, hannW);

  auto aMineosStf0 = applyGaussianStf(
      processfunctions::fulltime2freq(mineosTime, myff, 0.01), true);
  aMineosStf0 *= prem.TimeNorm() * prem.TimeNorm();
  auto vecFiltTMineos =
      processfunctions::filtfreq2time(aMineosStf0, myff, false);
  auto aFiltMineos =
      processfunctions::fulltime2freq(vecFiltTMineos, myff, hannW);

  // SpecNM already provides time-domain traces, so it only needs the shared
  // filtering step.
  std::string specnmPath = std::string(PROJECT_BUILD_DIR) +
                           "data/reference/specnm/"
                           "seismogram_mf_traces_semicolon.txt";
  Eigen::MatrixXd specnmTime =
      DSpecM::loadSpecnmTimeSeries(specnmPath, vecFiltT.cols());

  auto filteredSpecnm = DSpecM::applyFilter(specnmTime, myff, filterOptions);
  auto &vecFiltTSpecnm = filteredSpecnm.timeSeries;
  auto &aFiltSpecnm = filteredSpecnm.frequencySeries;

  // Write the final comparison products consumed by the paper plotting scripts.
  std::string pathToFreqFile =
      std::string(PROJECT_BUILD_DIR) + "../plotting/outputs/ex2_w.out";
  double nval = 1.0 / prem.TimeNorm();
  PaperExamples::writeFourWayFrequencyComparison(
      pathToFreqFile, vecW, myff, nval, aFilt, aFiltYSpec, aFiltMineos,
      aFiltSpecnm);

  std::string pathToTimeFile =
      std::string(PROJECT_BUILD_DIR) + "../plotting/outputs/ex2_t.out";
  PaperExamples::writeFourWayTimeComparison(
      pathToTimeFile, myff, prem.TimeNorm(), params.t_out(), vecFiltT,
      vecFiltTYSpec, vecFiltTMineos, vecFiltTSpecnm,
      -stTime * prem.TimeNorm());

  return 0;
}
