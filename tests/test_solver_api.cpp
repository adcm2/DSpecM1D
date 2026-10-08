#include <gtest/gtest.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <tuple>
#include <vector>
#include <Eigen/SparseLU>
#include <DSpecM1D/src/InputParametersNew.h>
#include <DSpecM1D/src/FullSpec.h>
#include <DSpecM1D/src/SpectraRunContext.h>
#include <DSpecM1D/src/SEM/SEM.h>
#include "test_utils.h"

namespace {

InputParametersNew
makeTinyPreferredParams() {
  DSpecMTest::TempDir temp;
  DSpecMTest::ParameterOptions options;
  options.type = 1;
  options.lmax = 2;
  options.f1 = 0.1;
  options.f2 = 0.2;
  options.f11 = 0.1;
  options.f12 = 0.12;
  options.f21 = 0.18;
  options.f22 = 0.2;
  options.tOutMinutes = 1.0;
  options.timeStepSec = 10.0;
  options.numReceivers = 1;
  options.receivers = {{45.0, 90.0}};
  options.relativeError = 1e-3;

  const auto path = DSpecMTest::writeFile(
      temp.path() / "tiny_params.txt",
      DSpecMTest::makeParameterText(DSpecMTest::modelPath().string(), options));

  InputParametersNew paramsNew(path.string(), 3, 2, 0.2, 1.0, 0.05, 0.0, 1);
  return paramsNew;
}

InputParametersNew
makeCowlingSolveParams(bool attenuation = false, int lmin = 1, int lmax = 2,
                       double toutMinutes = 1.0, int nq = 3,
                       double maxstep = 0.2) {
  DSpecMTest::TempDir temp;
  DSpecMTest::ParameterOptions options;
  options.type = 3;
  options.attenuation = attenuation;
  options.lmin = lmin;
  options.lmax = lmax;
  options.f1 = 5.0;
  options.f2 = 80.0;
  options.f11 = 5.0;
  options.f12 = 5.0;
  options.f21 = 80.0;
  options.f22 = 80.0;
  options.tOutMinutes = toutMinutes;
  options.timeStepSec = 5.0;
  options.numReceivers = 1;
  options.receivers = {{45.0, 90.0}};

  const auto model = DSpecMTest::repoRoot() / "data" / "models" /
                     (attenuation ? "prem.200.no.txt"
                                  : "prem.200.noatten.txt");
  const auto path = DSpecMTest::writeFile(
      temp.path() / "cowling_params.txt",
      DSpecMTest::makeParameterText(model.string(), options));
  return InputParametersNew(path.string(), nq, 2, maxstep, 1.0, 0.05, 0.0,
                            1);
}

InputParametersNew makeMultiCowlingSolveParams(double fmin, double fmax,
                                               bool attenuation,
                                               double toutMinutes = 1.0,
                                               int lmax = 1) {
  DSpecMTest::TempDir temp;
  DSpecMTest::ParameterOptions options;
  options.type = 3;
  options.attenuation = attenuation;
  options.lmin = 1;
  options.lmax = lmax;
  options.f1 = fmin;
  options.f2 = fmax;
  options.f11 = fmin;
  options.f12 = fmin;
  options.f21 = fmax;
  options.f22 = fmax;
  options.tOutMinutes = toutMinutes;
  options.timeStepSec = 5.0;
  options.numReceivers = 1;
  options.receivers = {{45.0, 90.0}};
  options.relativeError = 1.0;

  const auto model = DSpecMTest::repoRoot() / "data" / "models" /
                     (attenuation ? "prem.200.no.txt"
                                  : "prem.200.noatten.txt");
  const auto path = DSpecMTest::writeFile(
      temp.path() / "multi_cowling_params.txt",
      DSpecMTest::makeParameterText(model.string(), options));
  return InputParametersNew(path.string(), 3, 2, 0.2, 1.0, 0.05, 0.0, 1);
}

}   // namespace

TEST(PreferredSolverApiTests, SpectraRunContextExposesWorkflowObjects) {
  auto paramsNew = makeTinyPreferredParams();

  SPARSESPEC::SpectraRunContext request(paramsNew.freqFull(), paramsNew.cmt(),
                                        paramsNew.inputParameters(),
                                        paramsNew.tref(), 0);

  EXPECT_EQ(&request.freqFull(), &paramsNew.freqFull());
  EXPECT_EQ(&request.cmt(), &paramsNew.cmt());
  EXPECT_EQ(&request.params(), &paramsNew.inputParameters());
  EXPECT_DOUBLE_EQ(request.tref(), paramsNew.tref());
  EXPECT_EQ(request.nskip(), 1);
}

TEST(PreferredSolverApiTests, SemConstructorFromInputParametersNewUsesSettings) {
  auto paramsNew = makeTinyPreferredParams();
  paramsNew.setNq(3);
  paramsNew.setMaxstep(0.2);

  Full1D::SEM sem(paramsNew);
  Full1D::SEM directSem(paramsNew.earthModel(), paramsNew.maxstep(),
                        paramsNew.nq(), paramsNew.inputParameters().lmax());

  EXPECT_EQ(sem.mesh().NN(), paramsNew.nq());
  EXPECT_EQ(sem.mesh().NE(), directSem.mesh().NE());
  EXPECT_EQ(sem.mesh().PR(), directSem.mesh().PR());
  EXPECT_LE(sem.el(), sem.eu());
}

TEST(PreferredSolverApiTests, PreferredSolverOverloadsReturnStableShapes) {
  auto paramsNew = makeTinyPreferredParams();
  paramsNew.setNq(3);
  paramsNew.setNskip(2);
  paramsNew.setMaxstep(0.2);

  SPARSESPEC::SparseFSpec solver;
  Full1D::SEM sem(paramsNew);

  const auto withOwnedSem = solver.spectra(paramsNew);
  const auto withSharedSem = solver.spectra(paramsNew, sem);
  const auto withSharedSemReversed = solver.spectra(sem, paramsNew);

  const auto expectedRows = 3 * paramsNew.inputParameters().num_receivers();
  const auto expectedCols =
      static_cast<Eigen::Index>(paramsNew.freqFull().w().size());

  ASSERT_EQ(withOwnedSem.rows(), expectedRows);
  ASSERT_EQ(withOwnedSem.cols(), expectedCols);
  ASSERT_EQ(withSharedSem.rows(), expectedRows);
  ASSERT_EQ(withSharedSem.cols(), expectedCols);
  ASSERT_EQ(withSharedSemReversed.rows(), expectedRows);
  ASSERT_EQ(withSharedSemReversed.cols(), expectedCols);

  EXPECT_TRUE(withOwnedSem.real().array().isFinite().all());
  EXPECT_TRUE(withOwnedSem.imag().array().isFinite().all());
  EXPECT_TRUE(withSharedSem.real().array().isFinite().all());
  EXPECT_TRUE(withSharedSem.imag().array().isFinite().all());
  EXPECT_TRUE(withSharedSemReversed.real().array().isFinite().all());
  EXPECT_TRUE(withSharedSemReversed.imag().array().isFinite().all());

  EXPECT_EQ(withSharedSem.rows(), withSharedSemReversed.rows());
  EXPECT_EQ(withSharedSem.cols(), withSharedSemReversed.cols());
  EXPECT_TRUE(withSharedSem.isApprox(withSharedSemReversed, 1e-12));
}

TEST(PreferredSolverApiTests, StandaloneCowlingSpheroidalSolveOnPrem) {
  auto paramsNew = makeCowlingSolveParams();
  Full1D::SEM sem(paramsNew);
  std::cout << "Cowling PREM mesh elements=" << sem.mesh().NE()
            << " nodes_per_element=" << sem.mesh().NN() << "\n";
  SPARSESPEC::SparseFSpec solver;
  SPARSESPEC::SpectraRunContext request(paramsNew.freqFull(), paramsNew.cmt(),
                                        paramsNew.inputParameters(),
                                        paramsNew.tref(), 2);

  const auto defaultContext = solver.spectra(request, sem);
  const auto defaultApi = solver.spectra(paramsNew, sem);
  // forceCowling=true forces every spheroidal bin through Cowling;
  // forceCowling=false follows the configured cutoff, with zero leaving the
  // full-gravity solve unchanged.
  const auto full = solver.spectra(request, sem, false);
  const auto cowling = solver.spectra(request, sem, true);
  const int begin = paramsNew.freqFull().i1();
  const int end = paramsNew.freqFull().i2();
  ASSERT_GT(end, begin);
  ASSERT_EQ(full.rows(), 3);
  ASSERT_EQ(cowling.rows(), full.rows());
  ASSERT_EQ(cowling.cols(), full.cols());
  EXPECT_TRUE(defaultContext.isApprox(defaultApi, 0.0));
  EXPECT_TRUE(defaultContext.isApprox(full, 0.0));
  EXPECT_TRUE(full.real().array().isFinite().all());
  EXPECT_TRUE(full.imag().array().isFinite().all());
  EXPECT_TRUE(cowling.real().array().isFinite().all());
  EXPECT_TRUE(cowling.imag().array().isFinite().all());
  EXPECT_GT(cowling.middleCols(begin, end - begin).norm(), 0.0);

  const auto startIndices = SpectralTools::allIndicesSph(
      sem, 100, paramsNew.freqFull(), sem.sourceElement(paramsNew.cmt()), 1,
      true);
  EXPECT_GT(*std::max_element(startIndices.begin(), startIndices.end()), 0);

  std::vector<double> relativeDifferences;
  std::vector<double> selectedFrequencies;
  for (int idx = begin; idx < end; ++idx) {
    const double frequency = paramsNew.freqFull().f(idx) * 1000.0 /
                             paramsNew.timeNorm();
    if (frequency >= 5.0 && frequency <= 80.0 &&
        (selectedFrequencies.empty() ||
         frequency >= selectedFrequencies.back() + 4.0)) {
      const double fullNorm = full.col(idx).norm();
      const double diff = (full.col(idx) - cowling.col(idx)).norm();
      if (fullNorm > 0.0) {
        selectedFrequencies.push_back(frequency);
        relativeDifferences.push_back(diff / fullNorm);
      }
    }
  }
  ASSERT_GE(selectedFrequencies.size(), 2u);
  EXPECT_LT(relativeDifferences.back(), relativeDifferences.front());
  for (std::size_t i = 0; i < selectedFrequencies.size(); ++i)
    std::cout << "Cowling PREM comparison " << selectedFrequencies[i]
              << " mHz relative_difference=" << relativeDifferences[i]
              << "\n";

  auto attenParams = makeCowlingSolveParams(true);
  ASSERT_EQ(attenParams.inputParameters().attenuation(), 1);
  Full1D::SEM attenSem(attenParams);
  SPARSESPEC::SpectraRunContext attenRequest(
      attenParams.freqFull(), attenParams.cmt(), attenParams.inputParameters(),
      attenParams.tref(), 2);
  const auto attenCowling = solver.spectra(attenRequest, attenSem, true);
  EXPECT_TRUE(attenCowling.real().array().isFinite().all());
  EXPECT_TRUE(attenCowling.imag().array().isFinite().all());
  EXPECT_GT(attenCowling.middleCols(attenParams.freqFull().i1(),
                                    attenParams.freqFull().i2() -
                                        attenParams.freqFull().i1())
                .norm(),
            0.0);

  auto reducedParams = makeCowlingSolveParams(false, 100, 100);
  Full1D::SEM reducedSem(reducedParams);
  SPARSESPEC::SpectraRunContext reducedRequest(
      reducedParams.freqFull(), reducedParams.cmt(),
      reducedParams.inputParameters(), reducedParams.tref(), 2);
  const auto reducedStarts = SpectralTools::allIndicesSph(
      reducedSem, 100, reducedParams.freqFull(),
      reducedSem.sourceElement(reducedParams.cmt()), 1, true);
  ASSERT_GT(*std::max_element(reducedStarts.begin(), reducedStarts.end()), 0);
  const auto reducedCowling = solver.spectra(reducedRequest, reducedSem, true);
  EXPECT_TRUE(reducedCowling.real().array().isFinite().all());
  EXPECT_TRUE(reducedCowling.imag().array().isFinite().all());
  EXPECT_GT(reducedCowling.middleCols(reducedParams.freqFull().i1(),
                                      reducedParams.freqFull().i2() -
                                          reducedParams.freqFull().i1())
                .norm(),
            0.0);
}

TEST(PreferredSolverApiTests, DISABLED_CowlingValidationCampaignOnPrem) {
  const std::vector<double> targetMhz{5.0, 10.0, 20.0, 40.0, 80.0};
  SPARSESPEC::SparseFSpec solver;

  std::cout << "COWLING_VALIDATION_CSV_BEGIN\n"
            << "target_mhz,actual_mhz,l,component,full_real,full_imag,"
               "cowling_real,cowling_imag,relative_complex,relative_amplitude,"
               "phase_difference_deg\n";
  std::vector<std::tuple<int, Eigen::MatrixXcd, Eigen::MatrixXcd>>
      campaignSpectra;
  auto l1Params = makeCowlingSolveParams(false, 1, 1, 5.0, 5, 0.005);
  for (const int degree : {1, 20, 100}) {
    auto paramsNew = makeCowlingSolveParams(false, degree, degree, 5.0, 5,
                                            0.005);
    Full1D::SEM sem(paramsNew);
    SPARSESPEC::SpectraRunContext request(
        paramsNew.freqFull(), paramsNew.cmt(), paramsNew.inputParameters(),
        paramsNew.tref(), 1);
    const auto full = solver.spectra(request, sem, false);
    const auto cowling = solver.spectra(request, sem, true);
    EXPECT_TRUE(full.real().array().isFinite().all());
    EXPECT_TRUE(full.imag().array().isFinite().all());
    EXPECT_TRUE(cowling.real().array().isFinite().all());
    EXPECT_TRUE(cowling.imag().array().isFinite().all());
    EXPECT_GT(full.middleCols(paramsNew.freqFull().i1(),
                              paramsNew.freqFull().i2() -
                                  paramsNew.freqFull().i1())
                  .norm(),
              0.0);
    EXPECT_GT(cowling.middleCols(paramsNew.freqFull().i1(),
                                 paramsNew.freqFull().i2() -
                                     paramsNew.freqFull().i1())
                  .norm(),
              0.0);
    campaignSpectra.emplace_back(degree, full, cowling);

    const auto &freq = paramsNew.freqFull();
    for (const double target : targetMhz) {
      int sample = freq.i1();
      double bestDistance = std::numeric_limits<double>::infinity();
      for (int idx = freq.i1(); idx < freq.i2(); ++idx) {
        const double actual = freq.f(idx) * 1000.0 / freq.timeNorm();
        const double distance = std::abs(actual - target);
        if (distance < bestDistance) {
          bestDistance = distance;
          sample = idx;
        }
      }
      const double actual = freq.f(sample) * 1000.0 / freq.timeNorm();
      for (int component = 0; component < full.rows(); ++component) {
        const auto f = full(component, sample);
        const auto c = cowling(component, sample);
        const double rowScale = std::max(
            full.row(component).cwiseAbs().maxCoeff(),
            cowling.row(component).cwiseAbs().maxCoeff());
        const double floor = 1e-12 * rowScale;
        const double denominator = std::max(std::abs(f), floor);
        const double relComplex = std::abs(c - f) / denominator;
        const double relAmplitude = (std::abs(c) - std::abs(f)) / denominator;
        const double phase = (std::abs(f) > floor && std::abs(c) > floor)
                                 ? std::arg(c * std::conj(f)) * 180.0 / M_PI
                                 : std::numeric_limits<double>::quiet_NaN();
        std::cout << std::setprecision(17) << target << ',' << actual << ','
                  << degree << ',' << component << ',' << f.real() << ','
                  << f.imag() << ',' << c.real() << ',' << c.imag() << ','
                  << relComplex << ',' << relAmplitude << ',';
        if (std::isfinite(phase))
          std::cout << phase;
        else
          std::cout << "NA";
        std::cout << '\n';
      }
    }
  }

  std::cout << "COWLING_RESOLUTION_CSV_BEGIN\n"
            << "target_mhz,coarse_mhz,fine_mhz,l,component,formulation,"
               "coarse_real,coarse_imag,fine_real,fine_imag,"
               "relative_complex_difference\n";
  for (const int degree : {1, 100}) {
    auto coarseParams = makeCowlingSolveParams(false, degree, degree, 5.0, 5,
                                               0.005);
    auto fineParams = makeCowlingSolveParams(false, degree, degree, 5.0, 5,
                                             0.0025);
    Full1D::SEM coarseSem(coarseParams);
    Full1D::SEM fineSem(fineParams);
    SPARSESPEC::SpectraRunContext coarseRequest(
        coarseParams.freqFull(), coarseParams.cmt(),
        coarseParams.inputParameters(), coarseParams.tref(), 1);
    SPARSESPEC::SpectraRunContext fineRequest(
        fineParams.freqFull(), fineParams.cmt(), fineParams.inputParameters(),
        fineParams.tref(), 1);
    const auto dofs = [](const Full1D::SEM &sem) {
      const int ne = sem.mesh().NE();
      const int nq = sem.mesh().NN();
      return std::pair{sem.ltgS(2, ne - 1, nq - 1) + 1,
                       sem.ltgSC(1, ne - 1, nq - 1) + 1};
    };
    const auto coarseDofs = dofs(coarseSem);
    const auto fineDofs = dofs(fineSem);
    std::cout << "COWLING_RESOLUTION_META,l=" << degree
              << ",coarse_NE=" << coarseSem.mesh().NE()
              << ",coarse_full_dof=" << coarseDofs.first
              << ",coarse_cowling_dof=" << coarseDofs.second
              << ",fine_NE=" << fineSem.mesh().NE()
              << ",fine_full_dof=" << fineDofs.first
              << ",fine_cowling_dof=" << fineDofs.second << '\n';
    const auto coarseFull = solver.spectra(coarseRequest, coarseSem, false);
    const auto coarseCowling = solver.spectra(coarseRequest, coarseSem, true);
    const auto fineFull = solver.spectra(fineRequest, fineSem, false);
    const auto fineCowling = solver.spectra(fineRequest, fineSem, true);
    for (const double target : {20.0, 80.0}) {
      const auto nearest = [target](const SpectraSolver::FreqFull &freq) {
        int sample = freq.i1();
        double bestDistance = std::numeric_limits<double>::infinity();
        for (int idx = freq.i1(); idx < freq.i2(); ++idx) {
          const double actual = freq.f(idx) * 1000.0 / freq.timeNorm();
          if (std::abs(actual - target) < bestDistance) {
            bestDistance = std::abs(actual - target);
            sample = idx;
          }
        }
        return sample;
      };
      const int coarseIndex = nearest(coarseParams.freqFull());
      const int fineIndex = nearest(fineParams.freqFull());
      const double coarseMhz = coarseParams.freqFull().f(coarseIndex) * 1000.0 /
                               coarseParams.freqFull().timeNorm();
      const double fineMhz = fineParams.freqFull().f(fineIndex) * 1000.0 /
                             fineParams.freqFull().timeNorm();
      for (int component = 0; component < coarseFull.rows(); ++component) {
        for (const auto &item : {
                 std::tuple{"full", &coarseFull, &fineFull},
                 std::tuple{"cowling", &coarseCowling, &fineCowling}}) {
          const auto *coarse = std::get<1>(item);
          const auto *fine = std::get<2>(item);
          const auto c = (*coarse)(component, coarseIndex);
          const auto f = (*fine)(component, fineIndex);
          const double rowScale = std::max(
              coarse->row(component).cwiseAbs().maxCoeff(),
              fine->row(component).cwiseAbs().maxCoeff());
          const double floor = 1e-12 * rowScale;
          const double rel = std::abs(f - c) / std::max(std::abs(f), floor);
          std::cout << std::setprecision(17) << target << ',' << coarseMhz
                    << ',' << fineMhz << ',' << degree << ',' << component
                    << ',' << std::get<0>(item) << ',' << c.real() << ','
                    << c.imag() << ',' << f.real() << ',' << f.imag() << ','
                    << rel << '\n';
        }
      }
    }
  }
  std::cout << "COWLING_RESOLUTION_CSV_END\n";

  auto &freq = l1Params.freqFull();
  int cutoffIndex = freq.i1();
  double cutoffDistance = std::numeric_limits<double>::infinity();
  for (int idx = freq.i1(); idx + 1 < freq.i2(); ++idx) {
    const double actual = freq.f(idx) * 1000.0 / freq.timeNorm();
    const double distance = std::abs(actual - 20.0);
    if (distance < cutoffDistance) {
      cutoffDistance = distance;
      cutoffIndex = idx;
    }
  }
  const double cutoffMhz = freq.f(cutoffIndex) * 1000.0 / freq.timeNorm();
  ASSERT_GT(cutoffIndex, freq.i1());
  ASSERT_LT(cutoffIndex + 1, freq.i2());
  const int lower = cutoffIndex - 1;
  const int adjacent = cutoffIndex + 1;
  std::cout << "COWLING_CUTOFF_CSV_BEGIN\n"
            << "l,lower_mhz,cutoff_mhz,upper_mhz,component,full_lower_real,"
               "full_lower_imag,full_at_cutoff_real,full_at_cutoff_imag,"
               "cowling_at_cutoff_real,cowling_at_cutoff_imag,"
               "cowling_upper_real,cowling_upper_imag,mixed_lower_real,"
               "mixed_lower_imag,mixed_cutoff_real,mixed_cutoff_imag,"
               "mixed_upper_real,mixed_upper_imag,full_natural_complex_jump,"
               "mixed_complex_jump,switch_excess_relative,"
               "switch_amplitude_excess_relative,switch_phase_excess_deg\n";
  for (const auto &[degree, full, cowling] : campaignSpectra) {
    auto mixedParams = makeCowlingSolveParams(false, degree, degree, 5.0, 5,
                                              0.005);
    mixedParams.setCowlingFrequencyMhz(cutoffMhz);
    Full1D::SEM mixedSem(mixedParams);
    SPARSESPEC::SpectraRunContext mixedRequest(
        mixedParams.freqFull(), mixedParams.cmt(),
        mixedParams.inputParameters(), mixedParams.tref(), 1);
    const auto mixed = solver.spectra(mixedRequest, mixedSem, false);
    EXPECT_TRUE(mixed.col(lower).isApprox(full.col(lower), 0.0));
    EXPECT_TRUE(mixed.col(cutoffIndex).isApprox(cowling.col(cutoffIndex), 0.0));
    EXPECT_TRUE(mixed.col(adjacent).isApprox(cowling.col(adjacent), 0.0));
    EXPECT_TRUE(mixed.real().array().isFinite().all());
    EXPECT_TRUE(mixed.imag().array().isFinite().all());
    for (int component = 0; component < full.rows(); ++component) {
      const double rowScale = std::max(
          full.row(component).cwiseAbs().maxCoeff(),
          cowling.row(component).cwiseAbs().maxCoeff());
      const double floor = 1e-12 * rowScale;
      const auto fl = full(component, lower);
      const auto fcut = full(component, cutoffIndex);
      const auto ccut = cowling(component, cutoffIndex);
      const auto ml = mixed(component, lower);
      const auto mcut = mixed(component, cutoffIndex);
      const auto mu = mixed(component, adjacent);
      const double scale = std::max({std::abs(fl), std::abs(fcut), floor});
      const double naturalJump = std::abs(fcut - fl) / scale;
      const double mixedJump = std::abs(mcut - ml) / scale;
      const double switchExcess = std::abs(ccut - fcut) / scale;
      const double ampExcess = (std::abs(ccut) - std::abs(fcut)) / scale;
      const double phaseExcess = (std::abs(ccut) > floor &&
                                   std::abs(fcut) > floor)
                                      ? std::arg(ccut * std::conj(fcut)) *
                                            180.0 / M_PI
                                      : std::numeric_limits<double>::quiet_NaN();
      const double lowerMhz = freq.f(lower) * 1000.0 / freq.timeNorm();
      const double atCutoffMhz = freq.f(cutoffIndex) * 1000.0 / freq.timeNorm();
      const double upperMhz = freq.f(adjacent) * 1000.0 / freq.timeNorm();
      std::cout << std::setprecision(17) << degree << ',' << lowerMhz << ','
                << atCutoffMhz << ',' << upperMhz << ',' << component << ','
                << fl.real() << ',' << fl.imag() << ',' << fcut.real() << ','
                << fcut.imag() << ',' << ccut.real() << ',' << ccut.imag()
                << ',' << mu.real() << ',' << mu.imag() << ',' << ml.real()
                << ',' << ml.imag() << ',' << mcut.real() << ',' << mcut.imag()
                << ',' << mu.real() << ',' << mu.imag() << ',' << naturalJump
                << ',' << mixedJump << ',' << switchExcess << ',' << ampExcess
                << ',';
      if (std::isfinite(phaseExcess))
        std::cout << phaseExcess;
      else
        std::cout << "NA";
      std::cout << '\n';
    }
  }
  std::cout << "COWLING_CUTOFF_CSV_END\nCOWLING_VALIDATION_CSV_END\n";
}

TEST(PreferredSolverApiTests, DISABLED_CowlingPerformanceOnPrem) {
  using Clock = std::chrono::steady_clock;
  using Complex = std::complex<double>;
  using SparseMatrixC = Eigen::SparseMatrix<Complex>;
  using Solver = Eigen::SparseLU<SparseMatrixC, Eigen::COLAMDOrdering<int>>;
  constexpr int degree = 20;
  constexpr int rhsCount = 4;
  constexpr int repeats = 5;
  auto paramsNew = makeCowlingSolveParams(false, degree, degree, 1.0, 5,
                                          0.005);
  std::unique_ptr<Full1D::SEM> semStorage;
  std::vector<double> semBuildMs;
  for (int i = 0; i < repeats + 1; ++i) {
    const auto start = Clock::now();
    auto nextSem = std::make_unique<Full1D::SEM>(paramsNew);
    const auto stop = Clock::now();
    if (i > 0)
      semBuildMs.push_back(
          std::chrono::duration<double, std::milli>(stop - start).count());
    semStorage = std::move(nextSem);
  }
  std::sort(semBuildMs.begin(), semBuildMs.end());
  std::cout << "COWLING_PERF_TIMING,shared_sem,mesh_and_both_layouts_build_ms,"
            << semBuildMs[semBuildMs.size() / 2] << '\n';
  Full1D::SEM &sem = *semStorage;
  SPARSESPEC::SparseFSpec solver;
  SPARSESPEC::SpectraRunContext request(
      paramsNew.freqFull(), paramsNew.cmt(), paramsNew.inputParameters(),
      paramsNew.tref(), 1);
  auto &freq = paramsNew.freqFull();
  const double targetMhz = 20.0;
  int representative = freq.i1();
  for (int idx = freq.i1() + 1; idx < freq.i2(); ++idx) {
    if (std::abs(freq.f(idx) * 1000.0 / freq.timeNorm() - targetMhz) <
        std::abs(freq.f(representative) * 1000.0 / freq.timeNorm() -
                 targetMhz))
      representative = idx;
  }
  const double actualMhz =
      freq.f(representative) * 1000.0 / freq.timeNorm();
  const int sourceElement = sem.sourceElement(paramsNew.cmt());
  const auto allStarts = SpectralTools::allIndicesSph(
      sem, degree, freq, sourceElement, 1, false);
  const auto allCowlingStarts = SpectralTools::allIndicesSph(
      sem, degree, freq, sourceElement, 1, true);
  const int representativeOffset = representative - freq.i1();
  ASSERT_GE(representativeOffset, 0);
  ASSERT_LT(representativeOffset, static_cast<int>(allStarts.size()));
  const int startElement = SpectralTools::startElementSph(
      sem, degree, freq.w()[representative], sourceElement);
  EXPECT_EQ(allStarts[representativeOffset], sem.ltgS(0, startElement, 0));
  EXPECT_EQ(allCowlingStarts[representativeOffset],
            sem.ltgSC(0, startElement, 0));
  const auto fullH = sem.hS(degree);
  const auto fullP = sem.pS(degree);
  const auto cowlingH = sem.hSC(degree);
  const auto cowlingP = sem.pSC(degree);
  const int fullDofs = fullH.rows();
  const int cowlingDofs = cowlingH.rows();
  ASSERT_GT(fullDofs, cowlingDofs);
  ASSERT_EQ(fullH.cols(), fullDofs);
  ASSERT_EQ(cowlingH.cols(), cowlingDofs);

  const auto matrixFor = [&](bool useCowling) {
    const int first = useCowling ? allCowlingStarts[representativeOffset]
                                 : allStarts[representativeOffset];
    const int fullDimension = useCowling ? cowlingH.rows() : fullH.rows();
    if (first >= fullDimension) {
      ADD_FAILURE() << "Representative truncation starts beyond matrix";
      return std::make_pair(0, SparseMatrixC{});
    }
    const int reducedDimension = fullDimension - first;
    const auto &h = useCowling ? cowlingH : fullH;
    const auto &p = useCowling ? cowlingP : fullP;
    const Complex w = freq.w()[representative] + Complex(0.0, -freq.ep());
    SparseMatrixC matrix = h.block(first, first, reducedDimension,
                                   reducedDimension)
                               .cast<Complex>();
    matrix -= w * w * p.block(first, first, reducedDimension,
                              reducedDimension)
                             .cast<Complex>();
    matrix.makeCompressed();
    return std::make_pair(first, std::move(matrix));
  };
  const auto fullSystem = matrixFor(false);
  const auto cowlingSystem = matrixFor(true);
  auto fullForce = sem.calculateForceAll(paramsNew.cmt(), degree, false);
  auto cowlingForce = sem.calculateForceAll(paramsNew.cmt(), degree, true);
  const auto fullRhs = fullForce.bottomRows(fullDofs - fullSystem.first);
  const auto cowlingRhs = cowlingForce.bottomRows(cowlingDofs -
                                                   cowlingSystem.first);
  ASSERT_EQ(fullRhs.cols(), rhsCount);
  ASSERT_EQ(cowlingRhs.cols(), rhsCount);

  const auto matrixStats = [](const SparseMatrixC &matrix) {
    int lower = 0;
    int upper = 0;
    for (int col = 0; col < matrix.outerSize(); ++col)
      for (SparseMatrixC::InnerIterator it(matrix, col); it; ++it) {
        lower = std::max(lower, static_cast<int>(it.row() - it.col()));
        upper = std::max(upper, static_cast<int>(it.col() - it.row()));
      }
    return std::make_tuple(lower, upper, matrix.nonZeros());
  };
  const auto printStructure = [&](const char *name, int dofs,
                                  const std::pair<int, SparseMatrixC> &system) {
    const auto [lower, upper, nnz] = matrixStats(system.second);
    std::cout << "COWLING_PERF_STRUCTURE," << name << ',' << dofs << ','
              << system.first << ',' << system.second.rows() << ',' << lower
              << ',' << upper << ',' << nnz << '\n';
  };
  std::cout << std::setprecision(17)
            << "COWLING_PERF_META,degree," << degree << '\n'
            << "COWLING_PERF_META,target_mhz," << targetMhz << '\n'
            << "COWLING_PERF_META,actual_mhz," << actualMhz << '\n'
            << "COWLING_PERF_META,spectrum_bin_count," << freq.i2() - freq.i1()
            << '\n'
            << "COWLING_PERF_META,spectrum_first_mhz,"
            << freq.f(freq.i1()) * 1000.0 / freq.timeNorm() << '\n'
            << "COWLING_PERF_META,spectrum_last_mhz,"
            << freq.f(freq.i2() - 1) * 1000.0 / freq.timeNorm() << '\n'
            << "COWLING_PERF_META,full_to_cowling_dof_ratio,"
            << static_cast<double>(fullDofs) / cowlingDofs << '\n';
  printStructure("full", fullDofs, fullSystem);
  printStructure("cowling", cowlingDofs, cowlingSystem);

  const auto benchmarkSystem = [&](const char *name,
                                   const std::pair<int, SparseMatrixC> &system,
                                   const Eigen::MatrixXcd &rhs) {
    std::vector<double> factorMs;
    std::vector<double> solveMs;
    Eigen::MatrixXcd solution;
    for (int i = 0; i < repeats + 1; ++i) {
      Solver localSolver;
      auto start = Clock::now();
      localSolver.compute(system.second);  // includes symbolic analysis
      auto stop = Clock::now();
      ASSERT_EQ(localSolver.info(), Eigen::Success);
      if (i > 0)
        factorMs.push_back(
            std::chrono::duration<double, std::milli>(stop - start).count());
      start = Clock::now();
      solution = localSolver.solve(rhs);
      stop = Clock::now();
      ASSERT_EQ(localSolver.info(), Eigen::Success);
      if (i > 0)
        solveMs.push_back(
            std::chrono::duration<double, std::milli>(stop - start).count());
      EXPECT_TRUE(solution.real().array().isFinite().all());
      EXPECT_TRUE(solution.imag().array().isFinite().all());
      const double relativeResidual =
          (system.second * solution - rhs).norm() /
          std::max(rhs.norm(), std::numeric_limits<double>::min());
      EXPECT_LT(relativeResidual, 1e-8);
      if (i == 0)
        continue;
    }
    const auto median = [](std::vector<double> values) {
      std::sort(values.begin(), values.end());
      return values[values.size() / 2];
    };
    std::cout << "COWLING_PERF_TIMING," << name << ",factor_ms_including_"
                 "symbolic_analysis," << median(factorMs) << '\n'
              << "COWLING_PERF_TIMING," << name << ",solve_4_rhs_ms,"
              << median(solveMs) << '\n';
  };
  benchmarkSystem("full", fullSystem, fullRhs);
  benchmarkSystem("cowling", cowlingSystem, cowlingRhs);

  // Prepared-SEM spectrum timings include degree/frequency matrix, source,
  // and receiver setup in spectra(), but exclude SEM/model construction.
  for (const bool useCowling : {false, true}) {
    std::vector<double> spectrumMs;
    Eigen::MatrixXcd spectrum;
    for (int i = 0; i < repeats + 1; ++i) {
      const auto start = Clock::now();
      spectrum = solver.spectra(request, sem, useCowling);
      const auto stop = Clock::now();
      ASSERT_TRUE(spectrum.real().array().isFinite().all());
      ASSERT_TRUE(spectrum.imag().array().isFinite().all());
      if (i > 0)
        spectrumMs.push_back(std::chrono::duration<double, std::milli>(stop -
                                                                      start)
                                 .count());
    }
    std::sort(spectrumMs.begin(), spectrumMs.end());
    std::cout << "COWLING_PERF_TIMING,"
              << (useCowling ? "cowling" : "full")
              << ",prepared_sem_spectrum_5_80_mhz_ms,"
              << spectrumMs[spectrumMs.size() / 2] << '\n';
  }
}

TEST(PreferredSolverApiTests, SingleSemCowlingCutoffSplitsFrequencyRegions) {
  for (const bool attenuation : {false, true}) {
    auto paramsNew = makeCowlingSolveParams(attenuation, 1, 1);
    Full1D::SEM sem(paramsNew);
    SPARSESPEC::SparseFSpec solver;
    auto &freq = paramsNew.freqFull();
    auto &params = paramsNew.inputParameters();
    SPARSESPEC::SpectraRunContext request(freq, paramsNew.cmt(), params,
                                          paramsNew.tref(), 1);
    const int begin = freq.i1();
    const int end = freq.i2();
    ASSERT_GT(end - begin, 3);
    const auto frequencyMhz = [&](int idx) {
      return freq.f(idx) * 1000.0 / freq.timeNorm();
    };

    paramsNew.setCowlingFrequencyMhz(0.0);
    const auto disabled = solver.spectra(request, sem);
    const auto full = solver.spectra(request, sem, false);
    const auto cowling = solver.spectra(request, sem, true);
    EXPECT_TRUE(disabled.isApprox(full, 0.0));

    paramsNew.setCowlingFrequencyMhz(frequencyMhz(end - 1) + 1.0);
    const auto aboveRange = solver.spectra(request, sem);
    EXPECT_TRUE(aboveRange.isApprox(disabled, 0.0));

    paramsNew.setCowlingFrequencyMhz(frequencyMhz(begin) - 1.0);
    const auto belowRange = solver.spectra(request, sem);
    EXPECT_TRUE(belowRange.isApprox(cowling, 0.0));

    const int boundary = begin + (end - begin) / 2;
    const double cutoffMhz = frequencyMhz(boundary);
    paramsNew.setCowlingFrequencyMhz(cutoffMhz);
    const auto mixed = solver.spectra(request, sem);
    EXPECT_TRUE(mixed.real().array().isFinite().all());
    EXPECT_TRUE(mixed.imag().array().isFinite().all());
    EXPECT_GT(mixed.middleCols(begin, end - begin).norm(), 0.0);
    EXPECT_DOUBLE_EQ(frequencyMhz(boundary), cutoffMhz);
    EXPECT_TRUE(mixed.col(boundary - 1).isApprox(full.col(boundary - 1), 0.0));
    EXPECT_TRUE(mixed.col(boundary).isApprox(cowling.col(boundary), 0.0));
    if (boundary + 1 < end)
      EXPECT_TRUE(mixed.col(boundary + 1).isApprox(cowling.col(boundary + 1),
                                                  0.0));
  }
}

TEST(PreferredSolverApiTests, MixedCutoffRestartsHighDegreeTruncationCadence) {
  auto paramsNew = makeCowlingSolveParams(false, 100, 100);
  Full1D::SEM sem(paramsNew);
  SPARSESPEC::SparseFSpec solver;
  auto &freq = paramsNew.freqFull();
  const int begin = freq.i1();
  const int end = freq.i2();
  const int source = sem.sourceElement(paramsNew.cmt());
  auto allW = freq.w();
  SPARSESPEC::SpectraRunContext fullReferenceRequest(
      freq, paramsNew.cmt(), paramsNew.inputParameters(), paramsNew.tref(), 1);
  const auto fullReference = solver.spectra(fullReferenceRequest, sem, false);
  const std::vector<int> nskipValues{1, 2, 3, end - begin + 1};
  for (const int nskip : nskipValues) {
    auto wholeBandFullIndices = SpectralTools::allIndicesSph(
        sem, 100, freq, source, nskip, false);
    int boundary = begin + 1;
    ASSERT_LT(boundary, end);
    paramsNew.setCowlingFrequencyMhz(freq.f(boundary) * 1000.0 /
                                    freq.timeNorm());
    std::vector<double> fullW(allW.begin() + begin, allW.begin() + boundary);
    std::vector<double> cowlingW(allW.begin() + boundary, allW.begin() + end);
    auto fullIndices = SpectralTools::allIndicesSph(
        sem, 100, fullW, source, nskip, false);
    auto cowlingIndices = SpectralTools::allIndicesSph(
        sem, 100, cowlingW, source, nskip, true);
    ASSERT_FALSE(fullIndices.empty());
    ASSERT_FALSE(cowlingIndices.empty());
    EXPECT_GT(*std::max_element(fullIndices.begin(), fullIndices.end()), 0);
    const int fullStartElement = SpectralTools::startElementSph(
        sem, 100, fullW.back(), source);
    const int cowlingStartElement = SpectralTools::startElementSph(
        sem, 100, cowlingW.back(), source);
    EXPECT_EQ(fullIndices.back(), sem.ltgS(0, fullStartElement, 0));
    EXPECT_EQ(cowlingIndices.back(), sem.ltgSC(0, cowlingStartElement, 0));
    if (nskip > end - begin)
      EXPECT_NE(fullIndices.back(), wholeBandFullIndices[boundary - 1 - begin]);

    SPARSESPEC::SpectraRunContext request(
        freq, paramsNew.cmt(), paramsNew.inputParameters(), paramsNew.tref(),
        nskip);
    const auto mixed = solver.spectra(request, sem);
    EXPECT_TRUE(mixed.real().array().isFinite().all());
    EXPECT_TRUE(mixed.imag().array().isFinite().all());
    EXPECT_GT(mixed.middleCols(begin, end - begin).norm(), 0.0);
    if (nskip > end - begin)
      EXPECT_TRUE(mixed.col(boundary - 1).isApprox(
          fullReference.col(boundary - 1), 0.0));
  }
}

TEST(PreferredSolverApiTests, MultiSemCowlingCutoffUsesChunkRegions) {
  SPARSESPEC::SparseFSpec solver;
  for (const bool attenuation : {false, true}) {
    // A 2–8 mHz band forms one multi-SEM chunk. Both solver paths use the
    // same explicitly capped 0.05 mesh step for the numerical comparison.
    auto paramsNew = makeMultiCowlingSolveParams(2.0, 8.0, attenuation);
    auto &freq = paramsNew.freqFull();
    auto &params = paramsNew.inputParameters();
    Full1D::SEM singleSem(paramsNew.earthModel(), 0.05, paramsNew.nq(),
                          params.lmax());
    SPARSESPEC::SpectraRunContext request(freq, paramsNew.cmt(), params,
                                          paramsNew.tref(), 1);
    const auto frequencyMhz = [&](int idx) {
      return freq.f(idx) * 1000.0 / freq.timeNorm();
    };
    const int begin = freq.i1();
    const int end = freq.i2();
    ASSERT_GE(end - begin, 2);

    paramsNew.setCowlingFrequencyMhz(0.0);
    const auto multiDisabled = solver.spectra(
        freq, paramsNew.earthModel(), paramsNew.cmt(), params, paramsNew.nq(),
        paramsNew.srInfo(), params.relative_error());
    const auto singleFull = solver.spectra(request, singleSem, false);
    EXPECT_TRUE(multiDisabled.isApprox(singleFull, 1e-10));

    paramsNew.setCowlingFrequencyMhz(frequencyMhz(end - 1) + 1.0);
    const auto multiAbove = solver.spectra(
        freq, paramsNew.earthModel(), paramsNew.cmt(), params, paramsNew.nq(),
        paramsNew.srInfo(), params.relative_error());
    EXPECT_TRUE(multiAbove.isApprox(multiDisabled, 0.0));

    paramsNew.setCowlingFrequencyMhz(frequencyMhz(begin) - 1.0);
    const auto multiBelow = solver.spectra(
        freq, paramsNew.earthModel(), paramsNew.cmt(), params, paramsNew.nq(),
        paramsNew.srInfo(), params.relative_error());
    const auto singleCowling = solver.spectra(request, singleSem, true);
    EXPECT_TRUE(multiBelow.isApprox(singleCowling, 1e-10));

    const int boundary = begin + (end - begin) / 2;
    paramsNew.setCowlingFrequencyMhz(frequencyMhz(boundary));
    const auto multiMixed = solver.spectra(
        freq, paramsNew.earthModel(), paramsNew.cmt(), params, paramsNew.nq(),
        paramsNew.srInfo(), params.relative_error());
    paramsNew.setCowlingFrequencyMhz(frequencyMhz(boundary));
    const auto singleMixed = solver.spectra(request, singleSem, false);
    EXPECT_TRUE(multiMixed.isApprox(singleMixed, 1e-10));
    EXPECT_TRUE(multiMixed.col(boundary - 1).isApprox(
        singleFull.col(boundary - 1), 1e-10));
    EXPECT_TRUE(multiMixed.col(boundary).isApprox(
        singleCowling.col(boundary), 1e-10));
    if (boundary + 1 < end)
      EXPECT_TRUE(multiMixed.col(boundary + 1).isApprox(
          singleCowling.col(boundary + 1), 1e-10));
  }

  // The 5–35 mHz band forms three chunks: one stays full-gravity, one
  // crosses the cutoff, and the last is entirely Cowling.
  auto paramsNew = makeMultiCowlingSolveParams(5.0, 35.0, false, 60.0);
  auto &freq = paramsNew.freqFull();
  auto &params = paramsNew.inputParameters();
  Full1D::SEM singleSem(paramsNew.earthModel(), 0.05, paramsNew.nq(),
                        params.lmax());
  SPARSESPEC::SpectraRunContext request(freq, paramsNew.cmt(), params,
                                        paramsNew.tref(), 1);
  const auto frequencyMhz = [&](int idx) {
    return freq.f(idx) * 1000.0 / freq.timeNorm();
  };
  const int begin = freq.i1();
  const int end = freq.i2();
  const int derivedNskip = std::max(1, (end - begin) / 20);
  const int expectedChunks = std::max(
      1, static_cast<int>(std::floor((freq.f22() - freq.f11()) / 10.0) + 1.0));
  ASSERT_GT(derivedNskip, 1);
  ASSERT_EQ(expectedChunks, 3);
  std::cout << "Multi-SEM cutoff fixture bins=" << end - begin
            << " chunks=" << expectedChunks
            << " derived_nskip=" << derivedNskip << "\n";
  const int boundary = begin + (end - begin) / 2;
  const double cutoffMhz = frequencyMhz(boundary);
  paramsNew.setCowlingFrequencyMhz(cutoffMhz);
  const auto multiMixed = solver.spectra(
      freq, paramsNew.earthModel(), paramsNew.cmt(), params, paramsNew.nq(),
      paramsNew.srInfo(), params.relative_error());
  const auto singleMixed = solver.spectra(request, singleSem, false);
  EXPECT_TRUE(multiMixed.real().array().isFinite().all());
  EXPECT_TRUE(multiMixed.imag().array().isFinite().all());
  EXPECT_TRUE(multiMixed.isApprox(singleMixed, 1e-10));
}

TEST(PreferredSolverApiTests, MultiSemChunkAccumulationMatchesSingleSemAcrossDegrees) {
  SPARSESPEC::SparseFSpec solver;
  for (const bool attenuation : {false, true}) {
    // Four degrees exercise the parallel reduction when OMP_NUM_THREADS > 1.
    // The three chunks use the same capped mesh as this single-SEM reference.
    auto paramsNew = makeMultiCowlingSolveParams(5.0, 35.0, attenuation, 5.0, 4);
    auto &freq = paramsNew.freqFull();
    auto &params = paramsNew.inputParameters();
    Full1D::SEM singleSem(paramsNew.earthModel(), 0.05, paramsNew.nq(),
                         params.lmax());
    const int derivedNskip = std::max(1, (freq.i2() - freq.i1()) / 20);
    SPARSESPEC::SpectraRunContext request(freq, paramsNew.cmt(), params,
                                        paramsNew.tref(), derivedNskip);
    const auto frequencyMhz = [&](int idx) {
      return freq.f(idx) * 1000.0 / freq.timeNorm();
    };
    const std::vector<double> cutoffs{
        0.0, frequencyMhz(freq.i1()) / 2.0,
        frequencyMhz(freq.i1() + (freq.i2() - freq.i1()) / 2)};
    for (const double cutoff : cutoffs) {
      SCOPED_TRACE(::testing::Message() << "attenuation=" << attenuation
                                       << " cutoff_mhz=" << cutoff);
      paramsNew.setCowlingFrequencyMhz(cutoff);
      const auto multi = solver.spectra(paramsNew);
      const auto single = solver.spectra(request, singleSem);
      EXPECT_TRUE(multi.real().array().isFinite().all());
      EXPECT_TRUE(multi.imag().array().isFinite().all());
      EXPECT_GT(multi.norm(), 0.0);
      EXPECT_TRUE(multi.isApprox(single, 1e-10));
      EXPECT_TRUE(multi.leftCols(freq.i1()).isZero(0.0));
      EXPECT_TRUE(multi.rightCols(multi.cols() - freq.i2()).isZero(0.0));
    }
  }
}

TEST(PreferredSolverApiTests, LegacyMultiSemOverloadReturnsFiniteOutput) {
  auto paramsNew = makeTinyPreferredParams();
  paramsNew.setNq(3);
  paramsNew.setMaxstep(0.2);

  SPARSESPEC::SparseFSpec solver;
  auto &params = paramsNew.inputParameters();

  const auto legacy = solver.spectra(paramsNew.freqFull(), paramsNew.earthModel(),
                                     paramsNew.cmt(), params, paramsNew.nq(),
                                     paramsNew.srInfo(),
                                     params.relative_error());

  EXPECT_EQ(legacy.rows(), 3 * params.num_receivers());
  EXPECT_EQ(legacy.cols(),
            static_cast<Eigen::Index>(paramsNew.freqFull().w().size()));
  EXPECT_TRUE(legacy.real().array().isFinite().all());
  EXPECT_TRUE(legacy.imag().array().isFinite().all());
}

TEST(PreferredSolverApiTests, LegacySingleSemOverloadReturnsFiniteOutput) {
  auto paramsNew = makeTinyPreferredParams();
  paramsNew.setNq(3);
  paramsNew.setNskip(2);
  paramsNew.setMaxstep(0.2);

  SPARSESPEC::SparseFSpec solver;
  Full1D::SEM sem(paramsNew);

  const auto legacy = solver.spectra(paramsNew.freqFull(), sem,
                                     paramsNew.earthModel(), paramsNew.cmt(),
                                     paramsNew.inputParameters(),
                                     paramsNew.nskip());

  EXPECT_EQ(legacy.rows(), 3 * paramsNew.inputParameters().num_receivers());
  EXPECT_EQ(legacy.cols(),
            static_cast<Eigen::Index>(paramsNew.freqFull().w().size()));
  EXPECT_TRUE(legacy.real().array().isFinite().all());
  EXPECT_TRUE(legacy.imag().array().isFinite().all());
}
