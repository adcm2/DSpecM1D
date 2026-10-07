#include <gtest/gtest.h>
#include <algorithm>
#include <iostream>
#include <vector>
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
makeCowlingSolveParams(bool attenuation = false, int lmin = 1, int lmax = 2) {
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
  options.tOutMinutes = 1.0;
  options.timeStepSec = 5.0;
  options.numReceivers = 1;
  options.receivers = {{45.0, 90.0}};

  const auto model = DSpecMTest::repoRoot() / "data" / "models" /
                     (attenuation ? "prem.200.no.txt"
                                  : "prem.200.noatten.txt");
  const auto path = DSpecMTest::writeFile(
      temp.path() / "cowling_params.txt",
      DSpecMTest::makeParameterText(model.string(), options));
  return InputParametersNew(path.string(), 3, 2, 0.2, 1.0, 0.05, 0.0, 1);
}

InputParametersNew makeMultiCowlingSolveParams(double fmin, double fmax,
                                               bool attenuation,
                                               double toutMinutes = 1.0) {
  DSpecMTest::TempDir temp;
  DSpecMTest::ParameterOptions options;
  options.type = 3;
  options.attenuation = attenuation;
  options.lmin = 1;
  options.lmax = 1;
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
