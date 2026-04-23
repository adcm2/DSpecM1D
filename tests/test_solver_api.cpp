#include <gtest/gtest.h>
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
