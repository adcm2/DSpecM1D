#include <gtest/gtest.h>
#include <PlanetaryModel/All>
#include <DSpecM1D/src/SEM/SEM.h>
#include <DSpecM1D/src/SourceInfo.h>
#include <DSpecM1D/src/InputParametersNew.h>
#include "test_utils.h"

namespace {

InputParametersNew
makeTinySemParams() {
  DSpecMTest::TempDir temp;
  DSpecMTest::ParameterOptions options;
  options.type = 4;
  options.lmax = 4;
  options.f1 = 0.1;
  options.f2 = 0.3;
  options.f11 = 0.1;
  options.f12 = 0.12;
  options.f21 = 0.25;
  options.f22 = 0.3;
  options.tOutMinutes = 1.0;
  options.timeStepSec = 10.0;
  options.numReceivers = 1;
  options.receivers = {{45.0, 90.0}};
  options.relativeError = 1e-3;

  const auto path = DSpecMTest::writeFile(
      temp.path() / "tiny_sem_params.txt",
      DSpecMTest::makeParameterText(DSpecMTest::modelPath().string(), options));

  return InputParametersNew(path.string(), 4, 2, 0.2, 1.0, 0.05, 0.0, 1);
}

}   // namespace

TEST(SEMComponentTests, LocalToGlobalMapsAreMonotonicOnMinimalModel) {
  prem_norm<double> norm;
  auto model =
      EarthModels::ModelInput(DSpecMTest::modelPath().string(), norm, "true");
  Full1D::SEM sem(model, 0.05, 4, 4);

  EXPECT_LT(sem.ltgS(0, 0, 0), sem.ltgS(1, 0, 0));
  EXPECT_LT(sem.ltgS(1, 0, 0), sem.ltgS(2, 0, 0));
  EXPECT_LT(sem.ltgR(0, 0, 0), sem.ltgR(1, 0, 0));
  EXPECT_LT(sem.ltgS(2, 0, 0), sem.ltgS(0, 0, 1));
  EXPECT_LT(sem.ltgR(1, 0, 0), sem.ltgR(0, 0, 1));
  EXPECT_LE(sem.el(), sem.eu());
}

TEST(SEMComponentTests, ReceiverAndSourceElementsStayWithinMeshBounds) {
  auto paramsNew = makeTinySemParams();
  Full1D::SEM sem(paramsNew);

  auto &params = paramsNew.inputParameters();
  auto &cmt = paramsNew.cmt();

  const auto receiverElems = sem.receiverElements(params);
  ASSERT_FALSE(receiverElems.empty());
  for (int idx : receiverElems) {
    EXPECT_GE(idx, 0);
    EXPECT_LT(idx, sem.mesh().NE());
  }

  const auto sourceElem = sem.sourceElement(cmt);
  EXPECT_GE(sourceElem, 0);
  EXPECT_LT(sourceElem, sem.mesh().NE());
  EXPECT_LE(receiverElems.front(), receiverElems.back());
}

TEST(SEMComponentTests, SystemMatricesExposeConsistentShapes) {
  auto paramsNew = makeTinySemParams();
  Full1D::SEM sem(paramsNew);

  const int idxl = 2;
  const auto hS = sem.hS(idxl);
  const auto pS = sem.pS(idxl);
  const auto hTk = sem.hTk(idxl);
  const auto pTk = sem.pTk(idxl);
  const auto hR = sem.hR();
  const auto pR = sem.pR();

  ASSERT_EQ(hS.rows(), hS.cols());
  ASSERT_EQ(pS.rows(), pS.cols());
  ASSERT_EQ(hTk.rows(), hTk.cols());
  ASSERT_EQ(pTk.rows(), pTk.cols());
  ASSERT_EQ(hR.rows(), hR.cols());
  ASSERT_EQ(pR.rows(), pR.cols());

  EXPECT_EQ(hS.rows(), pS.rows());
  EXPECT_EQ(hTk.rows(), pTk.rows());
  EXPECT_EQ(hR.rows(), pR.rows());
  EXPECT_GT(hS.nonZeros(), 0);
  EXPECT_GT(hTk.nonZeros(), 0);
  EXPECT_GT(hR.nonZeros(), 0);
  EXPECT_TRUE(std::isfinite(sem.meshModel().Density(0, 0)));
  EXPECT_TRUE(std::isfinite(sem.meshModel().Gravity(0, 0)));
}

TEST(SEMComponentTests, ReceiverVectorsHaveExpectedDimensions) {
  auto paramsNew = makeTinySemParams();
  Full1D::SEM sem(paramsNew);

  auto &params = paramsNew.inputParameters();
  const int idxl = 2;
  const auto receiverElems = sem.receiverElements(params);

  const auto radial = sem.rvZR(params, 0);
  const auto reducedRadial = sem.rvRedZR(params);
  const auto baseFull = sem.rvBaseFull(params, idxl);
  const auto baseFullT = sem.rvBaseFullT(params, idxl);

  const auto fullRadialRows = sem.ltgR(1, sem.mesh().NE() - 1, sem.mesh().NN() - 1) + 1;
  const auto reducedRadialRows =
      sem.ltgR(1, receiverElems.back(), sem.mesh().NN() - 1) -
      sem.ltgR(0, receiverElems.front(), 0) + 1;
  const auto baseFullCols =
      sem.ltgS(1, receiverElems.back(), sem.mesh().NN() - 1) -
      sem.ltgS(0, receiverElems.front(), 0) + 1;
  const auto baseFullTCols =
      sem.ltgT(receiverElems.back(), sem.mesh().NN() - 1) -
      sem.ltgT(receiverElems.front(), 0) + 1;

  ASSERT_EQ(radial.rows(), fullRadialRows);
  ASSERT_EQ(radial.cols(), 1);
  ASSERT_EQ(reducedRadial.rows(), reducedRadialRows);
  ASSERT_EQ(reducedRadial.cols(), 1);
  ASSERT_EQ(baseFull.rows(), 3 * params.num_receivers());
  ASSERT_EQ(baseFull.cols(), baseFullCols);
  ASSERT_EQ(baseFullT.rows(), 3 * params.num_receivers());
  ASSERT_EQ(baseFullT.cols(), baseFullTCols);

  EXPECT_GT(radial.cwiseAbs().maxCoeff(), 0.0);
  EXPECT_GT(reducedRadial.cwiseAbs().maxCoeff(), 0.0);
  EXPECT_GT(baseFull.cwiseAbs().maxCoeff(), 0.0);
  EXPECT_GT(baseFullT.cwiseAbs().maxCoeff(), 0.0);
}
