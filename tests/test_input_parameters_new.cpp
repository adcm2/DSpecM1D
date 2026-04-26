#include <gtest/gtest.h>
#include <cmath>
#include <DSpecM1D/src/InputParametersNew.h>
#include "test_utils.h"

TEST(InputParametersNewTests, ResolvesAbsoluteModelPathAndBuildsContext) {
  DSpecMTest::TempDir temp;
  auto path = DSpecMTest::writeFile(
      temp.path() / "params.txt",
      DSpecMTest::makeParameterText(DSpecMTest::modelPath().string()));

  InputParametersNew paramsNew(path.string());
  EXPECT_EQ(paramsNew.earthModelPath(), DSpecMTest::modelPath().string());
  EXPECT_EQ(paramsNew.nq(), 5);
  EXPECT_EQ(paramsNew.nskip(), 10);
  EXPECT_EQ(paramsNew.inputParameters().num_receivers(), 2);
}

TEST(InputParametersNewTests, NormFactorMatchesOutputTypeSelection) {
  DSpecMTest::TempDir temp;

  auto dispPath = DSpecMTest::writeFile(
      temp.path() / "disp.txt",
      DSpecMTest::makeParameterText(DSpecMTest::modelPath().string(), 0));
  auto velPath = DSpecMTest::writeFile(
      temp.path() / "vel.txt",
      DSpecMTest::makeParameterText(DSpecMTest::modelPath().string(), 1));
  auto accPath = DSpecMTest::writeFile(
      temp.path() / "acc.txt",
      DSpecMTest::makeParameterText(DSpecMTest::modelPath().string(), 2));

  InputParametersNew disp(dispPath.string());
  InputParametersNew vel(velPath.string());
  InputParametersNew acc(accPath.string());

  const auto lengthNorm = disp.earthModel().LengthNorm();
  const auto timeNorm = disp.earthModel().TimeNorm();

  EXPECT_NEAR(disp.normFactor(), lengthNorm, 1e-12);
  EXPECT_NEAR(vel.normFactor(), lengthNorm / timeNorm, 1e-12);
  EXPECT_NEAR(acc.normFactor(), lengthNorm / (timeNorm * timeNorm), 1e-9);
}

TEST(InputParametersNewTests, SetterClampsToMinimumPositiveIntegers) {
  DSpecMTest::TempDir temp;
  auto path = DSpecMTest::writeFile(
      temp.path() / "params.txt",
      DSpecMTest::makeParameterText(DSpecMTest::modelPath().string()));

  InputParametersNew paramsNew(path.string());
  paramsNew.setNq(0);
  paramsNew.setNskip(-1);

  EXPECT_EQ(paramsNew.nq(), 1);
  EXPECT_EQ(paramsNew.nskip(), 1);
}

TEST(InputParametersNewTests, PreservesAnisotropicColumnsFromModelInputFile) {
  DSpecMTest::TempDir temp;
  auto path = DSpecMTest::writeFile(
      temp.path() / "params.txt",
      DSpecMTest::makeParameterText(DSpecMTest::modelPath().string()));

  InputParametersNew paramsNew(path.string());
  const auto &model = paramsNew.earthModel();

  bool foundAnisotropicLayer = false;
  for (int layer = 0; layer < model.NumberOfLayers(); ++layer) {
    const double r = 0.5 * (model.LowerRadius(layer) + model.UpperRadius(layer));
    const double vpDiff = std::abs(model.VPH(layer)(r) - model.VPV(layer)(r));
    const double vsDiff = std::abs(model.VSH(layer)(r) - model.VSV(layer)(r));
    const double etaDiff = std::abs(model.Eta(layer)(r) - 1.0);
    if (vpDiff > 1e-12 || vsDiff > 1e-12 || etaDiff > 1e-12) {
      foundAnisotropicLayer = true;
      break;
    }
  }

  EXPECT_TRUE(foundAnisotropicLayer);
}
