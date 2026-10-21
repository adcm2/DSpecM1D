#include <gtest/gtest.h>
#include <cmath>
#include <tuple>
#include <utility>
#include <array>
#include <DSpecM1D/src/InputParametersNew.h>
#include "test_utils.h"

template <typename T> class SmallModelNorm {
public:
  T LengthNorm() const { return 2.0; }
  T MassNorm() const { return 8.0; }
  T TimeNorm() const { return 4.0; }
};

static std::string isotropicRepeatedRadiusDeck(double qmu = 0.0) {
  std::ostringstream out;
  out << "isotropic repeated radius\n0 0 0\n4 0 0\n"
      << "0 2 10 0 100 0\n"
      << "2 4 12 0 100 0\n"
      << "2 6 14 3 100 " << qmu << "\n"
      << "4 8 16 4 100 " << qmu << "\n";
  return out.str();
}

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

TEST(InputParametersNewTests, ModelInterpolantsSurviveOwnerCopyAndMove) {
  using Model = EarthModels::ModelInput<double>;
  const auto path = DSpecMTest::modelPath().string();
  const int layer = 1;

  auto [interpolant, copiedModel, movedModel, radius, referenceValue,
        referenceDerivative] = [&] {
    Model original(path);
    const double r =
        0.5 * (original.LowerRadius(layer) + original.UpperRadius(layer));
    const auto spline = original.Density(layer);
    const double value = spline.Evaluate(r);
    const double derivative = spline.Evaluate<1>(r);

    Model copied = original;
    Model moved = std::move(copied);
    Model assigned;
    assigned = original;
    Model moveAssigned;
    moveAssigned = std::move(assigned);
    return std::tuple{original.Density(layer), std::move(moved),
                      std::move(moveAssigned), r, value, derivative};
  }();

  for (const auto &spline : {interpolant, copiedModel.Density(layer),
                             movedModel.Density(layer)}) {
    EXPECT_DOUBLE_EQ(spline.Evaluate(radius), referenceValue);
    EXPECT_DOUBLE_EQ(spline.Evaluate<1>(radius), referenceDerivative);
  }
}

TEST(InputParametersNewTests, PlanetmodelFacadePreservesIsotropicEndpointsAndOwnership) {
  DSpecMTest::TempDir temp;
  const auto path = DSpecMTest::writeFile(temp.path() / "iso.deck",
                                          isotropicRepeatedRadiusDeck());
  using Model = EarthModels::ModelInput<double>;
  const auto [density, retainedVs, lower, upper, interface] = [&] {
    Model model(path.string(), SmallModelNorm<double>{});
    EXPECT_EQ(model.NumberOfLayers(), 2);
    EXPECT_TRUE(model.IsIsotropic());
    EXPECT_FALSE(model.IsSolid(0));
    EXPECT_TRUE(model.IsSolid(1));
    EXPECT_DOUBLE_EQ(model.LowerRadius(1), 1.0);
    EXPECT_DOUBLE_EQ(model.UpperRadius(0), 1.0);
    EXPECT_DOUBLE_EQ(model.OuterRadius(), 2.0);
    const auto rho = model.Density(1);
    const auto vs = model.VS(1);
    const double lo = model.LowerRadius(1);
    const double hi = model.UpperRadius(1);
    const double boundary = model.UpperRadius(0);
    EXPECT_DOUBLE_EQ(model.VPH(1)(lo), model.VPV(1)(lo));
    EXPECT_DOUBLE_EQ(model.VSH(1)(lo), model.VSV(1)(lo));
    EXPECT_DOUBLE_EQ(model.Eta(1)(lo), 1.0);
    EXPECT_DOUBLE_EQ(model.Density(0)(boundary), 4.0);
    EXPECT_DOUBLE_EQ(model.Density(1)(boundary), 6.0);
    EXPECT_THROW(rho(std::nextafter(lo, -INFINITY)), std::out_of_range);
    EXPECT_THROW(rho(std::nextafter(hi, INFINITY)), std::out_of_range);
    EXPECT_DOUBLE_EQ(rho(lo), 6.0);
    EXPECT_DOUBLE_EQ(rho(hi), 8.0);
    Model copied = model;
    Model moved = std::move(copied);
    return std::tuple{rho, vs, lo, hi, boundary};
  }();
  EXPECT_DOUBLE_EQ(density(lower), 6.0);
  EXPECT_DOUBLE_EQ(density(upper), 8.0);
  EXPECT_TRUE(std::isfinite(retainedVs(interface)));
}

TEST(InputParametersNewTests, PlanetmodelFacadeRejectsInvalidQAndHeaderShape) {
  DSpecMTest::TempDir temp;
  for (const auto [label, qmu] :
       std::array<std::pair<const char*, const char*>, 3>{{
           {"negative", "-1"}, {"nan", "nan"}, {"infinite", "inf"}}}) {
    const auto path = DSpecMTest::writeFile(
        temp.path() / (std::string(label) + ".deck"),
        isotropicRepeatedRadiusDeck(std::stod(qmu)));
    EXPECT_THROW((EarthModels::ModelInput<double>(path.string())),
                 std::invalid_argument);
  }
  const auto mismatch = DSpecMTest::writeFile(
      temp.path() / "mismatch.deck", "bad header\n1 0 0\n4 0 0\n"
                                     "0 2 10 0 100 0\n2 4 12 0 100 0\n"
                                     "2 6 14 0 100 0\n4 8 16 0 100 0\n");
  EXPECT_THROW((EarthModels::ModelInput<double>(mismatch.string())),
               std::invalid_argument);
}
