#include <gtest/gtest.h>
#include <array>
#include <set>
#include <DSpecM1D/ModelInput>
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

template <class Model>
void
expectCowlingTopology(const Full1D::SEM &sem, const Model &model) {
  const int ne = sem.mesh().NE();
  const int nn = sem.mesh().NN();
  std::vector<std::vector<std::array<int, 2>>> expected(
      ne, std::vector<std::array<int, 2>>(nn));
  int next = 0;

  // Independently number the quotient mesh: U shares every element endpoint;
  // V shares an endpoint only when the material type is unchanged.
  for (int e = 0; e < ne; ++e) {
    const bool fluid = model.IsFluid(sem.mesh().LayerNumber(e));
    for (int n = 0; n < nn; ++n) {
      if (n == 0 && e > 0) {
        expected[e][n][0] = expected[e - 1][nn - 1][0];
        const bool previousFluid =
            model.IsFluid(sem.mesh().LayerNumber(e - 1));
        expected[e][n][1] = fluid == previousFluid
                                ? expected[e - 1][nn - 1][1]
                                : next++;
      } else {
        expected[e][n][0] = next++;
        expected[e][n][1] = next++;
      }
    }
  }

  const int expectedDofs = 2 * (ne * (nn - 1) + 1) +
                           static_cast<int>(sem.mesh().FS_Boundaries().size());
  EXPECT_EQ(next, expectedDofs);

  std::set<int> seen;
  for (int e = 0; e < ne; ++e) {
    for (int n = 0; n < nn; ++n) {
      for (int field = 0; field < 2; ++field) {
        const int actual = sem.ltgSC(field, e, n);
        EXPECT_EQ(actual, expected[e][n][field]);
        EXPECT_GE(actual, 0);
        EXPECT_LT(actual, expectedDofs);
        seen.insert(actual);
      }
    }
  }
  EXPECT_EQ(seen.size(), static_cast<std::size_t>(expectedDofs));
  for (int dof = 0; dof < expectedDofs; ++dof)
    EXPECT_TRUE(seen.count(dof));

  for (int e = 1; e < ne; ++e) {
    const bool fluid = model.IsFluid(sem.mesh().LayerNumber(e));
    const bool previousFluid = model.IsFluid(sem.mesh().LayerNumber(e - 1));
    EXPECT_EQ(sem.ltgSC(0, e, 0), sem.ltgSC(0, e - 1, nn - 1));
    if (fluid == previousFluid)
      EXPECT_EQ(sem.ltgSC(1, e, 0), sem.ltgSC(1, e - 1, nn - 1));
    else
      EXPECT_NE(sem.ltgSC(1, e, 0), sem.ltgSC(1, e - 1, nn - 1));
  }
}

}   // namespace

TEST(SEMComponentTests, LocalToGlobalMapsAreMonotonicOnMinimalModel) {
  prem_norm<double> norm;
  auto model = EarthModels::ModelInput(DSpecMTest::modelPath().string(), norm);
  Full1D::SEM sem(model, 0.05, 4, 4);

  EXPECT_LT(sem.ltgS(0, 0, 0), sem.ltgS(1, 0, 0));
  EXPECT_LT(sem.ltgS(1, 0, 0), sem.ltgS(2, 0, 0));
  EXPECT_LT(sem.ltgR(0, 0, 0), sem.ltgR(1, 0, 0));
  EXPECT_LT(sem.ltgS(2, 0, 0), sem.ltgS(0, 0, 1));
  EXPECT_LT(sem.ltgR(1, 0, 0), sem.ltgR(0, 0, 1));
  EXPECT_LE(sem.el(), sem.eu());
}

TEST(SEMComponentTests, CowlingMapHasCompactSolidFluidSolidTopology) {
  prem_norm<double> norm;
  const auto path = DSpecMTest::repoRoot() / "love_numbers" / "tests" / "data" /
                    "internal_fluid_sfs.txt";
  auto model = EarthModels::ModelInput(path.string(), norm);
  Full1D::SEM sem(model, 0.25, 4, 4);

  ASSERT_EQ(sem.mesh().FS_Boundaries().size(), 2u);
  expectCowlingTopology(sem, model);
}

TEST(SEMComponentTests, CowlingMapHasExpectedAllSolidTopology) {
  prem_norm<double> norm;
  const auto path = DSpecMTest::repoRoot() / "love_numbers" / "tests" / "data" /
                    "all_solid_isotropic.txt";
  auto model = EarthModels::ModelInput(path.string(), norm);
  Full1D::SEM sem(model, 0.3, 4, 4);

  ASSERT_TRUE(sem.mesh().FS_Boundaries().empty());
  expectCowlingTopology(sem, model);
}

TEST(SEMComponentTests, CowlingMapSharesSolidSolidMaterialBoundaries) {
  prem_norm<double> norm;
  auto model = EarthModels::ModelInput(DSpecMTest::modelPath().string(), norm);
  Full1D::SEM sem(model, 0.05, 4, 4);

  int solidSolidBoundaries = 0;
  for (int e = 1; e < sem.mesh().NE(); ++e) {
    const auto previousLayer = sem.mesh().LayerNumber(e - 1);
    const auto currentLayer = sem.mesh().LayerNumber(e);
    if (previousLayer != currentLayer && !model.IsFluid(previousLayer) &&
        !model.IsFluid(currentLayer))
      ++solidSolidBoundaries;
  }
  ASSERT_GT(solidSolidBoundaries, 0);
  expectCowlingTopology(sem, model);
}

TEST(SEMComponentTests, CowlingMatricesHaveCompactSymmetricMassAndK2Layout) {
  prem_norm<double> norm;
  const auto path = DSpecMTest::repoRoot() / "love_numbers" / "tests" / "data" /
                    "all_solid_isotropic.txt";
  auto model = EarthModels::ModelInput(path.string(), norm);
  Full1D::SEM sem(model, 0.3, 4, 4);

  const auto stiffness = sem.hSC(2);
  const auto attenuated = sem.hSCa(2);
  const auto mass = sem.pSC(2);
  const auto full = sem.hS(2);
  const auto fullMass = sem.pS(2);
  const auto nCowling = sem.ltgSC(1, sem.mesh().NE() - 1,
                                  sem.mesh().NN() - 1) + 1;

  EXPECT_EQ(stiffness.rows(), nCowling);
  EXPECT_EQ(stiffness.cols(), nCowling);
  EXPECT_EQ(mass.rows(), nCowling);
  EXPECT_EQ(attenuated.rows(), nCowling);
  EXPECT_LT(stiffness.rows(), full.rows());
  EXPECT_LT(mass.rows(), fullMass.rows());
  EXPECT_GT(stiffness.nonZeros(), 0);
  EXPECT_LT(stiffness.nonZeros(), full.nonZeros());
  EXPECT_EQ(attenuated.nonZeros(), stiffness.nonZeros());
  EXPECT_EQ(mass.nonZeros(), nCowling);

  const auto symmetric = [](const Eigen::SparseMatrix<double> &matrix) {
    return (matrix - Eigen::SparseMatrix<double>(matrix.transpose())).norm();
  };
  EXPECT_LT(symmetric(stiffness), 1e-12 * stiffness.norm());
  EXPECT_LT(symmetric(attenuated), 1e-12 * attenuated.norm());
  EXPECT_LT(symmetric(mass), 1e-12 * mass.norm());
  int positiveMassEntries = 0;
  for (int dof = 0; dof < mass.rows(); ++dof) {
    EXPECT_GE(mass.coeff(dof, dof), 0.0);
    positiveMassEntries += mass.coeff(dof, dof) > 0.0;
  }
  EXPECT_GT(positiveMassEntries, 0);
  EXPECT_EQ(mass.coeff(sem.ltgSC(0, 0, 0), sem.ltgSC(0, 0, 0)), 0.0);

  const auto h0 = sem.hSC(0);
  const auto h1 = sem.hSC(1);
  const auto h2 = sem.hSC(2);
  const auto h3 = sem.hSC(3);
  const Eigen::SparseMatrix<double> k2Prediction = 5.0 * h0 - 9.0 * h1 + 5.0 * h2;
  EXPECT_LT((h3 - k2Prediction).norm(), 1e-12 * h3.norm());
}

TEST(SEMComponentTests, CowlingStiffnessMatchesDirectWeakFormWithBackgroundGravity) {
  prem_norm<double> norm;
  const auto path = DSpecMTest::repoRoot() / "love_numbers" / "tests" / "data" /
                    "all_solid_isotropic.txt";
  auto model = EarthModels::ModelInput(path.string(), norm);
  Full1D::SEM sem(model, 0.3, 4, 4);

  const int element = sem.mesh().NE() / 2;
  const int node = 1;
  const int n = sem.mesh().NN();
  const auto q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(n);
  auto lagrange = Interpolation::LagrangePolynomial(q.Points().begin(),
                                                     q.Points().end());
  const double width = sem.mesh().EW(element);
  std::vector<double> testDerivative(n);
  for (int k = 0; k < n; ++k)
    testDerivative[k] = 2.0 / width * lagrange.Derivative(node, q.X(k));

  double expected = 0.0;
  const auto radius = sem.mesh().NodeRadius(element, node);
  const auto rho = sem.meshModel().Density(element, node);
  const auto gravity = sem.meshModel().Gravity(element, node);
  expected += 4.0 * width / 2.0 * q.W(node) *
              (sem.meshModel().A(element, node) -
               sem.meshModel().N(element, node) - rho * gravity * radius);
  expected += 2.0 * width * q.W(node) * sem.meshModel().F(element, node) *
              radius * testDerivative[node];
  for (int k = 0; k < n; ++k) {
    const auto rk = sem.mesh().NodeRadius(element, k);
    expected += width / 2.0 * q.W(k) * sem.meshModel().C(element, k) * rk * rk *
                testDerivative[k] * testDerivative[k];
  }

  const auto dof = sem.ltgSC(0, element, node);
  EXPECT_NEAR(sem.hSC(0).coeff(dof, dof), expected,
              1e-12 * std::max(1.0, std::abs(expected)));

  const double expectedMass = width / 2.0 * q.W(node) * rho * radius * radius;
  const auto vDof = sem.ltgSC(1, element, node);
  EXPECT_NEAR(sem.pSC(2).coeff(dof, dof), expectedMass,
              1e-12 * std::max(1.0, std::abs(expectedMass)));
  EXPECT_NEAR(sem.pSC(2).coeff(vDof, vDof), 6.0 * expectedMass,
              1e-12 * std::max(1.0, std::abs(6.0 * expectedMass)));

  double expectedAttenuation =
      4.0 * width / 2.0 * q.W(node) *
      (sem.meshModel().A_atten(element, node) -
       sem.meshModel().N_atten(element, node));
  expectedAttenuation += 2.0 * width * q.W(node) *
      sem.meshModel().F_atten(element, node) * radius * testDerivative[node];
  for (int k = 0; k < n; ++k) {
    const auto rk = sem.mesh().NodeRadius(element, k);
    expectedAttenuation += width / 2.0 * q.W(k) *
        sem.meshModel().C_atten(element, k) * rk * rk *
        testDerivative[k] * testDerivative[k];
  }
  ASSERT_GT(std::abs(expectedAttenuation), 0.0);
  EXPECT_NEAR(sem.hSCa(0).coeff(dof, dof), expectedAttenuation,
              1e-12 * std::max(1.0, std::abs(expectedAttenuation)));

  const double k2 = 2.0;
  const double expectedUV = k2 * (
      width / 2.0 * q.W(node) *
          (rho * gravity * radius - sem.meshModel().L(element, node) -
           2.0 * (sem.meshModel().A(element, node) -
                  sem.meshModel().N(element, node))) +
      width / 2.0 * q.W(node) * sem.meshModel().L(element, node) * radius *
          testDerivative[node] -
      width / 2.0 * q.W(node) * sem.meshModel().F(element, node) * radius *
          testDerivative[node]);
  ASSERT_GT(std::abs(rho * gravity * radius), 0.0);
  EXPECT_NEAR(sem.hSC(1).coeff(dof, vDof), expectedUV,
              1e-12 * std::max(1.0, std::abs(expectedUV)));

  const double selfGravity = 4.0 * width / 2.0 * q.W(node) *
      EIGEN_PI * 6.67230e-11 / model.GravitationalConstant() * rho * rho *
      radius * radius;
  const auto fullDof = sem.ltgS(0, element, node);
  EXPECT_NEAR(sem.hS(0).coeff(fullDof, fullDof) - expected, selfGravity,
              1e-12 * std::max(1.0, std::abs(selfGravity)));
  EXPECT_GT(std::abs(4.0 * width / 2.0 * q.W(node) * rho * gravity * radius),
            0.0);
}

TEST(SEMComponentTests, CowlingMatricesCoverLayeredFluidSolidMeshWithoutPotential) {
  prem_norm<double> norm;
  const auto path = DSpecMTest::repoRoot() / "love_numbers" / "tests" / "data" /
                    "internal_fluid_sfs.txt";
  auto model = EarthModels::ModelInput(path.string(), norm);
  Full1D::SEM sem(model, 0.25, 4, 4);
  const auto stiffness = sem.hSC(2);
  const auto mass = sem.pSC(2);
  const auto nCowling = sem.ltgSC(1, sem.mesh().NE() - 1,
                                  sem.mesh().NN() - 1) + 1;

  EXPECT_EQ(stiffness.rows(), nCowling);
  EXPECT_EQ(mass.rows(), nCowling);
  EXPECT_EQ(mass.nonZeros(), nCowling);
  EXPECT_LT(stiffness.rows(), sem.hS(2).rows());
  EXPECT_LT((stiffness - Eigen::SparseMatrix<double>(stiffness.transpose())).norm(),
            1e-12 * stiffness.norm());
  int positiveMassEntries = 0;
  for (int dof = 0; dof < mass.rows(); ++dof) {
    EXPECT_GE(mass.coeff(dof, dof), 0.0);
    positiveMassEntries += mass.coeff(dof, dof) > 0.0;
  }
  EXPECT_GT(positiveMassEntries, 0);
  EXPECT_EQ(mass.coeff(sem.ltgSC(0, 0, 0), sem.ltgSC(0, 0, 0)), 0.0);
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

TEST(SEMComponentTests, CowlingForcesAndReceiversReusePhysicalUVEntries) {
  auto paramsNew = makeTinySemParams();
  Full1D::SEM sem(paramsNew);
  auto &params = paramsNew.inputParameters();
  auto &cmt = paramsNew.cmt();
  const int idxl = 2;
  const int nq = sem.mesh().NN();

  const auto forceFull = sem.calculateForceAll(cmt, idxl);
  const auto forceCowling = sem.calculateForceAll(cmt, idxl, true);
  ASSERT_EQ(forceFull.cols(), forceCowling.cols());
  EXPECT_EQ(forceFull.rows(),
            sem.ltgS(2, sem.mesh().NE() - 1, nq - 1) + 1);
  EXPECT_EQ(forceCowling.rows(),
            sem.ltgSC(1, sem.mesh().NE() - 1, nq - 1) + 1);
  EXPECT_LT(forceCowling.rows(), forceFull.rows());
  for (int e = 0; e < sem.mesh().NE(); ++e) {
    for (int n = 0; n < nq; ++n) {
      for (int field = 0; field < 2; ++field) {
        const int fullIndex = sem.ltgS(field, e, n);
        const int cowlingIndex = sem.ltgSC(field, e, n);
        EXPECT_EQ((forceFull.row(fullIndex) -
                   forceCowling.row(cowlingIndex)).norm(),
                  0.0);
      }
    }
  }

  const auto receiverElems = sem.receiverElements(params);
  const int fullLow = sem.ltgS(0, receiverElems.front(), 0);
  const int cowlingLow = sem.ltgSC(0, receiverElems.front(), 0);
  const auto receiverFull = sem.rvBaseFull(params, idxl);
  const auto receiverCowling = sem.rvBaseFull(params, idxl, true);
  const int fullHigh = sem.ltgS(1, receiverElems.back(), nq - 1);
  const int cowlingHigh = sem.ltgSC(1, receiverElems.back(), nq - 1);
  ASSERT_EQ(receiverFull.cols(), fullHigh - fullLow + 1);
  ASSERT_EQ(receiverCowling.cols(), cowlingHigh - cowlingLow + 1);
  EXPECT_LT(receiverCowling.cols(), receiverFull.cols());

  Eigen::VectorXcd solutionFull = Eigen::VectorXcd::Constant(
      receiverFull.cols(), std::complex<double>(37.0, -4.0));
  Eigen::VectorXcd solutionCowling(receiverCowling.cols());
  for (int e = receiverElems.front(); e <= receiverElems.back(); ++e) {
    for (int n = 0; n < nq; ++n) {
      for (int field = 0; field < 2; ++field) {
        const std::complex<double> value(
            1.0 + sem.mesh().NodeRadius(e, n) /
                      sem.mesh().EUR(sem.mesh().NE() - 1),
            0.1 * field);
        const int fullIndex = sem.ltgS(field, e, n) - fullLow;
        const int cowlingIndex = sem.ltgSC(field, e, n) - cowlingLow;
        EXPECT_EQ((receiverFull.col(fullIndex) -
                   receiverCowling.col(cowlingIndex)).norm(),
                  0.0);
        solutionFull(fullIndex) = value;
        solutionCowling(cowlingIndex) = value;
      }
    }
  }
  EXPECT_TRUE((receiverFull * solutionFull)
                  .isApprox(receiverCowling * solutionCowling, 1e-14));
}
