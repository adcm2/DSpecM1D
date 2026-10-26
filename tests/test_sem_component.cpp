#include <gtest/gtest.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <set>
#include <vector>
#include <DSpecM1D/ModelInput>
#include <DSpecM1D/src/SEM/SEM.h>
#include <DSpecM1D/src/StartElement.h>
#include <DSpecM1D/src/SourceInfo.h>
#include <DSpecM1D/src/InputParametersNew.h>
#include <GSHTrans/Core>
#include <GaussQuad/All>
#include "test_utils.h"

namespace {

InputParametersNew
makeTinySemParams(double receiverDepth = 0.0,
                  const std::filesystem::path &earthModelPath =
                      DSpecMTest::modelPath()) {
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
  options.receiverDepth = receiverDepth;

  const auto path = DSpecMTest::writeFile(
      temp.path() / "tiny_sem_params.txt",
      DSpecMTest::makeParameterText(earthModelPath.string(), options));

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

TEST(DependencyApiTests, GaussLobattoIntegratesItsExactPolynomialRange) {
  const auto quadrature =
      GaussQuad::GaussLobattoLegendreQuadrature1D<double>(7);
  for (int degree = 0; degree <= 11; ++degree) {
    double integral = 0.0;
    for (int point = 0; point < quadrature.N(); ++point)
      integral += quadrature.W(point) * std::pow(quadrature.X(point), degree);
    const double exact = degree % 2 == 0 ? 2.0 / (degree + 1) : 0.0;
    EXPECT_NEAR(integral, exact, 2e-14) << "degree " << degree;
  }
}

TEST(DependencyApiTests, LagrangeBasisPreservesValuesAndDerivatives) {
  const std::vector<double> nodes{-1.0, -0.2, 0.5, 1.0};
  const auto basis = Interpolation::LagrangeBasis(nodes);
  const auto polynomial = [](double x) { return x * x * x - 0.5 * x + 0.25; };
  for (const double x : {-1.0, -0.4, -0.2, 0.1, 0.5, 0.8, 1.0}) {
    double value = 0.0;
    double derivative = 0.0;
    double partition = 0.0;
    double derivativeSum = 0.0;
    for (std::size_t i = 0; i < nodes.size(); ++i) {
      const double li = basis.Evaluate(i, x);
      const double dli = basis.Evaluate<1>(i, x);
      value += polynomial(nodes[i]) * li;
      derivative += polynomial(nodes[i]) * dli;
      partition += li;
      derivativeSum += dli;
    }
    EXPECT_NEAR(value, polynomial(x), 3e-14);
    EXPECT_NEAR(derivative, 3.0 * x * x - 0.5, 3e-14);
    EXPECT_NEAR(partition, 1.0, 3e-14);
    EXPECT_NEAR(derivativeSum, 0.0, 3e-14);
  }
}

TEST(DependencyApiTests, WignerKeepsSignedOrdersPolesAndHighDegreeFinite) {
  using Wigner =
      GSHTrans::Wigner<double, GSHTrans::All, GSHTrans::All,
                       GSHTrans::Multiple>;
  const std::array<double, 2> angles{0.0, 0.73};
  const Wigner table{750, 750, 1, angles};

  const double y00 = 1.0 / std::sqrt(4.0 * std::acos(-1.0));
  EXPECT_DOUBLE_EQ((table[0, 0][0, 0]), y00);
  const auto pole = table[0, 0];
  for (const int order : {-1, 1})
    EXPECT_DOUBLE_EQ((pole[1, order]), 0.0);
  EXPECT_DOUBLE_EQ((table[1, 0][1, 1]), std::sqrt(3.0) * y00);
  EXPECT_DOUBLE_EQ((table[1, 0][1, -1]), 0.0);

  const auto highDegree = table[1, 1];
  for (const int order : {-750, -1, 0, 1, 750})
    EXPECT_TRUE(std::isfinite(highDegree[750, order]));

  const GSHTrans::Wigner<double, GSHTrans::All, GSHTrans::All,
                         GSHTrans::Single>
      radialDegreeZero{0, 0, 0, 0.0};
  EXPECT_DOUBLE_EQ((radialDegreeZero[0][0, 0]), y00);
}

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
  auto lagrange = Interpolation::LagrangeBasis(q.Points());
  const double width = sem.mesh().EW(element);
  std::vector<double> testDerivative(n);
  for (int k = 0; k < n; ++k)
    testDerivative[k] = 2.0 / width * lagrange.Evaluate<1>(node, q.X(k));

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


TEST(SEMComponentTests, ExactOceanBottomReceiverUsesSolidSide) {
  auto paramsNew = makeTinySemParams(3.0);
  Full1D::SEM sem(paramsNew);
  auto &params = paramsNew.inputParameters();
  const double radius = sem.mesh().PR() - 1000.0 * params.receiver_depth() /
                                             paramsNew.earthModel().LengthNorm();

  std::vector<int> expectedSolid;
  std::vector<int> expectedFluid;
  bool foundSolid = false;
  bool foundFluid = false;
  for (int e = 0; e < sem.mesh().NE(); ++e) {
    if (sem.mesh().ELR(e) <= radius && sem.mesh().EUR(e) >= radius) {
      const bool fluid = paramsNew.earthModel().IsFluid(
          sem.mesh().LayerNumber(e));
      foundFluid = foundFluid || fluid;
      foundSolid = foundSolid || !fluid;
      if (fluid)
        expectedFluid.push_back(e);
      else
        expectedSolid.push_back(e);
    }
  }
  ASSERT_TRUE(foundSolid);
  ASSERT_TRUE(foundFluid);
  ASSERT_EQ(expectedFluid.size(), 1);
  const auto selected = sem.receiverElements(params);
  EXPECT_EQ(selected, expectedSolid);
  ASSERT_EQ(selected.size(), 1);

  const int solidElement = expectedSolid.front();
  const int fluidElement = expectedFluid.front();
  const auto interpolateElementField = [&](int element, double intercept,
                                           double slope) {
    std::vector<double> elementNodes(sem.mesh().NN());
    for (int n = 0; n < sem.mesh().NN(); ++n)
      elementNodes[n] = sem.mesh().NodeRadius(element, n);
    const auto elementBasis = Interpolation::LagrangeBasis(elementNodes);
    double result = 0.0;
    for (int n = 0; n < sem.mesh().NN(); ++n)
      result += elementBasis.Evaluate(n, radius) *
                (intercept + slope * elementNodes[n]);
    return result;
  };
  const double solidFieldAtReceiver =
      interpolateElementField(solidElement, 4.0, 2.0);
  const double fluidFieldAtReceiver =
      interpolateElementField(fluidElement, 80.0, -3.0);
  ASSERT_GT(std::abs(solidFieldAtReceiver - fluidFieldAtReceiver), 1.0);

  const int l = 2;
  const auto receiver = sem.rvBaseFull(params, l);
  const auto receiverCowling = sem.rvBaseFull(params, l, true);
  const int low = sem.ltgS(0, selected.front(), 0);
  const int cowlingLow = sem.ltgSC(0, selected.front(), 0);
  ASSERT_EQ(receiverCowling.cols(),
            sem.ltgSC(1, selected.back(), sem.mesh().NN() - 1) - cowlingLow + 1);
  ASSERT_EQ(receiver.cols(),
            sem.ltgS(1, selected.back(), sem.mesh().NN() - 1) - low + 1);
  for (int n = 0; n < sem.mesh().NN(); ++n) {
    const int fullV = sem.ltgS(1, selected.front(), n) - low;
    const int cowlingV = sem.ltgSC(1, selected.front(), n) - cowlingLow;
    for (int row = 1; row < 3; ++row)
      EXPECT_EQ(receiver(row, fullV), receiverCowling(row, cowlingV));
  }
  // The surface radial value is reproduced from an independently evaluated
  // solid-element Lagrange interpolant for a manufactured affine field.
  std::vector<double> nodes(sem.mesh().NN());
  for (int n = 0; n < sem.mesh().NN(); ++n)
    nodes[n] = sem.mesh().NodeRadius(selected.front(), n);
  auto solidBasis = Interpolation::LagrangeBasis(nodes);
  Eigen::VectorXcd manufacturedSolid = Eigen::VectorXcd::Zero(receiver.cols());
  for (int n = 0; n < sem.mesh().NN(); ++n) {
    const int col = sem.ltgS(0, selected.front(), n) - low;
    manufacturedSolid(col) = 4.0 + 2.0 * nodes[n];
  }
  const auto radialValue = receiver.row(0).dot(manufacturedSolid);
  double independentSolidValue = 0.0;
  for (int n = 0; n < sem.mesh().NN(); ++n)
    independentSolidValue += solidBasis.Evaluate(n, radius) *
                             (4.0 + 2.0 * nodes[n]);
  EXPECT_NEAR(radialValue.real(), independentSolidValue, 1e-12);
  EXPECT_GT(std::abs(radialValue.real() - fluidFieldAtReceiver), 1.0);
  EXPECT_NEAR(radialValue.imag(), 0.0, 1e-12);

  const auto toroidal = sem.rvThetaT(params, l, 0);
  EXPECT_EQ(toroidal.rows(), sem.ltgT(sem.eu() - 1,
                                       sem.mesh().NN() - 1) + 1);
  const auto baseToroidal = sem.rvBaseFullT(params, l);
  EXPECT_EQ(baseToroidal.cols(), sem.ltgT(selected.back(), sem.mesh().NN() - 1) -
                                     sem.ltgT(selected.front(), 0) + 1);
  EXPECT_GT(baseToroidal.cwiseAbs().maxCoeff(), 0.0);

  // The toroidal receiver vectors must interpolate the discontinuous solid
  // field at the interface, not the distinct fluid-side field.
  Eigen::VectorXcd toroidalSolidField = Eigen::VectorXcd::Zero(toroidal.rows());
  for (int e = sem.el(); e < sem.eu(); ++e) {
    for (int n = 0; n < sem.mesh().NN(); ++n)
      toroidalSolidField(sem.ltgT(e, n)) =
          4.0 + 2.0 * sem.mesh().NodeRadius(e, n);
  }
  const auto baseThetaT = sem.rvBaseThetaT(params, l, 0);
  const auto basePhiT = sem.rvBasePhiT(params, l, 0);
  const auto baseThetaValue = (baseThetaT * toroidalSolidField)(0, 0);
  const auto basePhiValue = (basePhiT * toroidalSolidField)(0, 0);
  EXPECT_NEAR(baseThetaValue.real(), 0.5 * solidFieldAtReceiver, 1e-12);
  EXPECT_NEAR(basePhiValue.real(), 0.5 * solidFieldAtReceiver, 1e-12);
  EXPECT_GT(std::abs(2.0 * baseThetaValue.real() -
                     fluidFieldAtReceiver), 1.0);

  const auto phiT = sem.rvPhiT(params, l, 0);
  const double theta = (90.0 - params.receivers()[0].first) * EIGEN_PI / 180.0;
  const double phi = params.receivers()[0].second * EIGEN_PI / 180.0;
  auto wigd = GSHTrans::Wigner<double, GSHTrans::All, GSHTrans::All,
                               GSHTrans::Single>(l, l, 1, theta);
  const std::complex<double> imaginary(0.0, 1.0);
  double maxThetaFluidDifference = 0.0;
  for (int m = -l; m <= l; ++m) {
    const auto phase = std::exp(std::complex<double>(0.0, m * phi));
    const auto ylm = wigd[-1][l, m] * phase;
    const auto ylp = wigd[1][l, m] * phase;
    const auto thetaActual =
        (toroidal.col(m + l).transpose() * toroidalSolidField)(0, 0);
    const auto phiActual =
        (phiT.col(m + l).transpose() * toroidalSolidField)(0, 0);
    const auto thetaExpected =
        -imaginary * (ylp + ylm) * solidFieldAtReceiver / 2.0;
    const auto phiExpected = (ylp - ylm) * solidFieldAtReceiver / 2.0;
    EXPECT_NEAR(std::abs(thetaActual - thetaExpected), 0.0, 1e-12);
    EXPECT_NEAR(std::abs(phiActual - phiExpected), 0.0, 1e-12);
    maxThetaFluidDifference = std::max(
        maxThetaFluidDifference,
        std::abs(thetaExpected -
                 (-imaginary * (ylp + ylm) * fluidFieldAtReceiver / 2.0)));
  }
  EXPECT_GT(maxThetaFluidDifference, 1.0);

  Eigen::VectorXcd localSolidField = Eigen::VectorXcd::Zero(baseToroidal.cols());
  const int localLow = sem.ltgT(solidElement, 0);
  for (int n = 0; n < sem.mesh().NN(); ++n)
    localSolidField(sem.ltgT(solidElement, n) - localLow) =
        4.0 + 2.0 * sem.mesh().NodeRadius(solidElement, n);
  const auto fullToroidalValue = baseToroidal.block(1, 0, 1, baseToroidal.cols()) *
                                 localSolidField;
  EXPECT_NEAR(fullToroidalValue(0, 0).real(),
              std::sqrt(l * (l + 1.0)) / 2.0 * solidFieldAtReceiver, 1e-12);
  EXPECT_GT(std::abs(fullToroidalValue(0, 0).real() * 2.0 /
                         std::sqrt(l * (l + 1.0)) -
                     fluidFieldAtReceiver),
            1.0);
}

TEST(SEMComponentTests, SourceCutInsideSolidBoundaryElementPreservesTopology) {
  prem_norm<double> norm;
  auto model = EarthModels::ModelInput(DSpecMTest::modelPath().string(), norm);
  EarthMesh::RadialMesh mesh(model, 4, 1.0, 0.2, false);
  ASSERT_FALSE(mesh.FS_Boundaries().empty());
  const int boundaryElement = mesh.FS_Boundaries().front();
  ASSERT_TRUE(mesh.IsSolid(boundaryElement));
  ASSERT_LT(boundaryElement + 1, mesh.NE());
  ASSERT_FALSE(mesh.IsSolid(boundaryElement + 1));
  const double interfaceRadius = mesh.EUR(boundaryElement);
  const double oneUlpBelowInterface =
      std::nextafter(interfaceRadius, 0.0);
  EXPECT_THROW(Full1D::SEM(model, 0.2, 4, 2, oneUlpBelowInterface),
               std::invalid_argument);
  const int originalElements = mesh.NE();
  const double sourceRadius =
      0.5 * (mesh.ELR(boundaryElement) + mesh.EUR(boundaryElement));
  const int layer = mesh.LayerNumber(boundaryElement);

  mesh.splitAtRadius(sourceRadius);

  EXPECT_EQ(mesh.NE(), originalElements + 1);
  ASSERT_EQ(mesh.FS_Boundaries().front(), boundaryElement + 1);
  EXPECT_EQ(mesh.LayerNumber(boundaryElement), layer);
  EXPECT_EQ(mesh.LayerNumber(boundaryElement + 1), layer);
  EXPECT_TRUE(mesh.IsSolid(boundaryElement));
  EXPECT_TRUE(mesh.IsSolid(boundaryElement + 1));
  EXPECT_DOUBLE_EQ(mesh.EUR(boundaryElement), sourceRadius);
  EXPECT_DOUBLE_EQ(mesh.ELR(boundaryElement + 1), sourceRadius);

  DSpecMTest::TempDir temp;
  const auto solidPath = DSpecMTest::writeFile(
      temp.path() / "homogeneous_solid.deck",
      "homogeneous solid\n"
      "0 1.0 1\n"
      "2 0 0\n"
      "0.0 5515.0 8000.0 4500.0 300.0 150.0\n"
      "6371000.0 5515.0 8000.0 4500.0 300.0 150.0\n");
  auto solidModel = EarthModels::ModelInput(solidPath.string(), norm);
  EarthMesh::RadialMesh alignedMesh(solidModel, 4, 1.0, 0.2, false);
  const int cutElement = 3;
  ASSERT_LT(cutElement + 1, alignedMesh.NE());
  const double exactNode = alignedMesh.EUR(cutElement);
  const double oneUlpBelow = std::nextafter(exactNode, 0.0);
  const int elementsBefore = alignedMesh.NE();
  alignedMesh.splitAtRadius(oneUlpBelow);
  EXPECT_EQ(alignedMesh.NE(), elementsBefore);
  EXPECT_DOUBLE_EQ(alignedMesh.SourceCutRadius(), exactNode);

  Full1D::SEM nearNodeSem(solidModel, 0.2, 4, 2, oneUlpBelow);
  EXPECT_EQ(nearNodeSem.mesh().NE(), elementsBefore);
  EXPECT_DOUBLE_EQ(nearNodeSem.sourceRadius(), exactNode);
  const double depthKm = (solidModel.OuterRadius() - exactNode) *
                         solidModel.LengthNorm() / 1000.0;
  SourceInfo::EarthquakeCMT nearSource(0.0, 0.0, depthKm,
                                       1.0, 0.0, 0.0, 1.0, 0.0, 1.0);
  EXPECT_EQ(nearNodeSem.sourceElement(nearSource), cutElement + 1);
}

TEST(SEMComponentTests, SpheroidalSourceLiftMatchesAssembledOperatorAwayFromTrace) {
  auto paramsNew = makeTinySemParams();
  Full1D::SEM sem(paramsNew);
  auto &cmt = paramsNew.cmt();
  ASSERT_TRUE(sem.sourceMatches(cmt));
  SourceInfo::EarthquakeCMT nearCut(
      cmt.Latitude(), cmt.Longitude(), cmt.Depth() + 1e-11,
      cmt.Mrr(), cmt.Mrt(), cmt.Mrp(), cmt.Mtt(), cmt.Mtp(), cmt.Mpp());
  ASSERT_TRUE(sem.sourceMatches(nearCut));
  EXPECT_EQ(sem.sourceElement(nearCut), sem.sourceElement(cmt));
  constexpr int degree = 2;
  const std::complex<double> z{0.37, 0.015};
  const auto treatment =
      sem.sourceTreatedSpheroidal(cmt, degree, z, {}, false);
  const int element = treatment.sourceElement;
  const double radius = sem.sourceRadius();
  const double lambda = degree * (degree + 1.0);
  const double kd2 = std::sqrt(lambda / 2.0);
  const std::complex<double> expectedU =
      1.0 / (radius * radius * sem.meshModel().C(element, 0));
  const std::complex<double> expectedV =
      kd2 / (lambda * radius * radius * sem.meshModel().L(element, 0));
  EXPECT_NEAR(std::abs(treatment.lift(sem.ltgS(0, element, 0), 0) -
                       expectedU),
              0.0, 1e-14 * std::abs(expectedU));
  EXPECT_NEAR(std::abs(treatment.lift(sem.ltgS(1, element, 0), 2) -
                       expectedV),
              0.0, 1e-14 * std::abs(expectedV));

  Eigen::MatrixXcd original = sem.calculateForceAll(cmt, degree, false);
  original.col(0).setZero();
  original.col(2).setZero();
  const Eigen::MatrixXcd qOperator =
      sem.hS(degree).cast<std::complex<double>>() -
      z * z * sem.pS(degree).cast<std::complex<double>>();
  const auto assembledAction = qOperator * treatment.lift;
  const int uSource = sem.ltgS(0, element, 0);
  const int vSource = sem.ltgS(1, element, 0);
  const int phiSource = sem.ltgS(2, element, 0);
  std::set<int> sharedTraceRows{uSource, vSource, phiSource};
  if (element > 0)
    for (int field = 0; field < 3; ++field)
      for (int node = 0; node < sem.mesh().NN(); ++node)
        sharedTraceRows.insert(sem.ltgS(field, element - 1, node));

  for (const int column : {0, 2}) {
    Eigen::VectorXcd qj = Eigen::VectorXcd::Zero(treatment.force.rows());
    if (column == 0) {
      const double j = expectedU.real();
      qj[uSource] = 2.0 * radius * sem.meshModel().F(element, 0) * j;
      qj[vSource] = -lambda * radius * sem.meshModel().F(element, 0) * j;
    } else {
      const double j = expectedV.real();
      qj[uSource] = lambda * radius * sem.meshModel().L(element, 0) * j;
      qj[vSource] = -lambda * radius * sem.meshModel().L(element, 0) * j;
    }
    const Eigen::VectorXcd localAction =
        original.col(column) - qj - treatment.force.col(column);
    double error = 0.0, scale = 0.0;
    for (Eigen::Index row = 0; row < localAction.size(); ++row) {
      if (sharedTraceRows.count(static_cast<int>(row)) != 0) continue;
      error += std::norm(localAction[row] - assembledAction(row, column));
      scale += std::norm(assembledAction(row, column));
    }
    EXPECT_LE(std::sqrt(error) / std::max(std::sqrt(scale), 1e-30), 2e-12);
  }
}

TEST(SEMComponentTests, SpheroidalLocalWeakFormScattersToProductionOperator) {
  DSpecMTest::TempDir temp;
  const auto solidPath = DSpecMTest::writeFile(
      temp.path() / "homogeneous_solid.deck",
      "homogeneous solid\n"
      "0 1.0 1\n"
      "2 0 0\n"
      "0.0 5515.0 8000.0 4500.0 300.0 150.0\n"
      "6371000.0 5515.0 8000.0 4500.0 300.0 150.0\n");
  auto paramsNew = makeTinySemParams(0.0, solidPath);
  Full1D::SEM sem(paramsNew);
  constexpr int degree = 2;
  const std::complex<double> z{0.37, 0.015};
  const int nq = sem.mesh().NN(), dofs = sem.hS(degree).rows();
  Eigen::MatrixXcd scattered = Eigen::MatrixXcd::Zero(dofs, dofs);
  for (int e = 0; e < sem.mesh().NE(); ++e) {
    ASSERT_FALSE(sem.elementIsFluid(e));
    const auto local = Full1D::source_treatment_detail::spheroidalElementOperator(
        sem, e, degree, z, {}, false);
    for (int a = 0; a < 3; ++a)
      for (int i = 0; i < nq; ++i)
        for (int b = 0; b < 3; ++b)
          for (int j = 0; j < nq; ++j)
            scattered(sem.ltgS(a, e, i), sem.ltgS(b, e, j)) +=
                local(a * nq + i, b * nq + j);
  }
  const int surfacePhi = sem.ltgS(2, sem.mesh().NE() - 1, nq - 1);
  scattered(surfacePhi, surfacePhi) += sem.mesh().OR() * (degree + 1.0) /
      (4.0 * static_cast<double>(EIGEN_PI) * sem.gravityCoefficient());
  const Eigen::MatrixXcd production =
      sem.hS(degree).cast<std::complex<double>>() -
      z * z * sem.pS(degree).cast<std::complex<double>>();
  const double operatorScale = std::max(production.norm(), 1e-30);
  std::array<std::vector<int>, 3> fieldIndices;
  for (int field = 0; field < 3; ++field) {
    for (int e = 0; e < sem.mesh().NE(); ++e)
      for (int i = 0; i < nq; ++i)
        fieldIndices[field].push_back(sem.ltgS(field, e, i));
    auto &indices = fieldIndices[field];
    std::sort(indices.begin(), indices.end());
    indices.erase(std::unique(indices.begin(), indices.end()), indices.end());
  }
  for (int rowField = 0; rowField < 3; ++rowField)
    for (int colField = 0; colField < 3; ++colField) {
      const auto &rows = fieldIndices[rowField];
      const auto &cols = fieldIndices[colField];
      Eigen::MatrixXcd blockError(rows.size(), cols.size());
      Eigen::MatrixXcd blockScale(rows.size(), cols.size());
      for (std::size_t i = 0; i < rows.size(); ++i)
        for (std::size_t j = 0; j < cols.size(); ++j) {
          blockError(i, j) = scattered(rows[i], cols[j]) - production(rows[i], cols[j]);
          blockScale(i, j) = production(rows[i], cols[j]);
        }
      const double scale = std::max(blockScale.norm(), operatorScale * 1e-14);
      EXPECT_LE(blockError.norm(), 2e-12 * scale)
          << "field block " << rowField << ',' << colField;
    }

  const int radialDofs = sem.hR().rows();
  Eigen::MatrixXcd radialScattered = Eigen::MatrixXcd::Zero(radialDofs, radialDofs);
  for (int e = 0; e < sem.mesh().NE(); ++e) {
    const auto local = Full1D::source_treatment_detail::radialElementOperator(
        sem, e, z);
    for (int a = 0; a < 2; ++a)
      for (int i = 0; i < nq; ++i)
        for (int b = 0; b < 2; ++b)
          for (int j = 0; j < nq; ++j)
            radialScattered(sem.ltgR(a, e, i), sem.ltgR(b, e, j)) +=
                local(a * nq + i, b * nq + j);
  }
  const int surfacePhiR = sem.ltgR(1, sem.mesh().NE() - 1, nq - 1);
  radialScattered(surfacePhiR, surfacePhiR) += sem.mesh().OR() /
      (4.0 * static_cast<double>(EIGEN_PI) * sem.gravityCoefficient());
  const Eigen::MatrixXcd radialProduction = sem.hR().cast<std::complex<double>>() -
      z * z * sem.pR().cast<std::complex<double>>();
  ASSERT_EQ(radialScattered.rows(), radialProduction.rows());
  EXPECT_LE((radialScattered - radialProduction).norm(),
            2e-12 * radialProduction.norm()) << "radial full scatter";

  const Eigen::MatrixXcd toroidalProduction =
      sem.hTk(degree).cast<std::complex<double>>() -
      z * z * sem.pTk(degree).cast<std::complex<double>>();
  Eigen::MatrixXcd toroidalScattered = Eigen::MatrixXcd::Zero(
      toroidalProduction.rows(), toroidalProduction.cols());
  for (int e = sem.el(); e < sem.eu(); ++e) {
    const auto local = Full1D::source_treatment_detail::toroidalElementOperator(
        sem, e, degree, z);
    for (int i = 0; i < nq; ++i)
      for (int j = 0; j < nq; ++j)
        toroidalScattered(sem.ltgT(e, i), sem.ltgT(e, j)) += local(i, j);
  }
  EXPECT_LE((toroidalScattered - toroidalProduction).norm(),
            2e-12 * toroidalProduction.norm()) << "toroidal full scatter";
}

TEST(SEMComponentTests, SpheroidalSourceLiftRemainsFiniteAtDegreeThreeHundred) {
  auto paramsNew = makeTinySemParams();
  auto &cmt = paramsNew.cmt();
  Full1D::SEM sem(paramsNew.earthModel(), paramsNew.maxstep(), 4, 300, cmt);
  const auto treatment = sem.sourceTreatedSpheroidal(cmt, 300, {0.31, 0.01});
  EXPECT_TRUE(treatment.force.allFinite());
  EXPECT_TRUE(treatment.lift.allFinite());
}

TEST(SEMComponentTests, SurfaceOceanCapsToroidalRangeAboveSolidSphere) {
  DSpecMTest::TempDir temp;
  const auto oceanPath = DSpecMTest::writeFile(
      temp.path() / "homogeneous_ocean.deck",
      "benchmark homogeneous mantle with a 3 km ocean\n"
      "0 1.0 1\n"
      "4 0 0\n"
      "0.0 5515.0 8000.0 4500.0 300.0 150.0\n"
      "6368000.0 5515.0 8000.0 4500.0 300.0 150.0\n"
      "6368000.0 1025.0 1500.0 0.0 1000.0 0.0\n"
      "6371000.0 1025.0 1500.0 0.0 1000.0 0.0\n");
  auto paramsNew = makeTinySemParams(3.0, oceanPath);
  Full1D::SEM sem(paramsNew);
  const auto &model = paramsNew.earthModel();
  auto &params = paramsNew.inputParameters();

  int firstFluid = sem.mesh().NE();
  for (int e = 0; e < sem.mesh().NE(); ++e) {
    if (model.IsFluid(sem.mesh().LayerNumber(e))) {
      firstFluid = e;
      break;
    }
  }
  ASSERT_GT(firstFluid, 0);
  ASSERT_LT(firstFluid, sem.mesh().NE());
  for (int e = firstFluid; e < sem.mesh().NE(); ++e)
    ASSERT_TRUE(model.IsFluid(sem.mesh().LayerNumber(e)));
  EXPECT_EQ(sem.el(), 0);
  EXPECT_EQ(sem.eu(), firstFluid);

  const auto selected = sem.receiverElements(params, true);
  ASSERT_FALSE(selected.empty());
  for (const int e : selected) {
    EXPECT_LT(e, firstFluid);
    EXPECT_FALSE(model.IsFluid(sem.mesh().LayerNumber(e)));
  }
  EXPECT_NEAR(sem.mesh().EUR(selected.back()),
              sem.mesh().ELR(firstFluid), 1e-10);

  const auto toroidalReceiver = sem.rvThetaT(params, 2, 0);
  EXPECT_EQ(toroidalReceiver.rows(), sem.ltgT(sem.eu() - 1,
                                               sem.mesh().NN() - 1) + 1);
  EXPECT_TRUE(toroidalReceiver.allFinite());
  EXPECT_GT(toroidalReceiver.norm(), 0.0);
  const auto toroidalStiffness = sem.hTk(2);
  EXPECT_GT(toroidalStiffness.rows(), 0);
  EXPECT_GT(toroidalStiffness.norm(), 0.0);
}

TEST(SEMComponentTests, FluidAndOutOfMeshReceiversFailBeforeAssembly) {
  auto waterParams = makeTinySemParams(2.999);
  Full1D::SEM sem(waterParams);
  EXPECT_NO_THROW(sem.receiverElements(waterParams.inputParameters()));
  EXPECT_THROW(sem.receiverElements(waterParams.inputParameters(), true),
               std::invalid_argument);
  EXPECT_NO_THROW(sem.rvBaseFull(waterParams.inputParameters(), 2));

  auto outsideParams = makeTinySemParams(7000.0);
  Full1D::SEM outsideSem(outsideParams);
  EXPECT_THROW(outsideSem.receiverElements(outsideParams.inputParameters()),
               std::out_of_range);
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
  EXPECT_THROW(sem.rvBaseFullT(params, idxl), std::invalid_argument);

  const auto fullRadialRows = sem.ltgR(1, sem.mesh().NE() - 1, sem.mesh().NN() - 1) + 1;
  const auto reducedRadialRows =
      sem.ltgR(1, receiverElems.back(), sem.mesh().NN() - 1) -
      sem.ltgR(0, receiverElems.front(), 0) + 1;
  const auto baseFullCols =
      sem.ltgS(1, receiverElems.back(), sem.mesh().NN() - 1) -
      sem.ltgS(0, receiverElems.front(), 0) + 1;

  ASSERT_EQ(radial.rows(), fullRadialRows);
  ASSERT_EQ(radial.cols(), 1);
  ASSERT_EQ(reducedRadial.rows(), reducedRadialRows);
  ASSERT_EQ(reducedRadial.cols(), 1);
  ASSERT_EQ(baseFull.rows(), 3 * params.num_receivers());
  ASSERT_EQ(baseFull.cols(), baseFullCols);

  EXPECT_GT(radial.cwiseAbs().maxCoeff(), 0.0);
  EXPECT_GT(reducedRadial.cwiseAbs().maxCoeff(), 0.0);
  EXPECT_GT(baseFull.cwiseAbs().maxCoeff(), 0.0);
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

TEST(SEMComponentTests, CowlingStartIndicesPreserveTruncationAndCadence) {
  auto paramsNew = makeTinySemParams();
  Full1D::SEM sem(paramsNew);
  const int l = 100;
  const int source = sem.mesh().NE() - 1;
  const int nskip = 3;
  std::vector<double> frequencies{0.02, 0.025, 0.03, 0.035, 0.04, 0.045};

  auto expected = [&](const std::vector<double> &w, int begin, int end,
                      int sourceElement, int cadence, bool cowling) {
    std::vector<int> indices(end - begin);
    std::vector<int> elements(end - begin);
    for (int idx = end - 1; idx >= begin; --idx) {
      const int offset = idx - begin;
      if ((end - 1 - idx) % cadence == 0) {
        elements[offset] = sourceElement < 0
                               ? SpectralTools::startElementSph(sem, l, w[idx])
                               : SpectralTools::startElementSph(
                                     sem, l, w[idx], sourceElement);
        indices[offset] = cowling
                              ? sem.ltgSC(0, elements[offset], 0)
                              : sem.ltgS(0, elements[offset], 0);
      } else {
        elements[offset] = elements[offset + 1];
        indices[offset] = indices[offset + 1];
      }
    }
    return std::pair<std::vector<int>, std::vector<int>>(indices, elements);
  };

  const auto fullEvery = SpectralTools::allIndicesSph(
      sem, l, frequencies, source, 1, false);
  const auto cowlingEvery = SpectralTools::allIndicesSph(
      sem, l, frequencies, source, 1, true);
  const auto expectedFullEvery = expected(
      frequencies, 0, static_cast<int>(frequencies.size()), source, 1, false);
  const auto expectedCowlingEvery = expected(
      frequencies, 0, static_cast<int>(frequencies.size()), source, 1, true);
  EXPECT_EQ(fullEvery, expectedFullEvery.first);
  EXPECT_EQ(cowlingEvery, expectedCowlingEvery.first);
  EXPECT_EQ(expectedFullEvery.second, expectedCowlingEvery.second);
  EXPECT_GT(*std::max_element(expectedFullEvery.second.begin(),
                              expectedFullEvery.second.end()),
            0);

  const auto fullCadenced = SpectralTools::allIndicesSph(
      sem, l, frequencies, source, nskip, false);
  const auto cowlingCadenced = SpectralTools::allIndicesSph(
      sem, l, frequencies, source, nskip, true);
  const auto expectedFullCadenced = expected(
      frequencies, 0, static_cast<int>(frequencies.size()), source, nskip,
      false);
  const auto expectedCowlingCadenced = expected(
      frequencies, 0, static_cast<int>(frequencies.size()), source, nskip,
      true);
  EXPECT_EQ(fullCadenced, expectedFullCadenced.first);
  EXPECT_EQ(cowlingCadenced, expectedCowlingCadenced.first);
  EXPECT_EQ(expectedFullCadenced.second, expectedCowlingCadenced.second);
  EXPECT_EQ(fullCadenced[5], fullEvery[5]);
  EXPECT_EQ(fullCadenced[2], fullEvery[2]);
  EXPECT_EQ(cowlingCadenced[5], cowlingEvery[5]);
  EXPECT_EQ(cowlingCadenced[2], cowlingEvery[2]);

  auto &freqFull = paramsNew.freqFull();
  const auto fullFreqExpected = expected(
      freqFull.w(), freqFull.i1(), freqFull.i2(), -1, 2, false);
  const auto cowlingFreqExpected = expected(
      freqFull.w(), freqFull.i1(), freqFull.i2(), -1, 2, true);
  EXPECT_EQ(SpectralTools::allIndicesSph(sem, l, freqFull, 2),
            fullFreqExpected.first);
  EXPECT_EQ(SpectralTools::allIndicesSph(sem, l, freqFull, 2, true),
            cowlingFreqExpected.first);

  const auto fullFreqSourceExpected = expected(
      freqFull.w(), freqFull.i1(), freqFull.i2(), source, 2, false);
  const auto cowlingFreqSourceExpected = expected(
      freqFull.w(), freqFull.i1(), freqFull.i2(), source, 2, true);
  EXPECT_EQ(SpectralTools::allIndicesSph(sem, l, freqFull, source, 2),
            fullFreqSourceExpected.first);
  EXPECT_EQ(SpectralTools::allIndicesSph(sem, l, freqFull, source, 2, true),
            cowlingFreqSourceExpected.first);

  std::vector<double> sourceAtZero{0.02, 0.03};
  const auto zeroSource = SpectralTools::allIndicesSph(
      sem, l, sourceAtZero, 0, 1, true);
  EXPECT_EQ(zeroSource,
            (std::vector<int>{sem.ltgSC(0, 0, 0), sem.ltgSC(0, 0, 0)}));
}
