#include <gtest/gtest.h>

#include <array>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <DSpecM1D/SensitivityKernels>
#include <DSpecM1D/Windowed3DFields>
#include <DSpecM1D/Windowed3DKernels>
#include <DSpecM1D/ModelInput>

#include "test_utils.h"

namespace {

Eigen::MatrixXd directPerturbationGravity(const Full1D::SEM &sem,
                                          const Eigen::MatrixXd &deltaDensity);

TEST(WindowedDisplacementMapTests, ForwardAndAdjointMatchRealInnerProduct) {
  SpectraSolver::FreqFull frequency(0.1, 1.0, 0.1, 0.2, 0.8, 1.0,
                                   10.0, 0.3, 0.1, 0.5, 0.0, 0.2,
                                   2, 3600.0);
  Eigen::VectorXd window(frequency.nt());
  for (int i = 0; i < window.size(); ++i)
    window[i] = 0.2 + 0.7 * std::cos(0.31 * i);
  SPARSESPEC::WindowedDisplacementMap map(frequency, 6.371e6, window);

  Eigen::MatrixXcd raw = Eigen::MatrixXcd::Random(3, map.frequencyCount());
  Eigen::MatrixXcd q = Eigen::MatrixXcd::Random(3, map.frequencyCount());
  q(0, 0) = {0.4, -0.9};
  q(1, map.frequencyCount() - 1) = {-0.7, 0.3};
  const Eigen::MatrixXcd mapped = map.forward(raw);
  const Eigen::MatrixXcd adjoint = map.adjoint(q);
  const double lhs = q.conjugate().cwiseProduct(mapped).sum().real();
  const double rhs = adjoint.conjugate().cwiseProduct(raw).sum().real();
  const double relative = std::abs(lhs - rhs) /
                          std::max({std::abs(lhs), std::abs(rhs), 1e-30});
  EXPECT_LE(relative, 1e-12);
  EXPECT_TRUE(mapped.allFinite());
  EXPECT_TRUE(adjoint.allFinite());
}

TEST(FixedMeshAdjointTests, SolvesFullHermitianSystemForCompactColumns) {
  prem_norm<double> norm;
  auto model = EarthModels::ModelInput(DSpecMTest::modelPath().string(), norm);
  Full1D::SEM sem(model, 0.1, 4, 3);
  const int degree = 2;
  const int dofs = static_cast<int>(sem.hS(degree).rows());
  const Eigen::MatrixXcd rhs = Eigen::MatrixXcd::Random(dofs, 3);
  const Eigen::MatrixXcd receiver = Eigen::MatrixXcd::Random(4, dofs);
  const Eigen::MatrixXcd weights = Eigen::MatrixXcd::Random(4, 3);
  const std::complex<double> z{0.31, -0.02};

  const auto fields = SPARSESPEC::FixedMeshAdjointSolver::solveFrequency(
      sem, SPARSESPEC::SensitivitySolveFamily::Spheroidal, degree, z,
      rhs, receiver, weights);
  const Eigen::SparseMatrix<std::complex<double>> h =
      sem.hS(degree).cast<std::complex<double>>() -
      z * z * sem.pS(degree).cast<std::complex<double>>();
  const Eigen::MatrixXcd displacementReceiver =
      (-std::complex<double>(0.0, 1.0) / z) * receiver;
  const double forwardBackwardError = (h * fields.forward - rhs).norm() /
      (h.norm() * fields.forward.norm() + rhs.norm());
  const double adjointBackwardError = (h.adjoint() * fields.adjoint -
      displacementReceiver.adjoint() * weights).norm() /
      (h.norm() * fields.adjoint.norm() +
       (displacementReceiver.adjoint() * weights).norm());
  EXPECT_LT(forwardBackwardError, 1e-12);
  EXPECT_LT(adjointBackwardError, 1e-12);
  const double objective = weights.conjugate().cwiseProduct(
      displacementReceiver * fields.forward).sum().real();
  EXPECT_NEAR(fields.objective, objective,
              1e-13 * std::max(std::abs(objective), 1e-30));
  const double adjointObjective = fields.adjoint.conjugate().cwiseProduct(rhs).sum().real();
  EXPECT_NEAR(adjointObjective, fields.objective,
              1e-11 * std::max(std::abs(fields.objective), 1e-30));
}

} // namespace

namespace {

InputParameters makeSensitivityParameters(int type, int outputType = 0,
                                         double receiverDepthKm = 3.0,
                                         std::pair<double,double> receiver={30.0,65.0}) {
  DSpecMTest::TempDir temp;
  DSpecMTest::ParameterOptions options;
  options.type = type;
  options.outputType = outputType;
  options.lmin = 0;
  options.lmax = 2;
  options.f1 = 0.1; options.f2 = 0.8;
  options.f11 = 0.1; options.f12 = 0.2;
  options.f21 = 0.6; options.f22 = 0.8;
  options.tOutMinutes = 0.5;
  options.timeStepSec = 20.0;
  options.numReceivers = 1;
  options.receiverDepth = receiverDepthKm;
  options.receivers = {receiver};
  const auto file = DSpecMTest::writeFile(
      temp.path() / "sensitivity_params.txt",
      DSpecMTest::makeParameterText(DSpecMTest::modelPath().string(), options));
  return InputParameters(file.string());
}

InputParameters makeMultiReceiverSensitivityParameters(
    int type,double receiverDepthKm,
    std::vector<std::pair<double,double>> receivers) {
  DSpecMTest::TempDir temp;
  DSpecMTest::ParameterOptions options;
  options.type=type;options.outputType=0;options.lmin=0;options.lmax=2;
  options.f1=0.1;options.f2=0.8;options.f11=0.1;options.f12=0.2;
  options.f21=0.6;options.f22=0.8;options.tOutMinutes=0.5;
  options.timeStepSec=20.0;options.numReceivers=static_cast<int>(receivers.size());
  options.receiverDepth=receiverDepthKm;options.receivers=std::move(receivers);
  const auto file=DSpecMTest::writeFile(temp.path()/"multi_receiver_params.txt",
      DSpecMTest::makeParameterText(DSpecMTest::modelPath().string(),options));
  return InputParameters(file.string());
}

SpectraSolver::FreqFull makeSensitivityFrequency(double timeNorm) {
  return SpectraSolver::FreqFull(0.1, 0.8, 0.1, 0.2, 0.6, 0.8,
                                 20.0, 0.5, 0.1, 0.5, 0.0, 0.3,
                                 2, timeNorm);
}

SourceInfo::EarthquakeCMT makeSensitivitySource(double depthKm = 30.0,
                                                double latitude = 15.0,
                                                double longitude = 110.0) {
  return SourceInfo::EarthquakeCMT(latitude, longitude, depthKm,
                                   1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
}

} // namespace

TEST(SensitivityBoundaryTests, RejectUnconstructedSemAtPublicEntryPoints) {
  Full1D::SEM sem;
  const Eigen::MatrixXcd one = Eigen::MatrixXcd::Ones(1, 1);
  const Eigen::MatrixXcd rhs = Eigen::MatrixXcd::Ones(1, 1);
  EXPECT_THROW(SPARSESPEC::FixedMeshAdjointSolver::solveFrequency(
                   sem, SPARSESPEC::SensitivitySolveFamily::Radial, 0,
                   {0.3, -0.1}, rhs, one, one), std::invalid_argument);
  EXPECT_THROW(SPARSESPEC::FixedMeshAdjointSolver::solveForwardFrequency(
                   sem, SPARSESPEC::SensitivitySolveFamily::Radial, 0,
                   {0.3, -0.1}, rhs), std::invalid_argument);

  SPARSESPEC::FrequencyAdjointFields fields;
  fields.forward = one;
  fields.adjoint = one;
  EXPECT_THROW(SPARSESPEC::modulusKernelSamples(
                   sem, SPARSESPEC::SensitivitySolveFamily::Radial, 0, fields),
               std::invalid_argument);

  auto parameters = makeSensitivityParameters(1);
  SourceInfo::EarthquakeCMT source(15.0, 110.0, 30.0,
                                   1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
  const auto frequency = makeSensitivityFrequency(1.0);
  const Eigen::VectorXd window = Eigen::VectorXd::Ones(frequency.nt());
  const Eigen::MatrixXcd weights = Eigen::MatrixXcd::Zero(
      3, frequency.nt() / 2 + 1);
  EXPECT_THROW(SPARSESPEC::WindowedSensitivityEngine::evaluate(
                   sem, source, parameters, frequency, window, weights,
                   1, static_cast<int>(weights.cols())), std::invalid_argument);
}

TEST(WindowedSensitivityEngineTests, ToroidalProjectionHasExactVerticalZeroAndHorizontalSignal) {
  auto parameters = makeSensitivityParameters(2);
  prem_norm<double> norm;
  auto model = EarthModels::ModelInput(DSpecMTest::modelPath().string(), norm);
  Full1D::SEM sem(model, 0.25, 4, 2, makeSensitivitySource());
  SourceInfo::EarthquakeCMT source(15.0, 110.0, 30.0,
                                   1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
  const auto frequency = makeSensitivityFrequency(model.TimeNorm());
  Eigen::VectorXd window = Eigen::VectorXd::Ones(frequency.nt());
  Eigen::MatrixXcd qVertical = Eigen::MatrixXcd::Zero(3, frequency.nt() / 2 + 1);
  qVertical(0, 3) = {0.7, 0.2};
  const auto vertical = SPARSESPEC::WindowedSensitivityEngine::evaluate(
      sem, source, parameters, frequency, window, qVertical, 2, 4);
  EXPECT_EQ(vertical.rawDisplacementSpectrum.row(0).norm(), 0.0);
  EXPECT_EQ(vertical.objective, 0.0);
  EXPECT_EQ(vertical.kernels.L.norm(), 0.0);
  EXPECT_EQ(vertical.kernels.N.norm(), 0.0);

  Eigen::MatrixXcd qHorizontal = Eigen::MatrixXcd::Zero(3, frequency.nt() / 2 + 1);
  qHorizontal(1, 3) = {0.7, 0.2};
  const auto horizontal = SPARSESPEC::WindowedSensitivityEngine::evaluate(
      sem, source, parameters, frequency, window, qHorizontal, 2, 4);
  EXPECT_GT(horizontal.rawDisplacementSpectrum.row(1).norm(), 0.0);
  EXPECT_GT(horizontal.kernels.L.norm() + horizontal.kernels.N.norm(), 0.0);
}

TEST(Windowed3DFieldEngineTests, ReconstructsRadialAndCoupledFieldsAtReceiversAndPoles) {
  for (const int modeType : {1, 2, 3, 4}) {
    auto parameters = makeSensitivityParameters(modeType, 0, (modeType == 2 || modeType == 4) ? 1000.0 : 3000.0);
    prem_norm<double> norm;
    auto model = EarthModels::ModelInput(DSpecMTest::modelPath().string(), norm);
    auto source = makeSensitivitySource();
    Full1D::SEM sem(model, 0.25, 4, 2, source);
    const auto frequency = makeSensitivityFrequency(model.TimeNorm());
    Eigen::VectorXd window = Eigen::VectorXd::Ones(frequency.nt());
    Eigen::MatrixXcd weights = Eigen::MatrixXcd::Zero(
        3 * parameters.num_receivers(), frequency.nt() / 2 + 1);
    weights(0, 3) = {0.7, -0.2};
    if (modeType == 4) {
      weights(1, 3) = {-0.3, 0.5};
      weights(2, 3) = {0.2, 0.4};
    }
    const int bin = 3;
    const std::complex<double> z = frequency.w().at(bin) +
        std::complex<double>(0.0, -frequency.ep());
    SPARSESPEC::WindowedDisplacementMap map(
        frequency, sem.lengthScaleMetres(), window);
    const Eigen::VectorXcd adjointWeights = map.adjoint(weights).col(bin);
    const double receiverRadius = (sem.mesh().PR() -
        1000.0 * parameters.receiver_depth() / sem.lengthScaleMetres()) *
        sem.lengthScaleMetres();
    const auto rec = parameters.receivers().front();
    const double deg = EIGEN_PI / 180.0;
    const int sourceElement = sem.sourceElement(source);
    std::vector<SPARSESPEC::Field3DPoint> points{
        {receiverRadius, (90.0 - rec.first) * deg, rec.second * deg},
        {receiverRadius, 0.0, 0.0}, {receiverRadius, EIGEN_PI, 2.1},
        {sem.sourceRadius() * sem.lengthScaleMetres(), 0.9, 1.7,
         sourceElement},
        {sem.sourceRadius() * sem.lengthScaleMetres(), 0.9, 1.7,
         sourceElement - 1}};
    int fluidPoint = -1, innerCorePoint = -1;
    if (modeType == 4) {
      for (int e = 0; e < sem.mesh().NE(); ++e) {
        const double radius = 0.5*(sem.mesh().ELR(e)+sem.mesh().EUR(e)) *
                              sem.lengthScaleMetres();
        if (fluidPoint < 0 && sem.elementIsFluid(e)) {
          fluidPoint = static_cast<int>(points.size());
          points.push_back({radius, 0.8, 1.1, e});
        }
        if (innerCorePoint < 0 && e < sem.el() && !sem.elementIsFluid(e)) {
          innerCorePoint = static_cast<int>(points.size());
          points.push_back({radius, 0.8, 1.1, e});
        }
      }
      ASSERT_GE(fluidPoint, 0);
      ASSERT_GE(innerCorePoint, 0);
    }
    const auto fields = SPARSESPEC::Windowed3DFieldEngine::evaluateFrequency(
        sem, source, parameters, z,
        adjointWeights, points);
    if(modeType==4){
      auto workspace=SPARSESPEC::Windowed3DFieldEngine::prepareFrequencyWorkspace(
          sem,source,parameters,z,adjointWeights);
      const auto ownedCount=workspace.modalBlocks().size();
      const auto all=workspace.sample(points);
      const std::vector<SPARSESPEC::Field3DPoint> first(points.begin(),points.begin()+3);
      const std::vector<SPARSESPEC::Field3DPoint> second(points.begin()+3,points.end());
      const auto left=workspace.sample(first),right=workspace.sample(second);
      ASSERT_EQ(all.size(),fields.size());
      ASSERT_EQ(ownedCount,workspace.modalBlocks().size());
      for(std::size_t i=0;i<all.size();++i){
        const auto &chunk=i<left.size()?left[i]:right[i-left.size()];
        EXPECT_LT((Eigen::Map<const Eigen::VectorXcd>(all[i].forwardRawDisplacement.data(),3)-
                   Eigen::Map<const Eigen::VectorXcd>(chunk.forwardRawDisplacement.data(),3)).norm(),1e-12);
        EXPECT_LT((Eigen::Map<const Eigen::VectorXcd>(all[i].adjointDualDisplacement.data(),3)-
                   Eigen::Map<const Eigen::VectorXcd>(chunk.adjointDualDisplacement.data(),3)).norm(),1e-12);
        EXPECT_LT(std::abs(all[i].forwardRawStrain.radialPlus-chunk.forwardRawStrain.radialPlus),1e-12);
        EXPECT_LT(std::abs(all[i].forwardRawStrain.horizontalMinus-chunk.forwardRawStrain.horizontalMinus),1e-12);
      }
    }
    ASSERT_EQ(fields.size(), points.size());
    for (const auto &sample : fields) {
      EXPECT_TRUE(std::isfinite(sample.forwardRawStrain.rr.real()));
      EXPECT_TRUE(std::isfinite(sample.rScaledForwardStrainMetres.radialPlus.imag()));
      EXPECT_TRUE(std::isfinite(sample.forwardRawDisplacement[0].real()));
      EXPECT_TRUE(std::isfinite(sample.adjointDualDisplacement[1].imag()));
      EXPECT_GE(sample.radialElement, 0);
    }
    EXPECT_EQ(fields[3].radialElement, sourceElement);
    EXPECT_EQ(fields[4].radialElement, sourceElement - 1);
    double traceJump2 = 0.0;
    for (int c = 0; c < 3; ++c)
      traceJump2 += std::norm(fields[3].forwardRawDisplacement[c] -
                              fields[4].forwardRawDisplacement[c]);
    EXPECT_GT(traceJump2, 0.0) << "mode type " << modeType;

    const auto reference = SPARSESPEC::WindowedSensitivityEngine::evaluate(
        sem, source, parameters,
        frequency, window, weights, bin, bin + 1);
    Eigen::VectorXcd reconstructed(3);
    for (int c = 0; c < 3; ++c)
      reconstructed[c] =
          fields.front().observableForwardDisplacementSpectrumMetreSeconds[c];
    const Eigen::VectorXcd expected = sem.lengthScaleMetres() *
        sem.timeScaleSeconds() *
        reference.rawDisplacementSpectrum.col(bin);
    const double error = (reconstructed - expected).norm() /
        std::max({expected.norm(), reconstructed.norm(), 1e-60});
    EXPECT_LT(error, 1e-9) << "mode type " << modeType
                           << " reconstructed=" << reconstructed.transpose()
                           << " expected=" << expected.transpose();
    if (modeType == 4) {
      auto spheroidalParameters = makeSensitivityParameters(3, 0, 1000.0);
      const auto spheroidal = SPARSESPEC::Windowed3DFieldEngine::evaluateFrequency(
          sem, source, spheroidalParameters, z, adjointWeights, points);
      for (const int index : {fluidPoint, innerCorePoint}) {
        EXPECT_EQ(fields[index].fluid, index == fluidPoint);
        Eigen::Vector2cd combinedHorizontal, spheroidalHorizontal;
        for (int c = 0; c < 2; ++c) {
          combinedHorizontal[c] = fields[index].forwardRawDisplacement[c+1];
          spheroidalHorizontal[c] = spheroidal[index].forwardRawDisplacement[c+1];
        }
        const double toroidalLeak = (combinedHorizontal-spheroidalHorizontal).norm() /
            std::max({combinedHorizontal.norm(), spheroidalHorizontal.norm(), 1e-60});
        EXPECT_LT(toroidalLeak, 1e-12) << "out-of-range toroidal sampling at point " << index;
      }
    }
    EXPECT_THROW(SPARSESPEC::Windowed3DFieldEngine::evaluateFrequency(
        sem, source, parameters, z, adjointWeights, {{0.0, 0.7, 0.2}}),
        std::invalid_argument);
    const auto centerScaled = SPARSESPEC::Windowed3DFieldEngine::rScaledStrainFromModal(
        2, {1.0,0.2}, {0.3,-0.1}, {0.4,0.1}, {-0.2,0.3},
        {0.1,-0.4}, {0.2,0.2}, 0.0, {0.7,0.2}, {0.1,-0.3},
        {-0.2,0.1}, {0.4,0.2}, {-0.3,-0.1});
    EXPECT_TRUE(std::isfinite(centerScaled.rr.real()));
    EXPECT_TRUE(std::isfinite(centerScaled.horizontalMinus.imag()));

    if (modeType == 4) {
      const double theta = 0.71, phi = 1.23;
      const double hRModel = 1.0e-3;
      const double hR = hRModel * sem.lengthScaleMetres();
      const double hA = 2.0e-4;
      const double r = receiverRadius;
      std::vector<SPARSESPEC::Field3DPoint> fdPoints{{r, theta, phi}};
      for (const int k : {-2, -1, 1, 2})
        fdPoints.push_back({r + k*hR, theta, phi});
      for (const int k : {-2, -1, 1, 2})
        fdPoints.push_back({r, theta + k*hA, phi});
      for (const int k : {-2, -1, 1, 2})
        fdPoints.push_back({r, theta, phi + k*hA});
      const Eigen::VectorXcd zeroWeights = Eigen::VectorXcd::Zero(
          3 * parameters.num_receivers());
      const auto fd = SPARSESPEC::Windowed3DFieldEngine::evaluateFrequency(
          sem, source, parameters, z, zeroWeights, fdPoints);
      const auto derivative = [&](int base, double h, int c) {
        return (-fd[base+4].forwardRawDisplacement[c] +
                8.0*fd[base+3].forwardRawDisplacement[c] -
                8.0*fd[base+2].forwardRawDisplacement[c] +
                fd[base+1].forwardRawDisplacement[c])/(12.0*h);
      };
      // Each coordinate block is ordered [-2h,-h,+h,+2h] after its base.
      const auto dr = [&](int c) { return derivative(0, hRModel, c); };
      const auto dtheta = [&](int c) { return derivative(4, hA, c); };
      const auto dphi = [&](int c) { return derivative(8, hA, c); };
      const auto &u = fd.front().forwardRawDisplacement;
      const double rModel = r/sem.lengthScaleMetres();
      const double sinTheta = std::sin(theta), cotTheta = 1.0/std::tan(theta);
      // Reconstructed horizontal components are north/east; north is -e_theta.
      const auto eRNorth = 0.5*(dr(1)-u[1]/rModel-dtheta(0)/rModel);
      const auto eREast = 0.5*(dr(2)-u[2]/rModel+
                                dphi(0)/(rModel*sinTheta));
      const auto eNorthNorth = (-dtheta(1)+u[0])/rModel;
      const auto eEastEast = (dphi(2)/sinTheta+u[0]-cotTheta*u[1])/rModel;
      const auto eNorthEast = 0.5*(dphi(1)/sinTheta-dtheta(2)+cotTheta*u[2])/rModel;
      const auto &strain = fd.front().forwardRawStrain;
      const std::array<std::complex<double>,6> finiteDifference{
          dr(0), eNorthNorth+eEastEast, -eRNorth-std::complex<double>(0,1)*eREast,
          -eRNorth+std::complex<double>(0,1)*eREast,
          eNorthNorth-eEastEast+std::complex<double>(0,2)*eNorthEast,
          eNorthNorth-eEastEast-std::complex<double>(0,2)*eNorthEast};
      const std::array<std::complex<double>,6> reconstructed{
          strain.rr, strain.horizontalTrace, strain.radialPlus, strain.radialMinus,
          strain.horizontalPlus, strain.horizontalMinus};
      Eigen::VectorXcd fdVector(6), strainVector(6);
      for (int c = 0; c < 6; ++c) {
        fdVector[c] = finiteDifference[c];
        strainVector[c] = reconstructed[c];
      }
      const double strainError = (fdVector-strainVector).norm() /
          std::max(fdVector.norm(), strainVector.norm());
      EXPECT_LT(strainError, 1e-8) << "six-component independent strain relative error="
                                   << strainError << " FD=" << fdVector.transpose()
                                   << " reconstructed=" << strainVector.transpose();
    }
  }
}

TEST(Windowed3DFieldEngineTests, PerOrderAdjointMatchesIndependentHermitianSolves) {
  auto parameters = makeSensitivityParameters(4, 0, 1000.0);
  prem_norm<double> norm;
  auto model = EarthModels::ModelInput(DSpecMTest::modelPath().string(), norm);
  auto source = makeSensitivitySource();
  Full1D::SEM sem(model, 0.25, 4, 2, source);
  const auto frequency = makeSensitivityFrequency(model.TimeNorm());
  Eigen::VectorXd window = Eigen::VectorXd::Ones(frequency.nt());
  Eigen::MatrixXcd weights = Eigen::MatrixXcd::Zero(
      3 * parameters.num_receivers(), frequency.nt()/2 + 1);
  weights(0, 3) = {0.37, -0.61};
  weights(1, 3) = {-0.82, 0.29};
  weights(2, 3) = {0.43, 0.77};
  const int bin = 3;
  const std::complex<double> z = frequency.w().at(bin) +
      std::complex<double>(0.0, -frequency.ep());
  const std::complex<double> fz = -std::complex<double>(0.0, 1.0)/z;
  SPARSESPEC::WindowedDisplacementMap map(frequency, sem.lengthScaleMetres(), window);
  const Eigen::VectorXcd qhat = map.adjoint(weights).col(bin);
  const double deg = EIGEN_PI/180.0;
  const auto receiver = parameters.receivers().front();
  const double theta = (90.0-receiver.first)*deg, phi = receiver.second*deg;
  const double rSI = (sem.mesh().PR() - 1000.0*parameters.receiver_depth()/
      sem.lengthScaleMetres())*sem.lengthScaleMetres();
  const SPARSESPEC::Field3DPoint point{rSI, theta, phi};
  const auto reconstructed = SPARSESPEC::Windowed3DFieldEngine::evaluateFrequency(
      sem, source, parameters, z, qhat, {point});

  const auto makeOperator = [&](const Eigen::SparseMatrix<double> &stiffness,
                                const Eigen::SparseMatrix<double> &mass) {
    Eigen::SparseMatrix<std::complex<double>> h =
        stiffness.cast<std::complex<double>>() - z*z*mass.cast<std::complex<double>>();
    h.makeCompressed();
    return h;
  };
  const auto directAdjoint = [&](const Eigen::SparseMatrix<std::complex<double>> &h,
                                 const Eigen::MatrixXcd &receiverBase,
                                 const Eigen::MatrixXcd &angularRows,
                                 int m, int l) {
    Eigen::MatrixXcd rMode = receiverBase;
    for (int row = 0; row < rMode.rows(); ++row)
      rMode.row(row) *= angularRows(row, m+l);
    const Eigen::VectorXcd rhs = rMode.adjoint()*(std::conj(fz)*qhat);
    Eigen::SparseMatrix<std::complex<double>> hAdj = h.adjoint();
    hAdj.makeCompressed();
    Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>> solver;
    solver.analyzePattern(hAdj);
    solver.factorize(hAdj);
    EXPECT_EQ(solver.info(), Eigen::Success);
    const Eigen::VectorXcd lambda = solver.solve(rhs);
    EXPECT_EQ(solver.info(), Eigen::Success);
    return lambda;
  };
  const auto fullReceiver = [&](const Eigen::MatrixXcd &local, int low, int ndof) {
    Eigen::MatrixXcd full = Eigen::MatrixXcd::Zero(local.rows(), ndof);
    full.block(0, low, local.rows(), local.cols()) = local;
    return full;
  };
  Eigen::Vector3cd oracle = Eigen::Vector3cd::Zero();
  const Eigen::RowVectorXcd radialReceiver = sem.rvZR(parameters, 0).transpose();
  Eigen::MatrixXcd rR = Eigen::MatrixXcd::Zero(3, sem.hR().rows());
  rR.row(0) = radialReceiver;
  // The radial receiver has only scalar Y00; form it from the pinned Wigner table.
  GSHTrans::Wigner<double, GSHTrans::All, GSHTrans::All,
                   GSHTrans::Single> w0(0, 0, 0, theta);
  const auto y00 = w0[0][0, 0];
  Eigen::MatrixXcd yR = Eigen::MatrixXcd::Zero(3, 1);
  yR(0,0) = 1.0;
  const auto aR = directAdjoint(makeOperator(sem.hR(), sem.pR()), rR, yR, 0, 0);
  const auto eR = sem.receiverElements(parameters).front();
  std::vector<double> nodes(sem.mesh().NN());
  for (int q = 0; q < sem.mesh().NN(); ++q) nodes[q] = sem.mesh().NodeRadius(eR, q);
  const auto basis = Interpolation::LagrangeBasis(nodes);
  for (int q = 0; q < sem.mesh().NN(); ++q) {
    const auto i = sem.ltgR(0, eR, q);
    oracle[0] += basis.Evaluate(q, rSI/sem.lengthScaleMetres())*aR[i]*y00;
  }

  ParamInfo angular(parameters, parameters.lmax());
  for (int l = 1; l <= parameters.lmax(); ++l) {
    const double zeta = std::sqrt(l*(l+1.0));
    const auto sphElements = sem.receiverElements(parameters, false);
    const int lowS = sem.ltgS(0, sphElements.front(), 0);
    const Eigen::MatrixXcd rS = fullReceiver(sem.rvBaseFull(parameters, l, false),
        lowS, sem.hS(l).rows());
    const Eigen::MatrixXcd yS = angular.rvFullSph(l);
    const auto hS = makeOperator(sem.hS(l), sem.pS(l));
    const auto wignerS = GSHTrans::Wigner<double, GSHTrans::All, GSHTrans::All,
        GSHTrans::Single>(l, l, std::min(2,l), theta);
    for (int m = -l; m <= l; ++m) {
      const auto lambda = directAdjoint(hS, rS, yS, m, l);
      std::complex<double> U{}, V{};
      for (int q = 0; q < sem.mesh().NN(); ++q) {
        const auto iu = sem.ltgS(0, eR, q), iv = sem.ltgS(1, eR, q);
        const double shape = basis.Evaluate(q, rSI/sem.lengthScaleMetres());
        U += shape*lambda[iu]; V += shape*lambda[iv];
      }
      const auto phase = std::exp(std::complex<double>(0.0, m*phi));
      const auto x0 = wignerS[0][l,m]*phase;
      const auto xp = wignerS[1][l,m]*phase;
      const auto xm = wignerS[-1][l,m]*phase;
      oracle[0] += U*x0;
      oracle[1] += zeta/2.0*V*(xp-xm);
      oracle[2] += -std::complex<double>(0,1)*zeta/2.0*V*(xp+xm);
    }

    const auto torElements = sem.receiverElements(parameters, true);
    const int lowT = sem.ltgT(torElements.front(), 0);
    const Eigen::MatrixXcd rT = fullReceiver(sem.rvBaseFullT(parameters, l),
        lowT, sem.hTk(l).rows());
    const Eigen::MatrixXcd yT = angular.rvFullTor(l);
    const auto hT = makeOperator(sem.hTk(l), sem.pTk(l));
    const auto wignerT = GSHTrans::Wigner<double, GSHTrans::All, GSHTrans::All,
        GSHTrans::Single>(l, l, std::min(2,l), theta);
    for (int m = -l; m <= l; ++m) {
      const auto lambda = directAdjoint(hT, rT, yT, m, l);
      std::complex<double> W{};
      for (int q = 0; q < sem.mesh().NN(); ++q) {
        const auto iw = sem.ltgT(eR, q);
        W += basis.Evaluate(q, rSI/sem.lengthScaleMetres())*lambda[iw];
      }
      const auto phase = std::exp(std::complex<double>(0.0, m*phi));
      const auto xp = wignerT[1][l,m]*phase;
      const auto xm = wignerT[-1][l,m]*phase;
      oracle[1] += std::complex<double>(0,1)*zeta/2.0*W*(xp+xm);
      oracle[2] += zeta/2.0*W*(xp-xm);
    }
  }
  Eigen::Vector3cd actual;
  for (int c = 0; c < 3; ++c)
    actual[c] = reconstructed.front().adjointDualDisplacement[c];
  const double error = (actual-oracle).norm()/std::max(actual.norm(),oracle.norm());
  EXPECT_LT(error, 1e-9) << "per-order direct H^dagger reconstruction error=" << error;

  SourceInfo::EarthquakeCMT alternate(15.0, 110.0, 30.0,
                                       -3.0, 1.0, 4.0, 2.0, -5.0, 1.0);
  const auto alternateFields = SPARSESPEC::Windowed3DFieldEngine::evaluateFrequency(
      sem, alternate, parameters, z, qhat, {point});
  Eigen::Vector3cd alternateAdjoint, forwardDifference;
  for (int c = 0; c < 3; ++c) {
    alternateAdjoint[c] = alternateFields.front().adjointDualDisplacement[c];
    forwardDifference[c] = alternateFields.front().forwardRawDisplacement[c] -
                           reconstructed.front().forwardRawDisplacement[c];
  }
  EXPECT_LT((alternateAdjoint-actual).norm()/actual.norm(), 1e-12);
  EXPECT_GT(forwardDifference.norm(), 0.0);
}

TEST(Windowed3DFieldEngineTests, CanonicalizesSIEndpointRoundTripsAndRejectsExteriorPoints) {
  auto parameters = makeSensitivityParameters(1, 0, 1000.0);
  EarthModels::EarthConstants<double> norm;
  auto model = EarthModels::ModelInput(DSpecMTest::modelPath().string(), norm);
  auto source = makeSensitivitySource();
  Full1D::SEM sem(model, 0.01, 4, 2, source);
  const auto frequency = makeSensitivityFrequency(model.TimeNorm());
  const int bin = 3;
  const std::complex<double> z = frequency.w().at(bin) +
      std::complex<double>(0.0, -frequency.ep());
  const Eigen::VectorXcd qhat = Eigen::VectorXcd::Zero(
      3*parameters.num_receivers());
  const double scale = sem.lengthScaleMetres();
  const double theta = 0.63, phi = 1.17;
  std::vector<SPARSESPEC::Field3DPoint> points;
  std::vector<int> expectedElements;
  int endpointsRejectedByStrictBounds = 0;
  int defaultEndpointsRejectedByStrictSelection = 0;
  const auto strictDefaultElement = [&](double radiusModel) {
    int selected = -1;
    for (int e = 0; e < sem.mesh().NE(); ++e)
      if (sem.mesh().ELR(e) <= radiusModel &&
          (radiusModel < sem.mesh().EUR(e) ||
           (e == sem.mesh().NE()-1 && radiusModel <= sem.mesh().EUR(e))))
        selected = e;
    return selected;
  };
  for (int outer = 1; outer < sem.mesh().NE(); ++outer) {
    // SI conversion and round-trip deliberately exercise endpoint ulps.
    const double radiusSI = sem.mesh().ELR(outer)*scale;
    const double roundTripRadius = radiusSI/scale;
    const int innerSide = outer-1;
    const int outerSide = outer;
    endpointsRejectedByStrictBounds +=
        (roundTripRadius < sem.mesh().ELR(innerSide) ||
         roundTripRadius > sem.mesh().EUR(innerSide));
    endpointsRejectedByStrictBounds +=
        (roundTripRadius < sem.mesh().ELR(outerSide) ||
         roundTripRadius > sem.mesh().EUR(outerSide));
    points.push_back({radiusSI, theta, phi, innerSide});
    expectedElements.push_back(outer-1);
    points.push_back({radiusSI, theta, phi, outerSide});
    expectedElements.push_back(outer);
    const double innerUpperSI = sem.mesh().EUR(innerSide)*scale;
    const double innerUpperRoundTrip = innerUpperSI/scale;
    defaultEndpointsRejectedByStrictSelection +=
        strictDefaultElement(innerUpperRoundTrip) != outer;
    points.push_back({innerUpperSI, theta, phi});
    expectedElements.push_back(outer);
    const double outerLowerRoundTrip = radiusSI/scale;
    defaultEndpointsRejectedByStrictSelection +=
        strictDefaultElement(outerLowerRoundTrip) != outer;
    points.push_back({radiusSI, theta, phi});
    expectedElements.push_back(outer);
  }
  ASSERT_GT(sem.mesh().NE(), 1);
  const int defaultOuter = 1;
  points.push_back({sem.mesh().ELR(defaultOuter)*scale, theta, phi});
  expectedElements.push_back(defaultOuter);
  const int surfaceElement = sem.mesh().NE()-1;
  points.push_back({sem.mesh().EUR(surfaceElement)*scale, theta, phi,
                     surfaceElement});
  expectedElements.push_back(surfaceElement);
  points.push_back({sem.mesh().EUR(surfaceElement)*scale, theta, phi});
  expectedElements.push_back(surfaceElement);
  EXPECT_GT(endpointsRejectedByStrictBounds, 0)
      << "Fixture must exercise endpoint round trips rejected by the old strict bounds.";
  EXPECT_GT(defaultEndpointsRejectedByStrictSelection, 0)
      << "Fixture must exercise default endpoints rejected by strict selection.";
  std::vector<SPARSESPEC::Field3DSample> samples;
  samples.reserve(points.size());
  for (std::size_t i = 0; i < points.size(); ++i) {
    SCOPED_TRACE(::testing::Message() << "endpoint sample " << i
        << " selected element " << expectedElements[i]
        << " radius SI " << points[i].radiusMetres);
    std::vector<SPARSESPEC::Field3DSample> one;
    try {
      one = SPARSESPEC::Windowed3DFieldEngine::evaluateFrequency(
          sem, source, parameters, z, qhat, {points[i]});
    } catch (const std::exception &error) {
      const int e = expectedElements[i];
      ADD_FAILURE() << error.what() << "; expected element=" << e
                    << " requested element=" << points[i].radialElement
                    << " radiusModel=" << points[i].radiusMetres/scale
                    << " bounds=[" << sem.mesh().ELR(e) << ", "
                    << sem.mesh().EUR(e) << "]";
      return;
    }
    ASSERT_EQ(one.size(), 1u);
    EXPECT_EQ(one.front().radialElement, expectedElements[i]);
    samples.push_back(one.front());
  }

  const double outsideRadius = (sem.mesh().EUR(sem.mesh().NE()-1) + 1.0e-5)*scale;
  EXPECT_THROW(SPARSESPEC::Windowed3DFieldEngine::evaluateFrequency(
      sem, source, parameters, z, qhat, {{outsideRadius, theta, phi}}),
      std::invalid_argument);
}

namespace {

enum class TestModulus { A, C, F, L, N };

class PerturbedModulusModel : public EarthModels::ModelInput<double> {
  using Base = EarthModels::ModelInput<double>;
  template <class Field>
  auto perturbed(Field field, TestModulus modulus, int layer) const {
    const double amplitude = modulus == modulus_ && !Base::IsFluid(layer)
                                ? amplitude_ : 0.0;
    const double pressureScale = Base::StressNorm();
    const double radiusScale = Base::OuterRadius();
    return [field = std::move(field), amplitude, pressureScale, radiusScale](double radius) {
      return field(radius) + amplitude * (1e11 / pressureScale) *
             (1.0 + 0.3 * radius / radiusScale);
    };
  }
public:
  template <class Norm>
  PerturbedModulusModel(const std::string &path, const Norm &norm,
                        TestModulus modulus, double amplitude)
      : Base(path, norm), modulus_(modulus), amplitude_(amplitude) {}

  auto A(int layer) const { return perturbed(Base::A(layer), TestModulus::A, layer); }
  auto C(int layer) const { return perturbed(Base::C(layer), TestModulus::C, layer); }
  auto F(int layer) const { return perturbed(Base::F(layer), TestModulus::F, layer); }
  auto L(int layer) const { return perturbed(Base::L(layer), TestModulus::L, layer); }
  auto N(int layer) const { return perturbed(Base::N(layer), TestModulus::N, layer); }

private:
  TestModulus modulus_;
  double amplitude_;
};

class PerturbedDensityModel : public EarthModels::ModelInput<double> {
  using Base = EarthModels::ModelInput<double>;
public:
  template <class Norm>
  PerturbedDensityModel(const std::string &path, const Norm &norm, double amplitude,
                        double lowerRadius, double upperRadius, int targetLayer = -1)
      : Base(path, norm), amplitude_(amplitude), lower_(lowerRadius), upper_(upperRadius),
        layer_(targetLayer) {}

  auto Density(int layer) const {
    auto original = Base::Density(layer);
    const bool active = layer_ < 0 || layer_ == layer;
    const double increment = amplitude_ * (1000.0 / Base::DensityNorm());
    const double lower = lower_, upper = upper_;
    const double scale = Base::OuterRadius();
    return [original = std::move(original), active, increment, lower, upper, scale](double r) {
      if (!active || (lower >= 0.0 && !(r >= lower && r <= upper))) return original(r);
      const double profile = lower < 0.0 ? 0.7 + 0.3 * r / scale :
          4.0 * (r - lower) * (upper - r) / ((upper - lower) * (upper - lower));
      return original(r) + increment * profile;
    };
  }

private:
  double amplitude_, lower_, upper_;
  int layer_;
};

class PerturbedFluidBulkModel : public EarthModels::ModelInput<double> {
  using Base = EarthModels::ModelInput<double>;
  template <class Field>
  auto perturbed(Field field, int layer) const {
    const bool active = Base::IsFluid(layer) && layer == layer_;
    const double increment = active ? amplitude_ * (1e11 / Base::StressNorm()) : 0.0;
    const double radius = Base::OuterRadius();
    const double lower = lower_, upper = upper_;
    return [field = std::move(field), increment, radius, lower, upper](double r) {
      if (increment == 0.0 || r < lower || r > upper) return field(r);
      const double x = (2.0 * r - lower - upper) / (upper - lower);
      return field(r) + increment * (1.0 - x * x);
    };
  }
public:
  template <class Norm>
  PerturbedFluidBulkModel(const std::string &path, const Norm &norm, double amplitude,
                          double lower, double upper, int layer)
      : Base(path, norm), amplitude_(amplitude), lower_(lower), upper_(upper), layer_(layer) {}
  auto A(int layer) const { return perturbed(Base::A(layer), layer); }
  auto C(int layer) const { return perturbed(Base::C(layer), layer); }
  auto F(int layer) const { return perturbed(Base::F(layer), layer); }
private:
  double amplitude_, lower_, upper_;
  int layer_;
};

Eigen::MatrixXd independentWeakFormDelta(
    const Full1D::SEM &sem, SPARSESPEC::SensitivitySolveFamily family,
    int degree, TestModulus modulus, const Eigen::MatrixXd &delta,
    bool fluidDomain = false) {
  const auto &mesh = sem.mesh();
  const int nq = mesh.NN();
  int dofs = family == SPARSESPEC::SensitivitySolveFamily::Radial
                 ? sem.pR().rows()
                 : (family == SPARSESPEC::SensitivitySolveFamily::Spheroidal
                        ? sem.pS(degree).rows() : sem.pTk(degree).rows());
  Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(dofs, dofs);
  const auto gll = const_cast<EarthMesh::RadialMesh &>(mesh).GLL();
  const auto basis = Interpolation::LagrangeBasis(gll.Points());
  const double lambda = degree * (degree + 1.0);
  for (int e = 0; e < mesh.NE(); ++e) {
    if ((sem.elementIsFluid(e) != fluidDomain) ||
        (family == SPARSESPEC::SensitivitySolveFamily::Toroidal &&
         (e < sem.el() || e >= sem.eu())))
      continue;
    const double scale = 2.0 / mesh.EW(e);
    const double jacobian = 0.5 * mesh.EW(e);
    for (int q = 0; q < nq; ++q) {
      const double r = mesh.NodeRadius(e, q);
      const double weight = jacobian * gll.W(q) * delta(e, q);
      for (int i = 0; i < nq; ++i) {
        for (int j = 0; j < nq; ++j) {
          const double Li = basis.Evaluate(i, gll.X(q));
          const double Lj = basis.Evaluate(j, gll.X(q));
          const double Di = scale * basis.template Evaluate<1>(i, gll.X(q));
          const double Dj = scale * basis.template Evaluate<1>(j, gll.X(q));
          if (family == SPARSESPEC::SensitivitySolveFamily::Toroidal) {
            const int ti = sem.ltgT(e, i), tj = sem.ltgT(e, j);
            const double strainI = r * Di - Li;
            const double strainJ = r * Dj - Lj;
            double value = 0.0;
            if (modulus == TestModulus::L)
              value = lambda * strainI * strainJ;
            if (modulus == TestModulus::N)
              value = lambda * (lambda - 2.0) * Li * Lj;
            matrix(ti, tj) += weight * value;
          } else {
            const bool radial = family == SPARSESPEC::SensitivitySolveFamily::Radial;
            const std::array<std::size_t, 2> testDof{
                static_cast<std::size_t>(radial ? sem.ltgR(0, e, i) : sem.ltgS(0, e, i)),
                radial ? 0U : static_cast<std::size_t>(sem.ltgS(1, e, i))};
            const std::array<std::size_t, 2> trialDof{
                static_cast<std::size_t>(radial ? sem.ltgR(0, e, j) : sem.ltgS(0, e, j)),
                radial ? 0U : static_cast<std::size_t>(sem.ltgS(1, e, j))};
            for (int a = 0; a < (radial ? 1 : 2); ++a) {
              for (int b = 0; b < (radial ? 1 : 2); ++b) {
                const double ut = a == 0 ? Li : 0.0;
                const double vt = a == 1 ? Li : 0.0;
                const double dut = a == 0 ? Di : 0.0;
                const double dvt = a == 1 ? Di : 0.0;
                const double uu = b == 0 ? Lj : 0.0;
                const double vu = b == 1 ? Lj : 0.0;
                const double duu = b == 0 ? Dj : 0.0;
                const double dvu = b == 1 ? Dj : 0.0;
                const double Ttest = 2.0 * ut - lambda * vt;
                const double Ttrial = 2.0 * uu - lambda * vu;
                const double Stest = r * dvt + ut - vt;
                const double Strial = r * dvu + uu - vu;
                double value = 0.0;
                if (modulus == TestModulus::C) value = r * r * dut * duu;
                if (modulus == TestModulus::F)
                  value = r * (duu * Ttest + Ttrial * dut);
                if (modulus == TestModulus::A) value = Ttest * Ttrial;
                if (modulus == TestModulus::N)
                  value = -2.0 * (2.0 * uu * ut - lambda * (uu * vt + vu * ut) +
                                  lambda * vu * vt);
                if (modulus == TestModulus::L) value = lambda * Strial * Stest;
                matrix(testDof[a], trialDof[b]) += weight * value;
              }
            }
          }
        }
      }
    }
  }
  return matrix;
}

} // namespace

TEST(ModulusKernelTests, MatchIndependentWeakFormMatricesForAllFiveModuli) {
  prem_norm<double> norm;
  auto model = EarthModels::ModelInput(DSpecMTest::modelPath().string(), norm);
  Full1D::SEM sem(model, 0.4, 4, 2, makeSensitivitySource());
  Eigen::MatrixXd delta(sem.mesh().NE(), sem.mesh().NN());
  for (int e = 0; e < delta.rows(); ++e)
    for (int i = 0; i < delta.cols(); ++i)
      delta(e, i) = 0.07 * (1.0 + 0.3 * sem.mesh().NodeRadius(e, i));
  const std::array<TestModulus, 5> materials{
      TestModulus::A, TestModulus::C, TestModulus::F,
      TestModulus::L, TestModulus::N};
  for (const auto family : {SPARSESPEC::SensitivitySolveFamily::Radial,
                            SPARSESPEC::SensitivitySolveFamily::Spheroidal,
                            SPARSESPEC::SensitivitySolveFamily::Toroidal}) {
    const int degree = family == SPARSESPEC::SensitivitySolveFamily::Radial ? 0 : 2;
    const int dofs = family == SPARSESPEC::SensitivitySolveFamily::Radial
                         ? sem.pR().rows()
                         : family == SPARSESPEC::SensitivitySolveFamily::Spheroidal
                               ? sem.pS(degree).rows() : sem.pTk(degree).rows();
    SPARSESPEC::FrequencyAdjointFields fields;
    fields.forward = Eigen::MatrixXcd::Random(dofs, 2);
    fields.adjoint = Eigen::MatrixXcd::Random(dofs, 2);
    const auto samples = SPARSESPEC::modulusKernelSamples(sem, family, degree, fields);
    for (const auto material : materials) {
      const auto matrix = independentWeakFormDelta(sem, family, degree, material, delta);
      const double matrixDerivative = -
          (fields.adjoint.adjoint() * matrix.cast<std::complex<double>>() *
           fields.forward).trace().real();
      const Eigen::MatrixXd *kernel = nullptr;
      switch (material) {
      case TestModulus::A: kernel = &samples.A; break;
      case TestModulus::C: kernel = &samples.C; break;
      case TestModulus::F: kernel = &samples.F; break;
      case TestModulus::L: kernel = &samples.L; break;
      case TestModulus::N: kernel = &samples.N; break;
      }
      double integratedKernel = 0.0;
      for (int e = 0; e < delta.rows(); ++e)
        for (int i = 0; i < delta.cols(); ++i)
          integratedKernel += (*kernel)(e, i) * delta(e, i) *
              sem.pressureScalePascals() * samples.drWeightsMetres(e, i);
      const double relative = std::abs(matrixDerivative - integratedKernel) /
          std::max({std::abs(matrixDerivative), std::abs(integratedKernel), 1e-30});
      EXPECT_LE(relative, 1e-12) << static_cast<int>(family) << " "
                                 << static_cast<int>(material);
    }
  }
}

TEST(ModulusKernelTests, ProductionOperatorDeltasMatchIndependentWeakForm) {
  prem_norm<double> norm;
  const std::string path = DSpecMTest::modelPath().string();
  auto baselineModel = EarthModels::ModelInput<double>(path, norm);
  Full1D::SEM baseline(baselineModel, 0.4, 4, 2, makeSensitivitySource());
  Eigen::MatrixXd delta(baseline.mesh().NE(), baseline.mesh().NN());
  for (int e = 0; e < delta.rows(); ++e)
    for (int i = 0; i < delta.cols(); ++i)
      delta(e, i) = (1e11 / baselineModel.StressNorm()) *
                    (1.0 + 0.3 * baseline.mesh().NodeRadius(e, i) /
                               baseline.mesh().PR());
  const std::array<TestModulus, 5> materials{
      TestModulus::A, TestModulus::C, TestModulus::F,
      TestModulus::L, TestModulus::N};
  // The operator is linear in the modulus samples; the large test-only step
  // keeps sparse-matrix subtraction well above roundoff. No solve uses it.
  // Large algebra-only step stabilizes subtraction of the linear stiffness delta.
  constexpr double amplitude = 1e6;
  for (const auto material : materials) {
    PerturbedModulusModel changedModel(path, norm, material, amplitude);
    Full1D::SEM changed(changedModel, 0.4, 4, 2, makeSensitivitySource());
    for (const auto family : {SPARSESPEC::SensitivitySolveFamily::Radial,
                              SPARSESPEC::SensitivitySolveFamily::Spheroidal,
                              SPARSESPEC::SensitivitySolveFamily::Toroidal}) {
      const int degree = family == SPARSESPEC::SensitivitySolveFamily::Radial ? 0 : 2;
      Eigen::MatrixXd actual;
      if (family == SPARSESPEC::SensitivitySolveFamily::Radial)
        actual = (changed.hR() - baseline.hR()) / amplitude;
      else if (family == SPARSESPEC::SensitivitySolveFamily::Spheroidal)
        actual = (changed.hS(degree) - baseline.hS(degree)) / amplitude;
      else
        actual = (changed.hTk(degree) - baseline.hTk(degree)) / amplitude;
      const Eigen::MatrixXd expected = independentWeakFormDelta(
          baseline, family, degree, material, delta);
      const double relative = (actual - expected).norm() /
          std::max({actual.norm(), expected.norm(), 1e-30});
      EXPECT_LE(relative, 1e-12) << static_cast<int>(family) << " "
                                << static_cast<int>(material) << " actual="
                                << actual.norm() << " expected=" << expected.norm()
                                << " projection=" << actual.cwiseProduct(expected).sum() /
                                       std::max(expected.squaredNorm(), 1e-30);
    }
  }
}

TEST(MaterialKernelTests, DensityGradientMatchesProductionOperatorForAllFamilies) {
  EarthModels::EarthConstants<double> norm;
  const std::string path = DSpecMTest::modelPath().string();
  auto baselineModel = EarthModels::ModelInput<double>(path, norm);
  // Deliberately large for matrix-subtraction precision; these algebra-only
  // centered checks do not solve with the unphysical perturbed density models.
  constexpr double amplitude = 1e4;
  const std::complex<double> z{0.37, -0.015};
  double maxRelative = 0.0;
  for (int nq : {2, 3, 4}) {
   Full1D::SEM baseline(baselineModel, 0.4, nq, 2);
   for (const auto family : {SPARSESPEC::SensitivitySolveFamily::Radial,
                             SPARSESPEC::SensitivitySolveFamily::Spheroidal,
                             SPARSESPEC::SensitivitySolveFamily::Toroidal}) {
    const int degree = family == SPARSESPEC::SensitivitySolveFamily::Radial ? 0 : 2;
    const int dofs = family == SPARSESPEC::SensitivitySolveFamily::Radial
        ? baseline.pR().rows() : family == SPARSESPEC::SensitivitySolveFamily::Spheroidal
            ? baseline.pS(degree).rows() : baseline.pTk(degree).rows();
    SPARSESPEC::FrequencyAdjointFields fields;
    fields.forward = Eigen::MatrixXcd::Random(dofs, 2);
    fields.adjoint = Eigen::MatrixXcd::Random(dofs, 2);
    const auto samples = SPARSESPEC::materialKernelSamples(
        baseline, family, degree, z, fields);
    EXPECT_TRUE(samples.densityKgPerM3.allFinite());
    EXPECT_TRUE(samples.bulkModulusPa.allFinite());
    EXPECT_EQ(samples.layer.size(), static_cast<std::size_t>(baseline.mesh().NE()));
    for (int e = 0; e < baseline.mesh().NE(); ++e)
      for (int i = 0; i < baseline.mesh().NN(); ++i)
        if (!samples.fluid(e, i)) EXPECT_EQ(samples.bulkModulusPa(e, i), 0.0);
    Eigen::MatrixXd delta(baseline.mesh().NE(), baseline.mesh().NN());
    for (int e = 0; e < delta.rows(); ++e)
      for (int i = 0; i < delta.cols(); ++i)
        delta(e, i) = 1000.0 * (0.7 + 0.3 * baseline.mesh().NodeRadius(e, i) /
                                         baseline.mesh().PR());
    double kernel = 0.0;
    for (int e = 0; e < delta.rows(); ++e)
      for (int i = 0; i < delta.cols(); ++i)
        kernel += samples.densityKgPerM3(e, i) * delta(e, i) * samples.drWeightsMetres(e, i);

    PerturbedDensityModel plusModel(path, norm, amplitude, -1.0, -1.0);
    Full1D::SEM plus(plusModel, 0.4, nq, 2);
    PerturbedDensityModel minusModel(path, norm, -amplitude, -1.0, -1.0);
    Full1D::SEM minus(minusModel, 0.4, nq, 2);
    const auto matrixFor = [&](const Full1D::SEM &sem) {
      Eigen::SparseMatrix<std::complex<double>> h;
      if (family == SPARSESPEC::SensitivitySolveFamily::Radial)
        h = sem.hR().cast<std::complex<double>>() - z * z * sem.pR().cast<std::complex<double>>();
      else if (family == SPARSESPEC::SensitivitySolveFamily::Spheroidal)
        h = sem.hS(degree).cast<std::complex<double>>() - z * z * sem.pS(degree).cast<std::complex<double>>();
      else
        h = sem.hTk(degree).cast<std::complex<double>>() - z * z * sem.pTk(degree).cast<std::complex<double>>();
      return h;
    };
    const Eigen::MatrixXcd dH = (matrixFor(plus) - matrixFor(minus)) / (2.0 * amplitude);
    const double matrixDerivative = -
        (fields.adjoint.conjugate().cwiseProduct(dH * fields.forward)).sum().real();
    const double relative = std::abs(matrixDerivative - kernel) /
        std::max({std::abs(matrixDerivative), std::abs(kernel), 1e-30});
    maxRelative = std::max(maxRelative, relative);
    EXPECT_LE(relative, 1e-8) << "nq=" << nq << " family=" << static_cast<int>(family)
        << " matrix=" << matrixDerivative << " kernel=" << kernel;
   }
  }
  ::testing::Test::RecordProperty("max_density_operator_contraction_relative_error", maxRelative);
}

TEST(MaterialKernelTests, FluidBulkKernelMatchesAPlusCPlusFProductionDelta) {
  EarthModels::EarthConstants<double> norm;
  const std::string path = DSpecMTest::modelPath().string();
  auto baselineModel = EarthModels::ModelInput<double>(path, norm);
  Full1D::SEM baseline(baselineModel, 0.4, 4, 2, makeSensitivitySource());
  int fluidElement = -1;
  for (int e = 0; e < baseline.mesh().NE(); ++e)
    if (baseline.elementIsFluid(e)) { fluidElement = e; break; }
  ASSERT_GE(fluidElement, 0);
  const double lower = baseline.mesh().ELR(fluidElement), upper = baseline.mesh().EUR(fluidElement);
  const int layer = baseline.mesh().LayerNumber(fluidElement);
  // Large algebra-only step stabilizes subtraction; no wave solve uses it.
  constexpr double amplitude = 1e6;
  PerturbedFluidBulkModel changedModel(path, norm, amplitude, lower, upper, layer);
  Full1D::SEM changed(changedModel, 0.4, 4, 2, makeSensitivitySource());
  Eigen::MatrixXd delta(baseline.mesh().NE(), baseline.mesh().NN());
  delta.setZero();
  for (int i = 0; i < delta.cols(); ++i) {
    const double r = baseline.mesh().NodeRadius(fluidElement, i);
    const double x = (2.0 * r - lower - upper) / (upper - lower);
    delta(fluidElement, i) = (1e11 / baselineModel.StressNorm()) * (1.0 - x * x);
  }
  const Eigen::MatrixXd expected =
      independentWeakFormDelta(baseline, SPARSESPEC::SensitivitySolveFamily::Spheroidal,
                               2, TestModulus::A, delta, true) +
      independentWeakFormDelta(baseline, SPARSESPEC::SensitivitySolveFamily::Spheroidal,
                               2, TestModulus::C, delta, true) +
      independentWeakFormDelta(baseline, SPARSESPEC::SensitivitySolveFamily::Spheroidal,
                               2, TestModulus::F, delta, true);
  const Eigen::MatrixXd actual = (changed.hS(2) - baseline.hS(2)) / amplitude;
  const double relative = (actual - expected).norm() /
      std::max({actual.norm(), expected.norm(), 1e-30});
  EXPECT_LE(relative, 1e-12);
  EXPECT_GT(expected.norm(), 0.0);
  SPARSESPEC::FrequencyAdjointFields fields;
  fields.forward = Eigen::MatrixXcd::Random(baseline.hS(2).rows(), 2);
  fields.adjoint = Eigen::MatrixXcd::Random(baseline.hS(2).rows(), 2);
  const auto samples = SPARSESPEC::materialKernelSamples(
      baseline, SPARSESPEC::SensitivitySolveFamily::Spheroidal, 2, {0.3, -0.01}, fields);
  double integrated = 0.0;
  for (int e = 0; e < delta.rows(); ++e)
    for (int i = 0; i < delta.cols(); ++i)
      if (e == fluidElement)
        integrated += samples.bulkModulusPa(e, i) *
            1e11 * (1.0 - std::pow((2.0 * baseline.mesh().NodeRadius(e, i) - lower - upper) /
                                   (upper - lower), 2)) *
            samples.drWeightsMetres(e, i);
  const Eigen::MatrixXcd dH =
      ((changed.hS(2) - baseline.hS(2)) / amplitude).cast<std::complex<double>>();
  const double contraction = -
      (fields.adjoint.conjugate().cwiseProduct(dH * fields.forward)).sum().real();
  const double kernelRelative = std::abs(integrated - contraction) /
      std::max({std::abs(integrated), std::abs(contraction), 1e-30});
  ::testing::Test::RecordProperty("fluid_bulk_kernel_contraction_relative_error", kernelRelative);
  EXPECT_LE(kernelRelative, 1e-12) << "bulk=" << integrated << " matrix=" << contraction;
  for (int e = 0; e < baseline.mesh().NE(); ++e)
    if (!baseline.elementIsFluid(e)) EXPECT_EQ(samples.bulkModulusPa.row(e).norm(), 0.0);
}

TEST(MaterialKernelTests, OneSidedPolynomialTransposeAndShellGravityResponse) {
  EarthModels::EarthConstants<double> norm;
  auto model = EarthModels::ModelInput<double>(DSpecMTest::modelPath().string(), norm);
  Full1D::SEM sem(model, 0.4, 4, 2, makeSensitivitySource());
  constexpr int degree = 2, nq = 4;
  const std::complex<double> z{0.33, -0.017};
  const int dofs = static_cast<int>(sem.hS(degree).rows());
  SPARSESPEC::FrequencyAdjointFields fields;
  fields.forward = Eigen::MatrixXcd::Random(dofs, 2);
  fields.adjoint = Eigen::MatrixXcd::Random(dofs, 2);
  const auto kernels = SPARSESPEC::materialKernelSamples(
      sem, SPARSESPEC::SensitivitySolveFamily::Spheroidal, degree, z, fields);
  Eigen::MatrixXd delta(sem.mesh().NE(), nq);
  for (int e = 0; e < delta.rows(); ++e)
    for (int i = 0; i < nq; ++i)
      delta(e, i) = (i % 2 ? -1.0 : 1.0) * (0.2 + 0.013 * e) + 0.017 * i * i * i;
  const Eigen::MatrixXd deltaGravity = directPerturbationGravity(sem, delta);
  const auto gll = const_cast<EarthMesh::RadialMesh &>(sem.mesh()).GLL();
  const auto basis = Interpolation::LagrangeBasis(gll.Points());
  Eigen::MatrixXd D(nq, nq);
  for (int i = 0; i < nq; ++i)
    for (int j = 0; j < nq; ++j)
      D(i, j) = basis.template Evaluate<1>(j, gll.X(i));
  const double lambda = degree * (degree + 1.0);
  double direct = 0.0, predicted = 0.0;
  const double rhoScale = sem.pressureScalePascals() * std::pow(sem.timeScaleSeconds(), 2) /
                          std::pow(sem.lengthScaleMetres(), 2);
  for (int e = 0; e < sem.mesh().NE(); ++e) {
    Eigen::MatrixXcd u(nq, 2), ua(nq, 2), v(nq, 2), va(nq, 2), phi(nq, 2), phia(nq, 2);
    for (int i = 0; i < nq; ++i) {
      const auto iu = sem.ltgS(0, e, i), iv = sem.ltgS(1, e, i), ip = sem.ltgS(2, e, i);
      u.row(i) = fields.forward.row(iu); ua.row(i) = fields.adjoint.row(iu);
      v.row(i) = fields.forward.row(iv); va.row(i) = fields.adjoint.row(iv);
      phi.row(i) = fields.forward.row(ip); phia.row(i) = fields.adjoint.row(ip);
    }
    const Eigen::MatrixXcd du = (2.0 / sem.mesh().EW(e)) * D.cast<std::complex<double>>() * u;
    const Eigen::MatrixXcd dphi = (2.0 / sem.mesh().EW(e)) * D.cast<std::complex<double>>() * phi;
    const Eigen::MatrixXcd dphia = (2.0 / sem.mesh().EW(e)) * D.cast<std::complex<double>>() * phia;
    for (int i = 0; i < nq; ++i) {
      const double r = sem.mesh().NodeRadius(e, i);
      const double qw = 0.5 * sem.mesh().EW(e) * gll.W(i);
      const auto uu = SPARSESPEC::sensitivity_detail::inner(ua.row(i), u.row(i));
      const auto uv = SPARSESPEC::sensitivity_detail::inner(ua.row(i), v.row(i)) +
                      SPARSESPEC::sensitivity_detail::inner(va.row(i), u.row(i));
      const auto vv = SPARSESPEC::sensitivity_detail::inner(va.row(i), v.row(i));
      const double rho = sem.meshModel().Density(e, i), g = sem.meshModel().Gravity(e, i);
      const auto local = (8.0 * static_cast<double>(EIGEN_PI) * sem.gravityCoefficient() *
          rho * r * r - 4.0 * g * r) * uu + lambda * g * r * uv +
          r * r * (SPARSESPEC::sensitivity_detail::inner(dphia.row(i), u.row(i)) +
                   SPARSESPEC::sensitivity_detail::inner(ua.row(i), dphi.row(i))) +
          lambda * r * (SPARSESPEC::sensitivity_detail::inner(phia.row(i), v.row(i)) +
                        SPARSESPEC::sensitivity_detail::inner(va.row(i), phi.row(i))) -
          z * z * r * r * (uu + lambda * vv);
      const auto gravityTerm = -4.0 * rho * r * uu + lambda * rho * r * uv;
      direct -= qw * (local * delta(e, i) + gravityTerm * deltaGravity(e, i)).real();
      predicted += kernels.densityKgPerM3(e, i) * delta(e, i) * rhoScale *
                   kernels.drWeightsMetres(e, i);
    }
  }
  const double transposeError = std::abs(direct - predicted) /
      std::max({std::abs(direct), std::abs(predicted), 1e-30});
  ::testing::Test::RecordProperty("one_sided_polynomial_gravity_transpose_relative_error", transposeError);
  EXPECT_LE(transposeError, 1e-12) << "direct=" << direct << " transpose=" << predicted;

  int shell = -1;
  for (int e = 1; e < sem.mesh().NE(); ++e)
    if (sem.mesh().ELR(e) > 0.2 * sem.mesh().PR() && sem.mesh().EUR(e) < 0.8 * sem.mesh().PR()) {
      shell = e; break;
    }
  ASSERT_GE(shell, 0);
  Eigen::MatrixXd shellDelta = Eigen::MatrixXd::Zero(sem.mesh().NE(), nq);
  shellDelta.row(shell).setOnes();
  const auto shellGravity = directPerturbationGravity(sem, shellDelta);
  const double a = sem.mesh().ELR(shell), b = sem.mesh().EUR(shell);
  double maxShellError = 0.0;
  for (int e = 0; e < sem.mesh().NE(); ++e) for (int i = 0; i < nq; ++i) {
    const double r = sem.mesh().NodeRadius(e, i);
    const double upper = std::min(r, b);
    const double mass = r > a ? (upper > a ? (std::pow(upper, 3) - std::pow(a, 3)) / 3.0 : 0.0) : 0.0;
    const double expected = r > 0.0 ? 4.0 * static_cast<double>(EIGEN_PI) *
        sem.gravityCoefficient() * mass / (r * r) : 0.0;
    maxShellError = std::max(maxShellError, std::abs(shellGravity(e, i) - expected) /
        std::max({std::abs(shellGravity(e, i)), std::abs(expected), 1e-30}));
  }
  ::testing::Test::RecordProperty("constant_shell_gravity_relative_error", maxShellError);
  EXPECT_LE(maxShellError, 1e-12);
  EXPECT_EQ(shellGravity(0, 0), 0.0);
}

namespace {

const Eigen::MatrixXd &kernelFor(const SPARSESPEC::ModulusKernelSamples &samples,
                                 TestModulus modulus) {
  switch (modulus) {
  case TestModulus::A: return samples.A;
  case TestModulus::C: return samples.C;
  case TestModulus::F: return samples.F;
  case TestModulus::L: return samples.L;
  case TestModulus::N: return samples.N;
  }
  throw std::logic_error("Unknown test modulus");
}

double predictedModulusDerivative(const SPARSESPEC::WindowedSensitivityResult &result,
                                  const Full1D::SEM &sem, TestModulus modulus) {
  const auto &samples = result.kernels;
  const auto &kernel = kernelFor(samples, modulus);
  double derivative = 0.0;
  for (int e = 0; e < kernel.rows(); ++e)
    for (int i = 0; i < kernel.cols(); ++i) {
      const double profilePascals = 1e11 *
          (1.0 + 0.3 * sem.mesh().NodeRadius(e, i) / sem.mesh().PR());
      derivative += kernel(e, i) * profilePascals * samples.drWeightsMetres(e, i);
    }
  const double sourceProfilePascals = 1e11 *
      (1.0 + 0.3 * sem.sourceRadius() / sem.mesh().PR());
  switch (modulus) {
  case TestModulus::A: derivative += result.sourcePointKernels.A * sourceProfilePascals; break;
  case TestModulus::C: derivative += result.sourcePointKernels.C * sourceProfilePascals; break;
  case TestModulus::F: derivative += result.sourcePointKernels.F * sourceProfilePascals; break;
  case TestModulus::L: derivative += result.sourcePointKernels.L * sourceProfilePascals; break;
  case TestModulus::N: derivative += result.sourcePointKernels.N * sourceProfilePascals; break;
  }
  return derivative;
}

double predictedDensityDerivative(const SPARSESPEC::MaterialKernelSamples &samples,
                                  const Full1D::SEM &sem,
                                  bool localized = false, double lower = 0.0,
                                  double upper = 0.0) {
  double derivative = 0.0;
  for (int e = 0; e < samples.densityKgPerM3.rows(); ++e)
    for (int i = 0; i < samples.densityKgPerM3.cols(); ++i) {
      const double r = sem.mesh().NodeRadius(e, i);
      double direction = 1000.0 * (0.7 + 0.3 * r / sem.mesh().PR());
      if (localized) {
        if (r < lower || r > upper) direction = 0.0;
        else direction = 1000.0 * 4.0 * (r - lower) * (upper - r) /
                          ((upper - lower) * (upper - lower));
      }
      derivative += samples.densityKgPerM3(e, i) * direction *
                    samples.drWeightsMetres(e, i);
    }
  return derivative;
}

Eigen::MatrixXd directPerturbationGravity(const Full1D::SEM &sem,
                                          const Eigen::MatrixXd &deltaDensity) {
  const auto &mesh = sem.mesh();
  const int ne = mesh.NE(), nq = mesh.NN();
  const auto gll = const_cast<EarthMesh::RadialMesh &>(mesh).GLL();
  const auto basis = Interpolation::LagrangeBasis(gll.Points());
  Eigen::MatrixXd change = Eigen::MatrixXd::Zero(ne, nq);
  for (int e = 0; e < ne; ++e) for (int i = 0; i < nq; ++i) {
    const double upper = mesh.NodeRadius(e, i);
    double mass = 0.0;
    for (int source = 0; source <= e; ++source) {
      const double lo = mesh.ELR(source), hi = mesh.EUR(source);
      const auto &knots = sem.meshModel().DensityKnots(source);
      for (int seg = 0; seg + 1 < nq; ++seg) {
        double left = mesh.NodeRadius(source, seg);
        const double rightNode = mesh.NodeRadius(source, seg + 1);
        if (left >= upper) break;
        std::vector<double> cuts;
        for (double knot : knots) if (knot > left && knot < rightNode && knot < upper) cuts.push_back(knot);
        if (rightNode < upper) cuts.push_back(rightNode);
        else cuts.push_back(upper);
        for (double right : cuts) {
          const double half = 0.5 * (right - left);
          for (int q = 0; q < nq; ++q) {
            const double s = std::lerp(left, right, 0.5 * (gll.X(q) + 1.0));
            const double x = 2.0 * (s - lo) / (hi - lo) - 1.0;
            double value = 0.0;
            for (int j = 0; j < nq; ++j) value += deltaDensity(source, j) * basis.Evaluate(j, x);
            mass += half * gll.W(q) * value * s * s;
          }
          left = right;
        }
        if (rightNode >= upper) break;
      }
    }
    const double r = mesh.NodeRadius(e, i);
    if (r > 0.0) change(e, i) = 4.0 * static_cast<double>(EIGEN_PI) *
        sem.gravityCoefficient() * mass / (r * r);
  }
  return change;
}

Eigen::MatrixXcd derivativeWeights(const SpectraSolver::FreqFull &frequency,
                                   int type) {
  Eigen::MatrixXcd weights = Eigen::MatrixXcd::Zero(3, frequency.nt() / 2 + 1);
  if (type == 1) {
    weights(0, 3) = {0.7, 0.2};
  } else if (type == 2) {
    weights(1, 3) = {0.7, 0.2};
    weights(2, 3) = {-0.13, 0.19};
  } else {
    weights(0, 3) = {0.7, 0.2};
    weights(1, 3) = {-0.13, 0.19};
  }
  return weights;
}

Eigen::VectorXd sensitivityWindow(const SpectraSolver::FreqFull &frequency) {
  Eigen::VectorXd window(frequency.nt());
  for (int i = 0; i < window.size(); ++i)
    window[i] = 0.5 - 0.5 * std::cos(2.0 * PI * i / (window.size() - 1));
  return window;
}

} // namespace

TEST(WindowedSensitivityEngineTests, FixedMeshDirectionalDerivativesForIndependentModuli) {
  EarthModels::EarthConstants<double> norm;
  auto baseModel = EarthModels::ModelInput<double>(DSpecMTest::modelPath().string(), norm);
  Full1D::SEM baseSem(baseModel, 0.01, 4, 2, makeSensitivitySource());
  SourceInfo::EarthquakeCMT source(15.0, 110.0, 30.0,
                                   1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
  const auto frequency = makeSensitivityFrequency(baseModel.TimeNorm());
  const auto window = sensitivityWindow(frequency);
  constexpr int firstBin = 2, endBin = 4;
  constexpr double step = 0.02;
  const std::array<std::pair<int, TestModulus>, 11> directions{{
      {1, TestModulus::A}, {1, TestModulus::C}, {1, TestModulus::F}, {1, TestModulus::N},
      {3, TestModulus::A}, {3, TestModulus::C}, {3, TestModulus::F},
      {3, TestModulus::L}, {3, TestModulus::N},
      {2, TestModulus::L}, {2, TestModulus::N}}};
  std::map<int, SPARSESPEC::WindowedSensitivityResult> baselineResults;
  std::map<int, InputParameters> parameters;
  for (const int type : {1, 2, 3}) {
    parameters.emplace(type, makeSensitivityParameters(type));
    const auto weights = derivativeWeights(frequency, type);
    baselineResults.emplace(type, SPARSESPEC::WindowedSensitivityEngine::evaluate(
        baseSem, source, parameters.at(type), frequency, window, weights,
        firstBin, endBin));
  }

  for (const auto &[type, modulus] : directions) {
    const auto weights = derivativeWeights(frequency, type);
    const auto &baseline = baselineResults.at(type);
    const double predicted = predictedModulusDerivative(baseline, baseSem, modulus);
    ASSERT_TRUE(std::isfinite(predicted));
    ASSERT_GT(std::abs(predicted), 1e-100)
        << "family " << type << " modulus " << static_cast<int>(modulus);

    auto centeredDifference = [&](double h) {
      PerturbedModulusModel plusModel(DSpecMTest::modelPath().string(), norm,
                                      modulus, h);
      Full1D::SEM plusSem(plusModel, 0.01, 4, 2, makeSensitivitySource());
      const double plus = SPARSESPEC::WindowedSensitivityEngine::evaluate(
          plusSem, source, parameters.at(type), frequency, window, weights,
          firstBin, endBin).objective;
      PerturbedModulusModel minusModel(DSpecMTest::modelPath().string(), norm,
                                       modulus, -h);
      Full1D::SEM minusSem(minusModel, 0.01, 4, 2, makeSensitivitySource());
      const double minus = SPARSESPEC::WindowedSensitivityEngine::evaluate(
          minusSem, source, parameters.at(type), frequency, window, weights,
          firstBin, endBin).objective;
      return (plus - minus) / (2.0 * h);
    };
    const double d0 = centeredDifference(step);
    const double d1 = centeredDifference(step / 2.0);
    const double d2 = centeredDifference(step / 4.0);
    const double richardsonCoarse = (4.0 * d1 - d0) / 3.0;
    const double finiteDifference = (4.0 * d2 - d1) / 3.0;
    const double relative = std::abs(finiteDifference - predicted) /
        std::max({std::abs(finiteDifference), std::abs(predicted), 1e-100});
    EXPECT_LE(relative, 1e-6) << "family " << type << " modulus "
                             << static_cast<int>(modulus)
                             << " predicted=" << predicted
                             << " finiteDifference=" << finiteDifference;
  }
}

TEST(WindowedSensitivityEngineTests, SourcePointDerivativeIncludesReceiverInsideLiftPatch) {
  EarthModels::EarthConstants<double> norm;
  const std::string path = DSpecMTest::modelPath().string();
  auto baseModel = EarthModels::ModelInput<double>(path, norm);
  auto source = makeSensitivitySource();
  Full1D::SEM sem(baseModel, 0.01, 4, 2, source);
  auto parameters = makeSensitivityParameters(3, 0, 29.9);
  ASSERT_GT(sem.receiverElements(parameters).size(), 0u);
  const int sourceElement = sem.sourceElement(source);
  const double receiverRadius = (1.0 - 29.9 / 6371.0) * sem.mesh().PR();
  ASSERT_GE(receiverRadius, sem.sourceRadius());
  ASSERT_LE(receiverRadius, sem.mesh().EUR(sourceElement));

  const auto frequency = makeSensitivityFrequency(baseModel.TimeNorm());
  const auto window = sensitivityWindow(frequency);
  const auto weights = derivativeWeights(frequency, 3);
  const auto baseline = SPARSESPEC::WindowedSensitivityEngine::evaluate(
      sem, source, parameters, frequency, window, weights, 2, 4);
  double predicted = baseline.sourcePointKernels.C *
      1e11 * (1.0 + 0.3 * sem.sourceRadius() / sem.mesh().PR());
  for (int e = 0; e < sem.mesh().NE(); ++e)
    if (!sem.elementIsFluid(e))
      for (int i = 0; i < sem.mesh().NN(); ++i) {
        const double profile = 1.0 + 0.3 * sem.mesh().NodeRadius(e, i) /
                                           sem.mesh().PR();
        predicted += baseline.kernels.C(e, i) * 1e11 * profile *
                     baseline.kernels.drWeightsMetres(e, i);
      }
  ASSERT_TRUE(std::isfinite(predicted));
  ASSERT_GT(std::abs(predicted), 1e-100);

  auto objective = [&](double amplitude) {
    PerturbedModulusModel changed(path, norm, TestModulus::C, amplitude);
    Full1D::SEM changedSem(changed, 0.01, 4, 2, source);
    return SPARSESPEC::WindowedSensitivityEngine::evaluate(
        changedSem, source, parameters, frequency, window, weights, 2, 4).objective;
  };
  const double d0 = (objective(0.02) - objective(-0.02)) / 0.04;
  const double d1 = (objective(0.01) - objective(-0.01)) / 0.02;
  const double d2 = (objective(0.005) - objective(-0.005)) / 0.01;
  const double richardson = (4.0 * d2 - d1) / 3.0;
  const double relative = std::abs(richardson - predicted) /
      std::max({std::abs(richardson), std::abs(predicted), 1e-100});
  EXPECT_LE(relative, 1e-6) << "predicted=" << predicted
                            << " centered=" << richardson
                            << " coarsest=" << d0;
}

TEST(WindowedSensitivityEngineTests, DensityDirectionsIncludeGravityForAllModeFamilies) {
  EarthModels::EarthConstants<double> norm;
  const std::string path = DSpecMTest::modelPath().string();
  auto baseModel = EarthModels::ModelInput<double>(path, norm);
  Full1D::SEM baseSem(baseModel, 0.4, 4, 2, makeSensitivitySource());
  SourceInfo::EarthquakeCMT source(15.0, 110.0, 30.0,
                                   1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
  const auto frequency = makeSensitivityFrequency(baseModel.TimeNorm());
  const auto window = sensitivityWindow(frequency);
  constexpr int firstBin = 2, endBin = 4;
  const std::array<int, 3> types{1, 2, 3};
  double maxRelative = 0.0;
  for (const int type : types) {
    auto parameters = makeSensitivityParameters(type);
    const auto weights = derivativeWeights(frequency, type);
    const auto baseline = SPARSESPEC::WindowedSensitivityEngine::evaluate(
        baseSem, source, parameters, frequency, window, weights, firstBin, endBin);
    const double predicted = predictedDensityDerivative(baseline.materialKernels, baseSem);
    ASSERT_TRUE(std::isfinite(predicted));
    ASSERT_GT(std::abs(predicted), 1e-100) << "type=" << type;
    auto centered = [&](double h) {
      PerturbedDensityModel plusModel(path, norm, h, -1.0, -1.0);
      Full1D::SEM plusSem(plusModel, 0.4, 4, 2, makeSensitivitySource());
      const double plus = SPARSESPEC::WindowedSensitivityEngine::evaluate(
          plusSem, source, parameters, frequency, window, weights, firstBin, endBin).objective;
      PerturbedDensityModel minusModel(path, norm, -h, -1.0, -1.0);
      Full1D::SEM minusSem(minusModel, 0.4, 4, 2, makeSensitivitySource());
      const double minus = SPARSESPEC::WindowedSensitivityEngine::evaluate(
          minusSem, source, parameters, frequency, window, weights, firstBin, endBin).objective;
      return (plus - minus) / (2.0 * h);
    };
    const double d0 = centered(0.02), d1 = centered(0.01), d2 = centered(0.005);
    const double richardson = (4.0 * d2 - d1) / 3.0;
    const double relative = std::abs(richardson - predicted) /
        std::max({std::abs(richardson), std::abs(predicted), 1e-100});
    maxRelative = std::max(maxRelative, relative);
    EXPECT_LE(relative, 1e-6) << "type=" << type << " predicted=" << predicted
                              << " centered=" << richardson << " coarse=" << d0;
  }
  ::testing::Test::RecordProperty("max_density_end_to_end_relative_error", maxRelative);
}

TEST(WindowedSensitivityEngineTests, LocalizedSolidDensityDirectionIncludesOutsideGravityEffect) {
  EarthModels::EarthConstants<double> norm;
  const std::string path = DSpecMTest::modelPath().string();
  auto model = EarthModels::ModelInput<double>(path, norm);
  Full1D::SEM sem(model, 0.4, 4, 2, makeSensitivitySource());
  int element = -1;
  for (int e = 0; e < sem.mesh().NE(); ++e)
    if (!sem.elementIsFluid(e) && sem.mesh().ELR(e) > 0.2 * sem.mesh().PR()) {
      element = e; break;
    }
  ASSERT_GE(element, 0);
  const double lower = sem.mesh().ELR(element), upper = sem.mesh().EUR(element);
  const int layer = sem.mesh().LayerNumber(element);
  SourceInfo::EarthquakeCMT source(15.0, 110.0, 30.0,
                                   1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
  auto parameters = makeSensitivityParameters(3);
  const auto frequency = makeSensitivityFrequency(model.TimeNorm());
  const auto window = sensitivityWindow(frequency);
  const auto weights = derivativeWeights(frequency, 3);
  const auto baseline = SPARSESPEC::WindowedSensitivityEngine::evaluate(
      sem, source, parameters, frequency, window, weights, 2, 4);
  const double predicted = predictedDensityDerivative(
      baseline.materialKernels, sem, true, lower, upper);
  ASSERT_TRUE(std::isfinite(predicted));
  ASSERT_GT(std::abs(predicted), 1e-100);
  constexpr double h = 0.01;
  auto objective = [&](double amplitude) {
    PerturbedDensityModel changed(path, norm, amplitude, lower, upper, layer);
    Full1D::SEM changedSem(changed, 0.4, 4, 2, makeSensitivitySource());
    return SPARSESPEC::WindowedSensitivityEngine::evaluate(
        changedSem, source, parameters, frequency, window, weights, 2, 4).objective;
  };
  const double plus = objective(h), minus = objective(-h);
  const double centered = (plus - minus) / (2.0 * h);
  const double relative = std::abs(centered - predicted) /
      std::max({std::abs(centered), std::abs(predicted), 1e-100});
  ::testing::Test::RecordProperty("localized_density_relative_error", relative);
  EXPECT_LE(relative, 1e-6) << "predicted=" << predicted << " centered=" << centered;
}

TEST(WindowedSensitivityEngineTests, FluidBulkDirectionsMatchEndToEndMeasurements) {
  EarthModels::EarthConstants<double> norm;
  const std::string path = DSpecMTest::modelPath().string();
  auto model = EarthModels::ModelInput<double>(path, norm);
  Full1D::SEM sem(model, 0.25, 4, 2, makeSensitivitySource(20.0, 33.0, 10.0));
  int fluidElement = -1;
  double best = std::numeric_limits<double>::infinity();
  for (int e = 0; e < sem.mesh().NE(); ++e) if (sem.elementIsFluid(e)) {
    const double centre = 0.5 * (sem.mesh().ELR(e) + sem.mesh().EUR(e));
    if (std::abs(centre - 0.3 * sem.mesh().PR()) < best) {
      best = std::abs(centre - 0.3 * sem.mesh().PR()); fluidElement = e;
    }
  }
  ASSERT_GE(fluidElement, 0);
  const double lower = sem.mesh().ELR(fluidElement), upper = sem.mesh().EUR(fluidElement);
  const int fluidLayer = sem.mesh().LayerNumber(fluidElement);
  SourceInfo::EarthquakeCMT source(33.0, 10.0, 20.0,
                                   1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
  const auto frequency = makeSensitivityFrequency(model.TimeNorm());
  const auto window = sensitivityWindow(frequency);
  constexpr int firstBin = 2, endBin = 4;
  double maxRelative = 0.0;
  for (int type : {1, 3}) {
    auto parameters = makeSensitivityParameters(type);
    Eigen::MatrixXcd weights = Eigen::MatrixXcd::Zero(3, frequency.nt() / 2 + 1);
    weights(0, 2) = {0.7, 0.2}; weights(0, 3) = {-0.25, 0.45};
    const auto baseline = SPARSESPEC::WindowedSensitivityEngine::evaluate(
        sem, source, parameters, frequency, window, weights, firstBin, endBin);
    double predicted = 0.0;
    for (int e = 0; e < baseline.materialKernels.bulkModulusPa.rows(); ++e)
      for (int i = 0; i < baseline.materialKernels.bulkModulusPa.cols(); ++i)
        if (e == fluidElement) {
          const double r = sem.mesh().NodeRadius(e, i);
          const double x = (2.0 * r - lower - upper) / (upper - lower);
          predicted += baseline.materialKernels.bulkModulusPa(e, i) *
            1e11 * (1.0 - x * x) *
            baseline.materialKernels.drWeightsMetres(e, i);
        }
    ASSERT_GT(std::abs(predicted), 1e-100);
    auto centered = [&](double h) {
      PerturbedFluidBulkModel plusModel(path, norm, h, lower, upper, fluidLayer);
      Full1D::SEM plus(plusModel, 0.25, 4, 2, makeSensitivitySource(20.0, 33.0, 10.0));
      const double yp = SPARSESPEC::WindowedSensitivityEngine::evaluate(
          plus, source, parameters, frequency, window, weights, firstBin, endBin).objective;
      PerturbedFluidBulkModel minusModel(path, norm, -h, lower, upper, fluidLayer);
      Full1D::SEM minus(minusModel, 0.25, 4, 2, makeSensitivitySource(20.0, 33.0, 10.0));
      const double ym = SPARSESPEC::WindowedSensitivityEngine::evaluate(
          minus, source, parameters, frequency, window, weights, firstBin, endBin).objective;
      return (yp - ym) / (2.0 * h);
    };
    const double d0 = centered(0.16), d1 = centered(0.08);
    const double d2 = centered(0.04), d3 = centered(0.02);
    const double richardson1 = (4.0 * d2 - d1) / 3.0;
    const double richardson2 = (4.0 * d3 - d2) / 3.0;
    const double richardson = (16.0 * richardson2 - richardson1) / 15.0;
    const double relative = std::abs(richardson - predicted) /
        std::max({std::abs(richardson), std::abs(predicted), 1e-100});
    maxRelative = std::max(maxRelative, relative);
    EXPECT_LE(relative, 1e-6) << "type=" << type << " predicted=" << predicted
                              << " centered=" << richardson << " steps0.16-0.02="
                              << d0 << "," << d1 << "," << d2 << "," << d3;
  }
  ::testing::Test::RecordProperty("max_fluid_bulk_end_to_end_relative_error", maxRelative);
}

TEST(WindowedSensitivityEngineTests, DensityTaylorRemainderIsSecondOrder) {
  EarthModels::EarthConstants<double> norm;
  const std::string path = DSpecMTest::modelPath().string();
  auto model = EarthModels::ModelInput<double>(path, norm);
  Full1D::SEM sem(model, 0.4, 4, 2, makeSensitivitySource());
  SourceInfo::EarthquakeCMT source(15.0, 110.0, 30.0,
                                   1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
  auto parameters = makeSensitivityParameters(3);
  const auto frequency = makeSensitivityFrequency(model.TimeNorm());
  const auto window = sensitivityWindow(frequency);
  const auto weights = derivativeWeights(frequency, 3);
  const auto baseline = SPARSESPEC::WindowedSensitivityEngine::evaluate(
      sem, source, parameters, frequency, window, weights, 2, 4);
  const double derivative = predictedDensityDerivative(baseline.materialKernels, sem);
  ASSERT_GT(std::abs(derivative), 1e-100);
  const std::array<double, 4> steps{0.08, 0.04, 0.02, 0.01};
  std::array<double, 4> remainders{};
  for (std::size_t i = 0; i < steps.size(); ++i) {
    PerturbedDensityModel changed(path, norm, steps[i], -1.0, -1.0);
    Full1D::SEM changedSem(changed, 0.4, 4, 2, makeSensitivitySource());
    const double value = SPARSESPEC::WindowedSensitivityEngine::evaluate(
        changedSem, source, parameters, frequency, window, weights, 2, 4).objective;
    remainders[i] = std::abs(value - baseline.objective - steps[i] * derivative);
    ASSERT_GT(remainders[i], 1e-100);
  }
  for (std::size_t i = 0; i + 1 < steps.size(); ++i) {
    const double slope = std::log(remainders[i] / remainders[i + 1]) /
                         std::log(steps[i] / steps[i + 1]);
    EXPECT_GE(slope, 1.8) << "steps=" << steps[i] << " slope=" << slope;
  }
  ::testing::Test::RecordProperty("density_taylor_minimum_slope", std::min(
      std::log(remainders[0] / remainders[1]) / std::log(steps[0] / steps[1]),
      std::min(std::log(remainders[1] / remainders[2]) / std::log(steps[1] / steps[2]),
               std::log(remainders[2] / remainders[3]) / std::log(steps[2] / steps[3]))));
}

TEST(WindowedSensitivityEngineTests, SecondOrderTaylorRemainderForFixedMesh) {
  EarthModels::EarthConstants<double> norm;
  auto baseModel = EarthModels::ModelInput<double>(DSpecMTest::modelPath().string(), norm);
  Full1D::SEM baseSem(baseModel, 0.01, 4, 2, makeSensitivitySource());
  SourceInfo::EarthquakeCMT source(15.0, 110.0, 30.0,
                                   1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
  auto parameters = makeSensitivityParameters(3);
  const auto frequency = makeSensitivityFrequency(baseModel.TimeNorm());
  const auto window = sensitivityWindow(frequency);
  const auto weights = derivativeWeights(frequency, 3);
  const auto base = SPARSESPEC::WindowedSensitivityEngine::evaluate(
      baseSem, source, parameters, frequency, window, weights, 2, 4);
  const double derivative = predictedModulusDerivative(
      base, baseSem, TestModulus::C);
  ASSERT_GT(std::abs(derivative), 1e-100);
  const std::array<double, 4> steps{0.4, 0.2, 0.1, 0.05};
  std::array<double, 4> remainders{};
  for (std::size_t i = 0; i < steps.size(); ++i) {
    PerturbedModulusModel changedModel(DSpecMTest::modelPath().string(), norm,
                                       TestModulus::C, steps[i]);
    Full1D::SEM changedSem(changedModel, 0.01, 4, 2, makeSensitivitySource());
    const double changed = SPARSESPEC::WindowedSensitivityEngine::evaluate(
        changedSem, source, parameters, frequency, window, weights, 2, 4).objective;
    remainders[i] = std::abs(changed - base.objective - steps[i] * derivative);
    ASSERT_GT(remainders[i], 1e-100) << "Taylor remainder reached roundoff";
  }
  for (std::size_t i = 0; i + 1 < steps.size(); ++i) {
    const double slope = std::log(remainders[i] / remainders[i + 1]) /
                         std::log(steps[i] / steps[i + 1]);
    EXPECT_GE(slope, 1.8) << "steps " << steps[i] << " and " << steps[i + 1];
  }
}


template <typename T> class DoubledLengthScales {
public:
  DoubledLengthScales(T length, T mass, T time)
      : length_(2 * length), mass_(8 * mass), time_(2 * time) {}
  T LengthNorm() const { return length_; }
  T MassNorm() const { return mass_; }
  T TimeNorm() const { return time_; }
private:
  T length_, mass_, time_;
};

TEST(WindowedSensitivityEngineTests, SIOutputsAndKernelsAreInvariantToComputationalRescaling) {
  EarthModels::EarthConstants<double> baseScales;
  DoubledLengthScales<double> scaledScales(baseScales.LengthNorm(),
                                            baseScales.MassNorm(),
                                            baseScales.TimeNorm());
  const std::string path = DSpecMTest::modelPath().string();
  auto baseModel = EarthModels::ModelInput<double>(path, baseScales);
  auto scaledModel = EarthModels::ModelInput<double>(path, scaledScales);
  Full1D::SEM baseSem(baseModel, 0.01, 4, 2, makeSensitivitySource());
  Full1D::SEM scaledSem(scaledModel, 0.01, 4, 2, makeSensitivitySource());
  ASSERT_EQ(baseSem.mesh().NE(), scaledSem.mesh().NE());
  EXPECT_NEAR(baseSem.lengthScaleMetres() * 2.0,
              scaledSem.lengthScaleMetres(), 1e-12 * scaledSem.lengthScaleMetres());
  EXPECT_NEAR(baseSem.pressureScalePascals(), scaledSem.pressureScalePascals(),
              1e-12 * scaledSem.pressureScalePascals());
  EXPECT_NEAR(baseSem.timeScaleSeconds() * 2.0,
              scaledSem.timeScaleSeconds(), 1e-12 * scaledSem.timeScaleSeconds());

  SourceInfo::EarthquakeCMT source(15.0, 110.0, 30.0,
                                   1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
  auto baseParameters = makeSensitivityParameters(4);
  auto scaledParameters = makeSensitivityParameters(4);
  const auto baseFrequency = makeSensitivityFrequency(baseModel.TimeNorm());
  const auto scaledFrequency = makeSensitivityFrequency(scaledModel.TimeNorm());
  ASSERT_EQ(baseFrequency.nt(), scaledFrequency.nt());
  const auto window = sensitivityWindow(baseFrequency);
  Eigen::MatrixXcd weights = Eigen::MatrixXcd::Zero(3, baseFrequency.nt() / 2 + 1);
  weights(0, 3) = {0.7, 0.2};
  weights(1, 3) = {-0.13, 0.19};
  const auto base = SPARSESPEC::WindowedSensitivityEngine::evaluate(
      baseSem, source, baseParameters, baseFrequency, window, weights, 2, 4);
  const auto scaled = SPARSESPEC::WindowedSensitivityEngine::evaluate(
      scaledSem, source, scaledParameters, scaledFrequency, window, weights, 2, 4);

  auto relativeMatrixError = [](const auto &left, const auto &right) {
    return (left - right).norm() /
           std::max({left.norm(), right.norm(), 1e-300});
  };
  EXPECT_LE(relativeMatrixError(
                baseSem.lengthScaleMetres() * baseSem.timeScaleSeconds() *
                    base.rawDisplacementSpectrum,
                scaledSem.lengthScaleMetres() * scaledSem.timeScaleSeconds() *
                    scaled.rawDisplacementSpectrum),
            1e-9);
  EXPECT_LE(relativeMatrixError(base.processedSISpectrum,
                                scaled.processedSISpectrum), 1e-9);
  EXPECT_LE(std::abs(base.objective - scaled.objective) /
                std::max({std::abs(base.objective), std::abs(scaled.objective), 1e-300}),
            1e-9);
  EXPECT_LE(relativeMatrixError(base.kernels.A, scaled.kernels.A), 1e-9);
  EXPECT_LE(relativeMatrixError(base.kernels.C, scaled.kernels.C), 1e-9);
  EXPECT_LE(relativeMatrixError(base.kernels.F, scaled.kernels.F), 1e-9);
  EXPECT_LE(relativeMatrixError(base.kernels.L, scaled.kernels.L), 1e-9);
  EXPECT_LE(relativeMatrixError(base.kernels.N, scaled.kernels.N), 1e-9);
  EXPECT_LE(relativeMatrixError(base.materialKernels.bulkModulusPa,
                                scaled.materialKernels.bulkModulusPa), 1e-9);
  EXPECT_LE(relativeMatrixError(base.materialKernels.densityKgPerM3,
                                scaled.materialKernels.densityKgPerM3), 1e-9);
  EXPECT_LE(relativeMatrixError(base.kernels.radiusMetres,
                                scaled.kernels.radiusMetres), 1e-12);
  EXPECT_LE(relativeMatrixError(base.kernels.drWeightsMetres,
                                scaled.kernels.drWeightsMetres), 1e-12);
  EXPECT_GT(std::abs(base.objective), 1e-300);

  const int bin = 3;
  auto baseFieldParameters = makeSensitivityParameters(4, 0, 1000.0);
  auto scaledFieldParameters = makeSensitivityParameters(4, 0, 1000.0);
  const auto receiver = baseFieldParameters.receivers().front();
  const double degreeToRadians = EIGEN_PI/180.0;
  const double theta = (90.0-receiver.first)*degreeToRadians;
  const double phi = receiver.second*degreeToRadians;
  const double baseRadiusSI = (baseSem.mesh().PR() -
      1000.0*baseFieldParameters.receiver_depth()/baseSem.lengthScaleMetres()) *
      baseSem.lengthScaleMetres();
  const double scaledRadiusSI = (scaledSem.mesh().PR() -
      1000.0*scaledFieldParameters.receiver_depth()/scaledSem.lengthScaleMetres()) *
      scaledSem.lengthScaleMetres();
  EXPECT_NEAR(baseRadiusSI, scaledRadiusSI, 1e-12*baseRadiusSI);
  const std::vector<SPARSESPEC::Field3DPoint> basePoint{
      {baseRadiusSI, theta, phi}};
  const std::vector<SPARSESPEC::Field3DPoint> scaledPoint{
      {scaledRadiusSI, theta, phi}};
  SPARSESPEC::WindowedDisplacementMap baseMap(
      baseFrequency, baseSem.lengthScaleMetres(), window);
  SPARSESPEC::WindowedDisplacementMap scaledMap(
      scaledFrequency, scaledSem.lengthScaleMetres(), window);
  const Eigen::MatrixXcd baseAdjointWeights = baseMap.adjoint(weights);
  const Eigen::MatrixXcd scaledAdjointWeights = scaledMap.adjoint(weights);
  const Eigen::VectorXcd baseQhat = baseAdjointWeights.col(bin);
  const Eigen::VectorXcd scaledQhat = scaledAdjointWeights.col(bin);
  const std::complex<double> baseZ = baseFrequency.w().at(bin) +
      std::complex<double>(0.0, -baseFrequency.ep());
  const std::complex<double> scaledZ = scaledFrequency.w().at(bin) +
      std::complex<double>(0.0, -scaledFrequency.ep());
  const auto baseField = SPARSESPEC::Windowed3DFieldEngine::evaluateFrequency(
      baseSem, source, baseFieldParameters, baseZ, baseQhat, basePoint).front();
  const auto scaledField = SPARSESPEC::Windowed3DFieldEngine::evaluateFrequency(
      scaledSem, source, scaledFieldParameters, scaledZ, scaledQhat, scaledPoint).front();
  Eigen::Vector3cd baseSpectrum, scaledSpectrum;
  for (int c = 0; c < 3; ++c) {
    baseSpectrum[c] = baseField.observableForwardDisplacementSpectrumMetreSeconds[c];
    scaledSpectrum[c] =
        scaledField.observableForwardDisplacementSpectrumMetreSeconds[c];
  }
  const auto baseFieldReference = SPARSESPEC::WindowedSensitivityEngine::evaluate(
      baseSem, source, baseFieldParameters, baseFrequency, window, weights, bin, bin+1);
  const auto scaledFieldReference = SPARSESPEC::WindowedSensitivityEngine::evaluate(
      scaledSem, source, scaledFieldParameters, scaledFrequency, window, weights,
      bin, bin+1);
  const Eigen::VectorXcd baseExpected = baseSem.lengthScaleMetres() *
      baseSem.timeScaleSeconds() * baseFieldReference.rawDisplacementSpectrum.col(bin);
  const Eigen::VectorXcd scaledExpected = scaledSem.lengthScaleMetres() *
      scaledSem.timeScaleSeconds() * scaledFieldReference.rawDisplacementSpectrum.col(bin);
  const double baseFieldError = (baseSpectrum-baseExpected).norm() /
      std::max(baseSpectrum.norm(), baseExpected.norm());
  const double scaledFieldError = (scaledSpectrum-scaledExpected).norm() /
      std::max(scaledSpectrum.norm(), scaledExpected.norm());
  EXPECT_LT(baseFieldError, 1e-9) << "base 3D spectrum error=" << baseFieldError
      << " actual=" << baseSpectrum.transpose()
      << " expected=" << baseExpected.transpose()
      << " T0=" << baseSem.timeScaleSeconds();
  EXPECT_LT(scaledFieldError, 1e-9) << "scaled 3D spectrum error=" << scaledFieldError
      << " actual=" << scaledSpectrum.transpose()
      << " expected=" << scaledExpected.transpose()
      << " T0=" << scaledSem.timeScaleSeconds();
  const double rescalingError = (baseSpectrum-scaledSpectrum).norm() /
      std::max(baseSpectrum.norm(), scaledSpectrum.norm());
  EXPECT_LT(rescalingError, 1e-9)
      << "3D SI spectrum rescaling error=" << rescalingError;
  const auto baseKernel=SPARSESPEC::Windowed3DKernelEngine::evaluateFrequency(
      baseSem,source,baseFieldParameters,baseZ,baseQhat,basePoint).front();
  const auto scaledKernel=SPARSESPEC::Windowed3DKernelEngine::evaluateFrequency(
      scaledSem,source,scaledFieldParameters,scaledZ,scaledQhat,scaledPoint).front();
  Eigen::Matrix<double,5,1> base3DK,scaled3DK;
  base3DK<<baseKernel.A,baseKernel.C,baseKernel.F,baseKernel.L,baseKernel.N;
  scaled3DK<<scaledKernel.A,scaledKernel.C,scaledKernel.F,scaledKernel.L,scaledKernel.N;
  EXPECT_LT((base3DK-scaled3DK).norm()/std::max(base3DK.norm(),1e-100),1e-9)
      << "3D SI kernel rescaling error="<<(base3DK-scaled3DK).norm()/
          std::max(base3DK.norm(),1e-100);

  const std::vector<std::array<double,2>> angles{{0.8,1.1},{1.7,4.2}};
  const auto baseCombined=SPARSESPEC::Windowed3DKernelEngine::evaluate(
      baseSem,source,baseFieldParameters,baseFrequency,window,weights,bin,bin+1,
      basePoint,angles);
  const auto scaledCombined=SPARSESPEC::Windowed3DKernelEngine::evaluate(
      scaledSem,source,scaledFieldParameters,scaledFrequency,window,weights,bin,bin+1,
      scaledPoint,angles);
  const auto relative=[](const auto &x,const auto &y){
    return (x-y).norm()/std::max({x.norm(),y.norm(),1e-100});
  };
  const auto &bc=baseCombined.sourceContact;
  const auto &sc=scaledCombined.sourceContact;
  ASSERT_EQ(bc.size(),sc.size());
  for(std::size_t j=0;j<bc.size();++j){
    Eigen::Vector3d b{bc[j].C,bc[j].F,bc[j].L};
    Eigen::Vector3d s{sc[j].C,sc[j].F,sc[j].L};
    EXPECT_LT(relative(b,s),1e-9)<<"3D source contact rescaling angle="<<j;
  }
  const double radialWeight=1000.0,angularWeight=0.2;
  const std::vector<SPARSESPEC::WeightedField3DPoint> baseRule{
      {basePoint.front(),radialWeight,angularWeight}};
  const std::vector<SPARSESPEC::WeightedField3DPoint> scaledRule{
      {scaledPoint.front(),radialWeight,angularWeight}};
  const auto baseVolume=SPARSESPEC::Windowed3DKernelEngine::
      integrateWeightedVolumeFrequency(baseSem,source,baseFieldParameters,baseZ,
                                         baseQhat,baseRule);
  const auto scaledVolume=SPARSESPEC::Windowed3DKernelEngine::
      integrateWeightedVolumeFrequency(scaledSem,source,scaledFieldParameters,scaledZ,
                                         scaledQhat,scaledRule);
  Eigen::Matrix<double,6,1> baseV,scaledV;
  baseV<<baseVolume.A,baseVolume.C,baseVolume.F,baseVolume.L,baseVolume.N,baseVolume.kappa;
  scaledV<<scaledVolume.A,scaledVolume.C,scaledVolume.F,scaledVolume.L,
           scaledVolume.N,scaledVolume.kappa;
  EXPECT_LT(relative(baseV,scaledV),1e-9)<<"3D weighted volume SI rescaling";
}

TEST(Windowed3DKernelEngineTests, SpinCoefficientsMatchIndependentCartesianTIContractions) {
  EarthModels::EarthConstants<double> norm;
  auto model = EarthModels::ModelInput(DSpecMTest::modelPath().string(), norm);
  Full1D::SEM sem(model, 0.1, 4, 2);
  const auto polar = [](std::complex<double> rr, std::complex<double> tt,
                        std::complex<double> pp, std::complex<double> rt,
                        std::complex<double> rp, std::complex<double> tp) {
    const std::complex<double> i(0.0, 1.0);
    const auto h = tt + pp, d = tt - pp;
    return SPARSESPEC::SpinStrainField{
        rr, h, rt-i*rp, rt+i*rp,
        d+2.0*i*tp, d-2.0*i*tp};
  };
  const auto f = polar({0.7,-0.2}, {0.4,0.1}, {-0.3,0.6},
                       {0.2,0.8}, {-0.5,0.3}, {0.9,-0.1});
  const auto a = polar({-0.2,0.4}, {0.8,-0.3}, {0.1,0.5},
                       {-0.7,0.2}, {0.6,0.4}, {-0.1,0.9});
  SPARSESPEC::Field3DSample sample;
  sample.forwardRawStrain = f;
  sample.adjointDualStrain = a;
  const auto got = SPARSESPEC::Windowed3DKernelEngine::pointKernel(sem, sample);
  const double scale = sem.pressureScalePascals() *
      std::pow(sem.lengthScaleMetres(), 3);
  // Independent Cartesian TI tensor contractions in the (r,theta,phi) frame.
  const double expectedA=-std::real(std::conj(a.horizontalTrace)*f.horizontalTrace);
  const double expectedC=-std::real(std::conj(a.rr)*f.rr);
  const double expectedF=-std::real(std::conj(a.rr)*f.horizontalTrace+
                                    std::conj(a.horizontalTrace)*f.rr);
  // Direct component form of the TI stiffness derivative, independently of
  // the production spin-component conversion.
  const std::complex<double> arT(-0.7,0.2), arP(0.6,0.4);
  const std::complex<double> frT(0.2,0.8), frP(-0.5,0.3);
  const std::complex<double> atT(0.8,-0.3), atP(0.1,0.5);
  const std::complex<double> ftT(0.4,0.1), ftP(-0.3,0.6);
  const std::complex<double> atTP(-0.1,0.9), ftTP(0.9,-0.1);
  const double expectedL=-4.0*std::real(std::conj(arT)*frT+
                                         std::conj(arP)*frP);
  const double expectedN=std::real(2.0*std::conj(atT)*ftP+
      2.0*std::conj(atP)*ftT-4.0*std::conj(atTP)*ftTP);
  EXPECT_NEAR(got.A*scale,expectedA,1e-12*std::max(1.0,std::abs(expectedA)));
  EXPECT_NEAR(got.C*scale,expectedC,1e-12*std::max(1.0,std::abs(expectedC)));
  EXPECT_NEAR(got.F*scale,expectedF,1e-12*std::max(1.0,std::abs(expectedF)));
  EXPECT_NEAR(got.L*scale,expectedL,1e-12*std::max(1.0,std::abs(expectedL)));
  EXPECT_NEAR(got.N*scale,expectedN,1e-12*std::max(1.0,std::abs(expectedN)));

  sample.fluid = true;
  const auto fluid = SPARSESPEC::Windowed3DKernelEngine::pointKernel(sem, sample);
  const auto divF = f.rr + f.horizontalTrace;
  const auto divA = a.rr + a.horizontalTrace;
  EXPECT_DOUBLE_EQ(fluid.A, 0.0); EXPECT_DOUBLE_EQ(fluid.C, 0.0);
  EXPECT_DOUBLE_EQ(fluid.F, 0.0); EXPECT_DOUBLE_EQ(fluid.L, 0.0);
  EXPECT_DOUBLE_EQ(fluid.N, 0.0);
  const double expectedBulk=-std::real(std::conj(divA)*divF);
  EXPECT_NEAR(fluid.kappa*scale,expectedBulk,
              1e-12*std::max(1.0,std::abs(expectedBulk)));
}

TEST(Windowed3DKernelEngineTests, WeightedCenterIsAlgebraicAndPointAPIRejectsIt) {
  auto parameters = makeSensitivityParameters(1);
  EarthModels::EarthConstants<double> norm;
  auto model = EarthModels::ModelInput(DSpecMTest::modelPath().string(), norm);
  auto source = makeSensitivitySource();
  Full1D::SEM sem(model, 0.1, 4, 2, source);
  const auto frequency = makeSensitivityFrequency(model.TimeNorm());
  constexpr int bin = 3;
  Eigen::MatrixXcd q = Eigen::MatrixXcd::Zero(3, frequency.nt()/2+1);
  q(0,bin) = {0.7,-0.2}; q(1,bin) = {-0.3,0.4}; q(2,bin) = {0.2,0.5};
  Eigen::VectorXd window = Eigen::VectorXd::Ones(frequency.nt());
  SPARSESPEC::WindowedDisplacementMap map(frequency, sem.lengthScaleMetres(), window);
  const Eigen::VectorXcd qhat = map.adjoint(q).col(bin);
  const auto z = frequency.w().at(bin) + std::complex<double>(0.0,-frequency.ep());
  const double centerSI = sem.mesh().ELR(0)*sem.lengthScaleMetres();
  ASSERT_NEAR(centerSI, 0.0, 1e-12*sem.lengthScaleMetres());
  const std::vector<SPARSESPEC::Field3DPoint> point{{0.0,0.7,1.1}};
  EXPECT_THROW(SPARSESPEC::Windowed3DKernelEngine::evaluateFrequency(
      sem, source, parameters, z, qhat, point), std::invalid_argument);
  const int sourceElement=sem.sourceElement(source);
  const double interiorRadius=0.5*(sem.mesh().ELR(sourceElement)+
      sem.mesh().EUR(sourceElement))*sem.lengthScaleMetres();
  const std::vector<SPARSESPEC::WeightedField3DPoint> quadrature{
      {{0.0,0.7,1.1}, 2.0e3, 0.25},
      {{interiorRadius,1.0,2.0,sourceElement}, 4.0e3, 0.5}};
  const auto integrated = SPARSESPEC::Windowed3DKernelEngine::
      integrateWeightedVolumeFrequency(sem, source, parameters, z, qhat, quadrature);
  EXPECT_TRUE(std::isfinite(integrated.A+integrated.C+integrated.F+
                            integrated.L+integrated.N));
  EXPECT_TRUE(std::isfinite(integrated.kappa));
  auto workspace=SPARSESPEC::Windowed3DFieldEngine::prepareFrequencyWorkspace(
      sem,source,parameters,z,qhat);
  const auto chunked=SPARSESPEC::Windowed3DKernelEngine::
      integrateWeightedVolumeFrequency(workspace,quadrature,1);
  const auto whole=SPARSESPEC::Windowed3DKernelEngine::
      integrateWeightedVolumeFrequency(workspace,quadrature,quadrature.size());
  Eigen::VectorXd chunkValues(6),wholeValues(6);
  chunkValues<<chunked.A,chunked.C,chunked.F,chunked.L,chunked.N,chunked.kappa;
  wholeValues<<whole.A,whole.C,whole.F,whole.L,whole.N,whole.kappa;
  EXPECT_LT((chunkValues-wholeValues).norm(),1e-12*std::max(wholeValues.norm(),1e-100));

  using Complex=std::complex<double>;
  constexpr double eps=1e-9;
  const Complex y00=1.0/std::sqrt(4.0*static_cast<double>(EIGEN_PI));
  const Complex uf(0.7,0.2),upf(0.3,-0.1),ua(-0.2,0.4),upa(0.1,0.25);
  const auto centerF=SPARSESPEC::Windowed3DFieldEngine::rScaledStrainFromModal(
      0,uf,upf,{},{},{},{},0.0,y00,{},{},{},{});
  const auto centerA=SPARSESPEC::Windowed3DFieldEngine::rScaledStrainFromModal(
      0,ua,upa,{},{},{},{},0.0,y00,{},{},{},{});
  const auto nearF=SPARSESPEC::Windowed3DFieldEngine::rScaledStrainFromModal(
      0,uf+upf*eps,upf,{},{},{},{},eps,y00,{},{},{},{});
  const auto nearA=SPARSESPEC::Windowed3DFieldEngine::rScaledStrainFromModal(
      0,ua+upa*eps,upa,{},{},{},{},eps,y00,{},{},{},{});
  const double centerAterm=-std::real(std::conj(centerA.horizontalTrace)*centerF.horizontalTrace);
  const double positiveAterm=-eps*eps*std::real(
      std::conj(nearA.horizontalTrace/eps)*(nearF.horizontalTrace/eps));
  EXPECT_LT(std::abs(centerAterm-positiveAterm)/std::max(std::abs(centerAterm),1e-100),1e-8);
  EXPECT_EQ(centerF.rr,Complex{});
  EXPECT_EQ(centerA.rr,Complex{});
}

TEST(Windowed3DKernelEngineTests, ForwardOnlyReceiverMatchesWindowedEngineOnBothTraces) {
  EarthModels::EarthConstants<double> norm;
  const auto modelPath = DSpecMTest::modelPath().string();
  auto model = EarthModels::ModelInput<double>(modelPath, norm);
  auto source = makeSensitivitySource();
  Full1D::SEM sem(model, 0.01, 4, 2, source);
  const auto frequency = makeSensitivityFrequency(model.TimeNorm());
  const auto window = sensitivityWindow(frequency);
  const auto weights = derivativeWeights(frequency, 3);
  constexpr int firstBin = 2, endBin = 4;
  for (const double depthKm : {3.0, 29.9}) {
    auto parameters = makeSensitivityParameters(3, 0, depthKm);
    const auto point = SPARSESPEC::Windowed3DKernelEngine::evaluateForwardObjective(
        sem, source, parameters, frequency, window, weights, firstBin, endBin);
    const auto reference = SPARSESPEC::WindowedSensitivityEngine::evaluate(
        sem, source, parameters, frequency, window, weights, firstBin, endBin);
    const double rawError = (point.rawDisplacementSpectrum-
        reference.rawDisplacementSpectrum).norm() /
        std::max(point.rawDisplacementSpectrum.norm(),
                 reference.rawDisplacementSpectrum.norm());
    const double processedError = (point.processedSISpectrum-
        reference.processedSISpectrum).norm() /
        std::max(point.processedSISpectrum.norm(),
                 reference.processedSISpectrum.norm());
    EXPECT_LT(rawError, 1e-9) << "depth_km=" << depthKm
      << " actual=" << point.rawDisplacementSpectrum.col(3).transpose()
      << " expected=" << reference.rawDisplacementSpectrum.col(3).transpose();
    EXPECT_LT(processedError, 1e-9) << "depth_km=" << depthKm;
    const auto recElem=sem.receiverElements(parameters,false).front();
    const auto receiver=parameters.receivers().front();
    const std::vector<SPARSESPEC::Field3DPoint> sample{{
        (sem.mesh().PR()-1000.0*depthKm/sem.lengthScaleMetres())*
            sem.lengthScaleMetres(),
        (90.0-receiver.first)*static_cast<double>(EIGEN_PI)/180.0,
        receiver.second*static_cast<double>(EIGEN_PI)/180.0,
        static_cast<int>(recElem)}};
    const auto highLevel=SPARSESPEC::Windowed3DKernelEngine::evaluate(
        sem,source,parameters,frequency,window,weights,firstBin,endBin,sample);
    EXPECT_LT((highLevel.processedSISpectrum-point.processedSISpectrum).norm()/
        std::max(point.processedSISpectrum.norm(),1e-100),1e-9);
    EXPECT_NEAR(highLevel.objective,point.objective,
        1e-9*std::max(std::abs(point.objective),1e-100));
  }
}

TEST(Windowed3DKernelEngineTests, MultipleReceiversShareSelectedRadialTrace) {
  EarthModels::EarthConstants<double> norm;
  const auto path=DSpecMTest::modelPath().string();
  auto model=EarthModels::ModelInput<double>(path,norm);
  auto source=makeSensitivitySource();
  Full1D::SEM sem(model,0.01,4,2,source);
  const auto frequency=makeSensitivityFrequency(model.TimeNorm());
  const auto window=sensitivityWindow(frequency);
  Eigen::MatrixXcd weights=Eigen::MatrixXcd::Zero(6,frequency.nt()/2+1);
  weights(0,3)={0.7,0.2};weights(1,3)={-0.13,0.19};
  weights(3,3)={-0.23,0.11};weights(5,3)={0.17,-0.09};
  constexpr int firstBin=2,endBin=4;
  const double radians=static_cast<double>(EIGEN_PI)/180.0;
  for(const double depthKm:{1000.0,3.0}){
    auto parameters=makeMultiReceiverSensitivityParameters(4,depthKm,
        {{30.0,65.0},{-10.0,120.0}});
    ASSERT_EQ(parameters.num_receivers(),2);
    const auto receiverElements=sem.receiverElements(parameters,false);
    ASSERT_FALSE(receiverElements.empty());
    const auto expected=SPARSESPEC::WindowedSensitivityEngine::evaluate(
        sem,source,parameters,frequency,window,weights,firstBin,endBin);
    const auto forward=SPARSESPEC::Windowed3DKernelEngine::evaluateForwardObjective(
        sem,source,parameters,frequency,window,weights,firstBin,endBin);
    const auto relative=[](const auto &a,const auto &b){
      return (a-b).norm()/std::max({a.norm(),b.norm(),1e-100});
    };
    EXPECT_LT(relative(forward.rawDisplacementSpectrum,
                       expected.rawDisplacementSpectrum),1e-9)<<"depth="<<depthKm;
    EXPECT_LT(relative(forward.processedSISpectrum,
                       expected.processedSISpectrum),1e-9)<<"depth="<<depthKm;
    EXPECT_LT(std::abs(forward.objective-expected.objective)/
        std::max({std::abs(forward.objective),std::abs(expected.objective),1e-100}),1e-9)
        <<"depth="<<depthKm;
    const double radius=(sem.mesh().PR()-1000.0*depthKm/
        sem.lengthScaleMetres())*sem.lengthScaleMetres();
    std::vector<SPARSESPEC::Field3DPoint> points{{radius,0.8,1.4}};
    const auto highLevel=SPARSESPEC::Windowed3DKernelEngine::evaluate(
        sem,source,parameters,frequency,window,weights,firstBin,endBin,points);
    EXPECT_LT(relative(highLevel.rawDisplacementSpectrum,
                       expected.rawDisplacementSpectrum),1e-9)<<"depth="<<depthKm;
    EXPECT_LT(relative(highLevel.processedSISpectrum,
                       expected.processedSISpectrum),1e-9)<<"depth="<<depthKm;
    EXPECT_LT(std::abs(highLevel.objective-expected.objective)/
        std::max({std::abs(highLevel.objective),std::abs(expected.objective),1e-100}),1e-9)
        <<"depth="<<depthKm;
    for(int rec=0;rec<parameters.num_receivers();++rec){
      const auto coordinate=parameters.receivers()[rec];
      const SPARSESPEC::Field3DPoint receiverPoint{
          radius,(90.0-coordinate.first)*radians,coordinate.second*radians,
          receiverElements.front()};
      points.push_back(receiverPoint);
    }
    const auto receiverSamples=SPARSESPEC::Windowed3DFieldEngine::evaluateFrequency(
        sem,source,parameters,frequency.w().at(3)+
            std::complex<double>(0.0,-frequency.ep()),
        SPARSESPEC::WindowedDisplacementMap(
            frequency,sem.lengthScaleMetres(),window).adjoint(weights).col(3),points);
    ASSERT_EQ(receiverSamples.size(),3u);
    EXPECT_EQ(receiverSamples[1].radialElement,receiverElements.front());
    EXPECT_EQ(receiverSamples[2].radialElement,receiverElements.front());
  }
}

TEST(Windowed3DKernelEngineTests, RejectsFrequencyWithMismatchedModelTimeScale) {
  EarthModels::EarthConstants<double> norm;
  auto model=EarthModels::ModelInput<double>(DSpecMTest::modelPath().string(),norm);
  auto source=makeSensitivitySource();
  Full1D::SEM sem(model,0.1,4,2,source);
  auto parameters=makeSensitivityParameters(4,0,1000.0);
  const auto frequency=makeSensitivityFrequency(2.0*model.TimeNorm());
  const auto window=sensitivityWindow(frequency);
  const auto weights=derivativeWeights(frequency,4);
  const double radius=(sem.mesh().PR()-1000.0/sem.lengthScaleMetres())*
      sem.lengthScaleMetres();
  const std::vector<SPARSESPEC::Field3DPoint> points{{radius,0.8,1.4}};
  EXPECT_THROW(SPARSESPEC::Windowed3DKernelEngine::evaluateForwardObjective(
      sem,source,parameters,frequency,window,weights,2,4),std::invalid_argument);
  EXPECT_THROW(SPARSESPEC::Windowed3DKernelEngine::evaluate(
      sem,source,parameters,frequency,window,weights,2,4,points),std::invalid_argument);
}

TEST(Windowed3DKernelEngineTests, CachedDualShellIntegralMatchesOneDimensionalContacts) {
  EarthModels::EarthConstants<double> norm;
  const auto path=DSpecMTest::modelPath().string();
  auto model=EarthModels::ModelInput<double>(path,norm);
  auto source=makeSensitivitySource();
  Full1D::SEM sem(model,0.01,4,2,source);
  auto parameters=makeSensitivityParameters(4,0,1000.0);
  const auto frequency=makeSensitivityFrequency(model.TimeNorm());
  const auto window=sensitivityWindow(frequency);
  const auto weights=derivativeWeights(frequency,4);
  const auto expected=SPARSESPEC::WindowedSensitivityEngine::evaluate(
      sem,source,parameters,frequency,window,weights,2,4);
  const auto qhat=SPARSESPEC::WindowedDisplacementMap(
      frequency,sem.lengthScaleMetres(),window).adjoint(weights);
  const auto polar=GaussQuad::GaussLegendreQuadrature1D<double>(5);
  constexpr int nPhi=9;
  std::vector<std::array<double,2>> angles;
  std::vector<double> omegaWeights;
  for(int it=0;it<polar.N();++it)
    for(int ip=0;ip<nPhi;++ip){
      angles.push_back({std::acos(polar.X(it)),
                        2.0*static_cast<double>(EIGEN_PI)*ip/nPhi});
      omegaWeights.push_back(polar.W(it)*2.0*static_cast<double>(EIGEN_PI)/nPhi);
    }
  std::vector<SPARSESPEC::SourceShellContact3D> density(angles.size());
  for(int bin=2;bin<4;++bin){
    const auto z=frequency.w().at(bin)+std::complex<double>(0.0,-frequency.ep());
    const std::vector<SPARSESPEC::Field3DPoint> point{{
        sem.sourceRadius()*sem.lengthScaleMetres(),0.6,1.2,sem.sourceElement(source)}};
    const auto combined=SPARSESPEC::Windowed3DKernelEngine::evaluateFrequencyWithContact(
        sem,source,parameters,z,qhat.col(bin),point,angles);
    ASSERT_EQ(combined.fields.size(),1u);
    const auto &one=combined.sourceContact;
    for(std::size_t i=0;i<density.size();++i){
      density[i].C+=one[i].C; density[i].F+=one[i].F;
      density[i].L+=one[i].L;
      EXPECT_DOUBLE_EQ(one[i].A,0.0); EXPECT_DOUBLE_EQ(one[i].N,0.0);
    }
  }
  double integralC=0.0,integralF=0.0,integralL=0.0;
  for(std::size_t i=0;i<density.size();++i){
    integralC+=omegaWeights[i]*density[i].C;
    integralF+=omegaWeights[i]*density[i].F;
    integralL+=omegaWeights[i]*density[i].L;
  }
  const auto relative=[](double a,double b){
    return std::abs(a-b)/std::max({std::abs(a),std::abs(b),1e-100});
  };
  EXPECT_GT(std::abs(expected.sourcePointKernels.C),1e-100);
  EXPECT_GT(std::abs(expected.sourcePointKernels.F),1e-100);
  EXPECT_GT(std::abs(expected.sourcePointKernels.L),1e-100);
  EXPECT_LT(relative(integralC,expected.sourcePointKernels.C),1e-9)
      << "C shell=" << integralC << " 1D=" << expected.sourcePointKernels.C;
  EXPECT_LT(relative(integralF,expected.sourcePointKernels.F),1e-9)
      << "F shell=" << integralF << " 1D=" << expected.sourcePointKernels.F;
  EXPECT_LT(relative(integralL,expected.sourcePointKernels.L),1e-9)
      << "L shell=" << integralL << " 1D=" << expected.sourcePointKernels.L;
}

TEST(Windowed3DKernelEngineTests, AngularPointKernelsBridgeSolidAndFluidRadialSamples) {
  EarthModels::EarthConstants<double> norm;
  const auto path=DSpecMTest::modelPath().string();
  auto model=EarthModels::ModelInput<double>(path,norm);
  auto source=makeSensitivitySource();
  Full1D::SEM sem(model,0.25,4,2,source);
  auto parameters=makeSensitivityParameters(4,0,1000.0);
  const auto frequency=makeSensitivityFrequency(model.TimeNorm());
  const auto window=sensitivityWindow(frequency);
  const auto weights=derivativeWeights(frequency,3);
  const auto baseline=SPARSESPEC::WindowedSensitivityEngine::evaluate(
      sem,source,parameters,frequency,window,weights,2,4);
  const auto map=SPARSESPEC::WindowedDisplacementMap(
      frequency,sem.lengthScaleMetres(),window);
  const auto polar=GaussQuad::GaussLegendreQuadrature1D<double>(5);
  constexpr int nPhi=9;
  std::vector<std::array<double,2>> angles;std::vector<double> omega;
  for(int it=0;it<polar.N();++it)for(int ip=0;ip<nPhi;++ip){
    angles.push_back({std::acos(polar.X(it)),2.0*static_cast<double>(EIGEN_PI)*ip/nPhi});
    omega.push_back(polar.W(it)*2.0*static_cast<double>(EIGEN_PI)/nPhi);
  }
  const int solid=sem.sourceElement(source);
  int fluid=-1;
  for(int e=0;e<sem.mesh().NE();++e)if(sem.elementIsFluid(e)){fluid=e;break;}
  ASSERT_GE(fluid,0);
  std::vector<Eigen::Matrix<double,5,1>> solidIntegrated(2,
      Eigen::Matrix<double,5,1>::Zero());
  double fluidIntegrated=0.0;
  for(const int bin:{2,3}){
    const auto z=frequency.w().at(bin)+std::complex<double>(0.0,-frequency.ep());
    const Eigen::VectorXcd qhat=map.adjoint(weights).col(bin);
    auto workspace=SPARSESPEC::Windowed3DFieldEngine::prepareFrequencyWorkspace(
        sem,source,parameters,z,qhat);
    int solidIndex=0;
    for(const int node:{0,sem.mesh().NN()-1}){
      const double radius=sem.mesh().NodeRadius(solid,node)*sem.lengthScaleMetres();
      std::vector<SPARSESPEC::Field3DPoint> points;
      for(const auto &a:angles)points.push_back({radius,a[0],a[1],solid});
      const auto kernels=SPARSESPEC::Windowed3DKernelEngine::evaluateFrequency(workspace,points);
      Eigen::Matrix<double,5,1> integrated=Eigen::Matrix<double,5,1>::Zero();
      for(std::size_t j=0;j<kernels.size();++j){
        integrated[0]+=omega[j]*kernels[j].A;integrated[1]+=omega[j]*kernels[j].C;
        integrated[2]+=omega[j]*kernels[j].F;integrated[3]+=omega[j]*kernels[j].L;
        integrated[4]+=omega[j]*kernels[j].N;
      }
      solidIntegrated[solidIndex++]+=radius*radius*integrated;
    }
    const int node=sem.mesh().NN()/2;
    const double radius=sem.mesh().NodeRadius(fluid,node)*sem.lengthScaleMetres();
    std::vector<SPARSESPEC::Field3DPoint> points;
    for(const auto &a:angles)points.push_back({radius,a[0],a[1],fluid});
    const auto kernels=SPARSESPEC::Windowed3DKernelEngine::evaluateFrequency(workspace,points);
    double integrated=0.0;
    for(std::size_t j=0;j<kernels.size();++j){
      EXPECT_TRUE(kernels[j].fluid);
      EXPECT_DOUBLE_EQ(kernels[j].A,0.0);EXPECT_DOUBLE_EQ(kernels[j].C,0.0);
      EXPECT_DOUBLE_EQ(kernels[j].F,0.0);EXPECT_DOUBLE_EQ(kernels[j].L,0.0);
      EXPECT_DOUBLE_EQ(kernels[j].N,0.0);
      integrated+=omega[j]*kernels[j].kappa;
    }
    fluidIntegrated+=radius*radius*integrated;
  }
  int solidIndex=0;
  for(const int node:{0,sem.mesh().NN()-1}){
    const Eigen::Matrix<double,5,1> expected{baseline.kernels.A(solid,node),
        baseline.kernels.C(solid,node),baseline.kernels.F(solid,node),
        baseline.kernels.L(solid,node),baseline.kernels.N(solid,node)};
    const auto actual=solidIntegrated[solidIndex++];
    const double componentScale=expected.cwiseAbs().maxCoeff();
    ASSERT_GT(componentScale,0.0);
    for(int c=0;c<5;++c){
      const double denominator=std::abs(expected[c])>1e-6*componentScale
          ?std::abs(expected[c]):componentScale;
      EXPECT_LT(std::abs(actual[c]-expected[c])/denominator,1e-9)
          << "solid e="<<solid<<" node="<<node<<" component="<<c
          << (std::abs(expected[c])>1e-6*componentScale?" informative":" weak-scaled")
          << " actual="<<actual[c]<<" expected="<<expected[c];
    }
  }
  const int fluidNode=sem.mesh().NN()/2;
  const double expectedFluid=baseline.materialKernels.bulkModulusPa(fluid,fluidNode);
  EXPECT_LT(std::abs(fluidIntegrated-expectedFluid)/std::max(std::abs(expectedFluid),1e-100),1e-9)
      << "fluid e="<<fluid<<" node="<<fluidNode<<" integrated="<<fluidIntegrated
      <<" expected="<<expectedFluid;
}

TEST(Windowed3DKernelEngineTests, SourceShellContactCoversEachModeFamilyAndL0) {
  EarthModels::EarthConstants<double> norm;
  const auto path=DSpecMTest::modelPath().string();
  auto model=EarthModels::ModelInput<double>(path,norm);
  auto source=makeSensitivitySource();Full1D::SEM sem(model,0.25,4,2,source);
  const auto frequency=makeSensitivityFrequency(model.TimeNorm());
  const auto window=sensitivityWindow(frequency);
  const auto polar=GaussQuad::GaussLegendreQuadrature1D<double>(5);
  constexpr int nPhi=9;
  std::vector<std::array<double,2>> angles;std::vector<double> omega;
  for(int it=0;it<polar.N();++it)for(int ip=0;ip<nPhi;++ip){
    angles.push_back({std::acos(polar.X(it)),2.0*static_cast<double>(EIGEN_PI)*ip/nPhi});
    omega.push_back(polar.W(it)*2.0*static_cast<double>(EIGEN_PI)/nPhi);
  }
  for(const int type:{1,2,3,4}){
    auto parameters=makeSensitivityParameters(type,0,1000.0);
    const auto weights=derivativeWeights(frequency,3);
    const auto baseline=SPARSESPEC::WindowedSensitivityEngine::evaluate(
        sem,source,parameters,frequency,window,weights,2,4);
    const auto qhat=SPARSESPEC::WindowedDisplacementMap(
        frequency,sem.lengthScaleMetres(),window).adjoint(weights);
    double c=0.0,f=0.0,l=0.0;
    for(int bin=2;bin<4;++bin){
      const auto z=frequency.w().at(bin)+std::complex<double>(0.0,-frequency.ep());
      const auto density=SPARSESPEC::Windowed3DKernelEngine::sourceShellContactFrequency(
          sem,source,parameters,z,qhat.col(bin),angles);
      for(std::size_t j=0;j<density.size();++j){
        c+=omega[j]*density[j].C;f+=omega[j]*density[j].F;l+=omega[j]*density[j].L;
        EXPECT_DOUBLE_EQ(density[j].A,0.0);EXPECT_DOUBLE_EQ(density[j].N,0.0);
      }
    }
    const std::array<double,3> actual{c,f,l};
    const std::array<double,3> expected{baseline.sourcePointKernels.C,
        baseline.sourcePointKernels.F,baseline.sourcePointKernels.L};
    const double scale=std::max({std::abs(expected[0]),std::abs(expected[1]),
                                 std::abs(expected[2]),1e-100});
    for(std::size_t j=0;j<actual.size();++j){
      const double denominator=std::abs(expected[j])>1e-6*scale
          ?std::abs(expected[j]):scale;
      EXPECT_LT(std::abs(actual[j]-expected[j])/denominator,1e-9)
          << "type="<<type<<" contact component="<<j
          << (std::abs(expected[j])>1e-6*scale?" informative":" weak-scaled")
          << " shell="<<actual[j]<<" 1D="<<expected[j];
    }
  }
}

TEST(Windowed3DKernelEngineTests, RotationPreservesPointVolumeAndShellDensities) {
  EarthModels::EarthConstants<double> norm;
  const auto path=DSpecMTest::modelPath().string();
  auto model=EarthModels::ModelInput<double>(path,norm);
  constexpr double rotationDegrees=37.0;
  auto source=makeSensitivitySource();
  auto rotatedSource=makeSensitivitySource(30.0,15.0,147.0);
  Full1D::SEM sem(model,0.1,4,2,source);
  Full1D::SEM rotatedSem(model,0.1,4,2,rotatedSource);
  auto parameters=makeSensitivityParameters(4,0,1000.0,{30.0,65.0});
  auto rotatedParameters=makeSensitivityParameters(4,0,1000.0,{30.0,102.0});
  const auto frequency=makeSensitivityFrequency(model.TimeNorm());
  const auto window=sensitivityWindow(frequency);
  const auto weights=derivativeWeights(frequency,3);
  const auto map=SPARSESPEC::WindowedDisplacementMap(
      frequency,sem.lengthScaleMetres(),window);
  const auto qhat=map.adjoint(weights);
  const auto rotatedQhat=SPARSESPEC::WindowedDisplacementMap(
      frequency,rotatedSem.lengthScaleMetres(),window).adjoint(weights);
  const double radians=static_cast<double>(EIGEN_PI)/180.0;
  const double radius=(sem.mesh().PR()-1000.0/sem.lengthScaleMetres())*
      sem.lengthScaleMetres();
  const SPARSESPEC::Field3DPoint point{radius,0.9,1.4};
  const SPARSESPEC::Field3DPoint rotatedPoint{
      radius,0.9,1.4+rotationDegrees*radians};
  const std::vector<SPARSESPEC::Field3DPoint> points{point},rotatedPoints{rotatedPoint};
  const std::vector<std::array<double,2>> angles{{0.7,0.4},{1.2,2.3},{2.1,5.4}};
  auto rotatedAngles=angles;
  for(auto &angle:rotatedAngles)angle[1]+=rotationDegrees*radians;
  const auto base=SPARSESPEC::Windowed3DKernelEngine::evaluate(
      sem,source,parameters,frequency,window,weights,2,4,points,angles);
  const auto rotated=SPARSESPEC::Windowed3DKernelEngine::evaluate(
      rotatedSem,rotatedSource,rotatedParameters,frequency,window,weights,2,4,
      rotatedPoints,rotatedAngles);
  const auto relative=[](const auto &x,const auto &y){
    return (x-y).norm()/std::max({x.norm(),y.norm(),1e-100});
  };
  EXPECT_LT(relative(base.processedSISpectrum,rotated.processedSISpectrum),1e-9);
  EXPECT_LT(relative(base.rawDisplacementSpectrum,rotated.rawDisplacementSpectrum),1e-9);
  EXPECT_LT(std::abs(base.objective-rotated.objective)/
      std::max({std::abs(base.objective),std::abs(rotated.objective),1e-100}),1e-9);
  ASSERT_EQ(base.pointKernels.size(),1u);
  ASSERT_EQ(rotated.pointKernels.size(),1u);
  Eigen::Matrix<double,5,1> baseK,rotatedK;
  const auto &bk=base.pointKernels.front();const auto &rk=rotated.pointKernels.front();
  baseK<<bk.A,bk.C,bk.F,bk.L,bk.N;
  rotatedK<<rk.A,rk.C,rk.F,rk.L,rk.N;
  EXPECT_LT(relative(baseK,rotatedK),1e-9)<<"rotated point kernels";
  ASSERT_EQ(base.sourceContact.size(),rotated.sourceContact.size());
  for(std::size_t j=0;j<base.sourceContact.size();++j){
    Eigen::Vector3d b{base.sourceContact[j].C,base.sourceContact[j].F,base.sourceContact[j].L};
    Eigen::Vector3d r{rotated.sourceContact[j].C,rotated.sourceContact[j].F,
                      rotated.sourceContact[j].L};
    EXPECT_LT(relative(b,r),1e-9)<<"rotated shell density="<<j;
  }
  const auto z=frequency.w().at(3)+std::complex<double>(0.0,-frequency.ep());
  const SPARSESPEC::WeightedField3DPoint baseWeighted{point,1200.0,0.3};
  const SPARSESPEC::WeightedField3DPoint rotatedWeighted{rotatedPoint,1200.0,0.3};
  const auto baseVolume=SPARSESPEC::Windowed3DKernelEngine::
      integrateWeightedVolumeFrequency(sem,source,parameters,z,qhat.col(3),{baseWeighted});
  const auto rotatedVolume=SPARSESPEC::Windowed3DKernelEngine::
      integrateWeightedVolumeFrequency(rotatedSem,rotatedSource,rotatedParameters,z,
          rotatedQhat.col(3),{rotatedWeighted});
  Eigen::Matrix<double,6,1> baseV,rotatedV;
  baseV<<baseVolume.A,baseVolume.C,baseVolume.F,baseVolume.L,baseVolume.N,baseVolume.kappa;
  rotatedV<<rotatedVolume.A,rotatedVolume.C,rotatedVolume.F,rotatedVolume.L,
            rotatedVolume.N,rotatedVolume.kappa;
  EXPECT_LT(relative(baseV,rotatedV),1e-9)<<"rotated weighted volume kernels";
}

TEST(Windowed3DKernelEngineTests, SourceProfileVolumeAndContactMatchFrozenFiniteDifferences) {
  EarthModels::EarthConstants<double> norm;
  const std::string path=DSpecMTest::modelPath().string();
  auto model=EarthModels::ModelInput<double>(path,norm);
  auto source=makeSensitivitySource();
  Full1D::SEM sem(model,0.01,4,2,source);
  constexpr int firstBin=2,endBin=4;
  constexpr double materialScale=1e11;
  constexpr std::array<double,4> steps{2e-3,1e-3,5e-4,2.5e-4};
  const auto frequency=makeSensitivityFrequency(model.TimeNorm());
  const auto window=sensitivityWindow(frequency);
  const double L0=sem.lengthScaleMetres(),P0=sem.pressureScalePascals(),
               T0=sem.timeScaleSeconds(),R=sem.mesh().PR()*L0;
  ASSERT_TRUE(std::isfinite(L0)&&L0>0.0&&std::isfinite(P0)&&P0>0.0&&
              std::isfinite(T0)&&T0>0.0&&std::isfinite(R)&&R>0.0);
  const auto polar=GaussQuad::GaussLegendreQuadrature1D<double>(5);
  constexpr int nPhi=9;
  std::vector<std::array<double,2>> angles;
  std::vector<double> omega;
  for(int it=0;it<polar.N();++it)for(int ip=0;ip<nPhi;++ip){
    angles.push_back({std::acos(polar.X(it)),2.0*static_cast<double>(EIGEN_PI)*ip/nPhi});
    omega.push_back(polar.W(it)*2.0*static_cast<double>(EIGEN_PI)/nPhi);
  }
  const int sourceElement=sem.sourceElement(source);
  const double sourceRadiusSI=sem.sourceRadius()*L0;
  const double sourceLowerSI=sem.mesh().ELR(sourceElement)*L0;
  const double sourceUpperSI=sem.mesh().EUR(sourceElement)*L0;
  ASSERT_GE(sourceElement,1);
  ASSERT_LE(sourceRadiusSI,sourceUpperSI);
  auto parametersC=makeSensitivityParameters(3,0,29.9);
  auto parametersF=makeSensitivityParameters(3,0,29.9);
  auto parametersL=makeSensitivityParameters(4,0,29.9);
  struct Direction { TestModulus modulus; int type; const char *name; };
  const std::array<Direction,3> directions{{
      {TestModulus::C,3,"C"},{TestModulus::F,3,"F"},{TestModulus::L,4,"L"}}};
  const std::array<InputParameters*,3> parameterSets{{
      &parametersC,&parametersF,&parametersL}};
  const auto modelC=SPARSESPEC::WindowedSensitivityEngine::evaluate(
      sem,source,parametersC,frequency,window,derivativeWeights(frequency,3),firstBin,endBin);
  const auto modelF=SPARSESPEC::WindowedSensitivityEngine::evaluate(
      sem,source,parametersF,frequency,window,derivativeWeights(frequency,3),firstBin,endBin);
  const auto modelL=SPARSESPEC::WindowedSensitivityEngine::evaluate(
      sem,source,parametersL,frequency,window,derivativeWeights(frequency,4),firstBin,endBin);
  std::cout<<std::setprecision(17)
      <<"C_SOURCE_PROFILE protocol=accepted source_element="<<sourceElement
      <<" source_radius_m="<<sourceRadiusSI<<" source_bounds_m=["<<sourceLowerSI
      <<","<<sourceUpperSI<<"] L0_m="<<L0<<" P0_Pa="<<P0<<" T0_s="<<T0
      <<" R_m="<<R<<" dr_weights_unit=m source_C_inner_outer="
      <<sem.meshModel().C(sourceElement-1,sem.mesh().NN()-1)<<","<<sem.meshModel().C(sourceElement,0)
      <<" source_F_inner_outer="<<sem.meshModel().F(sourceElement-1,sem.mesh().NN()-1)
      <<","<<sem.meshModel().F(sourceElement,0)<<" source_L_inner_outer="
      <<sem.meshModel().L(sourceElement-1,sem.mesh().NN()-1)<<","<<sem.meshModel().L(sourceElement,0)
      <<" receiver_element_rock="<<sem.receiverElements(parametersC,false).front()
      <<" receiver_depth_km=29.9 receiver_lat_lon_deg=30,65 bins=[2,4)"
      <<" frequency_dt_s=20 output_min=0.5 taper=0.5-0.5*cos(2*pi*n/(nt-1))"
      <<" q_bin3=(0.7+0.2i,-0.13+0.19i,0) angular_rule=5x9"
      <<" profile=1+0.3*r/R Sm_Pa="<<materialScale<<" steps=0.002,0.001,0.0005,0.00025\n";
  ASSERT_NEAR(sourceRadiusSI,sourceLowerSI,1e-9*R);

  for(std::size_t directionIndex=0;directionIndex<directions.size();++directionIndex){
    const auto direction=directions[directionIndex];
    auto &parameters=*parameterSets[directionIndex];
    const auto weights=derivativeWeights(frequency,direction.type);
    const auto baseline=direction.modulus==TestModulus::C?modelC:
        direction.modulus==TestModulus::F?modelF:modelL;
    const auto baselineC=SPARSESPEC::Windowed3DKernelEngine::evaluateForwardObjective(
        sem,source,parameters,frequency,window,weights,firstBin,endBin);
    EXPECT_NEAR(baselineC.objective,baseline.objective,
        1e-9*std::max(std::abs(baseline.objective),1e-100))<<direction.name;
    const double Sobj=weights.norm()*baseline.processedSISpectrum.middleCols(
        firstBin,endBin-firstBin).norm();
    ASSERT_TRUE(std::isfinite(Sobj)&&Sobj>0.0)<<direction.name;
    Eigen::MatrixXcd qhatAll=SPARSESPEC::WindowedDisplacementMap(
        frequency,L0,window).adjoint(weights);
    double volume=0.0,shell=0.0;
    const auto gll=const_cast<EarthMesh::RadialMesh&>(sem.mesh()).GLL();
    for(int bin=firstBin;bin<endBin;++bin){
      const auto z=frequency.w().at(bin)+std::complex<double>(0.0,-frequency.ep());
      const Eigen::VectorXcd qhat=qhatAll.col(bin);
      auto workspace=SPARSESPEC::Windowed3DFieldEngine::prepareFrequencyWorkspace(
          sem,source,parameters,z,qhat);
      const auto density=SPARSESPEC::Windowed3DKernelEngine::sourceShellContact(
          workspace,source,angles);
      for(std::size_t j=0;j<density.size();++j){
        const double value=direction.modulus==TestModulus::C?density[j].C:
            direction.modulus==TestModulus::F?density[j].F:density[j].L;
        shell+=omega[j]*value;
      }
      for(int e=0;e<sem.mesh().NE();++e){
        if(sem.elementIsFluid(e))continue;
        std::vector<SPARSESPEC::WeightedField3DPoint> quadrature;
        quadrature.reserve(static_cast<std::size_t>(sem.mesh().NN())*angles.size());
        for(int i=0;i<sem.mesh().NN();++i){
          const double radiusModel=sem.mesh().NodeRadius(e,i);
          const double radiusSI=radiusModel*L0;
          const double radialWeight=0.5*sem.mesh().EW(e)*gll.W(i)*L0*
              (1.0+0.3*radiusSI/R);
          for(std::size_t j=0;j<angles.size();++j)
            quadrature.push_back({{radiusSI,angles[j][0],angles[j][1],e},
                                  radialWeight,omega[j]});
        }
        const auto integrated=SPARSESPEC::Windowed3DKernelEngine::
            integrateWeightedVolumeFrequency(workspace,quadrature,64);
        volume+=direction.modulus==TestModulus::C?integrated.C:
            direction.modulus==TestModulus::F?integrated.F:integrated.L;
      }
    }
    const double profileAtSource=1.0+0.3*sourceRadiusSI/R;
    const double predictedVolume=materialScale*volume;
    const double predictedShell=materialScale*profileAtSource*shell;
    const double predicted=predictedVolume+predictedShell;
    const double shellReference=direction.modulus==TestModulus::C?
        baseline.sourcePointKernels.C:direction.modulus==TestModulus::F?
        baseline.sourcePointKernels.F:baseline.sourcePointKernels.L;
    const double shellError=std::abs(shell-shellReference)/
        std::max({std::abs(shell),std::abs(shellReference),1e-100});
    EXPECT_LT(shellError,1e-9)<<direction.name<<" 3D shell="<<shell
        <<" 1D sourcePoint="<<shellReference;
    auto objective=[&](double amplitude){
      PerturbedModulusModel changedModel(path,norm,direction.modulus,amplitude);
      Full1D::SEM changedSem(changedModel,0.01,4,2,source);
      return SPARSESPEC::Windowed3DKernelEngine::evaluateForwardObjective(
          changedSem,source,parameters,frequency,window,weights,firstBin,endBin).objective;
    };
    std::array<double,4> centered{},richardson{};
    std::array<double,4> plus{},minus{};
    for(std::size_t i=0;i<steps.size();++i){
      plus[i]=objective(steps[i]);minus[i]=objective(-steps[i]);
      centered[i]=(plus[i]-minus[i])/(2.0*steps[i]);
    }
    for(std::size_t i=0;i<3;++i)
      richardson[i]=(4.0*centered[i+1]-centered[i])/3.0;
    const double richardsonScale=std::max({std::abs(richardson[0]),
        std::abs(richardson[1]),std::abs(richardson[2])});
    const double spread=richardsonScale>0.0?std::max({
        std::abs(richardson[0]-richardson[1]),
        std::abs(richardson[0]-richardson[2]),
        std::abs(richardson[1]-richardson[2])})/richardsonScale:
        std::numeric_limits<double>::infinity();
    const double scaledAbsolute=std::abs(richardson[2]-predicted)/Sobj;
    const bool informative=std::abs(predicted)>1e-6*Sobj;
    std::cout<<"C_SOURCE_PROFILE property="<<direction.name<<" J0="<<baselineC.objective
        <<" Sobj="<<Sobj<<" volume="<<predictedVolume<<" shell="<<predictedShell
        <<" shell_angular="<<shell<<" shell_1D="<<shellReference<<" shell_rel="<<shellError
        <<" dJpred="<<predicted<<" informative="<<informative<<" D="
        <<centered[0]<<","<<centered[1]<<","<<centered[2]<<","<<centered[3]
        <<" R="<<richardson[0]<<","<<richardson[1]<<","<<richardson[2]
        <<" plateau_spread="<<spread<<" scaled_abs="<<scaledAbsolute<<" plus="
        <<plus[0]<<","<<plus[1]<<","<<plus[2]<<","<<plus[3]<<" minus="
        <<minus[0]<<","<<minus[1]<<","<<minus[2]<<","<<minus[3]<<"\n";
    EXPECT_TRUE(std::isfinite(predicted)&&std::isfinite(richardson[2]));
    EXPECT_LE(spread,1e-6)<<direction.name<<" centered-derivative plateau spread="<<spread;
    if(informative&&spread<=1e-6)
      EXPECT_LE(std::abs(richardson[2]-predicted)/std::abs(predicted),1e-6)
          <<direction.name<<" predicted="<<predicted<<" FD="<<richardson[2];
  }
}
