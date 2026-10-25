#include <gtest/gtest.h>

#include <array>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <DSpecM1D/SensitivityKernels>
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
                                         double receiverDepthKm = 3.0) {
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
  options.receivers = {{30.0, 65.0}};
  const auto file = DSpecMTest::writeFile(
      temp.path() / "sensitivity_params.txt",
      DSpecMTest::makeParameterText(DSpecMTest::modelPath().string(), options));
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
}
