#ifndef DSPECM1D_SENSITIVITY_KERNELS_H
#define DSPECM1D_SENSITIVITY_KERNELS_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <string>
#include <vector>

#include "SEM/SEM.h"
#include "EquilibratedSparseLU.h"
#include "ParamInfo.h"
#include "SpecHelpers.h"
#include "frequency_info/FreqFull.h"
#include "frequency_info/PostprocessFunctions.h"

namespace SPARSESPEC {

/**
 * @brief Real time-window measurement map from raw displacement spectra to SI.
 *
 * Input spectra are nondimensional displacement on the complete nonnegative
 * FFT grid. Damping is removed once, before the fixed time window; output is
 * the dt_SI-scaled real-to-complex transform.
 */
class WindowedDisplacementMap {
public:
  using ComplexMatrix = Eigen::MatrixXcd;

  WindowedDisplacementMap(const SpectraSolver::FreqFull &frequency,
                          double lengthScaleMetres,
                          const Eigen::VectorXd &timeWindow)
      : frequency_(frequency), lengthScale_(lengthScaleMetres),
        window_(timeWindow) {
    if (!(std::isfinite(lengthScale_) && lengthScale_ > 0.0))
      throw std::invalid_argument("Displacement length scale must be positive and finite.");
    if (frequency_.nt() < 2 || !(std::isfinite(frequency_.dt()) && frequency_.dt() > 0.0) ||
        !(std::isfinite(frequency_.df()) && frequency_.df() > 0.0) ||
        !(std::isfinite(frequency_.timeNorm()) && frequency_.timeNorm() > 0.0) ||
        !std::isfinite(frequency_.ep()))
      throw std::invalid_argument("FreqFull must contain a finite, initialized FFT grid.");
    if (window_.size() != frequency_.nt() || !window_.allFinite())
      throw std::invalid_argument("Time window must be finite and match FreqFull::nt().");
  }

  int frequencyCount() const { return frequency_.nt() / 2 + 1; }
  double lengthScaleMetres() const { return lengthScale_; }
  const Eigen::VectorXd &timeWindow() const { return window_; }

  ComplexMatrix forward(const ComplexMatrix &nondimensionalDisplacement) const {
    requireShape(nondimensionalDisplacement.cols() == frequencyCount(),
                 "Displacement spectrum has the wrong FFT-bin count.");
    requireShape(nondimensionalDisplacement.allFinite(),
                 "Displacement spectrum must be finite.");
    ComplexMatrix output(nondimensionalDisplacement.rows(), frequencyCount());
    for (Eigen::Index row = 0; row < nondimensionalDisplacement.rows(); ++row) {
      Eigen::MatrixXcd one(1, frequencyCount());
      one.row(0) = nondimensionalDisplacement.row(row);
      Eigen::MatrixXd time = processfunctions::freq2time(one, frequency_);
      time.row(0) *= lengthScale_;
      time.row(0).array() *= window_.transpose().array();
      output.row(row) = processfunctions::rawtime2freq(
          time, frequency_.nt(), frequency_.dt() * frequency_.timeNorm());
    }
    requireShape(output.allFinite(),
                 "Windowed displacement map produced non-finite values.");
    return output;
  }

  /** Apply P* under Re(x^dagger y), including DC and Nyquist weights. */
  ComplexMatrix adjoint(const ComplexMatrix &measurementWeights) const {
    requireShape(measurementWeights.cols() == frequencyCount(),
                 "Measurement weights have the wrong FFT-bin count.");
    requireShape(measurementWeights.allFinite(),
                 "Measurement weights must be finite.");
    const int nt = frequency_.nt();
    const double dtSI = frequency_.dt() * frequency_.timeNorm();
    Eigen::VectorXd damping(nt);
    for (int i = 0; i < nt; ++i)
      damping[i] = std::exp(frequency_.ep() * frequency_.dt() * i);
    Eigen::VectorXd eta = Eigen::VectorXd::Constant(frequencyCount(), 2.0);
    eta[0] = 1.0;
    if (nt % 2 == 0)
      eta[frequencyCount() - 1] = 1.0;

    ComplexMatrix result(measurementWeights.rows(), frequencyCount());
    for (Eigen::Index row = 0; row < measurementWeights.rows(); ++row) {
      Eigen::MatrixXcd q(1, frequencyCount());
      q.row(0) = measurementWeights.row(row).array() / eta.transpose().array();
      Eigen::MatrixXd time = processfunctions::rawfreq2time(q, nt);
      time.row(0) *= dtSI;
      time.row(0).array() *= damping.transpose().array() * window_.transpose().array();
      Eigen::MatrixXcd spectrum = processfunctions::rawtime2freq(time, nt, 1.0);
      result.row(row) = lengthScale_ * frequency_.df() *
                        (spectrum.row(0).array() * eta.transpose().array());
    }
    requireShape(result.allFinite(), "Measurement adjoint produced non-finite values.");
    return result;
  }

private:
  static void requireShape(bool condition, const char *message) {
    if (!condition)
      throw std::invalid_argument(message);
  }
  SpectraSolver::FreqFull frequency_;
  double lengthScale_;
  Eigen::VectorXd window_;
};

enum class SensitivitySolveFamily { Radial, Spheroidal, Toroidal };

namespace sensitivity_detail {
inline void requireConstructedSEM(const Full1D::SEM &sem) {
  if (sem.mesh().NE() == 0 || sem.mesh().NN() < 2)
    throw std::invalid_argument("SEM must be constructed before sensitivity operations.");
}
} // namespace sensitivity_detail

/** Full-mesh fields for a compact multi-column source and its Hermitian adjoint. */
struct FrequencyAdjointFields {
  Eigen::MatrixXcd forward;
  Eigen::MatrixXcd adjoint;
  Eigen::MatrixXcd output;
  double objective = 0.0;
};

/**
 * @brief Fixed-mesh elastic solve for a scalar weighted displacement objective.
 *
 * The receiver map has shape (output rows, SEM dofs), and weights has shape
 * (output rows, compact source columns). The same sparse factorization is used
 * for all compact RHS columns and its Hermitian transpose. Full SEM matrices
 * are used; no turning-point truncation or frequency-to-frequency reuse occurs.
 */
class FixedMeshAdjointSolver {
public:
  using Complex = std::complex<double>;
  using SparseComplex = Eigen::SparseMatrix<Complex>;
  using Solver = EquilibratedSparseLU<Complex>;

  static FrequencyAdjointFields solveFrequency(
      const Full1D::SEM &sem, SensitivitySolveFamily family, int degree,
      Complex z, const Eigen::MatrixXcd &compactRhs,
      const Eigen::MatrixXcd &receiverMap,
      const Eigen::MatrixXcd &measurementWeights) {
    sensitivity_detail::requireConstructedSEM(sem);
    if (compactRhs.rows() == 0 || compactRhs.cols() == 0 ||
        !compactRhs.allFinite())
      throw std::invalid_argument("Compact source RHS must be nonempty and finite.");
    if (receiverMap.rows() == 0 || receiverMap.cols() != compactRhs.rows() ||
        !receiverMap.allFinite())
      throw std::invalid_argument("Receiver map dimensions or values are invalid.");
    if (measurementWeights.rows() != receiverMap.rows() ||
        measurementWeights.cols() != compactRhs.cols() ||
        !measurementWeights.allFinite())
      throw std::invalid_argument("Measurement weight dimensions or values are invalid.");

    SparseComplex h = systemMatrix(sem, family, degree, z);
    if (h.rows() != compactRhs.rows() || h.cols() != compactRhs.rows())
      throw std::invalid_argument("Source columns do not match the selected SEM operator.");
    Eigen::MatrixXcd displacementReceiver =
        outputFactor(0, z, Complex(0.0, 1.0)) * receiverMap;

    Solver forwardSolver;
    forwardSolver.compute(h);
    if (forwardSolver.info() != Eigen::Success)
      throw std::runtime_error("Forward SEM factorization failed.");
    FrequencyAdjointFields result;
    result.forward = forwardSolver.solve(compactRhs);
    if (forwardSolver.info() != Eigen::Success || !result.forward.allFinite())
      throw std::runtime_error("Forward SEM solve failed or returned non-finite fields.");
    result.output = displacementReceiver * result.forward;
    if (!result.output.allFinite())
      throw std::runtime_error("Displacement receiver map returned non-finite values.");

    const Eigen::MatrixXcd receiverRhs = displacementReceiver.adjoint() * measurementWeights;
    result.adjoint = forwardSolver.solveAdjoint(receiverRhs);
    if (forwardSolver.info() != Eigen::Success || !result.adjoint.allFinite())
      throw std::runtime_error("Hermitian-adjoint SEM solve failed or returned non-finite fields.");
    result.objective = (measurementWeights.conjugate().cwiseProduct(result.output)).sum().real();
    if (!std::isfinite(result.objective))
      throw std::runtime_error("Displacement objective is non-finite.");
    return result;
  }

  static Eigen::MatrixXcd solveForwardFrequency(
      const Full1D::SEM &sem, SensitivitySolveFamily family, int degree,
      Complex z, const Eigen::MatrixXcd &compactRhs) {
    if (!compactRhs.allFinite() || compactRhs.rows() == 0 || compactRhs.cols() == 0)
      throw std::invalid_argument("Compact source RHS must be nonempty and finite.");
    SparseComplex h = systemMatrix(sem, family, degree, z);
    if (h.rows() != compactRhs.rows())
      throw std::invalid_argument("Source columns do not match the selected SEM operator.");
    Solver solver;
    solver.compute(h);
    Eigen::MatrixXcd fields = solver.solve(compactRhs);
    if (solver.info() != Eigen::Success || !fields.allFinite())
      throw std::runtime_error("Forward SEM solve failed or returned non-finite fields.");
    return fields;
  }

private:
  static SparseComplex systemMatrix(const Full1D::SEM &sem,
                                    SensitivitySolveFamily family, int degree,
                                    Complex z) {
    sensitivity_detail::requireConstructedSEM(sem);
    if (!std::isfinite(z.real()) || !std::isfinite(z.imag()) || z == Complex{})
      throw std::invalid_argument("Complex frequency must be finite and nonzero.");
    if (degree < 0 || degree > sem.maxDegree() ||
        (family != SensitivitySolveFamily::Radial && degree < 1))
      throw std::invalid_argument("Mode degree is invalid for the requested family or SEM.");
    Eigen::SparseMatrix<double> hReal, pReal;
    switch (family) {
    case SensitivitySolveFamily::Radial:
      if (degree != 0)
        throw std::invalid_argument("Radial modes require degree zero.");
      hReal = sem.hR(); pReal = sem.pR(); break;
    case SensitivitySolveFamily::Spheroidal:
      hReal = sem.hS(degree); pReal = sem.pS(degree); break;
    case SensitivitySolveFamily::Toroidal:
      hReal = sem.hTk(degree); pReal = sem.pTk(degree); break;
    }
    if (hReal.rows() != hReal.cols() || pReal.rows() != hReal.rows() ||
        pReal.cols() != hReal.cols())
      throw std::invalid_argument("Selected SEM stiffness and mass matrices are inconsistent.");
    SparseComplex h = hReal.cast<Complex>() - z * z * pReal.cast<Complex>();
    h.makeCompressed();
    for (int col = 0; col < h.outerSize(); ++col)
      for (SparseComplex::InnerIterator it(h, col); it; ++it)
        if (!std::isfinite(it.value().real()) || !std::isfinite(it.value().imag()))
          throw std::invalid_argument("Selected SEM operator contains non-finite values.");
    return h;
  }

};

/** Five absolute-SI elastic modulus integrands on the fixed SEM GLL mesh. */
struct ModulusKernelSamples {
  Eigen::MatrixXd A, C, F, L, N;
  Eigen::MatrixXd radiusMetres;
  Eigen::MatrixXd drWeightsMetres;
  Eigen::ArrayXX<bool> solid;
  std::vector<int> layer;
  std::string units = "objective per Pa per metre; multiply by drWeightsMetres and delta modulus in Pa";
};

/** Density and fluid bulk integrands on one-sided, element-local GLL samples. */
struct MaterialKernelSamples {
  Eigen::MatrixXd bulkModulusPa;
  Eigen::MatrixXd densityKgPerM3;
  Eigen::MatrixXd radiusMetres;
  Eigen::MatrixXd drWeightsMetres;
  Eigen::ArrayXX<bool> fluid;
  std::vector<int> layer;
  std::string bulkUnits = "objective per Pa per metre; fluid kappa perturbs A=C=F equally";
  std::string densityUnits = "objective per (kg/m^3) per metre; fixed moduli and geometry";
  std::string densityRepresentation =
      "independent one-sided element-local GLL samples; polynomial per element";
};

/** Integrated source-contact contribution, in objective per pascal. */
struct ModulusSourcePointTerms { double A = 0.0, C = 0.0, F = 0.0, L = 0.0, N = 0.0; };

namespace sensitivity_detail {
inline Eigen::MatrixXd derivativeMatrix(const Full1D::SEM &sem) {
  const auto &mesh = sem.mesh();
  // RadialMesh currently exposes its read-only quadrature through a
  // non-const accessor. This cast is limited to that const-preserving read.
  const auto gll = const_cast<EarthMesh::RadialMesh &>(mesh).GLL();
  auto basis = Interpolation::LagrangeBasis(gll.Points());
  Eigen::MatrixXd derivative(mesh.NN(), mesh.NN());
  for (int i = 0; i < mesh.NN(); ++i)
    for (int j = 0; j < mesh.NN(); ++j)
      derivative(i, j) = basis.template Evaluate<1>(j, gll.X(i));
  return derivative;
}

inline std::complex<double> inner(const Eigen::MatrixXcd &test,
                                  const Eigen::MatrixXcd &field) {
  return (test.conjugate().cwiseProduct(field)).sum();
}

inline void validateFields(const Full1D::SEM &sem,
                           const FrequencyAdjointFields &fields,
                           int expectedDofs) {
  if (fields.forward.rows() != expectedDofs ||
      fields.adjoint.rows() != expectedDofs ||
      fields.forward.cols() == 0 ||
      fields.forward.cols() != fields.adjoint.cols() ||
      !fields.forward.allFinite() || !fields.adjoint.allFinite())
    throw std::invalid_argument("Forward/adjoint fields do not match the selected SEM family.");
  if (!(std::isfinite(sem.lengthScaleMetres()) && sem.lengthScaleMetres() > 0.0 &&
        std::isfinite(sem.pressureScalePascals()) && sem.pressureScalePascals() > 0.0))
    throw std::invalid_argument("SEM SI length and pressure scales must be positive and finite.");
}
} // namespace sensitivity_detail

/**
 * Return one-sided element-local modulus integrands and dr quadrature weights.
 * Fluid samples are explicitly marked unsupported and carry structural zeros.
 */
inline ModulusKernelSamples modulusKernelSamples(
    const Full1D::SEM &sem, SensitivitySolveFamily family, int degree,
    const FrequencyAdjointFields &fields,
    const Eigen::MatrixXcd *elementLocalLift = nullptr,
    int liftElement = -1) {
  sensitivity_detail::requireConstructedSEM(sem);
  const auto &mesh = sem.mesh();
  const int ne = mesh.NE(), nq = mesh.NN();
  if (degree < 0 || degree > sem.maxDegree() ||
      (family == SensitivitySolveFamily::Radial && degree != 0) ||
      (family != SensitivitySolveFamily::Radial && degree < 1))
    throw std::invalid_argument("Degree does not match the requested sensitivity family or SEM.");
  int dofs = 0;
  if (family == SensitivitySolveFamily::Radial)
    dofs = static_cast<int>(sem.pR().rows());
  else if (family == SensitivitySolveFamily::Spheroidal)
    dofs = static_cast<int>(sem.pS(degree).rows());
  else
    dofs = static_cast<int>(sem.pTk(degree).rows());
  sensitivity_detail::validateFields(sem, fields, dofs);

  ModulusKernelSamples result;
  result.A = Eigen::MatrixXd::Zero(ne, nq);
  result.C = Eigen::MatrixXd::Zero(ne, nq);
  result.F = Eigen::MatrixXd::Zero(ne, nq);
  result.L = Eigen::MatrixXd::Zero(ne, nq);
  result.N = Eigen::MatrixXd::Zero(ne, nq);
  result.radiusMetres.resize(ne, nq);
  result.drWeightsMetres.resize(ne, nq);
  result.solid.resize(ne, nq);
  result.layer.resize(ne);
  const double lengthScale = sem.lengthScaleMetres();
  const double modulusScale = sem.pressureScalePascals();
  const double integrandScale = 1.0 / (lengthScale * modulusScale);
  const Eigen::MatrixXd D = sensitivity_detail::derivativeMatrix(sem);
  const double lambda = static_cast<double>(degree) * (degree + 1.0);

  for (int e = 0; e < ne; ++e) {
    result.layer[e] = mesh.LayerNumber(e);
    const double jacobian = 0.5 * mesh.EW(e);
    const bool isSolid = !sem.elementIsFluid(e);
    for (int i = 0; i < nq; ++i) {
      result.radiusMetres(e, i) = mesh.NodeRadius(e, i) * lengthScale;
      result.drWeightsMetres(e, i) = jacobian *
          const_cast<EarthMesh::RadialMesh &>(mesh).GLL().W(i) * lengthScale;
      result.solid(e, i) = isSolid;
    }
    if (!isSolid || (family == SensitivitySolveFamily::Toroidal &&
                     (e < sem.el() || e >= sem.eu())))
      continue;
    Eigen::MatrixXcd u(nq, fields.forward.cols()), ua(nq, fields.adjoint.cols());
    Eigen::MatrixXcd v = Eigen::MatrixXcd::Zero(nq, fields.forward.cols());
    Eigen::MatrixXcd va = Eigen::MatrixXcd::Zero(nq, fields.adjoint.cols());
    for (int i = 0; i < nq; ++i) {
      if (family == SensitivitySolveFamily::Radial) {
        auto index = sem.ltgR(0, e, i);
        u.row(i) = fields.forward.row(index);
        ua.row(i) = fields.adjoint.row(index);
        if (elementLocalLift && e == liftElement)
          u.row(i) += elementLocalLift->row(index);
      } else if (family == SensitivitySolveFamily::Spheroidal) {
        auto ui = sem.ltgS(0, e, i), vi = sem.ltgS(1, e, i);
        u.row(i) = fields.forward.row(ui);
        ua.row(i) = fields.adjoint.row(ui);
        v.row(i) = fields.forward.row(vi);
        va.row(i) = fields.adjoint.row(vi);
        if (elementLocalLift && e == liftElement) {
          u.row(i) += elementLocalLift->row(ui);
          v.row(i) += elementLocalLift->row(vi);
        }
      } else {
        auto wi = sem.ltgT(e, i);
        u.row(i) = fields.forward.row(wi);
        ua.row(i) = fields.adjoint.row(wi);
        if (elementLocalLift && e == liftElement)
          u.row(i) += elementLocalLift->row(wi);
      }
    }
    const Eigen::MatrixXcd du = (2.0 / mesh.EW(e)) * D.cast<std::complex<double>>() * u;
    const Eigen::MatrixXcd dua = (2.0 / mesh.EW(e)) * D.cast<std::complex<double>>() * ua;
    const Eigen::MatrixXcd dv = (2.0 / mesh.EW(e)) * D.cast<std::complex<double>>() * v;
    const Eigen::MatrixXcd dva = (2.0 / mesh.EW(e)) * D.cast<std::complex<double>>() * va;
    const Eigen::MatrixXcd T = 2.0 * u - lambda * v;
    const Eigen::MatrixXcd Ta = 2.0 * ua - lambda * va;
    const Eigen::MatrixXcd S = [&] {
      Eigen::MatrixXcd value(nq, fields.forward.cols());
      for (int i = 0; i < nq; ++i)
        value.row(i) = mesh.NodeRadius(e, i) * dv.row(i) + u.row(i) - v.row(i);
      return value;
    }();
    const Eigen::MatrixXcd Sa = [&] {
      Eigen::MatrixXcd value(nq, fields.adjoint.cols());
      for (int i = 0; i < nq; ++i)
        value.row(i) = mesh.NodeRadius(e, i) * dva.row(i) + ua.row(i) - va.row(i);
      return value;
    }();
    for (int i = 0; i < nq; ++i) {
      const double r = mesh.NodeRadius(e, i);

      if (family == SensitivitySolveFamily::Toroidal) {
        Eigen::MatrixXcd strain = r * du.row(i).eval() - u.row(i).eval();
        Eigen::MatrixXcd strainA = r * dua.row(i).eval() - ua.row(i).eval();
        result.L(e, i) = -integrandScale * lambda *
                         sensitivity_detail::inner(strainA, strain).real();
        result.N(e, i) = -integrandScale * lambda * (lambda - 2.0) *
                         sensitivity_detail::inner(ua.row(i), u.row(i)).real();
      } else {
        Eigen::MatrixXcd dUi = du.row(i), dUai = dua.row(i);
        Eigen::MatrixXcd Ti = T.row(i), Tai = Ta.row(i);
        Eigen::MatrixXcd Si = S.row(i), Sai = Sa.row(i);
        const auto cTerm = r * r * sensitivity_detail::inner(dUai, dUi);
        const auto fTerm = r * (sensitivity_detail::inner(Tai, dUi) +
                                sensitivity_detail::inner(dUai, Ti));
        const auto aTerm = sensitivity_detail::inner(Tai, Ti);
        const auto nTerm = -2.0 * (
            2.0 * sensitivity_detail::inner(ua.row(i), u.row(i)) -
            lambda * (sensitivity_detail::inner(ua.row(i), v.row(i)) +
                      sensitivity_detail::inner(va.row(i), u.row(i))) +
            lambda * sensitivity_detail::inner(va.row(i), v.row(i)));
        const auto lTerm = lambda * sensitivity_detail::inner(Sai, Si);
        result.C(e, i) = -integrandScale * cTerm.real();
        result.F(e, i) = -integrandScale * fTerm.real();
        result.A(e, i) = -integrandScale * aTerm.real();
        result.N(e, i) = -integrandScale * nTerm.real();
        result.L(e, i) = -integrandScale * lTerm.real();
      }
    }
  }
  if (!result.A.allFinite() || !result.C.allFinite() || !result.F.allFinite() ||
      !result.L.allFinite() || !result.N.allFinite())
    throw std::runtime_error("Modulus kernel integrands contain non-finite values.");
  return result;
}

/** Return fluid bulk and complete fixed-modulus density kernels. */
inline MaterialKernelSamples materialKernelSamples(
    const Full1D::SEM &sem, SensitivitySolveFamily family, int degree,
    std::complex<double> z, const FrequencyAdjointFields &fields,
    const Eigen::MatrixXcd *elementLocalLift = nullptr,
    int liftElement = -1) {
  sensitivity_detail::requireConstructedSEM(sem);
  const auto &mesh = sem.mesh();
  const int ne = mesh.NE(), nq = mesh.NN();
  if (degree < 0 || degree > sem.maxDegree() ||
      (family == SensitivitySolveFamily::Radial && degree != 0) ||
      (family != SensitivitySolveFamily::Radial && degree < 1) ||
      !std::isfinite(z.real()) || !std::isfinite(z.imag()) || z == std::complex<double>{})
    throw std::invalid_argument("Degree or frequency does not match material sensitivity family.");
  const int dofs = family == SensitivitySolveFamily::Radial ? sem.pR().rows() :
      family == SensitivitySolveFamily::Spheroidal ? sem.pS(degree).rows() : sem.pTk(degree).rows();
  sensitivity_detail::validateFields(sem, fields, dofs);
  const double length = sem.lengthScaleMetres(), pressure = sem.pressureScalePascals();
  const double time = sem.timeScaleSeconds();
  const double densityScale = pressure * time * time / (length * length);
  if (!(std::isfinite(densityScale) && densityScale > 0.0))
    throw std::invalid_argument("SEM SI density scale must be positive and finite.");
  const auto gll = const_cast<EarthMesh::RadialMesh &>(mesh).GLL();
  const auto basis = Interpolation::LagrangeBasis(gll.Points());
  const Eigen::MatrixXd D = sensitivity_detail::derivativeMatrix(sem);
  const double lambda = static_cast<double>(degree) * (degree + 1.0);
  const double grav = sem.gravityCoefficient();
  MaterialKernelSamples out;
  out.bulkModulusPa = Eigen::MatrixXd::Zero(ne, nq);
  out.densityKgPerM3 = Eigen::MatrixXd::Zero(ne, nq);
  out.radiusMetres.resize(ne, nq); out.drWeightsMetres.resize(ne, nq);
  out.fluid.resize(ne, nq); out.layer.resize(ne);
  Eigen::MatrixXd gravityLoads = Eigen::MatrixXd::Zero(ne, nq);
  Eigen::MatrixXd densityGradient = Eigen::MatrixXd::Zero(ne, nq);
  const double modelDensityPerSI = 1.0 / (length * densityScale);
  for (int e = 0; e < ne; ++e) {
    out.layer[e] = mesh.LayerNumber(e);
    const double jac = 0.5 * mesh.EW(e);
    const bool isFluid = sem.elementIsFluid(e);
    for (int i = 0; i < nq; ++i) {
      out.radiusMetres(e, i) = mesh.NodeRadius(e, i) * length;
      out.drWeightsMetres(e, i) = jac * gll.W(i) * length;
      out.fluid(e, i) = isFluid;
    }
    if (family == SensitivitySolveFamily::Toroidal && (e < sem.el() || e >= sem.eu())) continue;
    Eigen::MatrixXcd u(nq, fields.forward.cols()), ua(nq, fields.adjoint.cols());
    Eigen::MatrixXcd v = Eigen::MatrixXcd::Zero(nq, fields.forward.cols());
    Eigen::MatrixXcd va = Eigen::MatrixXcd::Zero(nq, fields.adjoint.cols());
    Eigen::MatrixXcd phi = Eigen::MatrixXcd::Zero(nq, fields.forward.cols());
    Eigen::MatrixXcd phia = Eigen::MatrixXcd::Zero(nq, fields.adjoint.cols());
    for (int i = 0; i < nq; ++i) {
      if (family == SensitivitySolveFamily::Radial) {
        auto iu = sem.ltgR(0, e, i), ip = sem.ltgR(1, e, i);
        u.row(i) = fields.forward.row(iu); ua.row(i) = fields.adjoint.row(iu);
        phi.row(i) = fields.forward.row(ip); phia.row(i) = fields.adjoint.row(ip);
        if (elementLocalLift && e == liftElement) {
          u.row(i) += elementLocalLift->row(iu);
          phi.row(i) += elementLocalLift->row(ip);
        }
      } else if (family == SensitivitySolveFamily::Spheroidal) {
        auto iu = sem.ltgS(0, e, i), iv = sem.ltgS(1, e, i), ip = sem.ltgS(2, e, i);
        u.row(i) = fields.forward.row(iu); ua.row(i) = fields.adjoint.row(iu);
        v.row(i) = fields.forward.row(iv); va.row(i) = fields.adjoint.row(iv);
        phi.row(i) = fields.forward.row(ip); phia.row(i) = fields.adjoint.row(ip);
        if (elementLocalLift && e == liftElement) {
          u.row(i) += elementLocalLift->row(iu);
          v.row(i) += elementLocalLift->row(iv);
          phi.row(i) += elementLocalLift->row(ip);
        }
      } else {
        auto iw = sem.ltgT(e, i); u.row(i) = fields.forward.row(iw); ua.row(i) = fields.adjoint.row(iw);
        if (elementLocalLift && e == liftElement)
          u.row(i) += elementLocalLift->row(iw);
      }
    }
    const Eigen::MatrixXcd du = (2.0 / mesh.EW(e)) * D.cast<std::complex<double>>() * u;
    const Eigen::MatrixXcd dua = (2.0 / mesh.EW(e)) * D.cast<std::complex<double>>() * ua;
    const Eigen::MatrixXcd dphi = (2.0 / mesh.EW(e)) * D.cast<std::complex<double>>() * phi;
    const Eigen::MatrixXcd dphia = (2.0 / mesh.EW(e)) * D.cast<std::complex<double>>() * phia;
    for (int i = 0; i < nq; ++i) {
      const double r = mesh.NodeRadius(e, i), qweight = jac * gll.W(i);
      const auto uu = sensitivity_detail::inner(ua.row(i), u.row(i));
      const auto uv = sensitivity_detail::inner(ua.row(i), v.row(i)) +
                      sensitivity_detail::inner(va.row(i), u.row(i));
      const auto vv = sensitivity_detail::inner(va.row(i), v.row(i));
      if (isFluid && family != SensitivitySolveFamily::Toroidal) {
        if (family == SensitivitySolveFamily::Radial) {
          const auto a = 4.0 * uu;
          const auto c = r * r * sensitivity_detail::inner(dua.row(i), du.row(i));
          const auto f = r * (sensitivity_detail::inner(2.0 * ua.row(i), du.row(i)) +
                              sensitivity_detail::inner(dua.row(i), 2.0 * u.row(i)));
          out.bulkModulusPa(e, i) = -(a + c + f).real() / (length * pressure);
        } else {
          const Eigen::MatrixXcd T = 2.0 * u.row(i).eval() - lambda * v.row(i).eval();
          const Eigen::MatrixXcd Ta = 2.0 * ua.row(i).eval() - lambda * va.row(i).eval();
          const auto kappa = sensitivity_detail::inner(Ta, T) +
              r * (sensitivity_detail::inner(Ta, du.row(i)) + sensitivity_detail::inner(dua.row(i), T)) +
              r * r * sensitivity_detail::inner(dua.row(i), du.row(i));
          out.bulkModulusPa(e, i) = -kappa.real() / (length * pressure);
        }
      }
      double local = 0.0;
      if (family == SensitivitySolveFamily::Toroidal) {
        local = (-z * z * lambda * r * r * uu).real();
      } else {
        const double rho = sem.meshModel().Density(e, i), g = sem.meshModel().Gravity(e, i);
        const double eightPiG = 8.0 * static_cast<double>(EIGEN_PI) * grav * rho * r * r;
        const auto explicitGravity = (eightPiG - 4.0 * g * r) * uu +
            lambda * g * r * uv + r * r * (sensitivity_detail::inner(dphia.row(i), u.row(i)) +
            sensitivity_detail::inner(ua.row(i), dphi.row(i))) + lambda * r *
            (sensitivity_detail::inner(phia.row(i), v.row(i)) + sensitivity_detail::inner(va.row(i), phi.row(i)));
        const auto inertia = -z * z * r * r * (uu + lambda * vv);
        local = (explicitGravity + inertia).real();
        gravityLoads(e, i) = qweight * (4.0 * rho * r * uu.real() - lambda * rho * r * uv.real());
      }
      densityGradient(e, i) = -qweight * local;
    }
  }
  // Transpose production gravity with a reverse element/node cumulative sum.
  // Its traversal follows MeshModel's cumulative element endpoints exactly.
  auto pointLoad = [&](int e, int i) {
    const double r = mesh.NodeRadius(e, i);
    return r > 0.0 ? gravityLoads(e, i) *
        4.0 * static_cast<double>(EIGEN_PI) * grav / (r * r) : 0.0;
  };
  double above = 0.0;
  for (int e = ne - 1; e >= 0; --e) {
    const double lo = mesh.ELR(e), hi = mesh.EUR(e);
    const auto &knots = sem.meshModel().DensityKnots(e);
    for (int seg = nq - 2; seg >= 0; --seg) {
      const double a = mesh.NodeRadius(e, seg), b = mesh.NodeRadius(e, seg + 1);
      above += pointLoad(e, seg + 1);
      double left = a;
      auto integrate = [&](double right) {
        const double half = 0.5 * (right - left);
        for (int q = 0; q < nq; ++q) {
          const double s = std::lerp(left, right, 0.5 * (gll.X(q) + 1.0));
          const double weight = half * gll.W(q) * above * s * s;
          const double x = 2.0 * (s - lo) / (hi - lo) - 1.0;
          for (int j = 0; j < nq; ++j)
            densityGradient(e, j) += weight * basis.Evaluate(j, x);
        }
        left = right;
      };
      for (double knot : knots) {
        if (knot > left && knot < b) integrate(knot);
      }
      integrate(b);
    }
    above += pointLoad(e, 0);
  }
  for (int e = 0; e < ne; ++e) for (int i = 0; i < nq; ++i) {
    const double nodalWeight = 0.5 * mesh.EW(e) * gll.W(i);
    out.densityKgPerM3(e, i) = densityGradient(e, i) / nodalWeight * modelDensityPerSI;
  }
  if (!out.bulkModulusPa.allFinite() || !out.densityKgPerM3.allFinite())
    throw std::runtime_error("Material kernel integrands contain non-finite values.");
  return out;
}

struct WindowedSensitivityResult {
  Eigen::MatrixXcd rawDisplacementSpectrum;
  Eigen::MatrixXcd processedSISpectrum;
  double objective = 0.0;
  int firstSolvedBin = 0;
  int endSolvedBin = 0;
  ModulusKernelSamples kernels;
  MaterialKernelSamples materialKernels;
  ModulusSourcePointTerms sourcePointKernels;
};

/**
 * @brief End-to-end fixed-source, fixed-receiver windowed displacement gradient.
 *
 * `weightsSI` is a complex matrix over (vertical,north,east) receiver rows and
 * the full RFFT grid. `firstBin:endBin` bounds all forward and adjoint solves;
 * all other raw bins are exactly zero. Each mode/bin receives a fresh sparse
 * factorization on the full elastic SEM mesh.
 */
class WindowedSensitivityEngine {
  struct Projection {
    SensitivitySolveFamily family;
    int degree;
    Eigen::MatrixXcd rhs;
    Eigen::MatrixXcd receiver;
    Eigen::MatrixXcd sourceCoefficients;
    Eigen::MatrixXcd sourceLift;
    Eigen::MatrixXcd directSourceOutput;
    int sourceElement = -1;
    int outputRows;
  };

public:
  static WindowedSensitivityResult evaluate(
      Full1D::SEM &sem, SourceInfo::EarthquakeCMT &source,
      InputParameters &parameters, const SpectraSolver::FreqFull &frequency,
      const Eigen::VectorXd &timeWindow, const Eigen::MatrixXcd &weightsSI,
      int firstBin, int endBin) {
    sensitivity_detail::requireConstructedSEM(sem);
    if (!sem.sourceMatches(source))
      throw std::invalid_argument(
          "Windowed sensitivities require an SEM aligned to the current source; construct it with the source-aware overload.");
    if (parameters.attenuation())
      throw std::invalid_argument("Sensitivity kernels do not support attenuation.");
    if (parameters.cowlingFrequencyMhz() != 0.0 || parameters.gravitation() != 2)
      throw std::invalid_argument("Sensitivity kernels require full-gravity elastic formulation.");
    if (parameters.output_type() != 0)
      throw std::invalid_argument("Sensitivity kernels support displacement output only.");
    if (parameters.num_receivers() <= 0)
      throw std::invalid_argument("Sensitivity kernels require at least one receiver.");
    if (parameters.lmin() < 0 || parameters.lmax() < parameters.lmin() ||
        parameters.lmax() > sem.maxDegree())
      throw std::invalid_argument("Sensitivity degree range is invalid for the constructed SEM.");
    const double timeScale = sem.timeScaleSeconds();
    if (!(std::isfinite(timeScale) && timeScale > 0.0))
      throw std::invalid_argument("SEM time scale must be positive and finite.");
    const double timeScaleError =
        std::abs(frequency.timeNorm() - timeScale) / std::abs(timeScale);
    if (!(std::isfinite(timeScaleError) && timeScaleError < 1e-12))
      throw std::invalid_argument("FreqFull and SEM must use the same model time scale.");
    if (weightsSI.rows() != 3 * parameters.num_receivers() ||
        weightsSI.cols() != frequency.nt() / 2 + 1 || !weightsSI.allFinite())
      throw std::invalid_argument("Sensitivity weights must be finite on all receiver components and FFT bins.");
    if (firstBin < 0 || endBin > weightsSI.cols() || firstBin >= endBin)
      throw std::invalid_argument("Sensitivity solved frequency interval is invalid.");
    if (!(std::isfinite(sem.lengthScaleMetres()) && sem.lengthScaleMetres() > 0.0 &&
          std::isfinite(sem.pressureScalePascals()) && sem.pressureScalePascals() > 0.0))
      throw std::invalid_argument("SEM SI length and stress scales must be positive and finite.");

    const int lmax = parameters.lmax();
    const auto flags = resolveModeFlags(parameters.type(), parameters.lmin(), lmax);
    if (!flags.inc_rad && !flags.inc_sph && !flags.inc_tor)
      throw std::invalid_argument("Sensitivity mode selection contains no supported mode family.");
    WindowedDisplacementMap measurement(frequency, sem.lengthScaleMetres(), timeWindow);

    WindowedSensitivityResult result;
    result.firstSolvedBin = firstBin;
    result.endSolvedBin = endBin;
    result.rawDisplacementSpectrum = Eigen::MatrixXcd::Zero(
        weightsSI.rows(), weightsSI.cols());
    result.kernels = emptySamples(sem);
    result.materialKernels = emptyMaterialSamples(sem);

    const auto w = frequency.w();
    const Complex ieps(0.0, -frequency.ep());
    forEachProjection(sem, source, parameters, lmax, flags,
                      [&](const Projection &projection) {
      for (int bin = firstBin; bin < endBin; ++bin) {
        const Complex z = w.at(bin) + ieps;
        Projection sourceProjection = projection;
        if (sourceProjection.family == SensitivitySolveFamily::Spheroidal) {
          const auto treatment = sem.sourceTreatedSpheroidal(
              source, sourceProjection.degree, z);
          sourceProjection.rhs = treatment.force;
          sourceProjection.sourceLift = treatment.lift;
          sourceProjection.sourceElement = treatment.sourceElement;
          sourceProjection.directSourceOutput = sem.sourceLiftReceiver(
              parameters, sourceProjection.degree, treatment.lift,
              treatment.sourceElement, false);
        } else if (sourceProjection.family == SensitivitySolveFamily::Radial) {
          const auto treatment = sem.sourceTreatedRadial(source, z);
          sourceProjection.rhs = treatment.force;
          sourceProjection.sourceLift = treatment.lift;
          sourceProjection.sourceElement = treatment.sourceElement;
          sourceProjection.directSourceOutput = sem.sourceLiftReceiverRadial(
              parameters, treatment.lift, treatment.sourceElement);
        } else {
          const auto treatment = sem.sourceTreatedToroidal(
              source, sourceProjection.degree, z);
          sourceProjection.rhs = treatment.force;
          sourceProjection.sourceLift = treatment.lift;
          sourceProjection.sourceElement = treatment.sourceElement;
          sourceProjection.directSourceOutput = sem.sourceLiftReceiverToroidal(
              parameters, sourceProjection.degree, treatment.lift,
              treatment.sourceElement);
        }
        const auto forward = FixedMeshAdjointSolver::solveForwardFrequency(
            sem, sourceProjection.family, sourceProjection.degree, z, sourceProjection.rhs);
        addProjectionOutput(result.rawDisplacementSpectrum, bin, sourceProjection,
                            forward, z);
      }
    });
    result.processedSISpectrum = measurement.forward(result.rawDisplacementSpectrum);
    result.objective = weightsSI.conjugate().cwiseProduct(
                           result.processedSISpectrum).sum().real();
    Eigen::MatrixXcd adjointWeights = measurement.adjoint(weightsSI);
    if (firstBin > 0)
      adjointWeights.leftCols(firstBin).setZero();
    if (endBin < adjointWeights.cols())
      adjointWeights.rightCols(adjointWeights.cols() - endBin).setZero();

    forEachProjection(sem, source, parameters, lmax, flags,
                      [&](const Projection &projection) {
      for (int bin = firstBin; bin < endBin; ++bin) {
        const Complex z = w.at(bin) + ieps;
        Projection sourceProjection = projection;
        if (sourceProjection.family == SensitivitySolveFamily::Spheroidal) {
          const auto treatment = sem.sourceTreatedSpheroidal(
              source, sourceProjection.degree, z);
          sourceProjection.rhs = treatment.force;
          sourceProjection.sourceLift = treatment.lift;
          sourceProjection.sourceElement = treatment.sourceElement;
          sourceProjection.directSourceOutput = sem.sourceLiftReceiver(
              parameters, sourceProjection.degree, treatment.lift,
              treatment.sourceElement, false);
        } else if (sourceProjection.family == SensitivitySolveFamily::Radial) {
          const auto treatment = sem.sourceTreatedRadial(source, z);
          sourceProjection.rhs = treatment.force;
          sourceProjection.sourceLift = treatment.lift;
          sourceProjection.sourceElement = treatment.sourceElement;
          sourceProjection.directSourceOutput = sem.sourceLiftReceiverRadial(
              parameters, treatment.lift, treatment.sourceElement);
        } else {
          const auto treatment = sem.sourceTreatedToroidal(
              source, sourceProjection.degree, z);
          sourceProjection.rhs = treatment.force;
          sourceProjection.sourceLift = treatment.lift;
          sourceProjection.sourceElement = treatment.sourceElement;
          sourceProjection.directSourceOutput = sem.sourceLiftReceiverToroidal(
              parameters, sourceProjection.degree, treatment.lift,
              treatment.sourceElement);
        }
        const auto modeWeights = compactWeights(sourceProjection, adjointWeights.col(bin));
        const auto fields = FixedMeshAdjointSolver::solveFrequency(
            sem, sourceProjection.family, sourceProjection.degree, z, sourceProjection.rhs,
            sourceProjection.receiver, modeWeights);
        if (sourceProjection.family == SensitivitySolveFamily::Spheroidal)
          addSpheroidalSourcePointKernels(
              result.sourcePointKernels, sem, parameters, sourceProjection,
              z, fields, modeWeights);
        else if (sourceProjection.family == SensitivitySolveFamily::Radial)
          addRadialSourcePointKernels(result.sourcePointKernels, sem, parameters,
                                      sourceProjection, z, fields, modeWeights);
        else
          addToroidalSourcePointKernels(result.sourcePointKernels, sem, parameters,
                                        sourceProjection, z, fields, modeWeights);
        addModulusKernelSamples(
            result.kernels,
            modulusKernelSamples(sem, sourceProjection.family, sourceProjection.degree,
                fields, sourceProjection.sourceLift.size() ?
                    &sourceProjection.sourceLift : nullptr,
                sourceProjection.sourceElement));
        addMaterialKernelSamples(result.materialKernels,
            materialKernelSamples(sem, sourceProjection.family, sourceProjection.degree,
                z, fields, sourceProjection.sourceLift.size() ?
                    &sourceProjection.sourceLift : nullptr,
                sourceProjection.sourceElement));
      }
    });
    return result;
  }

private:
  using Complex = std::complex<double>;

  static void addSpheroidalSourcePointKernels(
      ModulusSourcePointTerms &sum, Full1D::SEM &sem,
      InputParameters &parameters, const Projection &projection, Complex z,
      const FrequencyAdjointFields &fields,
      const Eigen::MatrixXcd &modeWeights) {
    const int degree = projection.degree, e = projection.sourceElement;
    const int nq = sem.mesh().NN(), u = sem.ltgS(0, e, 0), v = sem.ltgS(1, e, 0);
    const double r = sem.sourceRadius();
    const double lambda = static_cast<double>(degree) * (degree + 1.0);
    const auto &model = sem.meshModel();
    const double C = model.C(e, 0), F = model.F(e, 0), L = model.L(e, 0);
    const Complex jU = projection.sourceLift(u, 0);
    const Complex jV = projection.sourceLift(v, 2);
    Eigen::MatrixXcd deltaLiftC = Eigen::MatrixXcd::Zero(
        projection.sourceLift.rows(), projection.sourceLift.cols());
    Eigen::MatrixXcd deltaLiftL = Eigen::MatrixXcd::Zero(
        projection.sourceLift.rows(), projection.sourceLift.cols());
    Eigen::MatrixXcd deltaC = Eigen::MatrixXcd::Zero(
        projection.sourceLift.rows(), projection.sourceLift.cols());
    Eigen::MatrixXcd deltaF = deltaC;
    deltaLiftC(u, 0) = -jU / C;
    deltaC(u, 0) = 2.0 * r * F * jU / C;
    deltaC(v, 0) = -lambda * r * F * jU / C;
    deltaF(u, 0) = -2.0 * r * jU;
    deltaF(v, 0) = lambda * r * jU;
    deltaLiftL(v, 2) = -jV / L;
    Eigen::MatrixXcd localDeltaLift = Eigen::MatrixXcd::Zero(3 * nq, 4);
    Eigen::MatrixXcd localAdjoint = localDeltaLift;
    for (int field = 0; field < 3; ++field)
      for (int i = 0; i < nq; ++i) {
        const int index = sem.ltgS(field, e, i);
        localDeltaLift.row(field * nq + i) = deltaLiftC.row(index);
        localAdjoint.row(field * nq + i) = fields.adjoint.row(index);
      }
    const Eigen::MatrixXcd localOperator =
        Full1D::source_treatment_detail::spheroidalElementOperator(
            sem, e, degree, z, {}, false);
    const auto contraction = [&](const Eigen::MatrixXcd &pointLoad,
                                 const Eigen::MatrixXcd &globalLiftVariation) {
      Eigen::MatrixXcd localLiftVariation(3 * nq, globalLiftVariation.cols());
      for (int field = 0; field < 3; ++field)
        for (int i = 0; i < nq; ++i)
          localLiftVariation.row(field * nq + i) =
              globalLiftVariation.row(sem.ltgS(field, e, i));
      const auto residual = localOperator * localLiftVariation;
      const double weak = sensitivity_detail::inner(fields.adjoint, pointLoad).real() -
                          sensitivity_detail::inner(localAdjoint, residual).real();
      const auto direct = sem.sourceLiftReceiver(
          parameters, degree, globalLiftVariation, e, false);
      const Complex factor = outputFactor(0, z, Complex(0.0, 1.0));
      return weak + (modeWeights.conjugate().cwiseProduct(factor * direct)).sum().real();
    };
    const double scale = 1.0 / sem.pressureScalePascals();
    sum.C += scale * contraction(deltaC, deltaLiftC);
    sum.F += scale * contraction(deltaF, Eigen::MatrixXcd::Zero(
        deltaLiftC.rows(), deltaLiftC.cols()));
    sum.L += scale * contraction(Eigen::MatrixXcd::Zero(
        deltaLiftC.rows(), deltaLiftC.cols()), deltaLiftL);
  }

  static void addRadialSourcePointKernels(
      ModulusSourcePointTerms &sum, Full1D::SEM &sem,
      InputParameters &parameters, const Projection &projection, Complex z,
      const FrequencyAdjointFields &fields,
      const Eigen::MatrixXcd &modeWeights) {
    const int e = projection.sourceElement, nq = sem.mesh().NN();
    const int u = sem.ltgR(0, e, 0);
    const double r = sem.sourceRadius(), C = sem.meshModel().C(e, 0);
    const double F = sem.meshModel().F(e, 0);
    const Complex j = projection.sourceLift(u, 0);
    Eigen::MatrixXcd deltaLiftC = Eigen::MatrixXcd::Zero(
        projection.sourceLift.rows(), 1);
    Eigen::MatrixXcd deltaC = deltaLiftC, deltaF = deltaLiftC;
    deltaLiftC(u, 0) = -j / C;
    deltaC(u, 0) = 2.0 * r * F * j / C;
    deltaF(u, 0) = -2.0 * r * j;
    Eigen::MatrixXcd localAdjoint(2 * nq, 1), localDeltaLift(2 * nq, 1);
    for (int field = 0; field < 2; ++field)
      for (int i = 0; i < nq; ++i) {
        const int index = sem.ltgR(field, e, i);
        localAdjoint(field * nq + i, 0) = fields.adjoint(index, 0);
        localDeltaLift(field * nq + i, 0) = deltaLiftC(index, 0);
      }
    const auto localOperator = Full1D::source_treatment_detail::
        radialElementOperator(sem, e, z);
    const Complex factor = outputFactor(0, z, Complex(0.0, 1.0));
    const double scale = 1.0 / sem.pressureScalePascals();
    const auto contract = [&](const Eigen::MatrixXcd &pointLoad,
                              const Eigen::MatrixXcd &globalLift) {
      Eigen::MatrixXcd localLift(2 * nq, globalLift.cols());
      for (int field = 0; field < 2; ++field)
        for (int i = 0; i < nq; ++i)
          localLift.row(field * nq + i) =
              globalLift.row(sem.ltgR(field, e, i));
      const auto direct = sem.sourceLiftReceiverRadial(
          parameters, globalLift, e);
      return sensitivity_detail::inner(fields.adjoint, pointLoad).real() -
             sensitivity_detail::inner(localAdjoint, localOperator * localLift).real() +
             (modeWeights.conjugate().cwiseProduct(factor * direct)).sum().real();
    };
    sum.C += scale * contract(deltaC, deltaLiftC);
    sum.F += scale * contract(deltaF,
        Eigen::MatrixXcd::Zero(deltaLiftC.rows(), 1));
  }

  static void addToroidalSourcePointKernels(
      ModulusSourcePointTerms &sum, Full1D::SEM &sem,
      InputParameters &parameters, const Projection &projection, Complex z,
      const FrequencyAdjointFields &fields,
      const Eigen::MatrixXcd &modeWeights) {
    const int e = projection.sourceElement, nq = sem.mesh().NN();
    const int source = sem.ltgT(e, 0);
    const double L = sem.meshModel().L(e, 0);
    const Complex j = projection.sourceLift(source, 0);
    Eigen::MatrixXcd deltaLift = Eigen::MatrixXcd::Zero(
        projection.sourceLift.rows(), projection.sourceLift.cols());
    deltaLift(source, 0) = -j / L;
    Eigen::MatrixXcd localAdjoint(nq, fields.adjoint.cols());
    Eigen::MatrixXcd localLift(nq, deltaLift.cols());
    for (int i = 0; i < nq; ++i) {
      const int index = sem.ltgT(e, i);
      localAdjoint.row(i) = fields.adjoint.row(index);
      localLift.row(i) = deltaLift.row(index);
    }
    const auto localOperator = Full1D::source_treatment_detail::
        toroidalElementOperator(sem, e, projection.degree, z);
    const auto direct = sem.sourceLiftReceiverToroidal(
        parameters, projection.degree, deltaLift, e);
    const Complex factor = outputFactor(0, z, Complex(0.0, 1.0));
    const double point = -sensitivity_detail::inner(
        localAdjoint, localOperator * localLift).real() +
        (modeWeights.conjugate().cwiseProduct(factor * direct)).sum().real();
    sum.L += point / sem.pressureScalePascals();
  }

  static void addModulusKernelSamples(ModulusKernelSamples &sum,
                                      const ModulusKernelSamples &term) {
    if (sum.A.size() == 0) {
      sum = term;
      return;
    }
    if (sum.A.rows() != term.A.rows() || sum.A.cols() != term.A.cols() ||
        sum.layer != term.layer || (sum.solid != term.solid).any() ||
        !(sum.radiusMetres.array() == term.radiusMetres.array()).all() ||
        !(sum.drWeightsMetres.array() == term.drWeightsMetres.array()).all())
      throw std::invalid_argument("Cannot sum sensitivity samples from different SEM meshes.");
    sum.A += term.A; sum.C += term.C; sum.F += term.F;
    sum.L += term.L; sum.N += term.N;
  }

  static ModulusKernelSamples emptySamples(const Full1D::SEM &sem) {
    ModulusKernelSamples out;
    const auto &mesh = sem.mesh();
    out.A = Eigen::MatrixXd::Zero(mesh.NE(), mesh.NN());
    out.C = Eigen::MatrixXd::Zero(mesh.NE(), mesh.NN());
    out.F = Eigen::MatrixXd::Zero(mesh.NE(), mesh.NN());
    out.L = Eigen::MatrixXd::Zero(mesh.NE(), mesh.NN());
    out.N = Eigen::MatrixXd::Zero(mesh.NE(), mesh.NN());
    out.radiusMetres.resize(mesh.NE(), mesh.NN());
    out.drWeightsMetres.resize(mesh.NE(), mesh.NN());
    out.solid.resize(mesh.NE(), mesh.NN());
    out.layer.resize(mesh.NE());
    for (int e = 0; e < mesh.NE(); ++e) {
      out.layer[e] = mesh.LayerNumber(e);
      for (int i = 0; i < mesh.NN(); ++i) {
        out.radiusMetres(e, i) = mesh.NodeRadius(e, i) * sem.lengthScaleMetres();
        out.drWeightsMetres(e, i) = 0.5 * mesh.EW(e) *
            const_cast<EarthMesh::RadialMesh &>(mesh).GLL().W(i) * sem.lengthScaleMetres();
        out.solid(e, i) = !sem.elementIsFluid(e);
      }
    }
    return out;
  }

  static MaterialKernelSamples emptyMaterialSamples(const Full1D::SEM &sem) {
    MaterialKernelSamples out;
    const auto &mesh = sem.mesh();
    out.bulkModulusPa = Eigen::MatrixXd::Zero(mesh.NE(), mesh.NN());
    out.densityKgPerM3 = Eigen::MatrixXd::Zero(mesh.NE(), mesh.NN());
    out.radiusMetres.resize(mesh.NE(), mesh.NN());
    out.drWeightsMetres.resize(mesh.NE(), mesh.NN());
    out.fluid.resize(mesh.NE(), mesh.NN());
    out.layer.resize(mesh.NE());
    for (int e = 0; e < mesh.NE(); ++e) {
      out.layer[e] = mesh.LayerNumber(e);
      for (int i = 0; i < mesh.NN(); ++i) {
        out.radiusMetres(e, i) = mesh.NodeRadius(e, i) * sem.lengthScaleMetres();
        out.drWeightsMetres(e, i) = 0.5 * mesh.EW(e) *
            const_cast<EarthMesh::RadialMesh &>(mesh).GLL().W(i) * sem.lengthScaleMetres();
        out.fluid(e, i) = sem.elementIsFluid(e);
      }
    }
    return out;
  }

  static void addMaterialKernelSamples(MaterialKernelSamples &sum,
                                       const MaterialKernelSamples &term) {
    if (sum.bulkModulusPa.rows() != term.bulkModulusPa.rows() ||
        sum.bulkModulusPa.cols() != term.bulkModulusPa.cols() ||
        sum.layer != term.layer || (sum.fluid != term.fluid).any() ||
        !(sum.radiusMetres.array() == term.radiusMetres.array()).all() ||
        !(sum.drWeightsMetres.array() == term.drWeightsMetres.array()).all())
      throw std::invalid_argument("Cannot sum material samples from different SEM meshes.");
    sum.bulkModulusPa += term.bulkModulusPa;
    sum.densityKgPerM3 += term.densityKgPerM3;
  }

  static Eigen::MatrixXcd embedReceiver(const Full1D::SEM &sem,
                                       const Eigen::MatrixXcd &local,
                                       std::size_t low, int dofs) {
    if (low + local.cols() > dofs)
      throw std::invalid_argument("Existing receiver slice exceeds the SEM operator.");
    Eigen::MatrixXcd full = Eigen::MatrixXcd::Zero(local.rows(), dofs);
    full.block(0, low, local.rows(), local.cols()) = local;
    return full;
  }

  template <class Callback>
  static void forEachProjection(
      Full1D::SEM &sem, SourceInfo::EarthquakeCMT &source,
      InputParameters &parameters, int lmax, SPARSESPEC::ModeFlags flags,
      Callback &&callback) {
    const int nrec = parameters.num_receivers();
    const int outputRows = 3 * nrec;
    const int nq = sem.mesh().NN();
    const int firstDegree = std::max(1, parameters.lmin());
    if (flags.inc_rad) {
      Projection p;
      p.family = SensitivitySolveFamily::Radial; p.degree = 0;
      p.rhs = sem.calculateForceR(source);
      p.directSourceOutput = Eigen::MatrixXcd::Zero(outputRows, p.rhs.cols());
      p.outputRows = outputRows;
      p.sourceCoefficients = Eigen::MatrixXcd::Ones(outputRows, 1);
      p.receiver = Eigen::MatrixXcd::Zero(outputRows, p.rhs.rows());
      for (int r = 0; r < nrec; ++r)
        p.receiver.row(3 * r) = sem.rvZR(parameters, r).transpose();
      callback(p);
    }
    ParamInfo angular(parameters, lmax);
    for (int degree = firstDegree; degree <= parameters.lmax(); ++degree) {
      if (flags.inc_sph) {
        Projection p;
        p.family = SensitivitySolveFamily::Spheroidal; p.degree = degree;
        p.rhs = sem.calculateForceAll(source, degree, false);
        const Eigen::MatrixXcd sourceCoefficients =
            sem.calculateForceCoefficients(source, degree);
        const auto recElems = sem.receiverElements(parameters, false);
        const auto low = sem.ltgS(0, recElems.front(), 0);
        p.receiver = embedReceiver(sem, sem.rvBaseFull(parameters, degree, false),
                                   low, static_cast<int>(sem.hS(degree).rows()));
        p.sourceCoefficients = angular.rvFullSph(degree) * sourceCoefficients;
        p.outputRows = outputRows;
        callback(p);
      }
      if (flags.inc_tor) {
        Projection p;
        p.family = SensitivitySolveFamily::Toroidal; p.degree = degree;
        p.rhs = sem.calculateForceAllT(source, degree);
        p.sourceLift = Eigen::MatrixXcd::Zero(p.rhs.rows(), p.rhs.cols());
        p.directSourceOutput = Eigen::MatrixXcd::Zero(outputRows, p.rhs.cols());
        const Eigen::MatrixXcd sourceCoefficients =
            sem.calculateForceCoefficientsT(source, degree);
        const auto recElems = sem.receiverElements(parameters, true);
        const auto low = sem.ltgT(recElems.front(), 0);
        const auto local = sem.rvBaseFullT(parameters, degree);
        p.receiver = embedReceiver(sem, local, low,
                                   static_cast<int>(sem.hTk(degree).rows()));
        p.sourceCoefficients = angular.rvFullTor(degree) * sourceCoefficients;
        p.outputRows = outputRows;
        callback(p);
      }
    }
  }

  static Eigen::MatrixXcd compactWeights(const Projection &p,
                                         const Eigen::VectorXcd &q) {
    if (q.size() != p.outputRows || p.sourceCoefficients.rows() != p.outputRows)
      throw std::invalid_argument("Compact observation map has inconsistent rows.");
    Eigen::MatrixXcd weights(p.outputRows, p.rhs.cols());
    for (int row = 0; row < p.outputRows; ++row)
      for (Eigen::Index col = 0; col < p.rhs.cols(); ++col)
        weights(row, col) = q[row] * std::conj(p.sourceCoefficients(row, col));
    return weights;
  }

  static void addProjectionOutput(Eigen::MatrixXcd &spectrum, int bin,
                                  const Projection &p,
                                  const Eigen::MatrixXcd &forward,
                                  Complex z) {
    const Complex factor = outputFactor(0, z, Complex(0.0, 1.0));
    if (p.receiver.rows() != p.outputRows || p.receiver.cols() != forward.rows() ||
        p.sourceCoefficients.rows() != p.outputRows ||
        p.sourceCoefficients.cols() != forward.cols())
      throw std::invalid_argument("Compact forward/receiver maps have inconsistent dimensions.");
    Eigen::MatrixXcd base = p.receiver * forward;
    if (p.directSourceOutput.size() != 0)
      base += p.directSourceOutput;
    spectrum.col(bin) += factor *
        (p.sourceCoefficients.cwiseProduct(base)).rowwise().sum();
  }

};

} // namespace SPARSESPEC

#endif // DSPECM1D_SENSITIVITY_KERNELS_H
