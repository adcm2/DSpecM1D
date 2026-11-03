#pragma once

// Diagnostic-only algorithm prototypes. This header intentionally does not
// alter or replace the accepted production engine.
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <Eigen/SVD>
#include "profile_equilibrated_sparse_lu.hpp"
#include "profile.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <complex>
#include <cstdint>
#include <map>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <vector>

#ifdef DSFAST_HAVE_CBLAS
#include <cblas.h>
#endif

namespace dsfast {
using Complex = std::complex<double>;

struct Point {
  double radius = 0, theta = 0, phi = 0;
  int element = -1;
};

/** Semantic plan for an immutable physical tensor grid. Angle and radial
 * identity are independent; duplicate element-side traces remain distinct. */
class TensorGridPlan {
public:
  using AngleKey = std::pair<std::uint64_t, std::uint64_t>;
  using TraceKey = std::pair<int, std::uint64_t>;
  explicit TensorGridPlan(const std::vector<Point> &points) {
    pointAngle_.reserve(points.size());
    pointTrace_.reserve(points.size());
    for (const auto &p : points) {
      const AngleKey ak{bits(p.theta), bits(p.phi)};
      const TraceKey rk{p.element, bits(p.radius)};
      auto [ai, anew] = angleIds_.try_emplace(ak, angleIds_.size());
      if (anew) angles_.push_back({p.theta, p.phi});
      auto [ri, rnew] = traceIds_.try_emplace(rk, traceIds_.size());
      if (rnew) traces_.push_back({p.element, p.radius});
      pointAngle_.push_back(ai->second);
      pointTrace_.push_back(ri->second);
    }
  }
  static std::uint64_t bits(double x) { return std::bit_cast<std::uint64_t>(x); }
  const std::vector<std::array<double, 2>> &angles() const { return angles_; }
  const std::vector<std::pair<int, double>> &traces() const { return traces_; }
  const std::vector<std::size_t> &pointAngle() const { return pointAngle_; }
  const std::vector<std::size_t> &pointTrace() const { return pointTrace_; }
  std::size_t bytes() const {
    return angles_.capacity() * sizeof(angles_[0]) + traces_.capacity() * sizeof(traces_[0]) +
      (pointAngle_.capacity() + pointTrace_.capacity()) * sizeof(std::size_t);
  }
private:
  std::map<AngleKey, std::size_t> angleIds_;
  std::map<TraceKey, std::size_t> traceIds_;
  std::vector<std::array<double, 2>> angles_;
  std::vector<std::pair<int, double>> traces_;
  std::vector<std::size_t> pointAngle_, pointTrace_;
};

/** Real angular basis times complex coefficient vector, separate re/im sums. */
inline Complex splitDot(const double *basis, const Complex *coeff, std::size_t n) {
  double re = 0.0, im = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    re += basis[i] * coeff[i].real();
    im += basis[i] * coeff[i].imag();
  }
  return {re, im};
}

enum class DotBackend { Scalar, Eigen, OpenBlas };

/** Row-major real basis [point,q] by row-major complex coefficients [q,batch]. */
inline Eigen::MatrixXcd synthesize(const Eigen::MatrixXd &basis,
                                   const Eigen::MatrixXcd &coeff,
                                   DotBackend backend) {
  mathprofile::Scope total("tensor_angular_gemm_real_and_packing");
  mathprofile::count("tensor_real_gemm_m",static_cast<std::uint64_t>(basis.rows()));
  mathprofile::count("tensor_real_gemm_n",static_cast<std::uint64_t>(coeff.cols()));
  mathprofile::count("tensor_real_gemm_k",static_cast<std::uint64_t>(basis.cols()));
  mathprofile::count("tensor_real_gemm_batches");
  if (basis.cols() != coeff.rows()) throw std::invalid_argument("synthesis dimensions differ");
  Eigen::MatrixXcd out = Eigen::MatrixXcd::Zero(basis.rows(), coeff.cols());
  if (backend == DotBackend::Scalar) {
    for (Eigen::Index p = 0; p < basis.rows(); ++p)
      for (Eigen::Index b = 0; b < coeff.cols(); ++b) {
        double re = 0, im = 0;
        for (Eigen::Index q = 0; q < basis.cols(); ++q) {
          re += basis(p, q) * coeff(q, b).real();
          im += basis(p, q) * coeff(q, b).imag();
        }
        out(p, b) = {re, im};
      }
    return out;
  }
  if (backend == DotBackend::Eigen) return (basis.cast<Complex>() * coeff).eval();
#ifdef DSFAST_HAVE_CBLAS
  using RowMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  RowMatrix a = basis, cr = coeff.real(), ci = coeff.imag();
  RowMatrix rr = RowMatrix::Zero(basis.rows(), coeff.cols());
  RowMatrix ri = RowMatrix::Zero(basis.rows(), coeff.cols());
  const double alpha=1.0,beta=0.0;
  {mathprofile::Scope t("actual_DGEMM_real_part");mathprofile::count("dgemm_calls");
      mathprofile::count("tensor_actual_gemm_macs",static_cast<std::uint64_t>(basis.rows())*static_cast<std::uint64_t>(coeff.cols())*static_cast<std::uint64_t>(basis.cols()));cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,static_cast<int>(basis.rows()),
      static_cast<int>(coeff.cols()),static_cast<int>(basis.cols()),alpha,a.data(),
      static_cast<int>(basis.cols()),cr.data(),static_cast<int>(coeff.cols()),beta,rr.data(),
      static_cast<int>(coeff.cols()));}
  {mathprofile::Scope t("actual_DGEMM_imag_part");mathprofile::count("dgemm_calls");
      mathprofile::count("tensor_actual_gemm_macs",static_cast<std::uint64_t>(basis.rows())*static_cast<std::uint64_t>(coeff.cols())*static_cast<std::uint64_t>(basis.cols()));cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,static_cast<int>(basis.rows()),
      static_cast<int>(coeff.cols()),static_cast<int>(basis.cols()),alpha,a.data(),
      static_cast<int>(basis.cols()),ci.data(),static_cast<int>(coeff.cols()),beta,ri.data(),
      static_cast<int>(coeff.cols()));}
  return rr.cast<Complex>() + Complex(0,1)*ri.cast<Complex>();
#else
  throw std::runtime_error("OpenBLAS backend requested without DSFAST_HAVE_CBLAS");
#endif
}

inline Eigen::MatrixXcd synthesizeComplex(const Eigen::MatrixXcd &basis,
                                         const Eigen::MatrixXcd &coeff,
                                         DotBackend backend) {
  mathprofile::Scope total("tensor_angular_gemm_complex_and_packing");
  mathprofile::count("tensor_complex_gemm_m",static_cast<std::uint64_t>(basis.rows()));
  mathprofile::count("tensor_complex_gemm_n",static_cast<std::uint64_t>(coeff.cols()));
  mathprofile::count("tensor_complex_gemm_k",static_cast<std::uint64_t>(basis.cols()));
  mathprofile::count("tensor_complex_gemm_batches");
  if(basis.cols()!=coeff.rows())throw std::invalid_argument("complex synthesis dimensions differ");
  if(backend==DotBackend::Scalar){
    Eigen::MatrixXcd out=Eigen::MatrixXcd::Zero(basis.rows(),coeff.cols());
    for(Eigen::Index p=0;p<basis.rows();++p)for(Eigen::Index b=0;b<coeff.cols();++b){
      double re=0,im=0;
      for(Eigen::Index q=0;q<basis.cols();++q){
        const auto x=basis(p,q),y=coeff(q,b);
        re+=x.real()*y.real()-x.imag()*y.imag();
        im+=x.real()*y.imag()+x.imag()*y.real();
      }
      out(p,b)={re,im};
    }
    return out;
  }
  if(backend==DotBackend::Eigen)return (basis*coeff).eval();
#ifdef DSFAST_HAVE_CBLAS
  using RowComplex=Eigen::Matrix<Complex,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
  RowComplex a=basis,c=coeff,out=RowComplex::Zero(basis.rows(),coeff.cols());
  const Complex alpha{1,0},beta{0,0};
  {mathprofile::Scope t("actual_ZGEMM_complex");mathprofile::count("zgemm_calls");
      mathprofile::count("tensor_actual_gemm_macs",static_cast<std::uint64_t>(basis.rows())*static_cast<std::uint64_t>(coeff.cols())*static_cast<std::uint64_t>(basis.cols()));cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,static_cast<int>(basis.rows()),
      static_cast<int>(coeff.cols()),static_cast<int>(basis.cols()),&alpha,a.data(),
      static_cast<int>(basis.cols()),c.data(),static_cast<int>(coeff.cols()),&beta,out.data(),
      static_cast<int>(coeff.cols()));}
  return Eigen::MatrixXcd(out);
#else
  throw std::runtime_error("OpenBLAS backend requested without DSFAST_HAVE_CBLAS");
#endif
}

/** Same-spin local contraction using one relative rotation. */
inline Complex relativeSpinProduct(Complex adj, Complex fwd, int spin,
                                    Complex cF, Complex cA) {
  const Complex ratio = cF * std::conj(cA);
  Complex power{1.0, 0.0};
  for (int i = 0; i < std::abs(spin); ++i) power *= ratio;
  if (spin < 0) power = std::conj(power);
  return power * std::conj(adj) * fwd;
}

/** Exact frequency-first coefficient pair C = sum_k conj(a_k) f_k. */
inline Eigen::MatrixXcd pairCovariance(const Eigen::MatrixXcd &a,
                                       const Eigen::MatrixXcd &f) {
  if (a.rows() != f.rows()) throw std::invalid_argument("frequency dimensions differ");
  return a.adjoint() * f;
}

inline Eigen::MatrixXcd evaluatePair(const Eigen::MatrixXd &hA,
                                     const Eigen::MatrixXd &hF,
                                     const Eigen::MatrixXcd &cov) {
  if (hA.cols() != cov.rows() || hF.cols() != cov.cols() || hA.rows() != hF.rows())
    throw std::invalid_argument("pair evaluation dimensions differ");
  Eigen::MatrixXcd out(hA.rows(), 1);
  for (Eigen::Index p = 0; p < hA.rows(); ++p)
    out(p, 0) = hA.row(p).cast<Complex>().dot(cov * hF.row(p).transpose().cast<Complex>());
  return out;
}

struct LowRankPair {
  Eigen::MatrixXcd left, right;
  double relativeFrobeniusResidual = 0;
};
inline LowRankPair truncatePair(const Eigen::MatrixXcd &c, Eigen::Index rank) {
  Eigen::JacobiSVD<Eigen::MatrixXcd> svd(c, Eigen::ComputeThinU | Eigen::ComputeThinV);
  rank = std::clamp<Eigen::Index>(rank, 0, svd.singularValues().size());
  Eigen::MatrixXcd left = svd.matrixU().leftCols(rank) *
      svd.singularValues().head(rank).cast<Complex>().asDiagonal();
  Eigen::MatrixXcd right = svd.matrixV().leftCols(rank).adjoint();
  const double total = c.norm();
  const double residual = rank == svd.singularValues().size() ? 0.0 :
      std::sqrt(svd.singularValues().tail(svd.singularValues().size() - rank).squaredNorm());
  return {std::move(left), std::move(right), total == 0 ? 0 : residual / total};
}

/** Element-local differentiation and exact r*strain polynomial interpolation. */
inline Eigen::MatrixXcd differentiateNodes(const Eigen::MatrixXd &Dref,
                                           const Eigen::MatrixXcd &values,
                                           double elementJacobian) {
  if (Dref.rows() != values.rows() || Dref.cols() != values.rows() || elementJacobian == 0)
    throw std::invalid_argument("invalid nodal derivative inputs");
  return (Dref / elementJacobian).cast<Complex>() * values;
}
inline Eigen::MatrixXcd interpolateRScaledPolynomial(const Eigen::MatrixXd &shape,
                                                      const Eigen::MatrixXcd &nodalRStrain) {
  if (shape.cols() != nodalRStrain.rows()) throw std::invalid_argument("bridge dimensions differ");
  return shape.cast<Complex>() * nodalRStrain;
}

/** Sparse symbolic pattern identity; numeric values may differ, structure may not. */
template<class SparseMatrix>
inline bool samePattern(const SparseMatrix &a, const SparseMatrix &b) {
  if (a.rows() != b.rows() || a.cols() != b.cols() || a.nonZeros() != b.nonZeros()) return false;
  for (Eigen::Index r = 0; r < a.outerSize(); ++r) {
    auto ia = a.outerIndexPtr()[r], ib = b.outerIndexPtr()[r];
    if (ia != ib) return false;
  }
  for (Eigen::Index k = 0; k < a.nonZeros(); ++k)
    if (a.innerIndexPtr()[k] != b.innerIndexPtr()[k]) return false;
  return true;
}

/** One existing equilibrated factor serves forward and H^{-H} bases. */
template<class SparseMatrix>
class SharedSymmetricFactor {
public:
  void factorize(const SparseMatrix &h) {
    SparseMatrix difference=h-SparseMatrix(h.transpose());
    if(difference.norm()>1e-13*std::max(std::numeric_limits<double>::min(),h.norm()))
      throw std::invalid_argument("complex-symmetric shortcut requires symmetric H");
    solver_.compute(h);
    if (solver_.info() != Eigen::Success) throw std::runtime_error("sparse factorization failed");
  }
  Eigen::MatrixXcd forward(const Eigen::MatrixXcd &rhs) const {
    return solver_.solve(rhs);
  }
  Eigen::MatrixXcd adjoint(const Eigen::MatrixXcd &rhs) const {
    return solver_.solveAdjoint(rhs);
  }
private:
  mutable SPARSESPEC::EquilibratedSparseLU<typename SparseMatrix::Scalar> solver_;
};

/** Select a rank-minimal basis from receiver rows, retaining exact row map. */
inline std::vector<Eigen::Index> independentRows(const Eigen::MatrixXcd &a,
                                                  double relativeTolerance = 1e-13) {
  std::vector<Eigen::Index> rows;
  Eigen::MatrixXcd basis(a.cols(), 0);
  for (Eigen::Index r = 0; r < a.rows(); ++r) {
    Eigen::MatrixXcd candidate(a.cols(), basis.cols() + 1);
    if (basis.cols()) candidate.leftCols(basis.cols()) = basis;
    candidate.col(basis.cols()) = a.row(r).transpose();
    Eigen::ColPivHouseholderQR<Eigen::MatrixXcd> qr(candidate);
    qr.setThreshold(relativeTolerance);
    if (qr.rank() > basis.cols()) {
      rows.push_back(r);
      basis = std::move(candidate);
    }
  }
  return rows;
}

struct ReceiverRowBasis {
  Eigen::MatrixXcd rightHandSides; // dof x independent-row-count
  Eigen::MatrixXcd reconstruction; // original-row-count x independent-row-count
};

/** Build a rank-minimal exact basis and its row reconstruction map. */
inline ReceiverRowBasis exactReceiverRowBasis(const Eigen::MatrixXcd &rows,
                                               double relativeTolerance=1e-13) {
  mathprofile::Scope t("receiver_basis_and_QR_reconstruction");
  mathprofile::count("receiver_basis_input_rows",static_cast<std::uint64_t>(rows.rows()));
  mathprofile::count("receiver_basis_input_columns",static_cast<std::uint64_t>(rows.cols()));
  const auto selected=independentRows(rows,relativeTolerance);
  mathprofile::count("receiver_basis_selected_columns",selected.size());
  ReceiverRowBasis result;
  result.rightHandSides.resize(rows.cols(),static_cast<Eigen::Index>(selected.size()));
  for(std::size_t j=0;j<selected.size();++j)
    result.rightHandSides.col(static_cast<Eigen::Index>(j))=rows.row(selected[j]).adjoint();
  result.reconstruction.resize(rows.rows(),static_cast<Eigen::Index>(selected.size()));
  Eigen::ColPivHouseholderQR<Eigen::MatrixXcd> qr(result.rightHandSides);
  qr.setThreshold(relativeTolerance);
  for(Eigen::Index i=0;i<rows.rows();++i){
    const Eigen::VectorXcd coeff=qr.solve(rows.row(i).adjoint());
    if((result.rightHandSides*coeff-rows.row(i).adjoint()).norm()>
       relativeTolerance*std::max(1.0,rows.row(i).norm()))
      throw std::runtime_error("Receiver basis failed exact row reconstruction.");
    result.reconstruction.row(i)=coeff.adjoint();
  }
  return result;
}

} // namespace dsfast
