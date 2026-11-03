#pragma once

// Diagnostic-only exact pair and scalar-harmonic operations. This file does
// not alter or replace the accepted DSpecM1D implementation.
#include <Eigen/Core>
#include <Eigen/SVD>

#include <algorithm>
#include <complex>
#include <stdexcept>
#include <vector>

namespace dspectral {
using Complex = std::complex<double>;

// Input rows are frequencies; columns are compact angular coefficients.
// Keeping this product (rather than multiplying summed fields) preserves all
// same-frequency cross-degree/order terms and excludes k != k' products.
inline Eigen::MatrixXcd frequencyCovariance(const Eigen::MatrixXcd &adjoint,
                                           const Eigen::MatrixXcd &forward) {
  if (adjoint.rows() != forward.rows())
    throw std::invalid_argument("frequency dimensions differ");
  return adjoint.adjoint() * forward;
}

// `hA` and `hF` hold the frame-converted angular basis rows for each point.
// For row vectors the scalar is conj(hA[p])^T C hF[p].
inline Eigen::VectorXcd evaluateCovariance(const Eigen::MatrixXcd &hA,
                                          const Eigen::MatrixXcd &hF,
                                          const Eigen::MatrixXcd &covariance) {
  if (hA.rows() != hF.rows() || hA.cols() != covariance.rows() ||
      hF.cols() != covariance.cols())
    throw std::invalid_argument("angular covariance dimensions differ");
  Eigen::VectorXcd out(hA.rows());
  for (Eigen::Index p = 0; p < hA.rows(); ++p)
    out[p] = hA.row(p).conjugate() * covariance * hF.row(p).transpose();
  return out;
}

struct PairSvd {
  Eigen::MatrixXcd left;
  Eigen::MatrixXcd right;
  Eigen::VectorXd singular;
  double factorBytes = 0;
};

inline PairSvd factorCovariance(const Eigen::MatrixXcd &covariance) {
  Eigen::JacobiSVD<Eigen::MatrixXcd> svd(
      covariance, Eigen::ComputeThinU | Eigen::ComputeThinV);
  PairSvd out;
  out.left = std::move(svd.matrixU());
  out.right = std::move(svd.matrixV());
  out.singular = std::move(svd.singularValues());
  out.factorBytes = static_cast<double>(
      (out.left.size() + out.right.size()) * sizeof(Complex) +
      out.singular.size() * sizeof(double));
  return out;
}

inline Eigen::VectorXcd evaluateSvd(const Eigen::MatrixXcd &hA,
                                    const Eigen::MatrixXcd &hF,
                                    const PairSvd &svd,
                                    Eigen::Index rank) {
  if (rank < 0 || rank > svd.singular.size())
    throw std::invalid_argument("SVD rank out of range");
  Eigen::MatrixXcd left = hA.conjugate() * svd.left.leftCols(rank);
  Eigen::MatrixXcd right = hF * svd.right.leftCols(rank).conjugate();
  Eigen::VectorXcd out = Eigen::VectorXcd::Zero(hA.rows());
  for (Eigen::Index k = 0; k < rank; ++k)
    out.array() += svd.singular[k] * left.col(k).array() * right.col(k).array();
  // `left` is conj(hA) U and `right` is hF conj(V), which equals V^H hF
  // after interpreting each point's angular row as a column in the product.
  return out;
}

struct VerifiedRank {
  Eigen::Index rank = 0;
  double relativeL2 = 0;
  double scaledPointwise = 0;
  bool accepted = false;
};

// Scan every retained rank against the full covariance output. The rank is
// selected by the agreed output tolerances, never by a singular-value cutoff.
inline VerifiedRank verifyRank(const Eigen::VectorXcd &reference,
                               const Eigen::VectorXcd &candidate) {
  if (reference.size() != candidate.size())
    throw std::invalid_argument("rank comparison dimensions differ");
  VerifiedRank out;
  const Eigen::VectorXcd error = candidate - reference;
  out.relativeL2 = reference.norm() == 0 ? error.norm() : error.norm() / reference.norm();
  const double peak = reference.cwiseAbs().maxCoeff();
  const double floor = 0.01 * peak;
  out.scaledPointwise = 0;
  if (floor == 0.0) out.scaledPointwise = error.cwiseAbs().maxCoeff();
  else for (Eigen::Index i = 0; i < reference.size(); ++i)
    out.scaledPointwise = std::max(out.scaledPointwise,
        std::abs(error[i]) / std::max(std::abs(reference[i]), floor));
  out.accepted = out.relativeL2 <= 1e-12 && out.scaledPointwise <= 1e-10;
  return out;
}

inline std::pair<Eigen::Index, VerifiedRank> smallestVerifiedRank(
    const Eigen::MatrixXcd &hA, const Eigen::MatrixXcd &hF,
    const PairSvd &svd, const Eigen::VectorXcd &fullReference) {
  // Build each rank's contribution once and extend the candidate in place.
  // Re-evaluating all factors from scratch at every rank adds an unnecessary
  // O(P r^2) scan and overstates the cost of verified numerical rank.
  const Eigen::MatrixXcd left = hA.conjugate() * svd.left;
  const Eigen::MatrixXcd right = hF * svd.right.conjugate();
  Eigen::VectorXcd candidate = Eigen::VectorXcd::Zero(hA.rows());
  for (Eigen::Index rank = 0; rank <= svd.singular.size(); ++rank) {
    auto result = verifyRank(fullReference, candidate);
    if (result.accepted) {
      result.rank = rank;
      return {rank, result};
    }
    if (rank < svd.singular.size())
      candidate.array() += svd.singular[rank] * left.col(rank).array() * right.col(rank).array();
  }
  throw std::runtime_error("full SVD rank failed its own output tolerance");
}

} // namespace dspectral
