#ifndef DSPECM1D_EQUILIBRATED_SPARSE_LU_H
#define DSPECM1D_EQUILIBRATED_SPARSE_LU_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>

namespace SPARSESPEC {

/**
 * @brief SparseLU with symmetric positive diagonal equilibration.
 *
 * The factorized operator is D A D, where D_i is the inverse square root of
 * the largest magnitude in row i of A. Forward and Hermitian-adjoint solves
 * map their results back to the units of the original operator.
 */
template <typename Scalar> class EquilibratedSparseLU {
public:
  using SparseMatrix = Eigen::SparseMatrix<Scalar>;
  using DenseMatrix = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
  using Real = typename Eigen::NumTraits<Scalar>::Real;
  using RealVector = Eigen::Matrix<Real, Eigen::Dynamic, 1>;
  using Solver = Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<int>>;

  void compute(const SparseMatrix &matrix) {
    computed_ = false;
    SparseMatrix equilibrated = equilibrate(matrix);
    solver_.compute(equilibrated);
    info_ = solver_.info();
    computed_ = info_ == Eigen::Success;
  }

  void factorize(const SparseMatrix &matrix) {
    const bool reusePattern = computed_;
    SparseMatrix equilibrated = equilibrate(matrix);
    if (reusePattern)
      solver_.factorize(equilibrated);
    else
      solver_.compute(equilibrated);
    info_ = solver_.info();
    computed_ = info_ == Eigen::Success;
  }

  DenseMatrix solve(const DenseMatrix &rhs) {
    requireSolveInput(rhs);
    DenseMatrix scaledRhs = rhs;
    for (Eigen::Index row = 0; row < rhs.rows(); ++row)
      scaledRhs.row(row) *= scale_[row];
    DenseMatrix solution = solver_.solve(scaledRhs);
    info_ = solver_.info();
    if (info_ != Eigen::Success || !solution.allFinite())
      throw std::runtime_error(
          "EquilibratedSparseLU forward solve failed or returned non-finite values.");
    for (Eigen::Index row = 0; row < solution.rows(); ++row)
      solution.row(row) *= scale_[row];
    if (!solution.allFinite())
      throw std::runtime_error(
          "EquilibratedSparseLU forward solution is not finite in original units.");
    return solution;
  }

  DenseMatrix solveAdjoint(const DenseMatrix &rhs) {
    requireSolveInput(rhs);
    DenseMatrix scaledRhs = rhs;
    for (Eigen::Index row = 0; row < rhs.rows(); ++row)
      scaledRhs.row(row) *= scale_[row];
    DenseMatrix solution = solver_.adjoint().solve(scaledRhs);
    info_ = solver_.info();
    if (info_ != Eigen::Success || !solution.allFinite())
      throw std::runtime_error(
          "EquilibratedSparseLU adjoint solve failed or returned non-finite values.");
    for (Eigen::Index row = 0; row < solution.rows(); ++row)
      solution.row(row) *= scale_[row];
    if (!solution.allFinite())
      throw std::runtime_error(
          "EquilibratedSparseLU adjoint solution is not finite in original units.");
    return solution;
  }

  Eigen::ComputationInfo info() const { return info_; }

private:
  SparseMatrix equilibrate(const SparseMatrix &matrix) {
    computed_ = false;
    info_ = Eigen::InvalidInput;
    if (matrix.rows() <= 0 || matrix.rows() != matrix.cols())
      throw std::invalid_argument(
          "EquilibratedSparseLU requires a nonempty square matrix.");
    scale_.resize(matrix.rows());
    scale_.setZero();
    for (int column = 0; column < matrix.outerSize(); ++column)
      for (typename SparseMatrix::InnerIterator entry(matrix, column); entry; ++entry) {
        const Real magnitude = std::abs(entry.value());
        if (!std::isfinite(magnitude))
          throw std::invalid_argument(
              "EquilibratedSparseLU matrix must contain only finite values.");
        scale_[entry.row()] = std::max(scale_[entry.row()], magnitude);
      }
    for (Eigen::Index row = 0; row < scale_.size(); ++row) {
      if (!(scale_[row] > Real(0)))
        throw std::invalid_argument(
            "EquilibratedSparseLU matrix must not contain zero rows.");
      scale_[row] = Real(1) / std::sqrt(scale_[row]);
      if (!std::isfinite(scale_[row]))
        throw std::invalid_argument(
            "EquilibratedSparseLU row scaling is not finite.");
    }
    SparseMatrix result = matrix;
    for (int column = 0; column < result.outerSize(); ++column)
      for (typename SparseMatrix::InnerIterator entry(result, column); entry; ++entry) {
        entry.valueRef() *= scale_[entry.row()] * scale_[entry.col()];
        if (!std::isfinite(std::abs(entry.value())))
          throw std::invalid_argument(
              "EquilibratedSparseLU scaled matrix is not finite.");
      }
    return result;
  }

  void requireSolveInput(const DenseMatrix &rhs) const {
    if (!computed_)
      throw std::logic_error(
          "EquilibratedSparseLU must be successfully factorized before solving.");
    if (rhs.rows() != scale_.size() || rhs.cols() == 0 || !rhs.allFinite())
      throw std::invalid_argument(
          "EquilibratedSparseLU right-hand side dimensions or values are invalid.");
  }

  Solver solver_;
  RealVector scale_;
  Eigen::ComputationInfo info_ = Eigen::InvalidInput;
  bool computed_ = false;
};

} // namespace SPARSESPEC

#endif // DSPECM1D_EQUILIBRATED_SPARSE_LU_H
