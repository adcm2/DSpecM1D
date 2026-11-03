#ifndef DSPECM1D_EQUILIBRATED_SPARSE_LU_H
#define DSPECM1D_EQUILIBRATED_SPARSE_LU_H

#include "profile.hpp"

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <stdexcept>
#include <vector>

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
    mathprofile::Scope total("sparse_factorization_total");
    computed_ = false;
    mathprofile::count("sparse_factorization_calls");
    mathprofile::count("sparse_matrix_rows",static_cast<std::uint64_t>(matrix.rows()));
    mathprofile::count("sparse_matrix_nonzeros",static_cast<std::uint64_t>(matrix.nonZeros()));
    SparseMatrix equilibrated;
    {mathprofile::Scope t("sparse_equilibration_and_scaling");equilibrated = equilibrate(matrix);}
    // Eigen 3.4 SparseLU::compute is exactly analyzePattern(matrix), then
    // factorize(matrix), in that order. Keep those calls and their ordering.
    mathprofile::count("sparse_symbolic_analysis_calls");
    {mathprofile::Scope t("sparse_symbolic_analysis");solver_.analyzePattern(equilibrated);}
    savePattern(equilibrated);
    {mathprofile::Scope t("sparse_numerical_factorization");solver_.factorize(equilibrated);}
    info_ = solver_.info();
    computed_ = info_ == Eigen::Success;
  }

  void factorize(const SparseMatrix &matrix) {
    mathprofile::Scope total("sparse_factorization_total");
    const bool wasComputed = computed_;
    mathprofile::count("sparse_factorization_calls");
    mathprofile::count("sparse_matrix_rows",static_cast<std::uint64_t>(matrix.rows()));
    mathprofile::count("sparse_matrix_nonzeros",static_cast<std::uint64_t>(matrix.nonZeros()));
    SparseMatrix equilibrated;
    {mathprofile::Scope t("sparse_equilibration_and_scaling");equilibrated = equilibrate(matrix);}
    const bool reusePattern = wasComputed && samePattern(equilibrated);
    if (reusePattern) {
      mathprofile::count("sparse_symbolic_analysis_reuses");
      mathprofile::Scope t("sparse_numerical_factorization");solver_.factorize(equilibrated);
    } else {
      mathprofile::count("sparse_symbolic_analysis_calls");
      {mathprofile::Scope t("sparse_symbolic_analysis");solver_.analyzePattern(equilibrated);}
      savePattern(equilibrated);
      {mathprofile::Scope t("sparse_numerical_factorization");solver_.factorize(equilibrated);}
    }
    info_ = solver_.info();
    computed_ = info_ == Eigen::Success;
  }

  DenseMatrix solve(const DenseMatrix &rhs) {
    mathprofile::Scope total("sparse_forward_solve_including_rhs_scaling");
    mathprofile::count("forward_rhs_columns",static_cast<std::uint64_t>(rhs.cols()));
    {mathprofile::Scope t("sparse_forward_rhs_validation");requireSolveInput(rhs);}
    mathprofile::Scope rhsScaling("sparse_forward_rhs_scaling");
    DenseMatrix scaledRhs = rhs;
    for (Eigen::Index row = 0; row < rhs.rows(); ++row)
      scaledRhs.row(row) *= scale_[row];
    rhsScaling.stop();
    mathprofile::Scope eigenSolve("sparse_forward_eigen_solve");
    DenseMatrix solution = solver_.solve(scaledRhs);
    eigenSolve.stop();
    mathprofile::Scope resultValidation("sparse_forward_result_validation");
    info_ = solver_.info();
    if (info_ != Eigen::Success || !solution.allFinite())
      throw std::runtime_error(
          "EquilibratedSparseLU forward solve failed or returned non-finite values.");
    resultValidation.stop();
    mathprofile::Scope resultScaling("sparse_forward_result_unscaling_and_validation");
    for (Eigen::Index row = 0; row < solution.rows(); ++row)
      solution.row(row) *= scale_[row];
    if (!solution.allFinite())
      throw std::runtime_error(
          "EquilibratedSparseLU forward solution is not finite in original units.");
    return solution;
  }

  DenseMatrix solveAdjoint(const DenseMatrix &rhs) {
    mathprofile::Scope total("sparse_adjoint_solve_including_rhs_scaling");
    mathprofile::count("adjoint_rhs_columns",static_cast<std::uint64_t>(rhs.cols()));
    {mathprofile::Scope t("sparse_adjoint_rhs_validation");requireSolveInput(rhs);}
    mathprofile::Scope rhsScaling("sparse_adjoint_rhs_scaling");
    DenseMatrix scaledRhs = rhs;
    for (Eigen::Index row = 0; row < rhs.rows(); ++row)
      scaledRhs.row(row) *= scale_[row];
    rhsScaling.stop();
    mathprofile::Scope eigenSolve("sparse_adjoint_eigen_solve");
    DenseMatrix solution = solver_.adjoint().solve(scaledRhs);
    eigenSolve.stop();
    mathprofile::Scope resultValidation("sparse_adjoint_result_validation");
    info_ = solver_.info();
    if (info_ != Eigen::Success || !solution.allFinite())
      throw std::runtime_error(
          "EquilibratedSparseLU adjoint solve failed or returned non-finite values.");
    resultValidation.stop();
    mathprofile::Scope resultScaling("sparse_adjoint_result_unscaling_and_validation");
    for (Eigen::Index row = 0; row < solution.rows(); ++row)
      solution.row(row) *= scale_[row];
    if (!solution.allFinite())
      throw std::runtime_error(
          "EquilibratedSparseLU adjoint solution is not finite in original units.");
    return solution;
  }

  Eigen::ComputationInfo info() const { return info_; }

private:
  void savePattern(const SparseMatrix &matrix) {
    patternRows_ = matrix.rows();
    patternCols_ = matrix.cols();
    patternOuter_.assign(matrix.outerIndexPtr(),
                         matrix.outerIndexPtr() + matrix.outerSize() + 1);
    patternInner_.assign(matrix.innerIndexPtr(),
                         matrix.innerIndexPtr() + matrix.nonZeros());
  }

  bool samePattern(const SparseMatrix &matrix) const {
    if (matrix.rows() != patternRows_ || matrix.cols() != patternCols_ ||
        static_cast<std::size_t>(matrix.nonZeros()) != patternInner_.size() ||
        static_cast<std::size_t>(matrix.outerSize() + 1) != patternOuter_.size())
      return false;
    return std::equal(patternOuter_.begin(), patternOuter_.end(),
                      matrix.outerIndexPtr()) &&
           std::equal(patternInner_.begin(), patternInner_.end(),
                      matrix.innerIndexPtr());
  }

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
  Eigen::Index patternRows_ = 0, patternCols_ = 0;
  std::vector<int> patternOuter_, patternInner_;
  Eigen::ComputationInfo info_ = Eigen::InvalidInput;
  bool computed_ = false;
};

} // namespace SPARSESPEC

#endif // DSPECM1D_EQUILIBRATED_SPARSE_LU_H
