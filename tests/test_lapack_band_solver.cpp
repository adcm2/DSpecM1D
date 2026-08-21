#include <gtest/gtest.h>

#include <Eigen/SparseLU>

#include <DSpecM1D/src/LapackBandSolver.h>

namespace {

using Complex = std::complex<double>;
using SparseMatrix = Eigen::SparseMatrix<Complex>;
using DenseMatrix = Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic>;

SparseMatrix makeSystem() {
  SparseMatrix matrix(3, 3);
  matrix.insert(0, 0) = Complex(4.0, 1.0);
  matrix.insert(1, 0) = Complex(2.0, -1.0);
  matrix.insert(0, 1) = Complex(1.0, 0.5);
  matrix.insert(1, 1) = Complex(5.0, 0.0);
  matrix.insert(2, 1) = Complex(3.0, 0.5);
  matrix.insert(1, 2) = Complex(1.0, -0.5);
  matrix.insert(2, 2) = Complex(6.0, -1.0);
  return matrix;
}

double matrixInfinityNorm(const SparseMatrix &matrix) {
  Eigen::VectorXd rowSums = Eigen::VectorXd::Zero(matrix.rows());
  for (Eigen::Index column = 0; column < matrix.outerSize(); ++column)
    for (SparseMatrix::InnerIterator entry(matrix, column); entry; ++entry)
      rowSums(entry.row()) += std::abs(entry.value());
  return rowSums.maxCoeff();
}

double denseInfinityNorm(const DenseMatrix &matrix) {
  Eigen::VectorXd rowSums = Eigen::VectorXd::Zero(matrix.rows());
  for (Eigen::Index row = 0; row < matrix.rows(); ++row)
    for (Eigen::Index column = 0; column < matrix.cols(); ++column)
      rowSums(row) += std::abs(matrix(row, column));
  return rowSums.maxCoeff();
}

double relativeResidual(const SparseMatrix &matrix, const DenseMatrix &x,
                        const DenseMatrix &rhs) {
  const DenseMatrix residual = matrix * x - rhs;
  return denseInfinityNorm(residual) /
         (matrixInfinityNorm(matrix) * denseInfinityNorm(x) +
          denseInfinityNorm(rhs));
}

} // namespace

TEST(LapackBandSolverTests, SolvesSingleAndMultipleRhs) {
  const SparseMatrix matrix = makeSystem();
  DenseMatrix exact(3, 2);
  exact << Complex(1.0, -1.0), Complex(2.0, 0.5), Complex(-2.0, 0.25),
      Complex(1.0, 1.5), Complex(0.5, -0.75), Complex(-1.0, 2.0);
  const DenseMatrix rhs = matrix * exact;

  SPARSESPEC::detail::LapackBandSolver solver;
  solver.compute(matrix);
  ASSERT_EQ(solver.info(), 0);
  const DenseMatrix singleResult = solver.solve(rhs.leftCols(1));
  ASSERT_EQ(solver.info(), 0);
  EXPECT_LT(relativeResidual(matrix, singleResult, rhs.leftCols(1)), 1e-12);

  const DenseMatrix result = solver.solve(rhs);
  ASSERT_EQ(solver.info(), 0);
  EXPECT_LT(relativeResidual(matrix, result, rhs), 1e-12);

  Eigen::SparseLU<SparseMatrix> eigenSolver;
  eigenSolver.compute(matrix);
  ASSERT_EQ(eigenSolver.info(), Eigen::Success);
  const DenseMatrix eigenResult = eigenSolver.solve(rhs);
  EXPECT_LT((result - eigenResult).norm() / eigenResult.norm(), 1e-10);
}

TEST(LapackBandSolverTests, SingularMatrixReportsFactorizationFailure) {
  SparseMatrix matrix(2, 2);
  matrix.insert(0, 0) = Complex(1.0, 0.0);
  matrix.insert(1, 0) = Complex(2.0, 0.0);
  matrix.insert(0, 1) = Complex(2.0, 0.0);
  matrix.insert(1, 1) = Complex(4.0, 0.0);

  SPARSESPEC::detail::LapackBandSolver solver;
  solver.compute(matrix);
  EXPECT_NE(solver.info(), 0);
}

TEST(LapackBandSolverTests, FactorizeRejectsWiderEstablishedBand) {
  const SparseMatrix matrix = makeSystem();
  SparseMatrix wider = matrix;
  wider.insert(0, 2) = Complex(0.25, -0.5);
  wider.makeCompressed();

  SPARSESPEC::detail::LapackBandSolver solver;
  solver.compute(matrix);
  ASSERT_EQ(solver.info(), 0);
  solver.factorize(wider);
  EXPECT_EQ(solver.info(),
            SPARSESPEC::detail::LapackBandSolver::kBandStructureMismatch);
}
