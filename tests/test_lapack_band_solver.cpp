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

SparseMatrix makeTrailingSystem() {
  Eigen::MatrixXcd dense = Eigen::MatrixXcd::Zero(6, 6);
  for (Eigen::Index row = 0; row < dense.rows(); ++row) {
    dense(row, row) = Complex(8.0 + row, 0.25 * row);
    if (row > 0)
      dense(row, row - 1) = Complex(-1.0, 0.1);
    if (row > 1)
      dense(row, row - 2) = Complex(0.5, -0.2);
    if (row + 1 < dense.cols())
      dense(row, row + 1) = Complex(0.75, 0.3);
    if (row + 2 < dense.cols())
      dense(row, row + 2) = Complex(-0.25, 0.15);
  }
  return dense.sparseView();
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

TEST(LapackBandSolverTests, SolvesPackedFullAndTrailingSystems) {
  const SparseMatrix matrix = makeTrailingSystem();
  const auto band = SPARSESPEC::detail::packLapackBand(matrix);
  ASSERT_EQ(band.n, 6);
  ASSERT_EQ(band.kl, 2);
  ASSERT_EQ(band.ku, 2);
  ASSERT_EQ(band.ldab, 7);
  EXPECT_EQ(band.data.outerStride(), band.ldab);

  DenseMatrix exact(6, 2);
  exact.setRandom();
  const DenseMatrix rhs = matrix * exact;

  SPARSESPEC::detail::LapackBandSolver fullSolver;
  fullSolver.compute(band);
  ASSERT_EQ(fullSolver.info(), 0);
  const DenseMatrix fullResult = fullSolver.solve(rhs);
  ASSERT_EQ(fullSolver.info(), 0);
  EXPECT_LT(relativeResidual(matrix, fullResult, rhs), 1e-12);

  constexpr Eigen::Index ridx = 2;
  const SparseMatrix active = matrix.block(ridx, ridx, matrix.rows() - ridx,
                                           matrix.cols() - ridx);
  const DenseMatrix activeExact = exact.bottomRows(exact.rows() - ridx);
  const DenseMatrix activeRhs = active * activeExact;

  SPARSESPEC::detail::LapackBandSolver trailingSolver;
  trailingSolver.compute(band, ridx);
  ASSERT_EQ(trailingSolver.info(), 0);
  EXPECT_EQ(trailingSolver.activeOffset(), ridx);
  EXPECT_EQ(trailingSolver.activeSize(), matrix.rows() - ridx);
  const DenseMatrix trailingResult = trailingSolver.solve(activeRhs);
  ASSERT_EQ(trailingSolver.info(), 0);
  EXPECT_LT(relativeResidual(active, trailingResult, activeRhs), 1e-12);

  Eigen::SparseLU<SparseMatrix> reference;
  reference.compute(active);
  ASSERT_EQ(reference.info(), Eigen::Success);
  const DenseMatrix referenceResult = reference.solve(activeRhs);
  EXPECT_LT((trailingResult - referenceResult).norm() /
                referenceResult.norm(),
            1e-10);
}

TEST(LapackBandSolverTests, TrailingSolveIgnoresInactiveLeadingEntries) {
  const SparseMatrix matrix = makeTrailingSystem();
  auto band = SPARSESPEC::detail::packLapackBand(matrix);
  constexpr Eigen::Index ridx = 2;
  const SparseMatrix active = matrix.block(ridx, ridx, matrix.rows() - ridx,
                                           matrix.cols() - ridx);
  DenseMatrix exact(active.rows(), 1);
  exact.setRandom();
  const DenseMatrix rhs = active * exact;

  // These are entries from global rows before ridx in the first active
  // column.  They occupy the unused leading portion of the offset buffer.
  band.at(0, ridx) = Complex(1.0e100, -2.0e100);
  band.at(1, ridx) = Complex(-3.0e100, 4.0e100);

  SPARSESPEC::detail::LapackBandSolver solver;
  solver.compute(band, ridx);
  ASSERT_EQ(solver.info(), 0);
  const DenseMatrix result = solver.solve(rhs);
  ASSERT_EQ(solver.info(), 0);
  EXPECT_LT(relativeResidual(active, result, rhs), 1e-12);
}

TEST(LapackBandSolverTests, DirectWorkspaceFactorizesInPlaceAtTrailingOffset) {
  const SparseMatrix matrix = makeTrailingSystem();
  const auto packed = SPARSESPEC::detail::packLapackBand(matrix);

  SPARSESPEC::detail::LapackBandSolver solver;
  auto &workspace = solver.bandWorkspace(packed.n, packed.kl, packed.ku);
  auto *const dataAddress = workspace.data.data();
  const auto dataRows = workspace.data.rows();
  const auto dataColumns = workspace.data.cols();
  const auto dataLdab = workspace.ldab;

  auto solveActive = [&](Eigen::Index ridx, bool injectDiscarded = false) {
    const Eigen::Index activeSize = matrix.rows() - ridx;
    workspace.coefficients().middleCols(ridx, activeSize) =
        packed.coefficients().middleCols(ridx, activeSize);
    if (injectDiscarded) {
      workspace.at(1, ridx) = Complex(1.0e100, -2.0e100);
      workspace.at(2, ridx) = Complex(-3.0e100, 4.0e100);
    }
    const SparseMatrix active =
        matrix.block(ridx, ridx, activeSize, activeSize);
    DenseMatrix exact(activeSize, 2);
    exact.setRandom();
    const DenseMatrix rhs = active * exact;

    solver.factorize(ridx);
    ASSERT_EQ(solver.info(), 0);
    EXPECT_EQ(workspace.data.data(), dataAddress);
    EXPECT_EQ(workspace.data.rows(), dataRows);
    EXPECT_EQ(workspace.data.cols(), dataColumns);
    EXPECT_EQ(workspace.ldab, dataLdab);
    const DenseMatrix result = solver.solve(rhs);
    ASSERT_EQ(solver.info(), 0);
    EXPECT_LT(relativeResidual(active, result, rhs), 1e-12);

    Eigen::SparseLU<SparseMatrix> reference;
    reference.compute(active);
    ASSERT_EQ(reference.info(), Eigen::Success);
    const DenseMatrix referenceResult = reference.solve(rhs);
    EXPECT_LT((result - referenceResult).norm() / referenceResult.norm(),
              1e-10);
  };

  // The first run leaves factorization fill data throughout the workspace and
  // proves that discarded coefficients in its leading boundary are ignored.
  solveActive(3, true);

  // Grow the active system and refill only its physical coefficient rows,
  // including the newly active columns.  The top kl fill rows intentionally
  // retain stale LU values from the preceding factorization.
  constexpr Eigen::Index secondRidx = 1;
  solveActive(secondRidx);
}
