#include <algorithm>
#include <complex>

#include <gtest/gtest.h>
#include <Eigen/SparseLU>

#include <DSpecM1D/ModelInput>
#include <DSpecM1D/src/NormClass.h>
#include <DSpecM1D/src/SEM/SEM.h>
#include <DSpecM1D/src/SpecHelpers.h>
#include <DSpecM1D/src/LapackBandSolver.h>

#include "test_utils.h"

namespace {

using Complex = std::complex<double>;
using SparseMatrix = Eigen::SparseMatrix<Complex>;
using DenseMatrix = Eigen::MatrixXcd;

SparseMatrix frequencyMatrix(const Eigen::SparseMatrix<double> &h,
                             const Eigen::SparseMatrix<double> &p,
                             const Eigen::SparseMatrix<double> &hAtten,
                             double omega, double epsilon, double tref) {
  SPARSESPEC::SpecConstants constants(epsilon, tref);
  const Complex omegaHat = omega + constants.ieps;
  const Complex chi = SPARSESPEC::attenFactor(
      omega, constants.w0, constants.twodivpi, constants.myi);
  SparseMatrix result = h.cast<Complex>() - omegaHat * omegaHat * p.cast<Complex>();
  result += chi * hAtten.cast<Complex>();
  result.makeCompressed();
  return result;
}

DenseMatrix makeRhs(Eigen::Index rows, Eigen::Index columns) {
  DenseMatrix rhs(rows, columns);
  for (Eigen::Index column = 0; column < columns; ++column) {
    for (Eigen::Index row = 0; row < rows; ++row) {
      rhs(row, column) = Complex(0.7 + 0.001 * row + 0.13 * column,
                                 -0.2 + 0.002 * row - 0.07 * column);
    }
  }
  return rhs;
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

void compareDirectBand(const Eigen::SparseMatrix<double> &h,
                       const Eigen::SparseMatrix<double> &p,
                       const Eigen::SparseMatrix<double> &hAtten,
                       double omega, double epsilon, double tref,
                       Eigen::Index ridx, int rhsColumns) {
  const auto hComplex = h.cast<Complex>();
  const auto pComplex = p.cast<Complex>();
  const auto aComplex = hAtten.cast<Complex>();
  const auto hBandwidth = SPARSESPEC::detail::lapackBandBandwidth(hComplex);
  const auto pBandwidth = SPARSESPEC::detail::lapackBandBandwidth(pComplex);
  const auto aBandwidth = SPARSESPEC::detail::lapackBandBandwidth(aComplex);
  const Eigen::Index kl = std::max(
      {hBandwidth.first, pBandwidth.first, aBandwidth.first});
  const Eigen::Index ku = std::max(
      {hBandwidth.second, pBandwidth.second, aBandwidth.second});

  SPARSESPEC::detail::LapackBandMatrix hBand, pBand, aBand;
  SPARSESPEC::detail::packLapackBandInto(hComplex, hBand, kl, ku);
  SPARSESPEC::detail::packLapackBandInto(pComplex, pBand, kl, ku);
  SPARSESPEC::detail::packLapackBandInto(aComplex, aBand, kl, ku);

  SPARSESPEC::SpecConstants constants(epsilon, tref);
  const Complex omegaHat = omega + constants.ieps;
  const Complex chi = SPARSESPEC::attenFactor(
      omega, constants.w0, constants.twodivpi, constants.myi);
  const SparseMatrix matrix = frequencyMatrix(
      h, p, hAtten, omega, epsilon, tref);
  const Eigen::Index activeSize = matrix.rows() - ridx;

  SPARSESPEC::detail::LapackBandSolver solver;
  auto &workspace = solver.bandWorkspace(matrix.rows(), kl, ku);
  workspace.coefficients().middleCols(ridx, activeSize) =
      hBand.coefficients().middleCols(ridx, activeSize) -
      (omegaHat * omegaHat) * pBand.coefficients().middleCols(ridx, activeSize) +
      chi * aBand.coefficients().middleCols(ridx, activeSize);

  // Check the exact physical band entries before LAPACK overwrites them with
  // LU factors.  Entries outside the active principal subsystem are omitted.
  for (Eigen::Index column = ridx; column < matrix.cols(); ++column) {
    for (Eigen::Index row = ridx; row < matrix.rows(); ++row) {
      if (row - column > kl || column - row > ku)
        continue;
      const Complex expected = matrix.coeff(row, column);
      const Complex actual = workspace.at(row, column);
      EXPECT_NEAR(std::abs(actual - expected), 0.0, 1e-13)
          << "coefficient (" << row << "," << column << ")";
    }
  }

  solver.factorize(ridx);
  ASSERT_EQ(solver.info(), 0);
  const SparseMatrix active = matrix.block(ridx, ridx, activeSize, activeSize);
  const DenseMatrix rhs = makeRhs(activeSize, rhsColumns);
  const DenseMatrix lapackSolution = solver.solve(rhs);
  ASSERT_EQ(solver.info(), 0);

  Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<int>> reference;
  reference.compute(active);
  ASSERT_EQ(reference.info(), Eigen::Success);
  const DenseMatrix referenceSolution = reference.solve(rhs);
  EXPECT_LT(denseInfinityNorm(active * lapackSolution - rhs) /
                (matrixInfinityNorm(active) * denseInfinityNorm(lapackSolution) +
                 denseInfinityNorm(rhs)),
            1e-11);
  EXPECT_LT((lapackSolution - referenceSolution).norm() /
                referenceSolution.norm(),
            1e-9);
}

void compareSolvers(const SparseMatrix &matrix, int rhsColumns) {
  ASSERT_GT(matrix.rows(), 0);
  ASSERT_EQ(matrix.rows(), matrix.cols());

  const DenseMatrix rhs = makeRhs(matrix.rows(), rhsColumns);
  Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<int>> eigenSolver;
  eigenSolver.compute(matrix);
  ASSERT_EQ(eigenSolver.info(), Eigen::Success);
  const DenseMatrix eigenSolution = eigenSolver.solve(rhs);
  ASSERT_EQ(eigenSolver.info(), Eigen::Success);

  SPARSESPEC::detail::LapackBandSolver lapackSolver;
  lapackSolver.compute(matrix);
  ASSERT_EQ(lapackSolver.info(), 0);
  const DenseMatrix lapackSolution = lapackSolver.solve(rhs);
  ASSERT_EQ(lapackSolver.info(), 0);

  const double residual = denseInfinityNorm(matrix * lapackSolution - rhs) /
                          (matrixInfinityNorm(matrix) *
                               denseInfinityNorm(lapackSolution) +
                           denseInfinityNorm(rhs));
  const double agreement =
      (lapackSolution - eigenSolution).norm() / eigenSolution.norm();
  EXPECT_LT(residual, 1e-11);
  EXPECT_LT(agreement, 1e-9);

  lapackSolver.factorize(matrix);
  ASSERT_EQ(lapackSolver.info(), 0);
  const DenseMatrix refactoredSolution = lapackSolver.solve(rhs);
  ASSERT_EQ(lapackSolver.info(), 0);
  const double refactoredResidual =
      denseInfinityNorm(matrix * refactoredSolution - rhs) /
      (matrixInfinityNorm(matrix) * denseInfinityNorm(refactoredSolution) +
       denseInfinityNorm(rhs));
  EXPECT_LT(refactoredResidual, 1e-11);
}

} // namespace

TEST(LapackSEMComparisonTests, RealAttenuatedRadialToroidalAndSpheroidalSystems) {
  prem_norm<double> norm;
  auto model = EarthModels::ModelInput(
      (DSpecMTest::repoRoot() / "data/models/prem.200.no.txt").string(), norm);
  Full1D::SEM sem(model, 0.05, 5, 20);

  ASSERT_TRUE(sem.mesh().HasFluid());
  ASSERT_FALSE(sem.mesh().FS_Boundaries().empty());

  // Radial: full domain and one right-hand side at an interior frequency.
  const SparseMatrix radial = frequencyMatrix(
      sem.hR(), sem.pR(), sem.hRa(), 0.70, 0.01, model.TREF());
  ASSERT_EQ(radial.rows(), sem.hR().rows());
  compareDirectBand(sem.hR(), sem.pR(), sem.hRa(), 0.70, 0.01, model.TREF(),
                    0, 1);
  compareSolvers(radial, 1);

  // Toroidal: discard a genuine inner prefix of the mantle domain.
  const int toroidalDegree = 6;
  const auto toroidalFull = frequencyMatrix(
      sem.hTk(toroidalDegree), sem.pTk(toroidalDegree),
      sem.hTa(toroidalDegree), 0.85, 0.01, model.TREF());
  const int toroidalElement = sem.el() + 2;
  ASSERT_LT(toroidalElement, sem.eu());
  const Eigen::Index toroidalStart = sem.ltgT(toroidalElement, 0);
  ASSERT_GT(toroidalStart, 0);
  ASSERT_LT(toroidalStart, toroidalFull.rows());
  const auto toroidal = toroidalFull.block(
      toroidalStart, toroidalStart, toroidalFull.rows() - toroidalStart,
      toroidalFull.cols() - toroidalStart);
  compareDirectBand(sem.hTk(toroidalDegree), sem.pTk(toroidalDegree),
                    sem.hTa(toroidalDegree), 0.85, 0.01, model.TREF(),
                    toroidalStart, 2);
  compareSolvers(toroidal, 2);

  // Spheroidal: begin before a fluid-solid boundary so the retained domain
  // crosses that structure while still being a true truncation.
  const int fluidSolidBoundary = sem.mesh().FS_Boundaries().front();
  ASSERT_GT(fluidSolidBoundary, 1);
  const int spheroidalDegree = 8;
  const auto spheroidalFull = frequencyMatrix(
      sem.hS(spheroidalDegree), sem.pS(spheroidalDegree),
      sem.hSa(spheroidalDegree), 1.00, 0.01, model.TREF());
  const Eigen::Index spheroidalStart =
      sem.ltgS(0, fluidSolidBoundary - 1, 0);
  ASSERT_GT(spheroidalStart, 0);
  ASSERT_LT(spheroidalStart, spheroidalFull.rows());
  const auto spheroidal = spheroidalFull.block(
      spheroidalStart, spheroidalStart, spheroidalFull.rows() - spheroidalStart,
      spheroidalFull.cols() - spheroidalStart);
  compareDirectBand(sem.hS(spheroidalDegree), sem.pS(spheroidalDegree),
                    sem.hSa(spheroidalDegree), 1.00, 0.01, model.TREF(),
                    spheroidalStart, 4);
  compareSolvers(spheroidal, 4);
}
