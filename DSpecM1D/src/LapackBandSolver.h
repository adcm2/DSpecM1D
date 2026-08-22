#ifndef DSPECM1D_LAPACK_BAND_SOLVER_H
#define DSPECM1D_LAPACK_BAND_SOLVER_H

#include <algorithm>
#include <complex>
#include <type_traits>
#include <vector>

#include <Eigen/Dense>
#include <lapacke.h>

#include "LapackBandStorage.h"

namespace SPARSESPEC::detail {

/// Small internal LU wrapper for a complex general-band matrix.
class LapackBandSolver {
public:
  using Complex = std::complex<double>;
  using SparseMatrix = Eigen::SparseMatrix<Complex>;
  using DenseMatrix = Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic,
                                    Eigen::ColMajor>;

  // This is deliberately outside LAPACK's ordinary argument-error range.
  static constexpr lapack_int kBandStructureMismatch = -1001;
  static constexpr lapack_int kInvalidState = -1002;

  LapackBandSolver &compute(const SparseMatrix &matrix) {
    try {
      {
        profiling::Scope profilePack(profiling::Context::active(),
                                     profiling::Category::band_pack,
                                     profiling::Context::mode());
        packLapackBandInto(matrix, m_band);
        if (auto *profile = profiling::Context::active())
          profile->countBandPack();
      }
      resizePivotWorkspaceIfNeeded();
      m_factored = false;
      factorize();
    } catch (...) {
      clearFailure(kInvalidState);
    }
    return *this;
  }

  /// Copy and factorize a packed band matrix.  Direct assembly should use
  /// bandWorkspace() and factorize(ridx), which avoids this copy path.
  LapackBandSolver &compute(const LapackBandMatrix &band,
                            Eigen::Index ridx = 0) {
    try {
      m_band = band;
      resizePivotWorkspaceIfNeeded();
      m_factored = false;
      factorize(ridx);
    } catch (...) {
      clearFailure(kInvalidState);
    }
    return *this;
  }

  /// Repack and refactorize while retaining the established band dimensions.
  LapackBandSolver &factorize(const SparseMatrix &matrix) {
    if (m_band.n == 0) {
      clearFailure(kInvalidState);
      return *this;
    }
    if (matrix.rows() != m_band.n || matrix.cols() != m_band.n) {
      clearFailure(kBandStructureMismatch);
      return *this;
    }
    resizePivotWorkspaceIfNeeded();

    {
      profiling::Scope profilePack(profiling::Context::active(),
                                   profiling::Category::band_pack,
                                   profiling::Context::mode());
      Eigen::Index matrixKl = 0;
      Eigen::Index matrixKu = 0;
      for (Eigen::Index column = 0; column < matrix.outerSize(); ++column) {
        for (SparseMatrix::InnerIterator entry(matrix, column); entry; ++entry) {
          matrixKl = std::max(matrixKl, entry.row() - entry.col());
          matrixKu = std::max(matrixKu, entry.col() - entry.row());
        }
      }
      if (matrixKl > m_band.kl || matrixKu > m_band.ku) {
        clearFailure(kBandStructureMismatch);
        return *this;
      }

      m_band.data.setZero();
      for (Eigen::Index column = 0; column < matrix.outerSize(); ++column) {
        for (SparseMatrix::InnerIterator entry(matrix, column); entry; ++entry)
          m_band.at(entry.row(), entry.col()) = entry.value();
      }
      if (auto *profile = profiling::Context::active())
        profile->countBandPack();
    }
    m_factored = false;
    factorize();
    return *this;
  }

  /// Return the mutable solver-owned LAPACK workspace.
  LapackBandMatrix &bandWorkspace() { return m_band; }
  const LapackBandMatrix &bandWorkspace() const { return m_band; }

  /// Allocate the solver-owned workspace for direct band assembly.
  LapackBandMatrix &bandWorkspace(Eigen::Index n, Eigen::Index kl,
                                  Eigen::Index ku) {
    if (n <= 0 || kl < 0 || ku < 0)
      throw std::invalid_argument("invalid LAPACK band dimensions");
    m_band.n = n;
    m_band.kl = kl;
    m_band.ku = ku;
    m_band.ldab = 2 * kl + ku + 1;
    m_band.data = LapackBandMatrix::Storage::Zero(m_band.ldab, n);
    resizePivotWorkspaceIfNeeded();
    m_factored = false;
    m_info = kInvalidState;
    return m_band;
  }

  /// Reuse the existing coefficient/factor workspace when its shape is
  /// unchanged; direct callers fill the returned storage before factorizing.
  LapackBandMatrix &ensureBandWorkspace(Eigen::Index n, Eigen::Index kl,
                                        Eigen::Index ku) {
    if (m_band.n != n || m_band.kl != kl || m_band.ku != ku ||
        m_band.data.rows() != 2 * kl + ku + 1 ||
        m_band.data.cols() != n)
      return bandWorkspace(n, kl, ku);
    resizePivotWorkspaceIfNeeded();
    m_factored = false;
    m_info = kInvalidState;
    return m_band;
  }

  /// Factorize the current workspace in place, optionally at a trailing
  /// principal subsystem.  This is the direct-band hot path.
  LapackBandSolver &factorize(Eigen::Index ridx = 0) {
    m_factored = false;
    factorizePacked(ridx);
    return *this;
  }

  /// Copy and refactorize a packed band matrix; this is not the direct-band
  /// hot path.
  LapackBandSolver &factorize(const LapackBandMatrix &band,
                              Eigen::Index ridx = 0) {
    if (band.n <= 0 || band.kl != m_band.kl || band.ku != m_band.ku ||
        band.ldab != m_band.ldab || band.n != m_band.n) {
      clearFailure(kBandStructureMismatch);
      return *this;
    }
    m_band = band;
    resizePivotWorkspaceIfNeeded();
    m_factored = false;
    factorize(ridx);
    return *this;
  }

  DenseMatrix solve(const DenseMatrix &rhs) {
    if (!m_factored || rhs.rows() != m_activeN || rhs.cols() <= 0) {
      m_info = kInvalidState;
      return DenseMatrix{};
    }

    if (auto *profile = profiling::Context::active())
      profile->countSolve(static_cast<long>(rhs.cols()));

    DenseMatrix result = rhs;

    {
      profiling::Scope profileSolve(profiling::Context::active(),
                                    profiling::Category::solve,
                                    profiling::Context::mode());
      // Experimental post-zgbtrf solve: replay the exact TRANS='N' ZGBTRS
      // operations with Eigen expressions, avoiding millions of tiny BLAS
      // calls under outer OpenMP concurrency.
      solveWithEigen(result);
    }
    return result;
  }

  lapack_int info() const { return m_info; }
  const LapackBandMatrix &band() const { return m_band; }
  Eigen::Index activeOffset() const { return m_activeOffset; }
  Eigen::Index activeSize() const { return m_activeN; }

private:
  // Retained as the local LAPACKE reference implementation for validation of
  // the experimental Eigen/C++ solve path.  zgbtrf and its pivots are shared.
  lapack_int solveWithLapacke(DenseMatrix &result) {
    return LAPACKE_zgbtrs(
        LAPACK_COL_MAJOR, 'N', static_cast<lapack_int>(m_activeN),
        static_cast<lapack_int>(m_band.kl), static_cast<lapack_int>(m_band.ku),
        static_cast<lapack_int>(result.cols()),
        m_band.data.data() + m_band.ldab * m_activeOffset,
        static_cast<lapack_int>(m_band.ldab), m_pivots.data(), result.data(),
        static_cast<lapack_int>(m_activeN));
  }

  // Exact LAPACK 3.11.0 ZGBTRS TRANS='N' replay for all RHS at once.  The
  // factor buffer is a zgbtrf result and no BLAS/LAPACK call is made here.
  void solveWithEigen(DenseMatrix &result) {
    using FactorMap =
        Eigen::Map<const DenseMatrix, Eigen::Unaligned, Eigen::OuterStride<>>;
    const FactorMap factors(
        m_band.data.data() + m_band.ldab * m_activeOffset,
        m_band.ldab, m_activeN, Eigen::OuterStride<>(m_band.ldab));
    const Eigen::Index n = m_activeN;
    const Eigen::Index kd = m_band.kl + m_band.ku;

    for (Eigen::Index j = 0; j < n - 1; ++j) {
      const Eigen::Index lm = std::min(m_band.kl, n - j - 1);
      const Eigen::Index pivot = m_pivots[static_cast<std::size_t>(j)] - 1;
      if (pivot != j)
        result.row(pivot).swap(result.row(j));
      if (lm != 0) {
        const auto multipliers = factors.col(j).segment(kd + 1, lm);
        result.middleRows(j + 1, lm).noalias() -=
            multipliers * result.row(j);
      }
    }

    for (Eigen::Index j = n - 1; j >= 0; --j) {
      result.row(j) /= factors(kd, j);
      const Eigen::Index len = std::min(kd, j);
      if (len != 0) {
        const auto upperSegment = factors.col(j).segment(kd - len, len);
        result.middleRows(j - len, len).noalias() -=
            upperSegment * result.row(j);
      }
    }
    m_info = 0;
  }

  void factorizePacked(Eigen::Index ridx = 0) {
    static_assert(std::is_same_v<Complex, lapack_complex_double>,
                  "LAPACKE must use the C++ complex representation");
    if (m_band.n <= 0 || ridx < 0 || ridx >= m_band.n ||
        m_band.data.rows() != m_band.ldab ||
        m_band.data.cols() != m_band.n) {
      clearFailure(kInvalidState);
      return;
    }
    m_activeOffset = ridx;
    m_activeN = m_band.n - ridx;

    {
      profiling::Scope profileFactorize(profiling::Context::active(),
                                        profiling::Category::factorization,
                                        profiling::Context::mode());
      m_info = LAPACKE_zgbtrf(
          LAPACK_COL_MAJOR, static_cast<lapack_int>(m_activeN),
          static_cast<lapack_int>(m_activeN),
          static_cast<lapack_int>(m_band.kl),
          static_cast<lapack_int>(m_band.ku),
          m_band.data.data() + m_band.ldab * m_activeOffset,
          static_cast<lapack_int>(m_band.ldab), m_pivots.data());
    }
    m_factored = (m_info == 0);
    if (auto *profile = profiling::Context::active())
      profile->countFactorize(true);
  }

  void clearFailure(lapack_int status) {
    m_info = status;
    m_factored = false;
  }

  void resizePivotWorkspaceIfNeeded() {
    if (m_pivots.size() != static_cast<std::size_t>(m_band.n))
      m_pivots.resize(static_cast<std::size_t>(m_band.n));
  }

  LapackBandMatrix m_band;
  std::vector<lapack_int> m_pivots;
  lapack_int m_info = kInvalidState;
  bool m_factored = false;
  Eigen::Index m_activeOffset = 0;
  Eigen::Index m_activeN = 0;
};

} // namespace SPARSESPEC::detail

#endif
