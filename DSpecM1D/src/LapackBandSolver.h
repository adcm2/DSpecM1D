#ifndef DSPECM1D_LAPACK_BAND_SOLVER_H
#define DSPECM1D_LAPACK_BAND_SOLVER_H

#include <algorithm>
#include <complex>
#include <complex.h>
#include <vector>

#include <lapacke.h>

#include "LapackBandStorage.h"

namespace SPARSESPEC::detail {

/// Small internal LU wrapper for a complex general-band matrix.
class LapackBandSolver {
public:
  using Complex = std::complex<double>;
  using SparseMatrix = Eigen::SparseMatrix<Complex>;
  using DenseMatrix = Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic>;

  // This is deliberately outside LAPACK's ordinary argument-error range.
  static constexpr lapack_int kBandStructureMismatch = -1001;
  static constexpr lapack_int kInvalidState = -1002;

  LapackBandSolver &compute(const SparseMatrix &matrix) {
    try {
      m_band = packLapackBand(matrix);
      m_pivots.assign(static_cast<std::size_t>(m_band.n), 0);
      m_factored = false;
      factorizePacked();
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

      std::fill(m_band.data.begin(), m_band.data.end(), Complex{});
      for (Eigen::Index column = 0; column < matrix.outerSize(); ++column) {
        for (SparseMatrix::InnerIterator entry(matrix, column); entry; ++entry)
          m_band.at(entry.row(), entry.col()) = entry.value();
      }
      if (auto *profile = profiling::Context::active())
        profile->countBandPack();
    }
    m_factored = false;
    factorizePacked();
    return *this;
  }

  DenseMatrix solve(const DenseMatrix &rhs) {
    if (!m_factored || rhs.rows() != m_band.n || rhs.cols() <= 0) {
      m_info = kInvalidState;
      return DenseMatrix{};
    }

    if (auto *profile = profiling::Context::active())
      profile->countSolve(static_cast<long>(rhs.cols()));

    std::vector<lapack_complex_double> b(
        static_cast<std::size_t>(rhs.size()));
    for (Eigen::Index column = 0; column < rhs.cols(); ++column)
      for (Eigen::Index row = 0; row < rhs.rows(); ++row)
        b[static_cast<std::size_t>(row + rhs.rows() * column)] =
            toLapack(rhs(row, column));

    const lapack_int n = static_cast<lapack_int>(m_band.n);
    const lapack_int nrhs = static_cast<lapack_int>(rhs.cols());
    {
      profiling::Scope profileSolve(profiling::Context::active(),
                                    profiling::Category::solve,
                                    profiling::Context::mode());
      m_info = LAPACKE_zgbtrs(
          LAPACK_COL_MAJOR, 'N', n, static_cast<lapack_int>(m_band.kl),
          static_cast<lapack_int>(m_band.ku), nrhs, m_factors.data(),
          static_cast<lapack_int>(m_band.ldab), m_pivots.data(), b.data(), n);
    }

    DenseMatrix result(rhs.rows(), rhs.cols());
    for (Eigen::Index column = 0; column < result.cols(); ++column)
      for (Eigen::Index row = 0; row < result.rows(); ++row)
        result(row, column) = fromLapack(
            b[static_cast<std::size_t>(row + result.rows() * column)]);
    return result;
  }

  lapack_int info() const { return m_info; }

private:
  static lapack_complex_double toLapack(const Complex &value) {
    return lapack_make_complex_double(value.real(), value.imag());
  }

  static Complex fromLapack(const lapack_complex_double &value) {
    return {lapack_complex_double_real(value),
            lapack_complex_double_imag(value)};
  }

  void factorizePacked() {
    const lapack_int n = static_cast<lapack_int>(m_band.n);
    {
      profiling::Scope profilePrepare(profiling::Context::active(),
                                      profiling::Category::band_pack,
                                      profiling::Context::mode());
      m_factors.resize(m_band.data.size());
      for (std::size_t index = 0; index < m_band.data.size(); ++index)
        m_factors[index] = toLapack(m_band.data[index]);
    }

    {
      profiling::Scope profileFactorize(profiling::Context::active(),
                                        profiling::Category::factorization,
                                        profiling::Context::mode());
      m_info = LAPACKE_zgbtrf(
          LAPACK_COL_MAJOR, n, n, static_cast<lapack_int>(m_band.kl),
          static_cast<lapack_int>(m_band.ku), m_factors.data(),
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

  LapackBandMatrix m_band;
  std::vector<lapack_int> m_pivots;
  std::vector<lapack_complex_double> m_factors;
  lapack_int m_info = kInvalidState;
  bool m_factored = false;
};

} // namespace SPARSESPEC::detail

#endif
