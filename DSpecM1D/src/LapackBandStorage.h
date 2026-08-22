#ifndef DSPECM1D_LAPACK_BAND_STORAGE_H
#define DSPECM1D_LAPACK_BAND_STORAGE_H

#include <algorithm>
#include <cstddef>
#include <complex>
#include <stdexcept>
#include <vector>

#include <Eigen/SparseCore>

#include "Profiling.h"

namespace SPARSESPEC::detail {

/// Internal column-major storage for a square LAPACK general-band matrix.
struct LapackBandMatrix {
  using Complex = std::complex<double>;
  using Storage = Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic,
                                Eigen::ColMajor>;

  Eigen::Index n = 0;
  Eigen::Index kl = 0;
  Eigen::Index ku = 0;
  Eigen::Index ldab = 0;
  // LAPACK's general-band buffer is column-major with ldab rows.  Keeping an
  // Eigen dense matrix here makes the layout and the pointer passed to
  // LAPACKE identical, while retaining Eigen's normal indexing for packing.
  Storage data;

  Complex &at(Eigen::Index row, Eigen::Index column) {
    const Eigen::Index bandRow = kl + ku + row - column;
    return data(bandRow, column);
  }

  const Complex &at(Eigen::Index row, Eigen::Index column) const {
    const Eigen::Index bandRow = kl + ku + row - column;
    return data(bandRow, column);
  }

  /// The physical input coefficient rows in the LAPACK band buffer.
  auto coefficients() {
    return data.middleRows(kl, kl + ku + 1);
  }

  auto coefficients() const {
    return data.middleRows(kl, kl + ku + 1);
  }
};

inline void packLapackBandInto(
    const Eigen::SparseMatrix<std::complex<double>> &matrix,
    LapackBandMatrix &band) {
  if (matrix.rows() != matrix.cols())
    throw std::invalid_argument(
        "LAPACK general-band storage requires a square matrix");

  band.n = matrix.rows();
  band.kl = 0;
  band.ku = 0;
  for (Eigen::Index column = 0; column < matrix.outerSize(); ++column) {
    for (Eigen::SparseMatrix<std::complex<double>>::InnerIterator entry(
             matrix, column);
         entry; ++entry) {
      band.kl = std::max(band.kl, entry.row() - entry.col());
      band.ku = std::max(band.ku, entry.col() - entry.row());
    }
  }

  band.ldab = 2 * band.kl + band.ku + 1;
  band.data = LapackBandMatrix::Storage::Zero(band.ldab, band.n);
  for (Eigen::Index column = 0; column < matrix.outerSize(); ++column) {
    for (Eigen::SparseMatrix<std::complex<double>>::InnerIterator entry(
             matrix, column);
         entry; ++entry) {
      band.at(entry.row(), entry.col()) = entry.value();
    }
  }
}

inline LapackBandMatrix
packLapackBand(const Eigen::SparseMatrix<std::complex<double>> &matrix) {
  profiling::Scope profilePack(profiling::Context::active(),
                               profiling::Category::band_pack,
                               profiling::Context::mode());
  LapackBandMatrix band;
  packLapackBandInto(matrix, band);
  if (auto *profile = profiling::Context::active())
    profile->countBandPack();
  return band;
}

} // namespace SPARSESPEC::detail

#endif
