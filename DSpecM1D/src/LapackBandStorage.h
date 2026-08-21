#ifndef DSPECM1D_LAPACK_BAND_STORAGE_H
#define DSPECM1D_LAPACK_BAND_STORAGE_H

#include <algorithm>
#include <cstddef>
#include <complex>
#include <stdexcept>
#include <vector>

#include <Eigen/SparseCore>

namespace SPARSESPEC::detail {

/// Internal column-major storage for a square LAPACK general-band matrix.
struct LapackBandMatrix {
  using Complex = std::complex<double>;

  Eigen::Index n = 0;
  Eigen::Index kl = 0;
  Eigen::Index ku = 0;
  Eigen::Index ldab = 0;
  std::vector<Complex> data;

  Complex &at(Eigen::Index row, Eigen::Index column) {
    const Eigen::Index bandRow = kl + ku + row - column;
    return data[static_cast<std::size_t>(bandRow + ldab * column)];
  }

  const Complex &at(Eigen::Index row, Eigen::Index column) const {
    const Eigen::Index bandRow = kl + ku + row - column;
    return data[static_cast<std::size_t>(bandRow + ldab * column)];
  }
};

inline LapackBandMatrix
packLapackBand(const Eigen::SparseMatrix<std::complex<double>> &matrix) {
  if (matrix.rows() != matrix.cols())
    throw std::invalid_argument(
        "LAPACK general-band storage requires a square matrix");

  LapackBandMatrix band;
  band.n = matrix.rows();
  for (Eigen::Index column = 0; column < matrix.outerSize(); ++column) {
    for (Eigen::SparseMatrix<std::complex<double>>::InnerIterator entry(
             matrix, column);
         entry; ++entry) {
      band.kl = std::max(band.kl, entry.row() - entry.col());
      band.ku = std::max(band.ku, entry.col() - entry.row());
    }
  }

  band.ldab = 2 * band.kl + band.ku + 1;
  band.data.assign(static_cast<std::size_t>(band.ldab * band.n), {});
  for (Eigen::Index column = 0; column < matrix.outerSize(); ++column) {
    for (Eigen::SparseMatrix<std::complex<double>>::InnerIterator entry(
             matrix, column);
         entry; ++entry) {
      band.at(entry.row(), entry.col()) = entry.value();
    }
  }
  return band;
}

} // namespace SPARSESPEC::detail

#endif
