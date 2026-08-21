#include <complex>

#include <gtest/gtest.h>
#include <Eigen/Sparse>

#include <DSpecM1D/src/LapackBandStorage.h>

namespace {

using Complex = std::complex<double>;

Eigen::SparseMatrix<Complex>
makeBandMatrix() {
  Eigen::SparseMatrix<Complex> matrix(5, 5);
  matrix.insert(0, 0) = Complex(1.0, 1.0);
  matrix.insert(1, 0) = Complex(2.0, -1.0);
  matrix.insert(2, 0) = Complex(0.0, 0.0); // Stored zero fixes the lower band.
  matrix.insert(1, 1) = Complex(3.0, 2.0);
  matrix.insert(2, 1) = Complex(4.0, -2.0);
  matrix.insert(0, 2) = Complex(5.0, 3.0);
  matrix.insert(1, 3) = Complex(6.0, -3.0);
  matrix.insert(2, 4) = Complex(7.0, 4.0);
  return matrix;
}

} // namespace

TEST(LapackBandStorageTests, PacksStoredSparseEntriesInColumnMajorLayout) {
  const auto matrix = makeBandMatrix();
  ASSERT_EQ(matrix.nonZeros(), 8);
  const auto band = SPARSESPEC::detail::packLapackBand(matrix);

  EXPECT_EQ(band.n, 5);
  EXPECT_EQ(band.kl, 2);
  EXPECT_EQ(band.ku, 2);
  EXPECT_EQ(band.ldab, 7);
  EXPECT_EQ(band.data.size(), 35U);
  EXPECT_EQ(band.data[4], Complex(1.0, 1.0));
  EXPECT_EQ(band.data[16], Complex(5.0, 3.0));
  EXPECT_EQ(band.data[30], Complex(7.0, 4.0));

  for (Eigen::Index row = 0; row < band.n; ++row) {
    for (Eigen::Index column = 0; column < band.n; ++column) {
      if (row - column <= band.kl && column - row <= band.ku)
        EXPECT_EQ(band.at(row, column), matrix.coeff(row, column));
    }
  }
  EXPECT_EQ(band.at(2, 0), Complex(0.0, 0.0));
  EXPECT_EQ(band.at(0, 2), Complex(5.0, 3.0));
  EXPECT_EQ(band.at(2, 4), Complex(7.0, 4.0));
}
