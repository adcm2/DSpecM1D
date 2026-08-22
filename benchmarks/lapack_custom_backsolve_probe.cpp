#include <lapacke.h>

#include <Eigen/Dense>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include <omp.h>

namespace {

using Complex = std::complex<double>;
using Matrix = Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic,
                             Eigen::ColMajor>;
using Clock = std::chrono::steady_clock;

struct FactorSystem {
  int n = 0;
  int kl = 0;
  int ku = 0;
  int ldab = 0;
  std::vector<Complex> factors;
  std::vector<lapack_int> pivots;
  Matrix rhs;
};

// This is the exact TRANS='N' solve performed by LAPACK 3.11.0 ZGBTRS:
// replay the row interchanges and L updates, then solve U.  The input buffer
// is a zgbtrf result and no BLAS/LAPACK routine is called here.
void customSolve(const Complex *factorData, int n, int kl, int ku, int ldab,
                 const lapack_int *pivots, Matrix &result) {
  using FactorMap = Eigen::Map<const Matrix, Eigen::Unaligned,
                               Eigen::OuterStride<>>;
  const FactorMap factors(factorData, ldab, n, Eigen::OuterStride<>(ldab));
  const int kd = kl + ku;

  for (int j = 0; j < n - 1; ++j) {
    const int lm = std::min(kl, n - j - 1);
    const int pivot = static_cast<int>(pivots[j]) - 1;
    if (pivot != j)
      result.row(pivot).swap(result.row(j));
    if (lm != 0) {
      const auto multipliers = factors.col(j).segment(kd + 1, lm);
      result.middleRows(j + 1, lm).noalias() -=
          multipliers * result.row(j);
    }
  }

  for (int j = n - 1; j >= 0; --j) {
    result.row(j) /= factors(kd, j);
    const int len = std::min(kd, j);
    if (len != 0) {
      const auto upperSegment = factors.col(j).segment(kd - len, len);
      result.middleRows(j - len, len).noalias() -=
          upperSegment * result.row(j);
    }
  }
}

Matrix makeBandDense(int n, int kl, int ku, bool inducePivot,
                     std::mt19937_64 &rng) {
  std::uniform_real_distribution<double> random(-0.15, 0.15);
  Matrix matrix = Matrix::Zero(n, n);
  for (int col = 0; col < n; ++col) {
    for (int delta = -ku; delta <= kl; ++delta) {
      const int row = col + delta;
      if (row < 0 || row >= n)
        continue;
      if (delta == 0)
        matrix(row, col) = Complex(8.0 + 0.03 * col, 0.25);
      else
        matrix(row, col) = Complex(random(rng), random(rng));
    }
  }
  if (inducePivot && kl > 0) {
    const int first = std::min(2, n - 1);
    matrix(first, 0) = Complex(7.0, -1.0);
    matrix(0, 0) = Complex(1.0e-12, 0.0);
    if (n > 4 && kl > 1) {
      matrix(4, 2) = Complex(6.0, 0.5);
      matrix(2, 2) = Complex(1.0e-12, 0.0);
    }
  }
  return matrix;
}

std::vector<Complex> packBand(const Matrix &matrix, int kl, int ku) {
  const int n = static_cast<int>(matrix.rows());
  const int ldab = 2 * kl + ku + 1;
  std::vector<Complex> packed(static_cast<std::size_t>(ldab) * n,
                              Complex(0.0, 0.0));
  const int kd = kl + ku;
  for (int col = 0; col < n; ++col)
    for (int row = std::max(0, col - ku); row <= std::min(n - 1, col + kl);
         ++row)
      packed[static_cast<std::size_t>(kd + row - col) +
             static_cast<std::size_t>(ldab) * col] = matrix(row, col);
  return packed;
}

double infinityNorm(const Matrix &matrix) {
  double result = 0.0;
  for (int row = 0; row < matrix.rows(); ++row) {
    double rowSum = 0.0;
    for (int col = 0; col < matrix.cols(); ++col)
      rowSum += std::abs(matrix(row, col));
    result = std::max(result, rowSum);
  }
  return result;
}

double normalizedDifference(const Matrix &a, const Matrix &b) {
  return infinityNorm(a - b) / (1.0 + infinityNorm(a) + infinityNorm(b));
}

double normalizedResidual(const Matrix &matrix, const Matrix &solution,
                          const Matrix &rhs) {
  return infinityNorm(matrix * solution - rhs) /
         (1.0 + infinityNorm(matrix) * infinityNorm(solution) +
          infinityNorm(rhs));
}

struct ValidationCase {
  std::string name;
  int n;
  int kl;
  int ku;
  int nrhs;
  int ridx;
  bool inducePivot;
};

bool validateCase(const ValidationCase &config, std::mt19937_64 &rng,
                  int &pivotCases, double &maxDifference,
                  double &maxResidual) {
  Matrix original = makeBandDense(config.n, config.kl, config.ku,
                                  config.inducePivot, rng);
  if (config.inducePivot && config.ridx > 0) {
    const int pivotRow = std::min(config.n - 1, config.ridx + 1);
    original(pivotRow, config.ridx) = Complex(7.0, -1.0);
    original(config.ridx, config.ridx) = Complex(1.0e-12, 0.0);
  }
  const auto packed = packBand(original, config.kl, config.ku);
  const int activeN = config.n - config.ridx;
  const int ldab = 2 * config.kl + config.ku + 1;
  std::vector<Complex> factors = packed;
  std::vector<lapack_int> pivots(static_cast<std::size_t>(config.n), 0);
  const lapack_int factorInfo = LAPACKE_zgbtrf(
      LAPACK_COL_MAJOR, activeN, activeN, config.kl, config.ku,
      factors.data() + static_cast<std::size_t>(ldab) * config.ridx, ldab,
      pivots.data());
  if (factorInfo != 0)
    return false;

  bool pivoted = false;
  for (int j = 0; j < activeN; ++j)
    pivoted |= pivots[j] - 1 != j;
  if (pivoted)
    ++pivotCases;

  Matrix exact(activeN, config.nrhs);
  for (int col = 0; col < config.nrhs; ++col)
    for (int row = 0; row < activeN; ++row)
      exact(row, col) = Complex(0.2 + 0.01 * row + 0.02 * col,
                                -0.1 + 0.005 * row - 0.01 * col);
  const Matrix active = original.bottomRightCorner(activeN, activeN);
  const Matrix rhs = active * exact;
  Matrix lapackResult = rhs;
  const lapack_int solveInfo = LAPACKE_zgbtrs(
      LAPACK_COL_MAJOR, 'N', activeN, config.kl, config.ku, config.nrhs,
      factors.data() + static_cast<std::size_t>(ldab) * config.ridx, ldab,
      pivots.data(), lapackResult.data(), activeN);
  if (solveInfo != 0)
    return false;
  Matrix customResult = rhs;
  customSolve(factors.data() + static_cast<std::size_t>(ldab) * config.ridx,
              activeN, config.kl, config.ku, ldab, pivots.data(),
              customResult);
  maxDifference = std::max(maxDifference,
                           normalizedDifference(customResult, lapackResult));
  maxResidual = std::max(maxResidual,
                         normalizedResidual(active, customResult, rhs));
  return maxDifference < 1.0e-12 && maxResidual < 1.0e-12;
}

FactorSystem makeBenchmarkSystem(int n, int kl, int ku, int nrhs,
                                 int system, std::mt19937_64 &rng) {
  FactorSystem result;
  result.n = n;
  result.kl = kl;
  result.ku = ku;
  result.ldab = 2 * kl + ku + 1;
  Matrix original = makeBandDense(n, kl, ku, false, rng);
  result.factors = packBand(original, kl, ku);
  result.pivots.resize(static_cast<std::size_t>(n));
  const lapack_int factorInfo = LAPACKE_zgbtrf(
      LAPACK_COL_MAJOR, n, n, kl, ku, result.factors.data(), result.ldab,
      result.pivots.data());
  if (factorInfo != 0)
    throw std::runtime_error("benchmark zgbtrf failed");
  result.rhs.resize(n, nrhs);
  for (int col = 0; col < nrhs; ++col)
    for (int row = 0; row < n; ++row)
      result.rhs(row, col) = Complex(1.0 + 0.001 * row + 0.01 * col,
                                     -0.2 + 0.0003 * (row + system));
  return result;
}

struct BenchmarkConfig {
  const char *name;
  int n;
  int kl;
  int ku;
  int systems;
};

template <class Solve>
double measure(const std::vector<FactorSystem> &systems, Solve solve) {
  std::vector<Matrix> rhs;
  rhs.reserve(systems.size());
  for (const auto &system : systems)
    rhs.push_back(system.rhs);
  auto run = [&]() {
#pragma omp parallel for schedule(static)
    for (int index = 0; index < static_cast<int>(systems.size()); ++index)
      solve(systems[static_cast<std::size_t>(index)],
            rhs[static_cast<std::size_t>(index)]);
  };
  run();
  for (std::size_t index = 0; index < systems.size(); ++index)
    rhs[index] = systems[index].rhs;
  const auto start = Clock::now();
  run();
  return std::chrono::duration<double>(Clock::now() - start).count();
}

} // namespace

int main() {
  std::mt19937_64 rng(0x51a7e5eedULL);
  const std::vector<ValidationCase> validation = {
      {"small_single", 3, 1, 1, 1, 0, false},
      {"small_multiple", 7, 2, 1, 3, 0, false},
      {"pivoted", 9, 2, 2, 4, 0, true},
      {"random_wide", 17, 4, 4, 3, 0, false},
      {"trailing", 23, 3, 5, 2, 4, false},
      {"trailing_pivoted", 19, 4, 3, 3, 3, true},
  };
  int pivotCases = 0;
  double maxDifference = 0.0;
  double maxResidual = 0.0;
  int validationPassed = 0;
  for (const auto &config : validation)
    if (validateCase(config, rng, pivotCases, maxDifference, maxResidual))
      ++validationPassed;

  std::cout << "validation_cases\t" << validation.size() << '\n';
  std::cout << "validation_passed\t" << validationPassed << '\n';
  std::cout << "pivot_exercising_cases\t" << pivotCases << '\n';
  std::cout << std::setprecision(17);
  std::cout << "max_normalized_difference\t" << maxDifference << '\n';
  std::cout << "max_normalized_residual\t" << maxResidual << '\n';
  std::cout << "threads\tcase\tbackend\tseconds\tn\tkl\tku\tnrhs\tsystems\n";

  const std::vector<BenchmarkConfig> benchmark = {
      {"n17", 17, 4, 4, 4096},
      {"n256", 256, 8, 8, 1024},
      {"n1554", 1554, 15, 15, 128},
  };
  for (const auto &config : benchmark) {
    std::vector<FactorSystem> systems;
    systems.reserve(config.systems);
    for (int index = 0; index < config.systems; ++index)
      systems.push_back(makeBenchmarkSystem(config.n, config.kl, config.ku,
                                             3, index, rng));
    for (const int threads : {1, 10, 32, 50}) {
      omp_set_num_threads(threads);
      const double lapack = measure(
          systems, [](const FactorSystem &system, Matrix &rhs) {
            const lapack_int info = LAPACKE_zgbtrs(
                LAPACK_COL_MAJOR, 'N', system.n, system.kl, system.ku, 3,
                system.factors.data(), system.ldab,
                const_cast<lapack_int *>(system.pivots.data()), rhs.data(),
                system.n);
            if (info != 0)
              std::abort();
          });
      const double custom = measure(
          systems, [](const FactorSystem &system, Matrix &rhs) {
            customSolve(system.factors.data(), system.n, system.kl,
                         system.ku, system.ldab, system.pivots.data(), rhs);
          });
      double benchmarkDifference = 0.0;
      for (const auto &system : systems) {
        Matrix lapackResult = system.rhs;
        Matrix customResult = system.rhs;
        const lapack_int info = LAPACKE_zgbtrs(
            LAPACK_COL_MAJOR, 'N', system.n, system.kl, system.ku, 3,
            system.factors.data(), system.ldab,
            const_cast<lapack_int *>(system.pivots.data()), lapackResult.data(),
            system.n);
        if (info != 0)
          std::abort();
        customSolve(system.factors.data(), system.n, system.kl, system.ku,
                     system.ldab, system.pivots.data(), customResult);
        benchmarkDifference = std::max(
            benchmarkDifference,
            normalizedDifference(customResult, lapackResult));
      }
      std::cout << threads << '\t' << config.name << "\tlapacke_zgbtrs\t"
                << lapack << '\t' << config.n << '\t' << config.kl << '\t'
                << config.ku << "\t3\t" << config.systems << '\n';
      std::cout << threads << '\t' << config.name << "\teigen_custom\t"
                << custom << '\t' << config.n << '\t' << config.kl << '\t'
                << config.ku << "\t3\t" << config.systems << '\n';
      std::cout << "benchmark_max_normalized_difference\t" << config.name
                << '\t' << threads << '\t' << benchmarkDifference << '\n';
    }
  }
}
