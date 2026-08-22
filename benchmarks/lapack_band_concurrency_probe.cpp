#include <lapacke.h>

#include <chrono>
#include <complex>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <omp.h>

namespace {

using Complex = std::complex<double>;
using Clock = std::chrono::steady_clock;

// The 35 mHz spheroidal profile observed dimensions 17..1554 and bandwidths
// kl,ku 4..15.  Use the largest observed system to keep this probe on the
// production LAPACK general-band code path rather than a toy matrix.
constexpr int kN = 1554;
constexpr int kKl = 15;
constexpr int kKu = 15;
constexpr int kNrhs = 3;
constexpr int kSystems = 128;
constexpr int kLdab = 2 * kKl + kKu + 1;

struct Work {
  std::vector<Complex> base;
  std::vector<Complex> band;
  std::vector<lapack_int> pivots;
  std::vector<Complex> rhsBase;
  std::vector<Complex> rhsWork;
};

std::vector<Complex> makeBand(int system) {
  std::vector<Complex> result(static_cast<std::size_t>(kLdab) * kN,
                              Complex(0.0, 0.0));
  for (int column = 0; column < kN; ++column) {
    for (int delta = -kKu; delta <= kKl; ++delta) {
      const int row = column + delta;
      if (row < 0 || row >= kN)
        continue;
      const int bandRow = kKl + kKu + delta;
      if (delta == 0) {
        result[static_cast<std::size_t>(bandRow) +
               static_cast<std::size_t>(kLdab) * column] =
            Complex(20.0 + 0.01 * (system + 1), 0.1);
      } else {
        const double scale = 0.002 / (1.0 + std::abs(delta));
        result[static_cast<std::size_t>(bandRow) +
               static_cast<std::size_t>(kLdab) * column] =
            Complex(scale * (1.0 + (column + system) % 7),
                    -scale * (1.0 + (column + 3 * system) % 5));
      }
    }
  }
  return result;
}

void resetFactorInputs(std::vector<Work> &work) {
#pragma omp parallel for schedule(static)
  for (int index = 0; index < kSystems; ++index) {
    work[index].band = work[index].base;
    work[index].pivots.assign(kN, 0);
  }
}

void factorAll(std::vector<Work> &work) {
#pragma omp parallel for schedule(static)
  for (int index = 0; index < kSystems; ++index) {
    const lapack_int info = LAPACKE_zgbtrf(
        LAPACK_COL_MAJOR, kN, kN, kKl, kKu, work[index].band.data(), kLdab,
        work[index].pivots.data());
    if (info != 0) {
      std::cerr << "zgbtrf failed with info=" << info << "\n";
      std::abort();
    }
  }
}

void resetRhs(std::vector<Work> &work) {
#pragma omp parallel for schedule(static)
  for (int index = 0; index < kSystems; ++index)
    work[index].rhsWork = work[index].rhsBase;
}

void solveAll(std::vector<Work> &work, bool allocateResult) {
#pragma omp parallel for schedule(static)
  for (int index = 0; index < kSystems; ++index) {
    std::vector<Complex> result;
    if (allocateResult)
      result = work[index].rhsBase;
    else
      work[index].rhsWork = work[index].rhsBase;
    Complex *rhs = allocateResult ? result.data() : work[index].rhsWork.data();
    const lapack_int info = LAPACKE_zgbtrs(
        LAPACK_COL_MAJOR, 'N', kN, kKl, kKu, kNrhs,
        work[index].band.data(), kLdab, work[index].pivots.data(), rhs, kN);
    if (info != 0) {
      std::cerr << "zgbtrs failed with info=" << info << "\n";
      std::abort();
    }
  }
}

void combinedAll(std::vector<Work> &work) {
#pragma omp parallel for schedule(static)
  for (int index = 0; index < kSystems; ++index) {
    work[index].band = work[index].base;
    work[index].pivots.assign(kN, 0);
    const lapack_int factorInfo = LAPACKE_zgbtrf(
        LAPACK_COL_MAJOR, kN, kN, kKl, kKu, work[index].band.data(), kLdab,
        work[index].pivots.data());
    work[index].rhsWork = work[index].rhsBase;
    const lapack_int solveInfo = LAPACKE_zgbtrs(
        LAPACK_COL_MAJOR, 'N', kN, kKl, kKu, kNrhs,
        work[index].band.data(), kLdab, work[index].pivots.data(),
        work[index].rhsWork.data(), kN);
    if (factorInfo != 0 || solveInfo != 0) {
      std::cerr << "combined call failed: zgbtrf=" << factorInfo
                << " zgbtrs=" << solveInfo << "\n";
      std::abort();
    }
  }
}

template <class Prepare, class Run>
double measure(Prepare prepare, Run run, std::vector<Work> &work) {
  prepare(work);
  run(work); // untimed warm-up
  prepare(work);
  const auto start = Clock::now();
  run(work); // one measured run
  return std::chrono::duration<double>(Clock::now() - start).count();
}

} // namespace

int main() {
  std::vector<Work> work(kSystems);
  for (int index = 0; index < kSystems; ++index) {
    work[index].base = makeBand(index);
    work[index].band.resize(work[index].base.size());
    work[index].pivots.resize(kN);
    work[index].rhsBase.resize(static_cast<std::size_t>(kN) * kNrhs);
    for (int rhs = 0; rhs < kNrhs; ++rhs)
      for (int row = 0; row < kN; ++row)
        work[index].rhsBase[static_cast<std::size_t>(rhs) * kN + row] =
            Complex(1.0 + 0.001 * row, -0.25 + 0.0001 * (row + rhs));
    work[index].rhsWork.resize(work[index].rhsBase.size());
  }


  std::cout << "threads\tmetric\tseconds\tn\tkl\tku\tnrhs\tsystems\n";
  std::cout << std::setprecision(17);
  for (const int threads : {1, 10, 32, 50}) {
    omp_set_num_threads(threads);
    const double factor = measure(resetFactorInputs, factorAll, work);
    const double solveAllocate =
        measure(resetRhs, [&](auto &items) { solveAll(items, true); }, work);
    const double solveReused =
        measure(resetRhs, [&](auto &items) { solveAll(items, false); }, work);
    const double combined = measure(resetFactorInputs,
                                    combinedAll, work);
    std::cout << threads << "\tzgbtrf_only\t" << factor << '\t' << kN << '\t'
              << kKl << '\t' << kKu << '\t' << kNrhs << '\t' << kSystems
              << '\n';
    std::cout << threads << "\tzgbtrs_only_allocate\t" << solveAllocate << '\t'
              << kN << '\t' << kKl << '\t' << kKu << '\t' << kNrhs << '\t'
              << kSystems << '\n';
    std::cout << threads << "\tzgbtrs_only_reused\t" << solveReused << '\t'
              << kN << '\t' << kKl << '\t' << kKu << '\t' << kNrhs << '\t'
              << kSystems << '\n';
    std::cout << threads << "\tzgbtrf_plus_zgbtrs\t" << combined << '\t' << kN
              << '\t' << kKl << '\t' << kKu << '\t' << kNrhs << '\t'
              << kSystems << '\n';
  }
}
