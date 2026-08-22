#include <lapacke.h>

#include <chrono>
#include <complex>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>

#include <omp.h>

namespace {
using Complex = std::complex<double>;
using Clock = std::chrono::steady_clock;
constexpr int kN = 17;       // observed 35 mHz spheroidal minimum
constexpr int kKl = 4;       // observed minimum lower bandwidth
constexpr int kKu = 4;       // observed minimum upper bandwidth
constexpr int kNrhs = 3;
constexpr int kSystems = 128;
constexpr int kRepeats = 512;
constexpr int kLdab = 2 * kKl + kKu + 1;

struct Work {
  std::vector<Complex> band;
  std::vector<lapack_int> pivots;
  std::vector<Complex> rhs;
  std::vector<Complex> reused;
};

void solveMany(std::vector<Work> &work, bool allocate) {
#pragma omp parallel for schedule(static)
  for (int system = 0; system < kSystems; ++system) {
    for (int repeat = 0; repeat < kRepeats; ++repeat) {
      std::vector<Complex> result;
      if (allocate)
        result = work[system].rhs;
      else
        work[system].reused = work[system].rhs;
      Complex *target = allocate ? result.data() : work[system].reused.data();
      const lapack_int info = LAPACKE_zgbtrs(
          LAPACK_COL_MAJOR, 'N', kN, kKl, kKu, kNrhs,
          work[system].band.data(), kLdab, work[system].pivots.data(), target,
          kN);
      if (info != 0) {
        std::cerr << "zgbtrs failed with info=" << info << '\n';
        std::abort();
      }
    }
  }
}

double measure(std::vector<Work> &work, bool allocate) {
  solveMany(work, allocate); // untimed warm-up
  const auto start = Clock::now();
  solveMany(work, allocate); // one measured run
  return std::chrono::duration<double>(Clock::now() - start).count();
}
} // namespace

int main() {
  std::vector<Work> work(kSystems);
  for (int system = 0; system < kSystems; ++system) {
    work[system].band.assign(static_cast<std::size_t>(kLdab) * kN,
                             Complex(0.0, 0.0));
    for (int column = 0; column < kN; ++column) {
      for (int delta = -kKu; delta <= kKl; ++delta) {
        const int row = column + delta;
        if (row < 0 || row >= kN)
          continue;
        const int bandRow = kKl + kKu + delta;
        const double value = delta == 0 ? 20.0 : 0.002;
        work[system].band[static_cast<std::size_t>(bandRow) +
                          static_cast<std::size_t>(kLdab) * column] =
            Complex(value, delta == 0 ? 0.1 : -0.001);
      }
    }
    work[system].pivots.resize(kN);
    const lapack_int info = LAPACKE_zgbtrf(
        LAPACK_COL_MAJOR, kN, kN, kKl, kKu, work[system].band.data(), kLdab,
        work[system].pivots.data());
    if (info != 0) {
      std::cerr << "zgbtrf failed with info=" << info << '\n';
      return 1;
    }
    work[system].rhs.resize(static_cast<std::size_t>(kN) * kNrhs);
    work[system].reused.resize(work[system].rhs.size());
    for (std::size_t index = 0; index < work[system].rhs.size(); ++index)
      work[system].rhs[index] = Complex(1.0 + 0.01 * index, -0.1);
  }

  std::cout << "threads\tmetric\tseconds\tn\tkl\tku\tnrhs\tsystems\trepeats\n"
            << std::setprecision(17);
  for (const int threads : {1, 10, 32, 50}) {
    omp_set_num_threads(threads);
    std::cout << threads << "\tzgbtrs_allocate\t" << measure(work, true)
              << '\t' << kN << '\t' << kKl << '\t' << kKu << '\t' << kNrhs
              << '\t' << kSystems << '\t' << kRepeats << '\n';
    std::cout << threads << "\tzgbtrs_reused\t" << measure(work, false) << '\t'
              << kN << '\t' << kKl << '\t' << kKu << '\t' << kNrhs << '\t'
              << kSystems << '\t' << kRepeats << '\n';
  }
}
