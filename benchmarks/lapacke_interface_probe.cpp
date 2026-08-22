#include <lapacke.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>

#include <omp.h>

namespace {

using Complex = std::complex<double>;
using Clock = std::chrono::steady_clock;

constexpr int kNrhs = 3;
constexpr int kSystems = 128;

struct Configuration {
  int n;
  int kl;
  int ku;
  int solveRepeats;
};

struct System {
  std::vector<Complex> base;
  std::vector<Complex> factors;
  std::vector<lapack_int> pivots;
  std::vector<Complex> rhs;
  std::vector<Complex> highResult;
  std::vector<Complex> workResult;
};

int ldab(const Configuration &config) {
  return 2 * config.kl + config.ku + 1;
}

std::vector<Complex> makeBand(const Configuration &config, int system) {
  std::vector<Complex> band(
      static_cast<std::size_t>(ldab(config)) * config.n, Complex{});
  for (int column = 0; column < config.n; ++column) {
    for (int delta = -config.ku; delta <= config.kl; ++delta) {
      const int row = column + delta;
      if (row < 0 || row >= config.n)
        continue;
      const int bandRow = config.kl + config.ku + delta;
      const double scale = 0.002 / (1.0 + std::abs(delta));
      band[static_cast<std::size_t>(bandRow) +
           static_cast<std::size_t>(ldab(config)) * column] =
          delta == 0
              ? Complex(20.0 + 0.01 * (system + 1), 0.1)
              : Complex(scale * (1.0 + (column + system) % 7),
                        -scale * (1.0 + (column + 3 * system) % 5));
    }
  }
  return band;
}

std::vector<System> makeSystems(const Configuration &config) {
  std::vector<System> systems(kSystems);
  for (int index = 0; index < kSystems; ++index) {
    auto &system = systems[index];
    system.base = makeBand(config, index);
    system.factors = system.base;
    system.pivots.resize(config.n);
    const lapack_int info = LAPACKE_zgbtrf_work(
        LAPACK_COL_MAJOR, config.n, config.n, config.kl, config.ku,
        system.factors.data(), ldab(config), system.pivots.data());
    if (info != 0) {
      std::cerr << "initial zgbtrf_work failed with info=" << info << '\n';
      std::exit(1);
    }
    system.rhs.resize(static_cast<std::size_t>(config.n) * kNrhs);
    for (int rhs = 0; rhs < kNrhs; ++rhs)
      for (int row = 0; row < config.n; ++row)
        system.rhs[static_cast<std::size_t>(rhs) * config.n + row] =
            Complex(1.0 + 0.001 * row, -0.25 + 0.0001 * (row + rhs));
    system.highResult.resize(system.rhs.size());
    system.workResult.resize(system.rhs.size());
  }
  return systems;
}

void solveAll(std::vector<System> &systems, const Configuration &config,
              bool workInterface) {
#pragma omp parallel for schedule(static)
  for (int index = 0; index < kSystems; ++index) {
    auto &system = systems[index];
    auto &result = workInterface ? system.workResult : system.highResult;
    for (int repeat = 0; repeat < config.solveRepeats; ++repeat) {
      result = system.rhs;
      const lapack_int info = workInterface
                                  ? LAPACKE_zgbtrs_work(
                                        LAPACK_COL_MAJOR, 'N', config.n,
                                        config.kl, config.ku, kNrhs,
                                        system.factors.data(), ldab(config),
                                        system.pivots.data(), result.data(),
                                        config.n)
                                  : LAPACKE_zgbtrs(
                                        LAPACK_COL_MAJOR, 'N', config.n,
                                        config.kl, config.ku, kNrhs,
                                        system.factors.data(), ldab(config),
                                        system.pivots.data(), result.data(),
                                        config.n);
      if (info != 0) {
        std::cerr << "zgbtrs interface failed with info=" << info << '\n';
        std::abort();
      }
    }
  }
}

double measureSolve(std::vector<System> &systems, const Configuration &config,
                    bool workInterface) {
  solveAll(systems, config, workInterface); // one untimed warm-up
  const auto start = Clock::now();
  solveAll(systems, config, workInterface); // one measured run
  return std::chrono::duration<double>(Clock::now() - start).count();
}

double maximumSolveDifference(const std::vector<System> &systems) {
  double maximum = 0.0;
  for (const auto &system : systems)
    for (std::size_t index = 0; index < system.highResult.size(); ++index)
      maximum = std::max(
          maximum, std::abs(system.highResult[index] - system.workResult[index]));
  return maximum;
}

void factorAll(std::vector<System> &systems, const Configuration &config,
               bool workInterface) {
#pragma omp parallel for schedule(static)
  for (int index = 0; index < kSystems; ++index) {
    auto &system = systems[index];
    system.factors = system.base;
    const lapack_int info = workInterface
                                ? LAPACKE_zgbtrf_work(
                                      LAPACK_COL_MAJOR, config.n, config.n,
                                      config.kl, config.ku,
                                      system.factors.data(), ldab(config),
                                      system.pivots.data())
                                : LAPACKE_zgbtrf(
                                      LAPACK_COL_MAJOR, config.n, config.n,
                                      config.kl, config.ku,
                                      system.factors.data(), ldab(config),
                                      system.pivots.data());
    if (info != 0) {
      std::cerr << "zgbtrf interface failed with info=" << info << '\n';
      std::abort();
    }
  }
}

double measureFactor(std::vector<System> &systems,
                     const Configuration &config, bool workInterface) {
  factorAll(systems, config, workInterface); // one untimed warm-up
  const auto start = Clock::now();
  factorAll(systems, config, workInterface); // one measured run
  return std::chrono::duration<double>(Clock::now() - start).count();
}

void printPair(const char *operation, const Configuration &config, int threads,
               double highSeconds, double workSeconds, double difference,
               int repeats) {
  std::cout << threads << '\t' << operation << "\thigh_level\t"
            << highSeconds << '\t' << config.n << '\t' << config.kl << '\t'
            << config.ku << '\t' << kNrhs << '\t' << kSystems << '\t'
            << repeats << '\t' << difference << '\n';
  std::cout << threads << '\t' << operation << "\twork\t" << workSeconds
            << '\t' << config.n << '\t' << config.kl << '\t' << config.ku
            << '\t' << kNrhs << '\t' << kSystems << '\t' << repeats << '\t'
            << difference << '\n';
}

} // namespace

int main() {
  const std::vector<Configuration> configurations{
      {17, 4, 4, 512}, {256, 8, 8, 32}, {1554, 15, 15, 1}};
  LAPACKE_set_nancheck(1);

  std::cout << "threads\toperation\tinterface\tseconds\tn\tkl\tku\tnrhs\t"
               "systems\trepeats\tmax_abs_difference\n"
            << std::setprecision(17);
  for (std::size_t configIndex = 0; configIndex < configurations.size();
       ++configIndex) {
    const auto &config = configurations[configIndex];
    auto systems = makeSystems(config);
    for (const int threads : {1, 10, 32, 50}) {
      omp_set_num_threads(threads);
      double highSolve;
      double workSolve;
      if ((configIndex + static_cast<std::size_t>(threads)) % 2 == 0) {
        highSolve = measureSolve(systems, config, false);
        workSolve = measureSolve(systems, config, true);
      } else {
        workSolve = measureSolve(systems, config, true);
        highSolve = measureSolve(systems, config, false);
      }
      const double solveDifference = maximumSolveDifference(systems);
      if (!std::isfinite(solveDifference) || solveDifference > 1.0e-13) {
        std::cerr << "solve interfaces disagree: " << solveDifference << '\n';
        return 1;
      }
      printPair("zgbtrs", config, threads, highSolve, workSolve,
                solveDifference, config.solveRepeats);

      // Factorization is a secondary check: one factorization per independent
      // system, with the same one-warm-up/one-measured policy.
      const double highFactor = measureFactor(systems, config, false);
      std::vector<std::vector<Complex>> highFactors;
      std::vector<std::vector<lapack_int>> highPivots;
      highFactors.reserve(systems.size());
      highPivots.reserve(systems.size());
      for (const auto &system : systems) {
        highFactors.push_back(system.factors);
        highPivots.push_back(system.pivots);
      }
      const double workFactor = measureFactor(systems, config, true);
      double factorDifference = 0.0;
      for (std::size_t system = 0; system < systems.size(); ++system) {
        if (highPivots[system] != systems[system].pivots) {
          std::cerr << "factorization pivots disagree\n";
          return 1;
        }
        for (std::size_t index = 0; index < highFactors[system].size(); ++index)
          factorDifference = std::max(
              factorDifference,
              std::abs(highFactors[system][index] -
                       systems[system].factors[index]));
      }
      if (!std::isfinite(factorDifference) || factorDifference > 1.0e-13) {
        std::cerr << "factor interfaces disagree: " << factorDifference << '\n';
        return 1;
      }
      printPair("zgbtrf", config, threads, highFactor, workFactor,
                factorDifference, 1);
    }
  }
}
