// Focused cleanup check: compile unchanged against pre-cleanup/current headers.
// Uses the existing test input generator; no production benchmark API is added.
#include <DSpecM1D/src/FullSpec.h>
#include <DSpecM1D/src/InputParametersNew.h>
#include "test_utils.h"
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <vector>

int main(int argc, char **argv) {
  if (argc != 2)
    return 2;
  std::ofstream out(argv[1], std::ios::binary);
  if (!out)
    return 3;
  std::cout << std::setprecision(17);
  std::cout << "threads=" << omp_get_max_threads() << " dynamic="
            << omp_get_dynamic() << " degrees=1:24 nq=5 relative_error=0.001"
            << " requested_mhz=5:35 tout_minutes=60 receivers=2\n";
  for (const bool attenuation : {false, true}) {
    DSpecMTest::TempDir tmp;
    DSpecMTest::ParameterOptions options;
    options.type = 3;
    options.attenuation = attenuation;
    options.lmin = 1;
    options.lmax = 24;
    options.f1 = options.f11 = options.f12 = 5.0;
    options.f2 = options.f21 = options.f22 = 35.0;
    options.tOutMinutes = 60.0;
    options.timeStepSec = 5.0;
    options.relativeError = 1e-3;
    const auto model = DSpecMTest::repoRoot() / "data" / "models" /
                       (attenuation ? "prem.200.no.txt"
                                    : "prem.200.noatten.txt");
    const auto path = DSpecMTest::writeFile(
        tmp.path() / "params.txt",
        DSpecMTest::makeParameterText(model.string(), options));
    InputParametersNew params(path.string(), 5, 2, 0.05, 1.0, 0.05, 0.0, 1);
    const auto &freq = params.freqFull();
    const auto mhz = [&](int idx) {
      return freq.f(idx) * 1000.0 / freq.timeNorm();
    };
    const std::vector<double> cutoffs{
        0.0, mhz(freq.i1()) / 2.0,
        mhz(freq.i1() + (freq.i2() - freq.i1()) / 2)};
    std::cout << "bins=" << freq.i2() - freq.i1() << " actual_mhz="
              << mhz(freq.i1()) << ':' << mhz(freq.i2() - 1)
              << " nskip=" << std::max(1, (freq.i2() - freq.i1()) / 20)
              << '\n';
    SPARSESPEC::SparseFSpec solver;
    for (std::size_t caseIdx = 0; caseIdx < cutoffs.size(); ++caseIdx) {
      params.setCowlingFrequencyMhz(cutoffs[caseIdx]);
      Eigen::MatrixXcd spectrum = solver.spectra(params);
      // One benchmark case only. Warm-up above; time whole multi-SEM calls,
      // including mesh/both-layout construction, solves, rotation and scaling.
      if (!attenuation && caseIdx == 2) {
        std::vector<double> elapsed;
        for (int repeat = 0; repeat < 5; ++repeat) {
          const auto start = std::chrono::steady_clock::now();
          spectrum = solver.spectra(params);
          elapsed.push_back(std::chrono::duration<double>(
              std::chrono::steady_clock::now() - start).count());
        }
        for (const double seconds : elapsed)
          std::cout << "benchmark_seconds=" << seconds << '\n';
        std::sort(elapsed.begin(), elapsed.end());
        std::cout << "benchmark_median_seconds=" << elapsed[2] << '\n';
      }
      if (!spectrum.real().array().isFinite().all() ||
          !spectrum.imag().array().isFinite().all() || spectrum.norm() == 0.0)
        return 4;
      const long rows = spectrum.rows(), cols = spectrum.cols();
      out.write(reinterpret_cast<const char *>(&rows), sizeof(rows));
      out.write(reinterpret_cast<const char *>(&cols), sizeof(cols));
      for (long col = 0; col < cols; ++col) {
        for (long row = 0; row < rows; ++row) {
          const double real = spectrum(row, col).real();
          const double imag = spectrum(row, col).imag();
          out.write(reinterpret_cast<const char *>(&real), sizeof(real));
          out.write(reinterpret_cast<const char *>(&imag), sizeof(imag));
        }
      }
      std::cout << "case=" << caseIdx << " attenuation=" << attenuation
                << " cutoff_mhz=" << cutoffs[caseIdx]
                << " norm=" << spectrum.norm() << '\n';
    }
  }
  return out.good() ? 0 : 3;
}
