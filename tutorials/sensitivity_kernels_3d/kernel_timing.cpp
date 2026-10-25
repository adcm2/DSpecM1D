#include <DSpecM1D/All>
#include <DSpecM1D/ModelInput>
#include <DSpecM1D/Windowed3DKernels>
#include <DSpecM1D/src/InputParser.h>
#include <DSpecM1D/src/NormClass.h>
#include <DSpecM1D/src/SourceInfo.h>

#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/resource.h>
#include <vector>

namespace fs = std::filesystem;
using Complex = std::complex<double>;
constexpr int nt = 128;
constexpr int bins = nt / 2 + 1;
constexpr std::array<const char *, 5> moduli{"A", "C", "F", "L", "N"};

std::string paramsText(const fs::path &model, const fs::path &output, int lmax) {
  return "\"" + (output / "unused").string() + "\"\n\"" + model.string() +
      "\"\n4\n0\n2\n0\n0\n1e-5\n0\n" + std::to_string(lmax) +
      "\n0.1\n0.8\n30\n20\n0\n0\n0\n0\n50\n10\n20\n1e19\n2e18\n-3e18\n-5e18\n4e18\n-5e18\n3\n1\n45 90\n";
}

Eigen::VectorXd hannWindow() {
  Eigen::VectorXd window = Eigen::VectorXd::Zero(nt);
  for (int n = 0; n <= 90; ++n)
    window[n] = 0.5 - 0.5 * std::cos(2.0 * std::numbers::pi * n / 90.0);
  return window;
}

Eigen::MatrixXcd objectiveWeights(int objective) {
  Eigen::MatrixXcd q = Eigen::MatrixXcd::Zero(3, bins);
  q(objective, 1) = Complex(0.7, 0.2);
  q(objective, 2) = Complex(-0.25, 0.45);
  return q;
}

std::vector<SPARSESPEC::Field3DPoint> readPoints(const fs::path &path) {
  std::ifstream file(path);
  if (!file)
    throw std::runtime_error("cannot read shared tutorial samples CSV");
  std::string line;
  std::getline(file, line);
  std::vector<SPARSESPEC::Field3DPoint> points;
  while (std::getline(file, line)) {
    std::istringstream row(line);
    std::array<std::string, 11> values;
    for (auto &value : values)
      if (!std::getline(row, value, ','))
        throw std::runtime_error("invalid shared tutorial sample row");
    if (values[0] == "J_Z" && values[1] == "A")
      points.push_back({std::stod(values[2]), std::stod(values[3]),
                        std::stod(values[4]), std::stoi(values[5])});
  }
  if (points.empty())
    throw std::runtime_error("shared samples contain no A-modulus coordinates");
  return points;
}

long peakRssBytes() {
  rusage usage{};
  if (getrusage(RUSAGE_SELF, &usage) != 0)
    throw std::runtime_error("getrusage failed");
#if defined(__APPLE__)
  return usage.ru_maxrss;
#else
  return usage.ru_maxrss * 1024L;
#endif
}

int main(int argc, char **argv) try {
  if (argc != 7)
    throw std::runtime_error(
        "usage: dspec_kernel_timing MODEL SAMPLES OUTPUT LMAX MEASUREMENTS REPEATS");
  const auto processStart = std::chrono::steady_clock::now();
  const fs::path modelPath = fs::absolute(argv[1]);
  const fs::path samplesPath = fs::absolute(argv[2]);
  const fs::path output = fs::absolute(argv[3]);
  const int lmax = std::stoi(argv[4]);
  const int measurementCount = std::stoi(argv[5]);
  const int repeats = std::stoi(argv[6]);
  if ((lmax != 3 && lmax != 12) || (measurementCount != 1 && measurementCount != 2) ||
      repeats < 1 || repeats > 5)
    throw std::runtime_error("lmax must be 3 or 12, measurements 1 or 2, repeats 1..5");
  if (!fs::is_directory(output) || fs::directory_iterator(output) != fs::directory_iterator{})
    throw std::runtime_error("output must be an existing empty directory");

  const auto paramPath = output / "parameters.txt";
  {
    std::ofstream file(paramPath);
    if (!file)
      throw std::runtime_error("cannot create temporary parameter file in output");
    file << paramsText(modelPath, output, lmax);
  }
  InputParametersNew input(paramPath.string(), 8, 1, 0.05, 1.0, 0.1, 0.0, 1);
  auto &parameters = input.inputParameters();
  auto &frequency = input.freqFull();
  auto &source = input.cmt();
  auto &base = input.earthModel();
  Full1D::SEM sem(base, 0.05, 8, lmax, source);
  const auto points = readPoints(samplesPath);
  if (frequency.nt() != nt || frequency.i1() > 1 || frequency.i2() < 4 ||
      std::abs(frequency.dt() * frequency.timeNorm() - 20.0) > 1e-12 ||
      std::abs(frequency.ep() / frequency.timeNorm() - 1.0 / 360.0) > 1e-14)
    throw std::runtime_error("resolved DSpec grid differs from nt=128, dt=20, epsilon=1/360, bins 1..3");
  std::vector<std::array<double, 2>> noAngles;
  const auto window = hannWindow();
  const auto setupEnd = std::chrono::steady_clock::now();
  const auto evaluate = [&](int objective) {
    return SPARSESPEC::Windowed3DKernelEngine::evaluate(
        sem, source, parameters, frequency, window, objectiveWeights(objective),
        1, 4, points, noAngles);
  };

  std::vector<double> apiSeconds;
  std::vector<decltype(evaluate(0))> measured;
  const auto warmupStart = std::chrono::steady_clock::now();
  (void)evaluate(0);
  if (measurementCount == 2)
    (void)evaluate(1);
  const double warmupSeconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - warmupStart).count();
  for (int repeat = 0; repeat < repeats; ++repeat) {
    for (int objective = 0; objective < measurementCount; ++objective) {
      const auto start = std::chrono::steady_clock::now();
      measured.push_back(evaluate(objective));
      apiSeconds.push_back(std::chrono::duration<double>(
          std::chrono::steady_clock::now() - start).count());
    }
  }
  if (measured.empty())
    throw std::runtime_error("no timed kernel evaluations were produced");
  std::ofstream csv(output / "kernels.csv");
  if (!csv)
    throw std::runtime_error("cannot write actual timed kernels");
  csv << "objective,modulus,radius_m,theta_rad,phi_rad,element,kernel_per_Pa_m3\n"
      << std::setprecision(17);
  for (std::size_t repeat = 0; repeat < static_cast<std::size_t>(repeats); ++repeat) {
    for (int objective = 0; objective < measurementCount; ++objective) {
      const auto &result = measured[repeat * measurementCount + objective];
      if (result.pointKernels.size() != points.size())
        throw std::runtime_error("DSpec result count differs from requested sample count");
      for (std::size_t i = 0; i < points.size(); ++i) {
        const auto &p = points[i];
        const auto &k = result.pointKernels[i];
        const std::array<double, 5> values{k.A, k.C, k.F, k.L, k.N};
        for (double value : values)
          if (!std::isfinite(value))
            throw std::runtime_error("DSpec timed kernel contains a nonfinite value");
        for (int m = 0; m < 5; ++m)
          csv << (objective == 0 ? "J_Z" : "J_N") << ',' << moduli[m] << ','
              << p.radiusMetres << ',' << p.thetaRadians << ',' << p.phiRadians << ','
              << p.radialElement
              << ',' << values[m] << '\n';
      }
    }
  }
  csv.close();
  fs::remove(paramPath);
  const auto receiverCoordinate = parameters.receivers().at(0);
  std::cout << std::setprecision(17)
            << R"({"solver":"DSpecM1D","api":"Windowed3DKernelEngine::evaluate","lmax":)"
            << lmax << R"(,"measurements":)" << measurementCount
            << R"(,"repeats":)" << repeats << R"(,"sample_count":)" << points.size()
            << R"(,"warmup_batches":1,"warmup_api_calls":)" << measurementCount
            << R"(,"setup_s":)"
            << std::chrono::duration<double>(setupEnd - processStart).count()
            << R"(,"warmup_batch_s":)" << warmupSeconds
            << R"(,"api_s":[)";
  for (std::size_t i = 0; i < apiSeconds.size(); ++i)
    std::cout << (i ? "," : "") << apiSeconds[i];
  const double setupSeconds = std::chrono::duration<double>(setupEnd - processStart).count();
  std::cout << R"(],"first_workload_api_s":)" << warmupSeconds
            << R"(,"first_model_to_arrays_s":)" << setupSeconds + warmupSeconds
            << R"(,"calculation_s":)"
            << setupSeconds + std::accumulate(apiSeconds.begin(), apiSeconds.end(), 0.0)
            << R"(,"end_to_end_s":)"
            << std::chrono::duration<double>(std::chrono::steady_clock::now() - processStart).count()
            << R"(,"peak_rss_bytes":)" << peakRssBytes()
            << R"(,"objective_scale_S":{"J_Z":)"
            << objectiveWeights(0).norm() * measured[0].processedSISpectrum.norm()
            ;
  if (measurementCount == 2)
    std::cout << R"(,"J_N":)"
              << objectiveWeights(1).norm() * measured[1].processedSISpectrum.norm();
  std::cout << R"(},"config":{"frequency_grid":{"nt":)" << frequency.nt()
            << R"(,"dt_s":)" << frequency.dt() * frequency.timeNorm()
            << R"(,"df_hz":)" << frequency.df() / frequency.timeNorm()
            << R"(,"epsilon_s_inv":)" << frequency.ep() / frequency.timeNorm()
            << R"(,"solver_bins_zero_based":[1,2,3],"band_hz":[0.0001,0.0008])"
            << R"(,"degrees_inclusive":[0,)" << lmax << R"(]},)"
            << R"("physics":{"gravity":"full","attenuation":false,"physical_G":)"
            << sem.gravityCoefficient() * base.GravitationalConstant() << R"(},)"
            << R"("source":{"geocentric_latitude_deg":)" << source.Latitude()
            << R"(,"longitude_deg":)" << source.Longitude()
            << R"(,"depth_m":)" << source.Depth() * 1000.0
            << R"(,"moment_tensor_Nm":[)" << source.Mrr() << ',' << source.Mrt() << ','
            << source.Mrp() << ',' << source.Mtt() << ',' << source.Mtp() << ',' << source.Mpp()
            << R"(]},"receiver":{"geocentric_latitude_deg":)"
            << receiverCoordinate.first << R"(,"longitude_deg":)" << receiverCoordinate.second
            << R"(,"depth_m":)" << parameters.receiver_depth() * 1000.0 << R"(},)"
            << R"("window":{"type":"Hann","first_sample":0,"last_sample":90,"nt":128},)"
            << R"("objective":{"names":[{"name":"J_Z","component":"Z","q_nonzero":{"bin_1":[0.7,0.2],"bin_2":[-0.25,0.45]}})"
            << (measurementCount == 2 ? R"(,{"name":"J_N","component":"N","q_nonzero":{"bin_1":[0.7,0.2],"bin_2":[-0.25,0.45]}})" : "")
            << R"(]}})"
            << "}\n";
  return 0;
} catch (const std::exception &error) {
  std::cerr << "DSpec kernel timing: " << error.what() << '\n';
  return 1;
}
