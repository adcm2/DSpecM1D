#include <YSpec/YSpec.hpp>
#include <YSpec/Kernels.hpp>
#include <YSpec/Model.hpp>
#include <YSpec/Source.hpp>
#include <YSpec/Receiver.hpp>
#include <YSpec/Units.hpp>
#include "Spectra.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <numbers>
#include <numeric>
#include <stdexcept>
#include <string>
#include <sys/resource.h>
#include <vector>

namespace fs = std::filesystem;
using namespace YSpec;
constexpr int nt = 128;
constexpr int nSamples = 91;
constexpr double dt = 20.0;
constexpr double duration = 1800.0;
constexpr std::array<const char *, 5> moduli{"A", "C", "F", "L", "N"};

struct SamplePoint {
  double radius;
  double theta;
  double phi;
};

struct SharedSamples {
  std::vector<double> radii;
  std::vector<KernelPlace> places;
  std::vector<SamplePoint> points;
};

double geographicLatitude(double geocentricDegrees) {
  const double angle = geocentricDegrees * std::numbers::pi / 180.0;
  return std::atan(std::tan(angle) / 0.993277) * 180.0 / std::numbers::pi;
}

double geocentricLatitude(double geographicDegrees) {
  const double angle = geographicDegrees * std::numbers::pi / 180.0;
  return std::atan(0.993277 * std::tan(angle)) * 180.0 / std::numbers::pi;
}

SharedSamples readSamples(const fs::path &path, std::size_t anglesPerRadius) {
  std::ifstream file(path);
  if (!file)
    throw std::runtime_error("cannot read shared tutorial samples CSV");
  std::string line;
  std::getline(file, line);
  SharedSamples result;
  while (std::getline(file, line)) {
    std::istringstream row(line);
    std::array<std::string, 11> values;
    for (auto &value : values)
      if (!std::getline(row, value, ','))
        throw std::runtime_error("invalid shared tutorial sample row");
    if (values[0] != "J_Z" || values[1] != "A")
      continue;
    SamplePoint point{std::stod(values[2]), std::stod(values[3]), std::stod(values[4])};
    result.points.push_back(point);
  }
  if (result.points.empty() || result.points.size() % anglesPerRadius != 0)
    throw std::runtime_error("A rows do not contain complete radial blocks for this lmax");
  for (std::size_t i = 0; i < result.points.size(); ++i) {
    const auto &point = result.points[i];
    if (i % anglesPerRadius == 0)
      result.radii.push_back(point.radius);
    else if (point.radius != result.radii.back())
      throw std::runtime_error("shared radial block changes radius before all angles are read");
    if (i < anglesPerRadius)
      result.places.push_back({point.theta, point.phi});
    else {
      const auto &expected = result.points[i % anglesPerRadius];
      if (point.theta != expected.theta || point.phi != expected.phi)
        throw std::runtime_error("shared angular grid differs between radial blocks");
    }
  }
  return result;
}

Measurement makeMeasurement(const Receiver &receiver, int component) {
  std::vector<double> samples(nSamples);
  for (int n = 0; n < nSamples; ++n) {
    const double window = 0.5 - 0.5 * std::cos(2.0 * std::numbers::pi * n / (nSamples - 1));
    const auto term = [n](int k, double real, double imag) {
      const double angle = 2.0 * std::numbers::pi * k * n / nt;
      return real * std::cos(angle) - imag * std::sin(angle);
    };
    samples[n] = window * (term(1, 0.7, 0.2) + term(2, -0.25, 0.45));
  }
  return {receiver, component == 0 ? Direction{.up = 1} : Direction{.north = 1},
          std::move(samples)};
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
        "usage: yspec_kernel_timing MODEL SAMPLES OUTPUT LMAX MEASUREMENTS REPEATS");
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

  const Model model = Model::ReadDeck(modelPath);
  Calculation calculation;
  calculation.gravity = Gravity::Full;
  calculation.attenuation = false;
  calculation.quantity = Quantity::Displacement;
  calculation.degrees = {0, lmax};
  calculation.band = {0.1e-3, 0.8e-3};
  calculation.sampling = {duration, dt};
  calculation.solver.tolerance = 1e-9;
  calculation.solver.startDecay = 11;
  calculation.solver.stepFactor = 0.25;
  const Source source{{geographicLatitude(10.0), 20.0}, 50000.0,
                      {1e19, 2e18, -3e18, -5e18, 4e18, -5e18}, {}};
  const Receiver receiver{{geographicLatitude(45.0), 90.0}, 3000.0, {}};
  const std::size_t nTheta = std::max(8, 2 * lmax + 2);
  const std::size_t nPhi = std::max(24, 4 * lmax + 4);
  const auto shared = readSamples(samplesPath, nTheta * nPhi);
  KernelSampling sampling;
  for (double radius : shared.radii) {
    const int layer = model.LayerIndex(radius / model.GetScales().Length(), Side::Above);
    sampling.nodes.push_back({layer, radius});
  }
  sampling.places = shared.places;
  std::array<KernelSampling, 1> samplings{std::move(sampling)};
  std::vector<Measurement> measurements;
  measurements.push_back(makeMeasurement(receiver, 0));
  if (measurementCount == 2)
    measurements.push_back(makeMeasurement(receiver, 1));
  const auto grid = MakeGrid(calculation, model.GetScales());
  if (grid.Samples() != nt || grid.First() != 1 || grid.Last() != 3 ||
      std::abs(grid.Damping() / model.GetScales().Time() - 1.0 / 360.0) > 1e-14)
    throw std::runtime_error("resolved YSpec grid differs from bins 1..3 / epsilon 1/360");
  const auto setupEnd = std::chrono::steady_clock::now();
  const auto evaluate = [&]() {
    return ComputeSampledKernels(model, calculation, source, measurements,
                                 ParameterSet::Love, Scaling::Absolute, samplings);
  };

  const auto warmupStart = std::chrono::steady_clock::now();
  (void)evaluate();
  const double warmupSeconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - warmupStart).count();
  std::vector<double> apiSeconds;
  std::vector<decltype(evaluate())> measured;
  for (int repeat = 0; repeat < repeats; ++repeat) {
    const auto start = std::chrono::steady_clock::now();
    measured.push_back(evaluate());
    apiSeconds.push_back(std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count());
  }
  std::ofstream csv(output / "kernels.csv");
  if (!csv)
    throw std::runtime_error("cannot write actual timed kernels");
  csv << "objective,modulus,radius_m,theta_rad,phi_rad,kernel_per_Pa_m3\n"
      << std::setprecision(17);
  for (const auto &result : measured) {
    if (result.size() != static_cast<std::size_t>(measurementCount) || result[0].size() != 1)
      throw std::runtime_error("YSpec sampled result dimensions differ from requested measurements");
    for (int objective = 0; objective < measurementCount; ++objective) {
      const auto &sampled = result[objective][0];
      for (int m = 0; m < 5; ++m) {
        const auto found = std::find(sampled.Names().begin(), sampled.Names().end(), moduli[m]);
        if (found == sampled.Names().end())
          throw std::runtime_error("YSpec result is missing a Love modulus");
        const auto index = static_cast<std::size_t>(found - sampled.Names().begin());
        for (std::size_t r = 0; r < shared.radii.size(); ++r) {
          if (sampled.Values(index, r).size() != shared.places.size())
            throw std::runtime_error("YSpec angular count differs from shared grid");
          for (std::size_t a = 0; a < shared.places.size(); ++a) {
            const double value = sampled.Values(index, r)[a];
            if (!std::isfinite(value))
              throw std::runtime_error("YSpec timed kernel contains a nonfinite value");
            const auto &point = shared.points[r * shared.places.size() + a];
            csv << (objective == 0 ? "J_Z" : "J_N") << ',' << moduli[m] << ','
                << point.radius << ',' << point.theta << ',' << point.phi << ',' << value << '\n';
          }
        }
      }
    }
  }
  csv.close();
  std::cout << std::setprecision(17)
            << R"({"solver":"YSpecpp","api":"ComputeSampledKernels","lmax":)"
            << lmax << R"(,"measurements":)" << measurementCount
            << R"(,"repeats":)" << repeats << R"(,"point_count":)" << shared.points.size()
            << R"(,"warmup_batches":1,"warmup_api_calls":1,"setup_s":)"
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
            << R"(,"config":{"frequency_grid":{"nt":)" << grid.Samples()
            << R"(,"dt_s":)" << grid.TimeStep() * model.GetScales().Time()
            << R"(,"df_hz":)" << grid.Spacing() / model.GetScales().Time()
            << R"(,"epsilon_s_inv":)" << grid.Damping() / model.GetScales().Time()
            << R"(,"solver_bins_zero_based":[1,2,3],"band_hz":[)"
            << calculation.band.min << ',' << calculation.band.max << R"(],"degrees_inclusive":[0,)"
            << lmax << R"(]},)"
            << R"("physics":{"gravity":"full","attenuation":false,"physical_G":)"
            << GravitationalConstant << R"(},)"
            << R"("source":{"geocentric_latitude_deg":)"
            << geocentricLatitude(source.location.latitude)
            << R"(,"longitude_deg":)" << source.location.longitude
            << R"(,"depth_m":)" << source.depth
            << R"(,"moment_tensor_Nm":[)" << source.moment.rr << ',' << source.moment.rt << ','
            << source.moment.rp << ',' << source.moment.tt << ',' << source.moment.tp << ','
            << source.moment.pp << R"(]},)"
            << R"("receiver":{"geocentric_latitude_deg":)"
            << geocentricLatitude(receiver.location.latitude)
            << R"(,"longitude_deg":)" << receiver.location.longitude
            << R"(,"depth_m":)" << receiver.depth << R"(},)"
            << R"("window":{"type":"Hann","first_sample":0,"last_sample":90,"nt":128},)"
            << R"("objective":{"names":[{"name":"J_Z","component":"Z")"
            << R"(,"q_nonzero":{"bin_1":[0.7,0.2],"bin_2":[-0.25,0.45]}})"
            << (measurementCount == 2 ? R"(,{"name":"J_N","component":"N","q_nonzero":{"bin_1":[0.7,0.2],"bin_2":[-0.25,0.45]}})" : "")
            << R"(]},)"
            << R"("angular_node_count":)" << shared.places.size()
            << R"(,"radial_node_count":)" << shared.radii.size() << "}}\n";
  return 0;
} catch (const std::exception &error) {
  std::cerr << "YSpec kernel timing: " << error.what() << '\n';
  return 1;
}
