#include <YSpec/YSpec.hpp>
#include <YSpec/Calculation.hpp>
#include <YSpec/Kernels.hpp>
#include <YSpec/Model.hpp>
#include <YSpec/Receiver.hpp>
#include <YSpec/Source.hpp>
#include <YSpec/Units.hpp>

#include "Coefficients.hpp"
#include "PhysicalGrid.hpp"
#include "Spectra.hpp"
#include "TimeSeries.hpp"
#include "profile.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/resource.h>
#include <vector>
#include <omp.h>
#include <GSHTrans/GaussLegendreGrid.hpp>
#include <FFTWpp/Ranges>

namespace fs = std::filesystem;
namespace YSpec { bool PhysicalIncludeDC = false; }
namespace YSpec { unsigned PhysicalFftPadding = 1; }

namespace {
using Clock = std::chrono::steady_clock;
constexpr double kPi = std::numbers::pi;
constexpr double kDt = 5.0;
constexpr double kDuration = 3600.0;
constexpr double kFMax = 0.010;
constexpr std::array<const char*, 5> kModuli{"A", "C", "F", "L", "N"};

struct Options {
  std::string mode;
  fs::path model;
  fs::path output;
  fs::path baseline;
  fs::path samples;
  int lmax = 0;
  double tolerance = 1e-9;
  double stepFactor = 0.25;
  double kernelMemoryBudgetMiB = 6144;
  unsigned fftPadding = 1;
  int nodesPerElement = 0;
  bool includeDC = false;
};

struct SharedSamples {
  std::vector<double> radii;
  std::vector<int> layers;
  std::vector<YSpec::Side> sides;
  std::vector<YSpec::KernelPlace> places;
  std::size_t nTheta = 0;
  std::size_t nPhi = 0;
  std::size_t estimatedPeakBytes = 0;
};
struct NativeNode { int element, node; double radius, width; int fluid; };

double seconds(Clock::time_point start) {
  return std::chrono::duration<double>(Clock::now() - start).count();
}

void require(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}

double geocentricToGeographic(double degrees) {
  const double x = degrees * kPi / 180.0;
  return std::atan(std::tan(x) / 0.993277) * 180.0 / kPi;
}

double hann(double t) {
  if (t < 1000.0 || t > 1200.0) return 0.0;
  return 0.5 * (1.0 - std::cos(2.0 * kPi * (t - 1000.0) / 200.0));
}

double bandWeight(double f) {
  if (f < 0.008) return 1.0;
  if (f < 0.010) return 0.5 * (1.0 + std::cos(kPi * (f - 0.008) / 0.002));
  return 0.0;
}

std::string jsonString(const std::string& value) {
  std::ostringstream out;
  out << '"';
  for (char c : value) {
    if (c == '"' || c == '\\') out << '\\';
    if (c == '\n') out << "\\n";
    else if (c != '\r') out << c;
  }
  out << '"';
  return out.str();
}

Options parse(int argc, char** argv) {
  if (argc < 5) throw std::runtime_error(
      "usage: physical_yspec waveform MODEL OUTDIR LMAX [TOLERANCE STEPFAC] [--include-dc] [--fft-padding N]\n"
      "   or: physical_yspec kernel MODEL BASELINE_TRACE RADIAL_NODES.csv OUTDIR LMAX [TOLERANCE STEPFAC] [--include-dc] [--fft-padding N] [--nodes-per-element NN] [--kernel-memory-budget-mib N]");
  Options o;
  o.mode = argv[1];
  std::vector<std::string> positional;
  for (int i = 2; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--include-dc") { o.includeDC = true; continue; }
    if (arg == "--kernel-memory-budget-mib" || arg == "--max-output-mib") {
      if (++i >= argc) throw std::runtime_error("--kernel-memory-budget-mib needs a value");
      o.kernelMemoryBudgetMiB = std::stod(argv[i]);
      continue;
    }
    if (arg == "--fft-padding") {
      if (++i >= argc) throw std::runtime_error("--fft-padding needs a value");
      o.fftPadding = static_cast<unsigned>(std::stoul(argv[i]));
      continue;
    }
    if (arg == "--nodes-per-element") {
      if (++i >= argc) throw std::runtime_error("--nodes-per-element needs a value");
      o.nodesPerElement = std::stoi(argv[i]);
      continue;
    }
    positional.push_back(arg);
  }
  if (o.mode == "waveform") {
    if (positional.size() < 3 || positional.size() > 5)
      throw std::runtime_error("waveform needs MODEL OUTDIR LMAX and optional TOLERANCE STEPFAC");
    o.model = positional[0]; o.output = positional[1]; o.lmax = std::stoi(positional[2]);
    if (positional.size() > 3) o.tolerance = std::stod(positional[3]);
    if (positional.size() > 4) o.stepFactor = std::stod(positional[4]);
  } else if (o.mode == "kernel") {
    if (positional.size() < 5 || positional.size() > 7)
      throw std::runtime_error("kernel needs MODEL BASELINE_TRACE SAMPLES.csv OUTDIR LMAX and optional TOLERANCE STEPFAC");
    o.model = positional[0]; o.baseline = positional[1]; o.samples = positional[2];
    o.output = positional[3]; o.lmax = std::stoi(positional[4]);
    if (positional.size() > 5) o.tolerance = std::stod(positional[5]);
    if (positional.size() > 6) o.stepFactor = std::stod(positional[6]);
  } else throw std::runtime_error("mode must be waveform or kernel");
  require(o.lmax >= 1 && o.lmax <= 1536,
          "lmax must be in 1..1536 (diagnostic adapter ceiling; kernel memory guard is separate)");
  require(o.tolerance > 0 && std::isfinite(o.tolerance), "tolerance must be positive and finite");
  require(o.stepFactor > 0 && std::isfinite(o.stepFactor), "stepFactor must be positive and finite");
  require(o.kernelMemoryBudgetMiB > 0 && std::isfinite(o.kernelMemoryBudgetMiB),
          "kernel memory budget must be positive and finite");
  require(o.fftPadding == 1 || o.fftPadding == 2 || o.fftPadding == 4 || o.fftPadding == 8,
          "fft-padding must be 1, 2, 4, or 8");
  require(o.nodesPerElement == 0 || (o.nodesPerElement >= 2 && o.nodesPerElement <= 20),
          "nodes-per-element must be in 2..20 when supplied");
  return o;
}

YSpec::Calculation calculation(const Options& o) {
  YSpec::Calculation c;
  c.gravity = YSpec::Gravity::Full;
  c.attenuation = false;
  c.quantity = YSpec::Quantity::Displacement;
  c.degrees = {0, o.lmax};
  c.band = {0.0, kFMax};
  c.sampling = {kDuration, kDt};
  c.taper = YSpec::Taper{0.0, 0.0, 0.008, 0.010};
  c.solver.tolerance = o.tolerance;
  c.solver.stepFactor = o.stepFactor;
  return c;
}

YSpec::Source source() {
  return {{geocentricToGeographic(10.0), 20.0}, 50000.0,
          {1e19, 2e18, -3e18, -5e18, 4e18, -5e18}, {}};
}
YSpec::Receiver receiver() {
  return {{geocentricToGeographic(45.0), 75.0}, 3000.0, {}};
}

struct TraceData {
  std::vector<double> time, vertical, window;
  double energy = 0;
  double unwindowedEnergy = 0;
  double maxAbs = 0;
};

TraceData makeTrace(const YSpec::Model& model, const YSpec::Calculation& c,
                    const YSpec::FrequencyGrid& grid,
                    const YSpec::Spectra& spectra) {
  const auto taper = YSpec::MakeTaper(c, model.GetScales());
  const auto series = YSpec::InverseTransform(grid, taper, spectra.Values(), 3);
  const auto sampleCount = YSpec::SampleCount(grid, c.sampling);
  const double unit = YSpec::Unit(model.GetScales(), c.quantity);
  TraceData out;
  out.time.resize(sampleCount); out.vertical.resize(sampleCount); out.window.resize(sampleCount);
  for (std::size_t n = 0; n < sampleCount; ++n) {
    const double t = n * c.sampling.timeStep;
    const double u = series[n] * unit;
    const double w = hann(t);
    out.time[n] = t; out.vertical[n] = u; out.window[n] = w;
    out.energy += 0.5 * c.sampling.timeStep * w * u * u;
    if (t >= 1000.0 && t <= 1200.0)
      out.unwindowedEnergy += 0.5 * c.sampling.timeStep * u * u;
    out.maxAbs = std::max(out.maxAbs, std::abs(u));
  }
  require(out.energy > 0 && std::isfinite(out.energy), "baseline energy window is empty or non-finite");
  return out;
}

void writeTrace(const fs::path& path, const TraceData& trace) {
  std::ofstream out(path);
  require(bool(out), "cannot write trace.csv");
  out << "time_s,vertical_m,window_weight\n" << std::setprecision(17);
  for (std::size_t i = 0; i < trace.time.size(); ++i)
    out << trace.time[i] << ',' << trace.vertical[i] << ',' << trace.window[i] << '\n';
}

void writeSpectrum(const fs::path& path, const YSpec::Model& model,
                   const YSpec::FrequencyGrid& grid,
                   const YSpec::Spectra& spectra,
                   double& endpointImagMax) {
  std::ofstream out(path);
  require(bool(out), "cannot write spectrum.csv");
  out << "bin,frequency_hz,real,imag,filtered_real,filtered_imag\n" << std::setprecision(17);
  const auto& scales = model.GetScales();
  const double toSI = scales.Length() * scales.Time();
  const std::size_t last = static_cast<std::size_t>(std::floor(kFMax / (grid.Spacing() / scales.Time())));
  for (std::size_t bin = 0; bin <= last; ++bin) {
    std::complex<double> value{};
    if (bin >= grid.First() && bin <= grid.Last()) {
      const auto k = bin - grid.First();
      value = spectra(0, 0, YSpec::Vertical, k) * toSI;
    }
    const auto filtered = value * bandWeight(bin * grid.Spacing() / scales.Time());
    // The inverse transform projects endpoint imaginary parts to reals. Retain
    // the unprojected raw and filtered values here as diagnostics.
    const bool endpoint = bin == 0 || bin == grid.Samples() / 2;
    if (endpoint) endpointImagMax = std::max(endpointImagMax, std::abs(value.imag()));
    out << bin << ',' << bin * grid.Spacing() / scales.Time() << ','
        << value.real() << ',' << value.imag() << ','
        << filtered.real() << ',' << filtered.imag() << '\n';
  }
}

double number(const std::string& s, const char* field) {
  std::size_t used = 0;
  const double x = std::stod(s, &used);
  require(used == s.size() && std::isfinite(x), std::string("invalid ") + field);
  return x;
}

TraceData readTrace(const fs::path& path) {
  std::ifstream in(path);
  require(bool(in), "cannot open baseline trace");
  std::string line;
  require(bool(std::getline(in, line)) && line == "time_s,vertical_m,window_weight",
          "baseline trace header must be time_s,vertical_m,window_weight");
  TraceData trace;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    std::istringstream row(line);
    std::array<std::string, 3> fields;
    for (auto& f : fields) require(bool(std::getline(row, f, ',')), "malformed trace row");
    const double t = number(fields[0], "trace time");
    const double u = number(fields[1], "trace displacement");
    const double fileWindow = number(fields[2], "trace window");
    const double w = hann(t);
    require(std::abs(fileWindow - w) <= 1e-12, "baseline trace window differs from requested Hann window");
    trace.time.push_back(t); trace.vertical.push_back(u); trace.window.push_back(w);
    trace.energy += 0.5 * kDt * w * u * u;
    if (t >= 1000.0 && t <= 1200.0) trace.unwindowedEnergy += 0.5 * kDt * u * u;
    trace.maxAbs = std::max(trace.maxAbs, std::abs(u));
  }
  require(trace.time.size() == 721, "baseline trace must contain 721 samples at 5 s through 3600 s");
  for (std::size_t i = 0; i < trace.time.size(); ++i)
    require(std::abs(trace.time[i] - i * kDt) < 1e-10, "baseline trace has a mismatched time grid");
  require(trace.energy > 0 && std::isfinite(trace.energy), "baseline trace energy is zero or non-finite");
  return trace;
}

SharedSamples readSamples(const fs::path& path, const YSpec::Model& model,
                         int lmax, double maxMiB, int requestedNodesPerElement) {
  std::ifstream in(path);
  require(bool(in), "cannot open native radial_nodes.csv");
  std::string line;
  require(bool(std::getline(in, line)) &&
              line == "element,node,radius_m,element_width_m,fluid",
          "radial input must use element,node,radius_m,element_width_m,fluid header");
  std::vector<NativeNode> nodes;
  std::map<int, int> maxNode;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    std::istringstream row(line);
    std::array<std::string, 5> f;
    for (auto& v : f) require(bool(std::getline(row, v, ',')), "malformed radial_nodes.csv row");
    std::size_t usedElement = 0, usedNode = 0, usedFluid = 0;
    const int element = std::stoi(f[0], &usedElement);
    const int node = std::stoi(f[1], &usedNode);
    const int fluid = std::stoi(f[4], &usedFluid);
    require(usedElement == f[0].size() && usedNode == f[1].size() && usedFluid == f[4].size(),
            "element, node, and fluid columns must be integers");
    const double radius = number(f[2], "native radius");
    const double width = number(f[3], "native element width");
    require(element >= 0 && node >= 0 && radius > 0 && width > 0 && (fluid == 0 || fluid == 1),
            "native radial node has invalid element, node, radius, width, or fluid value");
    maxNode[element] = std::max(maxNode[element], node);
    nodes.push_back({element, node, radius, width, fluid});
  }
  require(!nodes.empty(), "radial_nodes.csv contains no radial nodes");
  int nodesPerElement = requestedNodesPerElement;
  if (nodesPerElement == 0) {
    nodesPerElement = maxNode.begin()->second + 1;
    for (const auto& [element, maximum] : maxNode) {
      (void)element;
      require(maximum + 1 == nodesPerElement,
              "native element node counts vary; provide --nodes-per-element from manifest metadata");
    }
  }
  for (const auto& node : nodes) {
    require(node.node < nodesPerElement,
            "radial_nodes.csv node index exceeds nodes-per-element");
    if (requestedNodesPerElement == 0)
      require(maxNode[node.element] == nodesPerElement - 1,
              "inferred element node extent is incomplete; provide --nodes-per-element for a subset");
  }

  SharedSamples shared;
  shared.radii.reserve(nodes.size());
  shared.layers.reserve(nodes.size());
  shared.sides.reserve(nodes.size());
  const double length = model.GetScales().Length();
  for (const auto& node : nodes) {
    const double radius = node.radius;
    const auto below = model.LayerIndex(radius / length, YSpec::Side::Below);
    const auto above = model.LayerIndex(radius / length, YSpec::Side::Above);
    YSpec::Side side = YSpec::Side::Above;
    if (below != above) {
      require(node.node == 0 || node.node == nodesPerElement - 1,
              "native radial node lies on a model interface away from an element endpoint");
      side = node.node == 0 ? YSpec::Side::Above : YSpec::Side::Below;
    }
    shared.radii.push_back(radius);
    shared.sides.push_back(side);
    shared.layers.push_back(model.LayerIndex(radius / length, side));
  }

  using Grid = GSHTrans::GaussLegendreGrid<double, GSHTrans::All, GSHTrans::All>;
  const auto angular = Grid::ForBand(static_cast<std::ptrdiff_t>(lmax), 2, 2.0,
                                     FFTWpp::Estimate);
  const auto thetas = angular.CoLatitudes();
  const auto phis = angular.Longitudes();
  shared.nTheta = thetas.size();
  shared.nPhi = phis.size();
  const long double angleCount = static_cast<long double>(shared.nTheta) * shared.nPhi;
  const long double radialCount = static_cast<long double>(shared.radii.size());
  // Conservative preflight for YSpec's simultaneous 8-array radial
  // accumulation and 5-array returned result, angular coordinates/geometry,
  // one 64-place temporary block, per-thread strain tables, and Legendre/Wigner
  // workspaces. This is a budget estimate for the complete API call, not just
  // the output file. Waveform mode does not use this kernel-only budget.
  const long double pointBytes = angleCount * radialCount * 24.0L * sizeof(double);
  const long double radialDegreeBytes = radialCount * (lmax + 1.0L) * 2048.0L;
  const long double angularWorkspaceBytes = (lmax + 1.0L) * 1024.0L * 1024.0L;
  constexpr long double fixedWorkspaceBytes = 64.0L * 1024 * 1024;
  const long double estimate = pointBytes + radialDegreeBytes +
                               angularWorkspaceBytes + fixedWorkspaceBytes;
  require(estimate <= static_cast<long double>(maxMiB) * 1024 * 1024,
          "conservative sampled-kernel peak-memory estimate exceeds --kernel-memory-budget-mib before angle/result allocation");
  shared.estimatedPeakBytes = static_cast<std::size_t>(std::ceil(estimate));
  shared.places.reserve(shared.nTheta * shared.nPhi);
  for (const double theta : thetas)
    for (const double phi : phis) shared.places.push_back({theta, phi});
  return shared;
}

std::string profileJson() {
  std::ostringstream out;
  out << std::setprecision(17) << '{';
  bool first = true;
  for (const auto& [name, entry] : mathprofile::entries) {
    out << (first ? "" : ",") << jsonString(name) << ":{\"inclusive_s\":"
        << entry.inclusive * 1e-9 << ",\"exclusive_s\":" << entry.exclusive * 1e-9
        << ",\"calls\":" << entry.calls << '}';
    first = false;
  }
  out << '}';
  return out.str();
}

long peakRssBytes() {
  rusage usage{};
  require(getrusage(RUSAGE_SELF, &usage) == 0, "getrusage failed");
#if defined(__APPLE__)
  return usage.ru_maxrss;
#else
  return usage.ru_maxrss * 1024L;
#endif
}

void prepareOutput(const fs::path& path) {
  require(!fs::exists(path), "refusing to overwrite an existing output directory");
  fs::create_directories(path);
}

void writeResult(const fs::path& path, const Options& o, const YSpec::Model& model,
                 const YSpec::FrequencyGrid& grid, const TraceData& trace,
                 double setupS, double apiS, double ioS, double endpointImag,
                 std::size_t pointCount, std::size_t radialCount,
                 std::size_t nTheta, std::size_t nPhi,
                 std::size_t kernelEstimatedPeakBytes, long rssBytes) {
  std::ofstream out(path);
  require(bool(out), "cannot write result.json");
  out << std::setprecision(17)
      << "{\n  \"solver\":\"YSpecpp\",\n  \"mode\":" << jsonString(o.mode)
      << ",\n  \"model\":" << jsonString(fs::absolute(o.model).string())
      << ",\n  \"lmax\":" << o.lmax << ",\n  \"degree_min\":0,\n  \"tolerance\":" << o.tolerance
      << ",\n  \"step_factor\":" << o.stepFactor
      << ",\n  \"gravity\":\"full\",\n  \"attenuation\":false,\n  \"quantity\":\"vertical_displacement\",\n"
      << "  \"source_geocentric_lat_lon_deg\":[10,20],\n  \"source_y_spec_geographic_lat_lon_deg\":["
      << geocentricToGeographic(10) << ",20],\n  \"source_depth_m\":50000,\n"
      << "  \"receiver_geocentric_lat_lon_deg\":[45,75],\n  \"receiver_y_spec_geographic_lat_lon_deg\":["
      << geocentricToGeographic(45) << ",75],\n  \"receiver_depth_m\":3000,\n"
      << "  \"dt_s\":" << kDt << ",\n  \"duration_s\":" << kDuration
      << ",\n  \"fft_samples\":" << grid.Samples() << ",\n  \"fft_period_s\":" << grid.Samples() * kDt
      << ",\n  \"fft_padding_factor\":" << o.fftPadding
      << ",\n  \"frequency_spacing_hz\":" << grid.Spacing() / model.GetScales().Time()
      << ",\n  \"epsilon_per_s\":" << grid.Damping() / model.GetScales().Time()
      << ",\n  \"solver_dc_included\":" << (o.includeDC ? "true" : "false")
      << ",\n  \"spectrum_endpoint_imaginary_max_before_projection_m_s\":" << endpointImag
      << ",\n  \"endpoint_real_projection\":true,\n  \"spectrum_csv_real_imag\":\"raw SI spectrum before 8-10 mHz filter\",\n"
      << "  \"spectrum_csv_filtered_real_imag\":\"SI spectrum after 8-10 mHz filter\",\n"
      << "  \"frequency_filter\":\"unity below 8 mHz; cosine 8-10 mHz; zero at and above 10 mHz\",\n"
      << "  \"energy_window_s\":[1000,1200],\n  \"windowed_energy_m2_s\":" << trace.energy
      << ",\n  \"linear_measurement_q_dot_p\":" << 2.0 * trace.energy
      << ",\n  \"unwindowed_energy_in_window_m2_s\":" << trace.unwindowedEnergy
      << ",\n  \"peak_vertical_m\":" << trace.maxAbs
      << ",\n  \"sampling_scope\":\"provided radial nodes and angular points only; no full-global claim\",\n"
      << "  \"kernel_point_count\":" << pointCount << ",\n  \"kernel_radial_blocks\":" << radialCount
      << ",\n  \"kernel_angles_per_radius\":" << (radialCount ? pointCount / radialCount : 0)
      << ",\n  \"kernel_n_theta\":" << nTheta
      << ",\n  \"kernel_n_phi\":" << nPhi
      << ",\n  \"kernel_angular_grid\":\"GSHTrans GaussLegendreGrid ForBand(lmax,2,2.0), theta-major then longitude\""
      << ",\n  \"kernel_interface_side_rule\":\"native q=0 maps Side::Above; q=NN-1 maps Side::Below when at a model interface\""
      << ",\n  \"kernel_memory_estimate_method\":\"24 doubles per radial-angular sample plus 2048 bytes per radial-degree pair, 1 MiB per degree of angular workspace, and 64 MiB fixed reserve\""
      << ",\n  \"kernel_memory_budget_mib\":" << o.kernelMemoryBudgetMiB
      << ",\n  \"kernel_estimated_peak_bytes\":" << kernelEstimatedPeakBytes
      << ",\n  \"kernel_output_file\":" << jsonString(o.mode == "kernel" ? "kernels.bin" : "")
      << ",\n  \"kernel_binary_layout\":\"float64 native endian; parameter-major A,C,F,L,N, then radial block, then angle\",\n"
      << "  \"kernel_baseline_trace\":" << jsonString(o.mode == "kernel" ? fs::absolute(o.baseline).string() : "")
      << ",\n  \"kernel_sampling_input\":" << jsonString(o.mode == "kernel" ? fs::absolute(o.samples).string() : "") << ",\n"
      << "  \"nodes_per_element_requested\":" << o.nodesPerElement << ",\n"
      << "  \"peak_rss_bytes\":" << rssBytes
      << ",\n  \"timings_s\":{\"setup\":" << setupS << ",\"numerical_api\":" << apiS
      << ",\"output_io\":" << ioS << "},\n  \"profile_stages\":" << profileJson() << "\n}\n";
}

void run(const Options& o) {
  omp_set_dynamic(0);
  omp_set_num_threads(1);
  require(omp_get_max_threads() == 1, "physical adapter requires exactly one OpenMP thread");
  const auto totalStart = Clock::now();
  const std::string profileName = "YSpecpp.physical_" + o.mode;
  mathprofile::Root profileRoot(profileName.c_str(), 0);
  auto physicalCase = std::make_unique<mathprofile::Scope>("physical_case_total");
  YSpec::PhysicalIncludeDC = o.includeDC;
  YSpec::PhysicalFftPadding = o.fftPadding;
  const auto setupStart = Clock::now();
  auto setupProfile = std::make_unique<mathprofile::Scope>("setup.model_calculation_source_receiver_grid");
  const auto model = YSpec::Model::ReadDeck(fs::absolute(o.model));
  const auto calc = calculation(o);
  const auto src = source();
  const auto rec = receiver();
  const auto grid = YSpec::MakeGrid(calc, model.GetScales());
  require(grid.Samples() == 1024 * o.fftPadding &&
          std::abs(grid.Spacing() / model.GetScales().Time() - 1.0 / (5120.0 * o.fftPadding)) < 1e-15,
          "resolved physical FFT grid differs from requested fixed-damping padding");
  require(std::abs(grid.Damping() / model.GetScales().Time() - 5.0 / 3600.0) < 1e-14,
          "resolved physical damping differs from eps=5/tout");
  setupProfile.reset();
  const double setupS = seconds(setupStart);
  double apiS = 0, ioS = 0, endpointImag = 0;
  TraceData trace;
  std::size_t pointCount = 0, radialCount = 0, nTheta = 0, nPhi = 0;
  std::size_t kernelEstimatedPeakBytes = 0;
  prepareOutput(o.output);

  if (o.mode == "waveform") {
    const auto apiStart = Clock::now();
    YSpec::Spectra spectra(1, 1, 3, grid.Count());
    {
      mathprofile::Scope api("api.SeismogramSpectra_and_trace_reconstruction");
      spectra = YSpec::SeismogramSpectra(model, calc, std::span{&src, 1}, std::span{&rec, 1});
      mathprofile::Scope reconstruction("reconstruction.inverse_transform_and_energy");
      trace = makeTrace(model, calc, grid, spectra);
    }
    apiS = seconds(apiStart);
    const auto ioStart = Clock::now();
    {
      mathprofile::Scope io("output.spectrum_and_trace_csv");
      writeSpectrum(o.output / "spectrum.csv", model, grid, spectra, endpointImag);
      writeTrace(o.output / "trace.csv", trace);
    }
    ioS = seconds(ioStart);
  } else {
    auto measurementSetup = std::make_unique<mathprofile::Scope>("setup.frozen_measurement_and_sample_grid");
    const auto samples = readSamples(fs::absolute(o.samples), model, o.lmax,
                                     o.kernelMemoryBudgetMiB, o.nodesPerElement);
    pointCount = samples.radii.size() * samples.places.size(); radialCount = samples.radii.size();
    nTheta = samples.nTheta; nPhi = samples.nPhi;
    kernelEstimatedPeakBytes = samples.estimatedPeakBytes;
    trace = readTrace(fs::absolute(o.baseline));
    YSpec::Measurement measurement{rec, {.up = 1.0}, {}};
    measurement.adjointSource.resize(trace.vertical.size());
    for (std::size_t i = 0; i < trace.vertical.size(); ++i)
      measurement.adjointSource[i] = trace.window[i] * trace.vertical[i];
    YSpec::KernelSampling sampling;
    for (std::size_t i = 0; i < samples.radii.size(); ++i) {
      const double radius = samples.radii[i];
      sampling.nodes.push_back({samples.layers[i], radius});
    }
    sampling.places = samples.places;
    std::array<YSpec::KernelSampling, 1> samplings{std::move(sampling)};
    std::array<YSpec::Measurement, 1> measurements{std::move(measurement)};
    measurementSetup.reset();
    const auto apiStart = Clock::now();
    std::vector<std::vector<YSpec::SampledKernels>> results;
    {
      mathprofile::Scope api("api.ComputeSampledKernels");
      results = YSpec::ComputeSampledKernels(model, calc, src, measurements,
          YSpec::ParameterSet::Love, YSpec::Scaling::Absolute, samplings);
    }
    apiS = seconds(apiStart);
    const auto ioStart = Clock::now();
    mathprofile::Scope io("output.sampled_kernel_binary");
    std::ofstream binary(o.output / "kernels.bin", std::ios::binary);
    require(bool(binary), "cannot write kernels.bin");
    require(results.size() == 1 && results[0].size() == 1, "YSpec sampled kernel result dimensions differ");
    const auto& result = results[0][0];
    for (int m = 0; m < 5; ++m) {
      const auto found = std::find(result.Names().begin(), result.Names().end(), kModuli[m]);
      require(found != result.Names().end(), "YSpec result is missing a Love modulus");
      const auto p = static_cast<std::size_t>(found - result.Names().begin());
      for (std::size_t r = 0; r < result.Nodes().size(); ++r) {
        const auto values = result.Values(p, r);
        for (std::size_t a = 0; a < values.size(); ++a) {
          require(std::isfinite(values[a]), "sampled kernel contains a non-finite value");
          binary.write(reinterpret_cast<const char*>(&values[a]), sizeof(double));
        }
      }
    }
    binary.close();
    require(bool(binary), "failed while writing kernels.bin");
    ioS = seconds(ioStart);
  }
  physicalCase.reset();
  writeResult(o.output / "result.json", o, model, grid, trace, setupS, apiS, ioS,
              endpointImag, pointCount, radialCount, nTheta, nPhi,
              kernelEstimatedPeakBytes, peakRssBytes());
  std::cout << std::setprecision(17) << "YSpec " << o.mode << " complete: E="
            << trace.energy << " q_dot_p=" << 2 * trace.energy
            << " output=" << fs::absolute(o.output) << " total_s=" << seconds(totalStart) << '\n';
}
}  // namespace

int main(int argc, char** argv) {
  try { run(parse(argc, argv)); return 0; }
  catch (const std::exception& error) { std::cerr << "FAIL: " << error.what() << '\n'; return 1; }
}
