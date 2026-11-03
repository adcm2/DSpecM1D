#include <GSHTrans/GaussLegendreGrid.hpp>

#include <bit>
#include <cerrno>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_set>
#include <vector>

namespace {

constexpr std::string_view kHeader =
    "objective,modulus,radius_m,theta_rad,phi_rad,element,fluid,"
    "kernel_objective_per_Pa_m3,deltaX_Pa,dr_m,dOmega";

struct Row {
  std::vector<std::string> fields;
  double radius{};
  double theta{};
  double phi{};
  double dr{};
  double dOmega{};
};

struct Trace {
  std::vector<std::string> metadata;
  std::uint64_t radiusBits{};
  double radius{};
  double dr{};
  std::vector<std::pair<double, double>> angles;
  std::vector<double> inputWeights;
};

double number(const std::string& text, const std::string& field) {
  char* end = nullptr;
  errno = 0;
  const double value = std::strtod(text.c_str(), &end);
  if (errno || end != text.c_str() + text.size() || !std::isfinite(value))
    throw std::runtime_error("invalid finite number in " + field + ": " + text);
  return value;
}

std::vector<std::string> split(const std::string& line) {
  std::vector<std::string> out;
  std::size_t start = 0;
  for (;;) {
    const auto comma = line.find(',', start);
    out.emplace_back(line.substr(start, comma == std::string::npos
                                            ? comma
                                            : comma - start));
    if (comma == std::string::npos) break;
    start = comma + 1;
  }
  return out;
}

std::string jsonQuote(std::string_view value) {
  std::string out = "\"";
  for (const char c : value) {
    if (c == '\\' || c == '"') out += '\\';
    if (c == '\n') out += "\\n";
    else if (c == '\r') out += "\\r";
    else out += c;
  }
  return out + '"';
}

void run(const std::filesystem::path& inputPath,
         const std::filesystem::path& outputPath, long long band) {
  if (band < 0 || band > std::numeric_limits<std::ptrdiff_t>::max() / 2)
    throw std::runtime_error("L must be a non-negative representable band");
  if (std::filesystem::exists(outputPath))
    throw std::runtime_error("refusing existing output: " + outputPath.string());

  std::ifstream input(inputPath);
  if (!input) throw std::runtime_error("cannot open input: " + inputPath.string());
  std::string line;
  if (!std::getline(input, line))
    throw std::runtime_error("input is empty");
  if (!line.empty() && line.back() == '\r') line.pop_back();
  if (line != kHeader)
    throw std::runtime_error("input header does not match the 11-column schema");

  std::vector<Trace> traces;
  std::map<std::pair<std::string, std::uint64_t>, std::size_t> traceByKey;
  std::vector<std::pair<double, double>> referenceAngles;
  std::unordered_set<std::string> seenAngles;
  std::size_t activeTrace = std::numeric_limits<std::size_t>::max();
  std::size_t rows = 0;
  std::size_t filteredRows = 0;
  while (std::getline(input, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty()) continue;
    const auto fields = split(line);
    if (fields.size() < 2)
      throw std::runtime_error("row does not contain objective and modulus columns");
    if (fields[0] != "J_Z" || fields[1] != "A") {
      ++filteredRows;
      continue;
    }
    if (fields.size() != 11) throw std::runtime_error("row does not have 11 columns");
    Row row{fields};
    row.radius = number(fields[2], "radius_m");
    row.theta = number(fields[3], "theta_rad");
    row.phi = number(fields[4], "phi_rad");
    row.dr = number(fields[9], "dr_m");
    row.dOmega = number(fields[10], "dOmega");
    if (!(row.radius > 0) || !(row.theta > 0 && row.theta < std::numbers::pi) ||
        !(row.phi >= 0 && row.phi < 2 * std::numbers::pi) || !(row.dr > 0) ||
        !(row.dOmega > 0))
      throw std::runtime_error("row has out-of-range radius, angle, spacing, or weight");

    const auto bits = std::bit_cast<std::uint64_t>(row.radius);
    const auto key = std::pair{fields[5], bits};
    auto [it, inserted] = traceByKey.emplace(key, traces.size());
    if (inserted) {
      Trace trace;
      trace.metadata = {fields[0], fields[1], fields[2], fields[5], fields[6], fields[9]};
      trace.radiusBits = bits;
      trace.radius = row.radius;
      trace.dr = row.dr;
      traces.push_back(std::move(trace));
    }
    auto& trace = traces[it->second];
    if (activeTrace != it->second) {
      if (activeTrace != std::numeric_limits<std::size_t>::max() && !inserted)
        throw std::runtime_error("radial trace rows are split into multiple blocks");
      activeTrace = it->second;
      seenAngles.clear();
    }
    if (fields[0] != trace.metadata[0] || fields[1] != trace.metadata[1] ||
        fields[2] != trace.metadata[2] || fields[6] != trace.metadata[4] ||
        fields[9] != trace.metadata[5] || row.dr != trace.dr)
      throw std::runtime_error("radial trace metadata or dr_m changes within a trace");

    const auto angleKey =
        std::to_string(std::bit_cast<std::uint64_t>(row.theta)) + "," +
        std::to_string(std::bit_cast<std::uint64_t>(row.phi));
    if (!seenAngles.insert(angleKey).second)
      throw std::runtime_error("duplicate angular sample within one radial trace block");
    trace.angles.emplace_back(row.theta, row.phi);
    trace.inputWeights.push_back(row.dOmega);
    if (it->second == 0) {
      // The first trace defines the ordered angular block. Later blocks are
      // compared against it after all rows are read.
      referenceAngles.emplace_back(row.theta, row.phi);
    }
    ++rows;

  }
  if (input.bad()) throw std::runtime_error("error reading input CSV");
  if (traces.empty() || referenceAngles.empty())
    throw std::runtime_error("input contains no radial trace samples");
  for (std::size_t t = 0; t < traces.size(); ++t) {
    const auto& trace = traces[t];
    if (trace.angles.size() != referenceAngles.size())
      throw std::runtime_error("radial traces have different angular block sizes");
    if (t != 0) {
      if (trace.angles != referenceAngles)
        throw std::runtime_error("radial traces have different ordered angular blocks");
      if (trace.inputWeights != traces.front().inputWeights)
        throw std::runtime_error("radial traces have inconsistent angular dOmega values");
    }
  }

  using Grid = GSHTrans::GaussLegendreGrid<double, GSHTrans::All,
                                            GSHTrans::All>;
  const auto grid = Grid::ForBand(static_cast<std::ptrdiff_t>(band), 2, 2.0,
                                  FFTWpp::Estimate);
  const auto thetas = grid.CoLatitudes();
  const auto phis = grid.Longitudes();
  const auto weights = grid.Weights();
  std::vector<double> thetaValues(thetas.begin(), thetas.end());
  std::vector<double> phiValues(phis.begin(), phis.end());
  std::vector<double> omegaValues(weights.begin(), weights.end());
  const auto expectedCount = thetaValues.size() * phiValues.size();
  if (omegaValues.size() != expectedCount)
    throw std::runtime_error("GSHTrans grid dimensions disagree with its weight count");

  std::ofstream output(outputPath, std::ios::out | std::ios::trunc);
  if (!output) throw std::runtime_error("cannot create output: " + outputPath.string());
  output << kHeader << '\n' << std::setprecision(17);
  double weightSum = 0;
  for (const auto weight : omegaValues) weightSum += weight;
  std::size_t outRows = 0;
  for (const auto& trace : traces) {
    for (std::size_t i = 0; i < thetaValues.size(); ++i) {
      for (std::size_t j = 0; j < phiValues.size(); ++j) {
        const auto weight = omegaValues[i * phiValues.size() + j];
        output << trace.metadata[0] << ',' << trace.metadata[1] << ','
               << trace.metadata[2] << ',' << thetaValues[i] << ',' << phiValues[j]
               << ',' << trace.metadata[3] << ',' << trace.metadata[4] << ",,0,"
               << trace.metadata[5] << ',' << weight << '\n';
        ++outRows;
      }
    }
  }
  output.close();
  if (!output) throw std::runtime_error("failed while writing output CSV");
  std::cout << std::setprecision(17)
            << "{\"grid_mode\":\"GSHTrans::GaussLegendreGrid::ForBand\","
            << "\"band_L\":" << band << ",\"resolved_band\":" << grid.MaxDegree()
            << ",\"n_theta\":" << thetaValues.size() << ",\"n_phi\":"
            << phiValues.size() << ",\"points_per_trace\":" << expectedCount
            << ",\"trace_count\":" << traces.size() << ",\"input_rows\":"
            << rows << ",\"filtered_rows\":" << filteredRows
            << ",\"output_rows\":" << outRows << ",\"weight_sum\":"
            << weightSum << ",\"input\":" << jsonQuote(inputPath.string())
            << ",\"output\":" << jsonQuote(outputPath.string()) << "}\n";
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 4) {
    std::cerr << "usage: kernel_grid_generator INPUT.csv OUTPUT.csv L\n";
    return 2;
  }
  try {
    std::size_t parsed = 0;
    const auto band = std::stoll(argv[3], &parsed);
    if (parsed != std::string(argv[3]).size())
      throw std::runtime_error("L must be an integer");
    run(argv[1], argv[2], band);
  } catch (const std::exception& error) {
    std::cerr << "kernel_grid_generator: " << error.what() << '\n';
    return 1;
  }
}
