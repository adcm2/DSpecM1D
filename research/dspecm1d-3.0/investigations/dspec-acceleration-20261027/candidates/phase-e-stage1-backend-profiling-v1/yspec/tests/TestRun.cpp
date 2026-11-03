#include <gtest/gtest.h>

#include <YSpec/Run.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <format>
#include <fstream>
#include <numbers>
#include <sstream>
#include <string>
#include <vector>

namespace {

using namespace YSpec;
namespace fs = std::filesystem;

TEST(Run, WritesASeismogramFilePerSourceAndReceiver) {
  const auto directory = fs::temp_directory_path() / "yspec-run-test";
  fs::remove_all(directory);
  fs::create_directories(directory);
  fs::copy_file(fs::path(YSPEC_EXAMPLES_DIR) / "prem.200",
                directory / "prem.200");
  std::ofstream(directory / "input.json") << R"({
  "model": "prem.200",
  "output": { "prefix": "out" },
  "degrees": { "max": 20 },
  "frequencies_mHz": { "min": 0.2, "max": 4 },
  "time_series": { "length_min": 30, "time_step_s": 60 },
  "sources": [
    { "name": "a", "depth_km": 30, "latitude": 0, "longitude": 0,
      "moment_tensor_Nm": { "rr": 1e18, "rt": 0, "rp": 0, "tt": 0, "tp": 0,
                            "pp": 0 } },
    { "name": "b", "depth_km": 100, "latitude": 10, "longitude": 10,
      "moment_tensor_Nm": { "rr": 0, "rt": 1e18, "rp": 0, "tt": 0, "tp": 0,
                            "pp": 0 } }
  ],
  "receivers": { "depth_km": 3, "locations": [[20, 30], [-20, 40]] }
})";
  const auto seismograms = RunSeismograms(directory / "input.json");
  ASSERT_EQ(seismograms.Sources(), 2u);
  ASSERT_EQ(seismograms.Receivers(), 2u);
  EXPECT_EQ(seismograms.Samples(), 31u);
  for (auto [s, name] : {std::pair{0, "a"}, {1, "b"}}) {
    for (auto k : {0, 1}) {
      const auto path = directory / std::format("out.{}.{}", name, k + 1);
      ASSERT_TRUE(fs::exists(path)) << path;
      auto file = std::ifstream(path);
      auto count = std::size_t{0};
      for (auto row = std::array<double, 4>{};
           file >> row[0] >> row[1] >> row[2] >> row[3]; ++count) {
        EXPECT_DOUBLE_EQ(row[0], 60.0 * count);
        for (auto c : {Vertical, North, East}) {
          const auto value = seismograms(s, k, c, count);
          EXPECT_NEAR(row[1 + c], value, 1e-6 * std::abs(value) + 1e-300);
        }
      }
      EXPECT_EQ(count, seismograms.Samples());
    }
  }
  // And the table of where the receivers are seen from the sources: from
  // source a, at the origin, the receiver at geographic latitude 20 and
  // longitude 30 is at the distance of its geocentric latitude.
  auto table = std::ifstream(directory / "out.receivers");
  ASSERT_TRUE(table);
  auto line = std::string{};
  std::getline(table, line);
  EXPECT_EQ(line.front(), '#');
  auto rows = 0;
  for (; std::getline(table, line); ++rows) {
    auto row = std::istringstream(line);
    auto file = std::string{};
    auto source = std::string{};
    auto values = std::array<double, 6>{};
    row >> file >> source;
    for (auto& value : values) row >> value;
    ASSERT_TRUE(row) << line;
    EXPECT_TRUE(fs::exists(directory / file)) << file;
    if (rows == 0) {
      EXPECT_EQ(file, "out.a.1");
      const auto degree = std::numbers::pi / 180;
      const auto latitude = std::atan(0.993277 * std::tan(20 * degree));
      EXPECT_NEAR(
          values[3],
          std::acos(std::cos(latitude) * std::cos(30 * degree)) / degree, 1e-3);
    }
  }
  EXPECT_EQ(rows, 4);
}

TEST(Run, WritesASnapshotFilePerTime) {
  const auto directory = fs::temp_directory_path() / "yspec-run-snapshot-test";
  fs::remove_all(directory);
  fs::create_directories(directory);
  fs::copy_file(fs::path(YSPEC_EXAMPLES_DIR) / "prem.200",
                directory / "prem.200");
  std::ofstream(directory / "input.json") << R"({
  "model": "prem.200",
  "output": { "prefix": "snap" },
  "degrees": { "max": 20 },
  "frequencies_mHz": { "min": 0.2, "max": 4 },
  "time_series": { "length_min": 30 },
  "sources": [
    { "depth_km": 30, "latitude": 10, "longitude": 20,
      "moment_tensor_Nm": { "rr": 1e18, "rt": 0, "rp": 0, "tt": 0, "tp": 0,
                            "pp": 0 } }
  ],
  "snapshots": { "times_s": [600, 1200], "depths_km": [100, 500],
                 "distances_deg": { "min": 0, "max": 90, "count": 4 },
                 "azimuths_deg": [0, 180] }
})";
  const auto snapshots = RunSnapshots(directory / "input.json");
  for (auto time : {0, 1}) {
    const auto path = directory / std::format("snap.{}", time + 1);
    ASSERT_TRUE(fs::exists(path)) << path;
    auto file = std::ifstream(path);
    auto line = std::string{};
    std::getline(file, line);
    EXPECT_EQ(line.front(), '#');
    std::getline(file, line);
    EXPECT_EQ(line.front(), '#');
    auto count = std::size_t{0};
    for (auto row = std::array<double, 8>{};
         file >> row[0] >> row[1] >> row[2] >> row[3] >> row[4] >> row[5] >>
         row[6] >> row[7];
         ++count) {
      // Depth by depth, azimuth by azimuth, distance by distance.
      const auto depth = count / 8;
      const auto azimuth = (count / 4) % 2;
      const auto distance = count % 4;
      EXPECT_DOUBLE_EQ(row[0], depth == 0 ? 100 : 500);
      EXPECT_DOUBLE_EQ(row[1], 30.0 * static_cast<double>(distance));
      EXPECT_DOUBLE_EQ(row[2], azimuth == 0 ? 0 : 180);
      for (auto c : {Up, Away, Transverse}) {
        const auto value = snapshots(0, time, depth, distance, azimuth, c);
        EXPECT_NEAR(row[5 + c], value, 1e-6 * std::abs(value) + 1e-300);
      }
    }
    EXPECT_EQ(count, 16u);
  }
  // At the epicentre the location is the source's.
  auto file = std::ifstream(directory / "snap.1");
  auto line = std::string{};
  std::getline(file, line);
  std::getline(file, line);
  auto row = std::array<double, 8>{};
  file >> row[0] >> row[1] >> row[2] >> row[3] >> row[4];
  EXPECT_DOUBLE_EQ(row[3], 10);
  EXPECT_DOUBLE_EQ(row[4], 20);
}

TEST(Run, WritesKernelsAndTheirMeasurement) {
  const auto directory = fs::temp_directory_path() / "yspec-run-kernel-test";
  fs::remove_all(directory);
  fs::create_directories(directory);
  fs::copy_file(fs::path(YSPEC_EXAMPLES_DIR) / "prem.200",
                directory / "prem.200");
  std::ofstream(directory / "input.json") << R"({
  "model": "prem.200",
  "output": { "prefix": "out", "quantity": "velocity" },
  "degrees": { "max": 20 },
  "frequencies_mHz": { "min": 0.2, "max": 4 },
  "time_series": { "length_min": 30, "time_step_s": 30 },
  "sources": [
    { "depth_km": 30, "latitude": 0, "longitude": 0,
      "moment_tensor_Nm": { "rr": 1e18, "rt": 0, "rp": 0, "tt": 0, "tp": 0,
                            "pp": 0 } }
  ],
  "receiver": { "depth_km": 3, "distance_deg": 30, "azimuth_deg": 45 },
  "measurements": [
    { "name": "Z", "component": "Z",
      "traveltime": { "window_s": [400, 1000], "taper_s": 100 } },
    { "name": "T", "component": "T",
      "traveltime": { "window_s": [400, 1000], "taper_s": 100 } } ],
  "kernels": { "parameters": "moduli", "scaling": "relative",
               "mesh": { "order": 3, "max_element_km": 200 },
               "sections": { "margin_deg": 5, "half_width_deg": 5,
                             "spacing_deg": 2, "radial_spacing_km": 500 },
               "maps": { "depths_km": 100, "margin_deg": 5,
                         "half_width_deg": 5, "spacing_deg": 2 } }
})";
  const auto result = RunKernels(directory / "input.json");
  ASSERT_EQ(result.measurements.size(), 2u);
  const auto& z = result.measurements[0];
  EXPECT_EQ(z.named.names,
            (std::vector<std::string>{"rho", "kappa", "mu", "qkappa", "qmu"}));
  // The kernels: the header, then a line for each radius.
  auto table = std::ifstream(directory / "out.Z.kernels");
  ASSERT_TRUE(table);
  auto line = std::string{};
  auto header = std::vector<std::string>{};
  auto rows = std::size_t{0};
  while (std::getline(table, line)) {
    if (line.starts_with('#')) {
      header.push_back(line);
      continue;
    }
    auto row = std::istringstream(line);
    auto values = std::array<double, 8>{};
    for (auto& value : values) row >> value;
    ASSERT_TRUE(row) << line;
    EXPECT_NEAR(values[0] * 1e3, z.kernels.Radii()[rows], 1e-2);
    EXPECT_NEAR(values[4], z.named.values[1][rows],
                1e-5 * std::abs(z.named.values[1][rows]) + 1e-300);
    ++rows;
  }
  EXPECT_EQ(rows, z.kernels.Size());
  ASSERT_GE(header.size(), 2u);
  EXPECT_EQ(header.back(),
            "# radius_km depth_km weight_m3 rho kappa mu qkappa qmu");
  EXPECT_TRUE(std::ranges::any_of(header, [](const std::string& h) {
    return h.starts_with("# at_source ");
  }));
  // The measurement: a line for each sample.
  auto measurement = std::ifstream(directory / "out.Z.measurement");
  ASSERT_TRUE(measurement);
  auto samples = std::size_t{0};
  while (std::getline(measurement, line)) {
    if (!line.starts_with('#')) ++samples;
  }
  EXPECT_EQ(samples, result.seismogram.Samples());
  // The sections and the maps of each, a line for each node and place.
  for (const auto& name : {"Z", "T"}) {
    const auto& kernels = result.measurements[name == std::string("Z") ? 0 : 1];
    ASSERT_TRUE(kernels.path && kernels.across && kernels.maps);
    auto count = [&](const std::string& suffix) {
      auto file =
          std::ifstream(directory / ("out." + std::string(name) + suffix));
      auto lines = std::size_t{0};
      for (auto line = std::string{}; std::getline(file, line);) {
        if (!line.starts_with('#')) ++lines;
      }
      return lines;
    };
    auto size = [](const SampledKernels& k) {
      return k.Nodes().size() * k.Places().size();
    };
    EXPECT_EQ(count(".sections"), size(*kernels.path) + size(*kernels.across));
    EXPECT_EQ(count(".maps"), size(*kernels.maps));
  }
}

}  // namespace
