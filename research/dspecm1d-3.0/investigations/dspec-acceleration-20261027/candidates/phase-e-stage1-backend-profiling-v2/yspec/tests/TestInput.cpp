#include <gtest/gtest.h>

#include <YSpec/Input.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

using namespace YSpec;
namespace fs = std::filesystem;

const auto kExample = fs::path(YSPEC_EXAMPLES_DIR) / "input.json";

// Writes a document to a file of its own and returns the path.
fs::path Write(const std::string& name, const std::string& text) {
  const auto directory = fs::temp_directory_path() / "yspec-input-tests";
  fs::create_directories(directory);
  const auto path = directory / name;
  std::ofstream(path) << text;
  return path;
}

// A minimal valid document, with members spliced in at the top and among the
// first source's.
std::string Document(const std::string& extra = "",
                     const std::string& source = "") {
  return R"({
  "model": "prem.200",)" +
         extra + R"(
  "output": { "prefix": "out" },
  "degrees": { "max": 10 },
  "frequencies_mHz": { "min": 0.2, "max": 10 },
  "time_series": { "length_min": 30, "time_step_s": 2 },
  "sources": [{)" +
         source + R"(
    "depth_km": 10, "latitude": 0, "longitude": 0,
    "moment_tensor_Nm": { "rr": 1, "rt": 0, "rp": 0, "tt": 0, "tp": 0, "pp": 0 }
  }],
  "receivers": { "depth_km": 0, "locations": [[10, 20]] }
})";
}

TEST(Input, ReadsTheExample) {
  const auto input = ReadSeismogramInput(kExample);
  const auto& c = input.calculation;
  EXPECT_EQ(input.model, kExample.parent_path() / "prem.200");
  EXPECT_EQ(input.prefix, kExample.parent_path() / "yspec.out");
  EXPECT_EQ(c.gravity, Gravity::Full);
  EXPECT_FALSE(c.attenuation);
  EXPECT_EQ(c.quantity, Quantity::Displacement);
  EXPECT_EQ(c.degrees.max, 300);
  // In SI units.
  EXPECT_DOUBLE_EQ(c.band.min, 1.5e-3);
  EXPECT_DOUBLE_EQ(c.band.max, 21e-3);
  EXPECT_DOUBLE_EQ(c.sampling.duration, 18000);
  EXPECT_DOUBLE_EQ(c.sampling.timeStep, 1);
  ASSERT_TRUE(c.taper.has_value());
  EXPECT_DOUBLE_EQ(c.taper->f1, 1.9e-3);
  EXPECT_DOUBLE_EQ(c.taper->f2, 2.0e-3);
  EXPECT_DOUBLE_EQ(c.taper->f3, 20.0e-3);
  EXPECT_DOUBLE_EQ(c.taper->f4, 20.1e-3);
  ASSERT_EQ(input.sources.size(), 1u);
  EXPECT_DOUBLE_EQ(input.sources[0].depth, 647.1e3);
  EXPECT_DOUBLE_EQ(input.sources[0].location.latitude, -13.909919806613164);
  EXPECT_DOUBLE_EQ(input.sources[0].location.longitude, -67.25);
  EXPECT_DOUBLE_EQ(input.sources[0].moment.rr, -7.59e20);
  EXPECT_DOUBLE_EQ(input.sources[0].moment.rt, -2.5e21);
  EXPECT_DOUBLE_EQ(input.sources[0].moment.rp, 4.2e19);
  EXPECT_DOUBLE_EQ(input.sources[0].moment.tt, 7.75e20);
  EXPECT_DOUBLE_EQ(input.sources[0].moment.tp, -2.48e20);
  EXPECT_DOUBLE_EQ(input.sources[0].moment.pp, -1.6e19);
  EXPECT_EQ(input.names, std::vector<std::string>{"1"});
  ASSERT_EQ(input.receivers.size(), 1u);
  EXPECT_DOUBLE_EQ(input.receivers[0].location.latitude, 80.06588632614357);
  EXPECT_DOUBLE_EQ(input.receivers[0].location.longitude, 0);
  EXPECT_DOUBLE_EQ(input.receivers[0].depth, 3e3);
}

TEST(Input, AppliesTheDefaults) {
  const auto input = ReadSeismogramInput(Write("defaults.json", Document()));
  const auto& c = input.calculation;
  EXPECT_EQ(c.gravity, Gravity::Full);
  EXPECT_TRUE(c.attenuation);
  EXPECT_EQ(c.quantity, Quantity::Displacement);
  EXPECT_FALSE(c.corrections);
  EXPECT_EQ(c.degrees.min, 0);
  EXPECT_FALSE(c.taper.has_value());
  EXPECT_EQ(c.solver.tolerance, SolverOptions{}.tolerance);
}

TEST(Input, ReadsSeveralNamedSources) {
  auto text = Document("", R"( "name": "a",)");
  const auto second = std::string(R"(,{ "name": "b",
    "depth_km": 600, "latitude": -10, "longitude": 100,
    "moment_tensor_Nm": { "rr": 0, "rt": 1, "rp": 0, "tt": 0, "tp": 0, "pp": 0 }
  }])");
  text.replace(text.find("}],"), 2, "}" + second);
  const auto input = ReadSeismogramInput(Write("sources.json", text));
  ASSERT_EQ(input.sources.size(), 2u);
  EXPECT_EQ(input.names, (std::vector<std::string>{"a", "b"}));
  EXPECT_DOUBLE_EQ(input.sources[1].depth, 600e3);
  EXPECT_DOUBLE_EQ(input.sources[1].moment.rt, 1);
}

TEST(Input, RefusesAnUnknownMember) {
  const auto path =
      Write("unknown.json", Document(R"( "gravitation": "none",)"));
  try {
    (void)ReadSeismogramInput(path);
    FAIL() << "no exception";
  } catch (const std::runtime_error& e) {
    EXPECT_NE(std::string(e.what()).find("gravitation"), std::string::npos)
        << e.what();
  }
}

TEST(Input, RefusesAnUnknownChoice) {
  const auto path = Write("choice.json", Document(R"( "gravity": "self",)"));
  try {
    (void)ReadSeismogramInput(path);
    FAIL() << "no exception";
  } catch (const std::runtime_error& e) {
    const auto message = std::string(e.what());
    EXPECT_NE(message.find("cowling"), std::string::npos) << message;
  }
}

TEST(Input, RefusesValuesOutOfRange) {
  for (const auto& extra :
       {R"( "solver": { "tolerance": -1e-6 },)",
        R"( "solver": { "start_decay": 0 },)",
        R"( "filter_mHz": [5, 5, 5, 5],)", R"( "filter_mHz": [1, 2, 4, 3],)"}) {
    EXPECT_THROW(
        (void)ReadSeismogramInput(Write("range.json", Document(extra))),
        std::runtime_error)
        << extra;
  }
  auto text = Document();
  text.replace(text.find("[[10, 20]]"), 10, "[[100, 20]]");
  EXPECT_THROW((void)ReadSeismogramInput(Write("latitude.json", text)),
               std::runtime_error);
  EXPECT_THROW((void)ReadSeismogramInput(
                   Write("name.json", Document("", R"( "name": "a/b",)"))),
               std::runtime_error);
}

TEST(Input, RefusesAMissingMember) {
  auto text = Document();
  text.replace(text.find("\"degrees\""), std::string("\"degrees\"").size(),
               "\"degree\"");
  EXPECT_THROW((void)ReadSeismogramInput(Write("missing.json", text)),
               std::runtime_error);
}

// A minimal document with its receivers replaced.
std::string WithReceivers(const std::string& receivers) {
  auto text = Document();
  const auto old =
      std::string(R"("receivers": { "depth_km": 0, "locations": [[10, 20]] })");
  text.replace(text.find(old), old.size(), "\"receivers\": " + receivers);
  return text;
}

TEST(Input, ReadsGroupsOfReceivers) {
  Write("station.pz", "* INPUT UNIT : M\nZEROS 1\nPOLES 1\n-1 0\n");
  const auto input = ReadSeismogramInput(Write("groups.json", WithReceivers(R"([
    { "depth_km": [0, 10], "locations": [[1, 2], [3, 4], [5, 6]] },
    { "depth_km": 3, "profile": { "from": [0, 0], "to": [0, 40],
                                  "count": 5 },
      "instrument": "station.pz" },
    { "depth_km": 0, "profile": { "from": [10, 20], "azimuth_deg": 90,
                                  "distances_deg": [10, 30], "count": 3 } }
  ])")));
  const auto& receivers = input.receivers;
  ASSERT_EQ(receivers.size(), 6u + 5u + 3u);
  // Each depth in turn, with every location at each.
  EXPECT_DOUBLE_EQ(receivers[2].depth, 0);
  EXPECT_DOUBLE_EQ(receivers[2].location.latitude, 5);
  EXPECT_DOUBLE_EQ(receivers[3].depth, 10e3);
  EXPECT_DOUBLE_EQ(receivers[3].location.latitude, 1);
  EXPECT_FALSE(receivers[0].instrument.has_value());
  // The profile between two points on the equator.
  EXPECT_NEAR(receivers[7].location.longitude, 10, 1e-9);
  EXPECT_NEAR(receivers[7].location.latitude, 0, 1e-9);
  ASSERT_TRUE(receivers[7].instrument.has_value());
  EXPECT_EQ(receivers[7].instrument->poles.size(), 1u);
  // The profile at an azimuth.
  const auto profile = GreatCircle({10, 20}, 90, 10, 30, 3);
  for (auto i : {0, 1, 2}) {
    EXPECT_EQ(receivers[11 + i].location.latitude, profile[i].latitude);
    EXPECT_EQ(receivers[11 + i].location.longitude, profile[i].longitude);
  }
}

TEST(Input, RefusesFaultyGroupsOfReceivers) {
  for (const auto& receivers :
       {R"({ "depth_km": 0 })",
        R"({ "depth_km": 0, "locations": [[0, 0]],
                "profile": { "from": [0, 0], "to": [0, 9], "count": 2 } })",
        R"({ "depth_km": [], "locations": [[0, 0]] })",
        R"({ "depth_km": [-1], "locations": [[0, 0]] })",
        R"({ "depth_km": 0, "locations": [] })", R"([])",
        R"({ "depth_km": 0,
                "profile": { "from": [0, 0], "to": [0, 180], "count": 2 } })",
        R"({ "depth_km": 0,
                "profile": { "from": [0, 0], "to": [0, 9], "count": 0 } })",
        R"({ "depth_km": 0,
                "profile": { "from": [0, 0], "to": [0, 9],
                             "azimuth_deg": 0, "count": 2 } })",
        R"({ "depth_km": 0,
                "profile": { "from": [0, 0], "azimuth_deg": 0,
                             "distances_deg": [0, 190], "count": 2 } })",
        R"({ "depth_km": 0, "locations": [[0, 0]],
                "instrument": "no-such-file.pz" })"}) {
    EXPECT_THROW((void)ReadSeismogramInput(
                     Write("receivers.json", WithReceivers(receivers))),
                 std::runtime_error)
        << receivers;
  }
}

TEST(Input, ReadsTheSnapshotExample) {
  const auto path = fs::path(YSPEC_EXAMPLES_DIR) / "snapshot.json";
  const auto input = ReadSnapshotInput(path);
  EXPECT_EQ(input.prefix, path.parent_path() / "snapshot");
  EXPECT_EQ(input.times, (std::vector<double>{300, 600, 900, 1200}));
  ASSERT_EQ(input.grid.depths.size(), 58u);
  EXPECT_DOUBLE_EQ(input.grid.depths.front(), 25e3);
  EXPECT_DOUBLE_EQ(input.grid.depths[1], 75e3);
  EXPECT_DOUBLE_EQ(input.grid.depths.back(), 2875e3);
  ASSERT_EQ(input.grid.distances.size(), 181u);
  EXPECT_DOUBLE_EQ(input.grid.distances[90], 90);
  EXPECT_EQ(input.grid.azimuths, (std::vector<double>{30, 210}));
  // Without a time step, a quarter of the shortest period.
  EXPECT_DOUBLE_EQ(input.calculation.sampling.timeStep, 0.25 / 15e-3);
}

TEST(Input, RefusesFaultySnapshots) {
  const auto good = std::string(R"({
  "model": "prem.200",
  "output": { "prefix": "out" },
  "degrees": { "max": 10 },
  "frequencies_mHz": { "min": 0.2, "max": 10 },
  "time_series": { "length_min": 30 },
  "sources": [{ "depth_km": 10, "latitude": 0, "longitude": 0,
    "moment_tensor_Nm": { "rr": 1, "rt": 0, "rp": 0, "tt": 0, "tp": 0,
                          "pp": 0 } }],
  "snapshots": { "times_s": 600, "depths_km": 100,
                 "distances_deg": { "min": 0, "max": 180, "count": 19 },
                 "azimuths_deg": [0] }
})");
  const auto input = ReadSnapshotInput(Write("snapshot.json", good));
  EXPECT_EQ(input.grid.distances.size(), 19u);
  for (const auto& [from, to] :
       std::initializer_list<std::pair<std::string_view, std::string_view>>{
           {R"("times_s": 600)", R"("times_s": 1900)"},
           {R"("times_s": 600)", R"("times_s": [])"},
           {R"("depths_km": 100)", R"("depths_km": -1)"},
           {R"("max": 180)", R"("max": 181)"},
           {R"("count": 19)", R"("count": 0)"},
           {R"("azimuths_deg": [0])", R"("azimuth_deg": [0])"},
           {R"("snapshots")", R"("receivers")"}}) {
    auto text = good;
    text.replace(text.find(from), from.size(), to);
    EXPECT_THROW((void)ReadSnapshotInput(Write("snapshot.json", text)),
                 std::runtime_error)
        << to;
  }
}

TEST(Input, ReadsForces) {
  auto text = Document();
  const auto moment = std::string(
      R"("moment_tensor_Nm": { "rr": 1, "rt": 0, "rp": 0, "tt": 0, "tp": 0, "pp": 0 })");
  text.replace(text.find(moment), moment.size(),
               R"("force_N": { "r": -1e12, "p": 2e12 })");
  const auto input = ReadSeismogramInput(Write("force.json", text));
  ASSERT_EQ(input.sources.size(), 1u);
  EXPECT_DOUBLE_EQ(input.sources[0].force.r, -1e12);
  EXPECT_DOUBLE_EQ(input.sources[0].force.t, 0);
  EXPECT_DOUBLE_EQ(input.sources[0].force.p, 2e12);
  EXPECT_DOUBLE_EQ(input.sources[0].moment.rr, 0);
  text.replace(text.find(R"("force_N")"), 9, R"("forces_N")");
  EXPECT_THROW((void)ReadSeismogramInput(Write("force.json", text)),
               std::runtime_error);
  auto neither = Document();
  neither.replace(neither.find(moment), moment.size(), R"("name": "x")");
  try {
    (void)ReadSeismogramInput(Write("neither.json", neither));
    FAIL() << "no exception";
  } catch (const std::runtime_error& e) {
    EXPECT_NE(std::string(e.what()).find("force_N"), std::string::npos)
        << e.what();
  }
}

TEST(Input, ReadsTheKernelExamples) {
  const auto directory = fs::path(YSPEC_EXAMPLES_DIR) / "kernels";
  const auto input = ReadKernelInput(directory / "traveltimes.json");
  EXPECT_EQ(input.prefix, directory / "traveltimes");
  ASSERT_EQ(input.measurements.size(), 2u);
  const auto& s = input.measurements[0];
  EXPECT_EQ(s.name, "S");
  EXPECT_EQ(s.component, "T");
  EXPECT_DOUBLE_EQ(s.window.start, 1000);
  EXPECT_DOUBLE_EQ(s.window.end, 1200);
  EXPECT_DOUBLE_EQ(s.window.taper, 40);
  // The transverse direction is horizontal and of unit length.
  EXPECT_DOUBLE_EQ(s.direction.up, 0);
  EXPECT_NEAR(std::hypot(s.direction.north, s.direction.east), 1, 1e-12);
  const auto& long_period = input.measurements[1];
  EXPECT_EQ(long_period.name, "long-period-Z-packet");
  EXPECT_EQ(long_period.component, "Z");
  EXPECT_DOUBLE_EQ(long_period.direction.up, 1);
  EXPECT_EQ(input.parameters, ParameterSet::Velocities);
  EXPECT_EQ(input.scaling, Scaling::Relative);
  EXPECT_EQ(input.calculation.quantity, Quantity::Velocity);
  EXPECT_FALSE(input.volume || input.sections || input.maps);

  const auto volume =
      ReadKernelInput(fs::path(YSPEC_EXAMPLES_DIR) / "kernels-3d/s-wave.json");
  ASSERT_TRUE(volume.volume);
  EXPECT_EQ(volume.volume->degree, 48);
}

TEST(Input, RefusesFaultyKernels) {
  const auto good = std::string(R"({
  "model": "prem.200",
  "output": { "prefix": "out" },
  "degrees": { "max": 10 },
  "frequencies_mHz": { "min": 0.2, "max": 10 },
  "time_series": { "length_min": 30 },
  "sources": [{ "depth_km": 10, "latitude": 0, "longitude": 0,
    "moment_tensor_Nm": { "rr": 1, "rt": 0, "rp": 0, "tt": 0, "tp": 0,
                          "pp": 0 } }],
  "receiver": { "depth_km": 3, "location": [0, 40] },
  "measurements": [
    { "component": "R",
      "traveltime": { "window_s": [600, 900], "taper_s": 50 } },
    { "name": "late", "component": "Z",
      "traveltime": { "window_s": [1200, 1500] } } ],
  "kernels": { "parameters": "anisotropic",
               "mesh": { "order": 4, "max_element_km": 20 },
               "sections": { "spacing_deg": 0.5, "radial_spacing_km": 50 },
               "maps": { "depths_km": [100, 600] } }
})");
  const auto input = ReadKernelInput(Write("kernels.json", good));
  EXPECT_EQ(input.parameters, ParameterSet::Anisotropic);
  EXPECT_EQ(input.scaling, Scaling::Absolute);
  EXPECT_EQ(input.mesh.order, 4);
  EXPECT_DOUBLE_EQ(input.mesh.maxElement, 20e3);
  ASSERT_EQ(input.measurements.size(), 2u);
  // Numbered where not named.
  EXPECT_EQ(input.measurements[0].name, "1");
  EXPECT_EQ(input.measurements[1].name, "late");
  // Radial, from a source on the equator to a receiver east along it.
  EXPECT_NEAR(input.measurements[0].direction.east, 1, 1e-9);
  EXPECT_NEAR(input.measurements[0].direction.north, 0, 1e-9);
  ASSERT_TRUE(input.sections && input.maps);
  EXPECT_DOUBLE_EQ(input.sections->spacing, 0.5);
  EXPECT_DOUBLE_EQ(input.sections->radialSpacing, 50e3);
  EXPECT_DOUBLE_EQ(input.sections->margin, SectionOptions{}.margin);
  EXPECT_EQ(input.maps->depths, (std::vector<double>{100e3, 600e3}));
  for (const auto& [from, to] :
       std::initializer_list<std::pair<std::string_view, std::string_view>>{
           {R"("component": "R")", R"("component": "X")"},
           {R"([600, 900])", R"([900, 600])"},
           {R"([600, 900])", R"([600, 1900])"},
           {R"("taper_s": 50)", R"("taper_s": 200)"},
           {R"("location": [0, 40])", R"("distance_deg": 40)"},
           {R"("name": "late")", R"("name": "1")"},
           {R"("name": "late")", R"("name": "a/b")"},
           {R"("spacing_deg": 0.5)", R"("spacing_deg": 0)"},
           {R"("depths_km": [100, 600])", R"("depths_km": [-1])"},
           {R"("anisotropic")", R"("elastic")"},
           {R"("order": 4)", R"("order": 0)"},
           {R"("traveltime")", R"("amplitude")"},
           {R"("pp": 0 } }],)",
            R"("pp": 0 } }, { "depth_km": 10, "latitude": 0, "longitude": 9,
               "force_N": { "r": 1 } }],)"}}) {
    auto text = good;
    text.replace(text.find(from), from.size(), to);
    EXPECT_THROW((void)ReadKernelInput(Write("kernels.json", text)),
                 std::runtime_error)
        << to;
  }
}

}  // namespace
