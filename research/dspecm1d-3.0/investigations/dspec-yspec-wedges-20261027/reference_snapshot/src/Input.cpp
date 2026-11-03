#include <YSpec/Input.hpp>

#include "InputSections.hpp"

namespace YSpec {

using namespace InputSections;

SeismogramInput ReadSeismogramInput(const std::filesystem::path& path) {
  const auto json = planetmodel::ReadJsonFile(path);
  auto members = SharedMembers();
  members.insert(members.end(), {"output", "receivers"});
  Allow(json, members);
  const auto directory = path.parent_path();

  auto shared = ReadShared(json, directory);
  auto input = SeismogramInput{.model = std::move(shared.model),
                               .calculation = shared.calculation,
                               .sources = std::move(shared.sources),
                               .names = std::move(shared.names),
                               .receivers = {},
                               .prefix = {}};

  const auto& output = json.At("output");
  Allow(output, {"prefix", "quantity", "potential_and_tilt_corrections"});
  input.prefix = directory / output.At("prefix").AsString();
  auto& calculation = input.calculation;
  calculation.quantity = ReadQuantity(output, calculation.quantity);
  if (const auto* corrections = output.Find("potential_and_tilt_corrections")) {
    calculation.corrections = corrections->AsBool();
  }

  input.receivers = ReadReceivers(json.At("receivers"), directory);
  return input;
}

SnapshotInput ReadSnapshotInput(const std::filesystem::path& path) {
  const auto json = planetmodel::ReadJsonFile(path);
  auto members = SharedMembers();
  members.insert(members.end(), {"output", "snapshots"});
  Allow(json, members);
  const auto directory = path.parent_path();

  auto shared = ReadShared(json, directory);
  auto input = SnapshotInput{.model = std::move(shared.model),
                             .calculation = shared.calculation,
                             .sources = std::move(shared.sources),
                             .names = std::move(shared.names),
                             .grid = {},
                             .times = {},
                             .prefix = {}};

  const auto& output = json.At("output");
  Allow(output, {"prefix", "quantity"});
  input.prefix = directory / output.At("prefix").AsString();
  input.calculation.quantity = ReadQuantity(output, input.calculation.quantity);

  const auto& snapshots = json.At("snapshots");
  Allow(snapshots, {"times_s", "depths_km", "distances_deg", "azimuths_deg"});
  const auto& times = snapshots.At("times_s");
  input.times = ReadValues(times);
  for (auto time : input.times) {
    if (!(time >= 0 && time <= input.calculation.sampling.duration)) {
      Fail(times, "should lie between zero and the length of the time series");
    }
  }
  const auto& depths = snapshots.At("depths_km");
  input.grid.depths = ReadValues(depths, 1e3);
  for (auto depth : input.grid.depths) {
    if (!(depth >= 0)) Fail(depths, "should not be negative");
  }
  const auto& distances = snapshots.At("distances_deg");
  input.grid.distances = ReadValues(distances);
  for (auto distance : input.grid.distances) {
    if (!(distance >= 0 && distance <= 180)) {
      Fail(distances, "should lie between 0 and 180");
    }
  }
  input.grid.azimuths = ReadValues(snapshots.At("azimuths_deg"));
  return input;
}

namespace {

// A measurement: a traveltime on a component.
KernelMeasurement ReadMeasurement(const Json& json, const Location& from,
                                  const Location& at, double duration) {
  Allow(json, {"name", "component", "traveltime"});
  auto measurement = KernelMeasurement{};
  const auto& component = json.At("component");
  measurement.component = component.AsString();
  measurement.direction =
      Choice<Direction>(component, {{"Z", Along(Vertical)},
                                    {"N", Along(North)},
                                    {"E", Along(East)},
                                    {"R", RadialDirection(from, at)},
                                    {"T", TransverseDirection(from, at)}});
  const auto& traveltime = json.At("traveltime");
  Allow(traveltime, {"window_s", "taper_s"});
  const auto& window = traveltime.At("window_s");
  if (window.Size() != 2) Fail(window, "should be [start, end]");
  measurement.window = {.start = window.At(0).AsNumber(),
                        .end = window.At(1).AsNumber(),
                        .taper = Number(traveltime, "taper_s", 0)};
  const auto& w = measurement.window;
  if (!(w.start >= 0 && w.start < w.end && w.end <= duration)) {
    Fail(window,
         "should satisfy 0 <= start < end <= the length of the time series");
  }
  if (!(w.taper >= 0 && 2 * w.taper <= w.end - w.start)) {
    Fail(traveltime.At("taper_s"),
         "should be at least zero and at most half the window");
  }
  return measurement;
}

// A positive number, in degrees or km (returned in m), with a default.
double Positive(const Json& json, std::string_view key, double fallback,
                double scale = 1) {
  const auto value = Number(json, key, fallback / scale) * scale;
  if (!(value > 0)) Fail(json.At(key), "should be positive");
  return value;
}

void ReadKernelOptions(const Json& kernels, KernelInput& input) {
  Allow(kernels,
        {"parameters", "scaling", "mesh", "volume", "sections", "maps"});
  if (const auto* parameters = kernels.Find("parameters")) {
    input.parameters = Choice<ParameterSet>(
        *parameters, {{"love", ParameterSet::Love},
                      {"moduli", ParameterSet::Moduli},
                      {"velocities", ParameterSet::Velocities},
                      {"anisotropic", ParameterSet::Anisotropic}});
  }
  if (const auto* scaling = kernels.Find("scaling")) {
    input.scaling = Choice<Scaling>(
        *scaling,
        {{"absolute", Scaling::Absolute}, {"relative", Scaling::Relative}});
  }
  if (const auto* mesh = kernels.Find("mesh")) {
    Allow(*mesh, {"order", "max_element_km"});
    if (const auto* order = mesh->Find("order")) {
      input.mesh.order = order->AsInt();
      if (input.mesh.order < 1) Fail(*order, "should be at least 1");
    }
    input.mesh.maxElement =
        Positive(*mesh, "max_element_km", input.mesh.maxElement, 1e3);
  }
  if (const auto* volume = kernels.Find("volume")) {
    Allow(*volume, {"degree", "max_element_km"});
    auto grid = VolumeGrid{};
    if (const auto* degree = volume->Find("degree")) {
      grid.degree = degree->AsInt();
      if (grid.degree < 1) Fail(*degree, "should be at least 1");
    }
    grid.maxElement = Positive(*volume, "max_element_km", grid.maxElement, 1e3);
    input.volume = grid;
  }
  if (const auto* sections = kernels.Find("sections")) {
    Allow(*sections,
          {"margin_deg", "half_width_deg", "spacing_deg", "radial_spacing_km"});
    auto options = SectionOptions{};
    options.margin = Number(*sections, "margin_deg", options.margin);
    if (!(options.margin >= 0)) {
      Fail(sections->At("margin_deg"), "should not be negative");
    }
    options.halfWidth =
        Positive(*sections, "half_width_deg", options.halfWidth);
    options.spacing = Positive(*sections, "spacing_deg", options.spacing);
    options.radialSpacing =
        Positive(*sections, "radial_spacing_km", options.radialSpacing, 1e3);
    input.sections = options;
  }
  if (const auto* maps = kernels.Find("maps")) {
    Allow(*maps, {"depths_km", "margin_deg", "half_width_deg", "spacing_deg"});
    auto options = MapOptions{};
    const auto& depths = maps->At("depths_km");
    options.depths = ReadValues(depths, 1e3);
    for (auto depth : options.depths) {
      if (!(depth >= 0)) Fail(depths, "should not be negative");
    }
    options.margin = Number(*maps, "margin_deg", options.margin);
    if (!(options.margin >= 0)) {
      Fail(maps->At("margin_deg"), "should not be negative");
    }
    options.halfWidth = Positive(*maps, "half_width_deg", options.halfWidth);
    options.spacing = Positive(*maps, "spacing_deg", options.spacing);
    input.maps = options;
  }
}

}  // namespace

KernelInput ReadKernelInput(const std::filesystem::path& path) {
  const auto json = planetmodel::ReadJsonFile(path);
  auto members = SharedMembers();
  members.insert(members.end(),
                 {"output", "receiver", "measurements", "kernels"});
  Allow(json, members);
  const auto directory = path.parent_path();

  auto shared = ReadShared(json, directory);
  if (shared.sources.size() != 1) {
    Fail(json.At("sources"), "should hold one source for kernels");
  }
  auto input = KernelInput{.model = std::move(shared.model),
                           .calculation = shared.calculation,
                           .source = shared.sources.front(),
                           .name = shared.names.front(),
                           .receiver = {},
                           .measurements = {},
                           .parameters = ParameterSet::Velocities,
                           .scaling = Scaling::Absolute,
                           .mesh = {},
                           .volume = {},
                           .sections = {},
                           .maps = {},
                           .prefix = {}};

  const auto& output = json.At("output");
  Allow(output, {"prefix", "quantity"});
  input.prefix = directory / output.At("prefix").AsString();
  input.calculation.quantity = ReadQuantity(output, input.calculation.quantity);

  input.receiver =
      ReadReceiver(json.At("receiver"), input.source.location, directory);
  const auto& measurements = json.At("measurements");
  for (const auto& measurement : measurements.AsArray()) {
    auto read = ReadMeasurement(measurement, input.source.location,
                                input.receiver.location,
                                input.calculation.sampling.duration);
    const auto* name = measurement.Find("name");
    read.name =
        name ? name->AsString() : std::to_string(input.measurements.size() + 1);
    if (read.name.empty() || read.name.find('/') != std::string::npos) {
      Fail(name ? *name : measurement, "should be a name without a '/'");
    }
    for (const auto& other : input.measurements) {
      if (other.name == read.name) {
        Fail(name ? *name : measurement, "should differ from the others'");
      }
    }
    input.measurements.push_back(std::move(read));
  }
  if (input.measurements.empty()) Fail(measurements, "is empty");

  if (const auto* kernels = json.Find("kernels")) {
    ReadKernelOptions(*kernels, input);
  }
  return input;
}

}  // namespace YSpec
