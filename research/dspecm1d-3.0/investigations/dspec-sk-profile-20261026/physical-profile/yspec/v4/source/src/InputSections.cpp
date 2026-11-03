#include "InputSections.hpp"

#include <YSpec/Instrument.hpp>
#include <algorithm>
#include <array>
#include <iterator>
#include <optional>
#include <ranges>
#include <set>
#include <stdexcept>

namespace YSpec::InputSections {

void Fail(const Json& json, std::string_view message) {
  throw planetmodel::Error(json.Source(), json.Path(), message);
}

void Allow(const Json& json, std::span<const std::string_view> names) {
  for (const auto& [name, value] : json.AsObject()) {
    if (std::ranges::find(names, name) == names.end()) {
      Fail(value, "is not a recognised member");
    }
  }
}

void Allow(const Json& json, std::initializer_list<std::string_view> names) {
  Allow(json, std::span{names.begin(), names.size()});
}

double Number(const Json& object, std::string_view key) {
  return object.At(key).AsNumber();
}

double Number(const Json& object, std::string_view key, double fallback) {
  const auto* found = object.Find(key);
  return found ? found->AsNumber() : fallback;
}

namespace {

double Latitude(const Json& json) {
  const auto latitude = json.AsNumber();
  if (!(latitude >= -90 && latitude <= 90)) {
    Fail(json, "should be a latitude, between -90 and 90 degrees");
  }
  return latitude;
}

// A depth in km, returned in m.
double Depth(const Json& json) {
  const auto depth = json.AsNumber();
  if (!(depth >= 0)) Fail(json, "should not be negative");
  return depth * 1e3;
}

void ReadPhysics(const Json& json, Calculation& calculation) {
  if (const auto* attenuation = json.Find("attenuation")) {
    calculation.attenuation = attenuation->AsBool();
  }
  if (const auto* gravity = json.Find("gravity")) {
    calculation.gravity =
        Choice<Gravity>(*gravity, {{"none", Gravity::None},
                                   {"cowling", Gravity::Cowling},
                                   {"full", Gravity::Full}});
  }
}

void ReadDegrees(const Json& degrees, Calculation& calculation) {
  Allow(degrees, {"min", "max"});
  if (const auto* min = degrees.Find("min")) {
    calculation.degrees.min = min->AsInt();
  }
  calculation.degrees.max = degrees.At("max").AsInt();
  if (calculation.degrees.min < 0 ||
      calculation.degrees.max < calculation.degrees.min) {
    Fail(degrees, "should satisfy 0 <= min <= max");
  }
}

void ReadBand(const Json& band, Calculation& calculation) {
  Allow(band, {"min", "max"});
  calculation.band = {Number(band, "min") * 1e-3, Number(band, "max") * 1e-3};
  if (!(calculation.band.min > 0 &&
        calculation.band.max > calculation.band.min)) {
    Fail(band, "should satisfy 0 < min < max");
  }
}

// After the band, whose highest frequency gives the default time step: half
// the interval at which it is the Nyquist frequency.
void ReadSampling(const Json& series, Calculation& calculation) {
  Allow(series, {"length_min", "time_step_s"});
  calculation.sampling = {
      Number(series, "length_min") * 60,
      Number(series, "time_step_s", 0.25 / calculation.band.max)};
  if (!(calculation.sampling.duration > 0 &&
        calculation.sampling.timeStep > 0)) {
    Fail(series, "should have a positive length and time step");
  }
}

Taper ReadTaper(const Json& filter) {
  if (filter.Size() != 4) Fail(filter, "should have four corners");
  auto corner = std::array<double, 4>{};
  for (auto i : {0, 1, 2, 3}) {
    corner[i] = filter.At(static_cast<std::size_t>(i)).AsNumber() * 1e-3;
  }
  if (!(corner[0] >= 0 && corner[0] <= corner[1] && corner[1] < corner[2] &&
        corner[2] <= corner[3])) {
    Fail(filter,
         "should be in increasing order, with a pass band between the second "
         "and third corners");
  }
  return {corner[0], corner[1], corner[2], corner[3]};
}

void ReadSolver(const Json& solver, Calculation& calculation) {
  Allow(solver, {"tolerance", "start_decay"});
  auto& options = calculation.solver;
  options.tolerance = Number(solver, "tolerance", options.tolerance);
  options.startDecay = Number(solver, "start_decay", options.startDecay);
  if (!(options.tolerance > 0 && options.tolerance < 1)) {
    Fail(solver.At("tolerance"), "should lie between 0 and 1");
  }
  if (!(options.startDecay > 0)) {
    Fail(solver.At("start_decay"), "should be positive");
  }
}

// A source: a moment tensor, a force, or both.
Source ReadSource(const Json& json) {
  Allow(json, {"name", "depth_km", "latitude", "longitude", "moment_tensor_Nm",
               "force_N"});
  auto source = Source{
      .location = {Latitude(json.At("latitude")), Number(json, "longitude")},
      .depth = Depth(json.At("depth_km")),
      .moment = {},
      .force = {}};
  const auto* moment = json.Find("moment_tensor_Nm");
  const auto* force = json.Find("force_N");
  if (!moment && !force) {
    Fail(json, "should have a \"moment_tensor_Nm\", a \"force_N\", or both");
  }
  if (moment) {
    Allow(*moment, {"rr", "rt", "rp", "tt", "tp", "pp"});
    source.moment = {Number(*moment, "rr"), Number(*moment, "rt"),
                     Number(*moment, "rp"), Number(*moment, "tt"),
                     Number(*moment, "tp"), Number(*moment, "pp")};
  }
  if (force) {
    Allow(*force, {"r", "t", "p"});
    source.force = {Number(*force, "r", 0), Number(*force, "t", 0),
                    Number(*force, "p", 0)};
  }
  return source;
}

void ReadSources(const Json& sources, Shared& shared) {
  for (const auto& source : sources.AsArray()) {
    shared.sources.push_back(ReadSource(source));
    const auto* name = source.Find("name");
    shared.names.push_back(name ? name->AsString()
                                : std::to_string(shared.sources.size()));
    if (shared.names.back().empty() ||
        shared.names.back().find('/') != std::string::npos) {
      Fail(name ? *name : source, "should be a name without a '/'");
    }
  }
  if (shared.sources.empty()) Fail(sources, "is empty");
  if (std::set(shared.names.begin(), shared.names.end()).size() !=
      shared.names.size()) {
    Fail(sources, "should have different names");
  }
}

}  // namespace

std::vector<double> ReadValues(const Json& json, double scale) {
  if (json.IsNumber()) return {json.AsNumber() * scale};
  auto values = std::vector<double>{};
  if (json.IsArray()) {
    for (const auto& value : json.AsArray()) {
      values.push_back(value.AsNumber() * scale);
    }
    if (values.empty()) Fail(json, "is empty");
    return values;
  }
  Allow(json, {"min", "max", "count"});
  const auto min = Number(json, "min");
  const auto max = Number(json, "max");
  const auto& count = json.At("count");
  if (count.AsInt() < 1) Fail(count, "should be at least 1");
  if (!(min <= max)) Fail(json, "should have min <= max");
  const auto n = count.AsInt();
  for (auto i : std::views::iota(0, n)) {
    const auto t = n == 1 ? 0.0 : static_cast<double>(i) / (n - 1);
    values.push_back(((1 - t) * min + t * max) * scale);
  }
  return values;
}

Location ReadLocation(const Json& json) {
  if (json.Size() != 2) Fail(json, "should be [latitude, longitude]");
  return {Latitude(json.At(0)), json.At(1).AsNumber()};
}

std::vector<std::string_view> SharedMembers() {
  return {"model",       "attenuation", "gravity", "degrees", "frequencies_mHz",
          "time_series", "filter_mHz",  "sources", "solver"};
}

Shared ReadShared(const Json& json, const std::filesystem::path& directory) {
  auto shared = Shared{};
  auto& calculation = shared.calculation;
  shared.model = directory / json.At("model").AsString();
  ReadPhysics(json, calculation);
  ReadDegrees(json.At("degrees"), calculation);
  ReadBand(json.At("frequencies_mHz"), calculation);
  ReadSampling(json.At("time_series"), calculation);
  if (const auto* filter = json.Find("filter_mHz")) {
    calculation.taper = ReadTaper(*filter);
  }
  ReadSources(json.At("sources"), shared);
  if (const auto* solver = json.Find("solver")) {
    ReadSolver(*solver, calculation);
  }
  return shared;
}

Quantity ReadQuantity(const Json& output, Quantity fallback) {
  const auto* quantity = output.Find("quantity");
  if (!quantity) return fallback;
  return Choice<Quantity>(*quantity,
                          {{"displacement", Quantity::Displacement},
                           {"velocity", Quantity::Velocity},
                           {"acceleration", Quantity::Acceleration}});
}

namespace {

// One depth or several, in km, returned in m.
std::vector<double> Depths(const Json& json) {
  auto depths = std::vector<double>{};
  if (!json.IsArray()) return {Depth(json)};
  for (const auto& depth : json.AsArray()) depths.push_back(Depth(depth));
  if (depths.empty()) Fail(json, "is empty");
  return depths;
}

std::vector<Location> ReadProfile(const Json& profile) {
  Allow(profile, {"from", "to", "azimuth_deg", "distances_deg", "count"});
  const auto from = ReadLocation(profile.At("from"));
  const auto& count = profile.At("count");
  if (count.AsInt() < 1) Fail(count, "should be at least 1");
  const auto points = static_cast<std::size_t>(count.AsInt());
  const auto* to = profile.Find("to");
  const auto* azimuth = profile.Find("azimuth_deg");
  const auto* distances = profile.Find("distances_deg");
  if (to && !azimuth && !distances) {
    try {
      return GreatCircle(from, ReadLocation(*to), points);
    } catch (const std::invalid_argument& error) {
      Fail(profile, error.what());
    }
  }
  if (!to && azimuth && distances) {
    if (distances->Size() != 2) Fail(*distances, "should be [min, max]");
    const auto min = distances->At(0).AsNumber();
    const auto max = distances->At(1).AsNumber();
    if (!(min >= 0 && min <= max && max <= 180)) {
      Fail(*distances, "should satisfy 0 <= min <= max <= 180");
    }
    return GreatCircle(from, azimuth->AsNumber(), min, max, points);
  }
  Fail(profile,
       "should have either \"to\", or \"azimuth_deg\" and \"distances_deg\"");
}

std::optional<Instrument> ReadInstrument(
    const Json& group, const std::filesystem::path& directory) {
  const auto* file = group.Find("instrument");
  if (!file) return std::nullopt;
  try {
    return ReadSacPoleZero(directory / file->AsString());
  } catch (const planetmodel::Error&) {
    throw;
  } catch (const std::runtime_error& error) {
    Fail(*file, error.what());
  }
}

std::vector<Receiver> ReadReceiverGroup(
    const Json& group, const std::filesystem::path& directory) {
  Allow(group, {"depth_km", "locations", "profile", "instrument"});
  const auto depths = Depths(group.At("depth_km"));
  const auto* list = group.Find("locations");
  const auto* profile = group.Find("profile");
  if ((list != nullptr) == (profile != nullptr)) {
    Fail(group, "should have either \"locations\" or \"profile\"");
  }
  auto locations = std::vector<Location>{};
  if (list) {
    for (const auto& location : list->AsArray()) {
      locations.push_back(ReadLocation(location));
    }
    if (locations.empty()) Fail(*list, "is empty");
  } else {
    locations = ReadProfile(*profile);
  }
  const auto instrument = ReadInstrument(group, directory);
  auto receivers = std::vector<Receiver>{};
  for (auto depth : depths) {
    for (const auto& location : locations) {
      receivers.push_back(
          {.location = location, .depth = depth, .instrument = instrument});
    }
  }
  return receivers;
}

}  // namespace

Receiver ReadReceiver(const Json& json, const Location& from,
                      const std::filesystem::path& directory) {
  Allow(json,
        {"depth_km", "location", "distance_deg", "azimuth_deg", "instrument"});
  auto receiver = Receiver{.location = {},
                           .depth = Depth(json.At("depth_km")),
                           .instrument = ReadInstrument(json, directory)};
  const auto* location = json.Find("location");
  const auto* distance = json.Find("distance_deg");
  const auto* azimuth = json.Find("azimuth_deg");
  if (location && !distance && !azimuth) {
    receiver.location = ReadLocation(*location);
  } else if (!location && distance && azimuth) {
    const auto d = distance->AsNumber();
    if (!(d >= 0 && d <= 180)) Fail(*distance, "should lie between 0 and 180");
    receiver.location = GreatCircle(from, azimuth->AsNumber(), d, d, 1).front();
  } else {
    Fail(json,
         "should have either \"location\", or \"distance_deg\" and "
         "\"azimuth_deg\"");
  }
  return receiver;
}

std::vector<Receiver> ReadReceivers(const Json& json,
                                    const std::filesystem::path& directory) {
  if (!json.IsArray()) return ReadReceiverGroup(json, directory);
  auto receivers = std::vector<Receiver>{};
  for (const auto& group : json.AsArray()) {
    std::ranges::move(ReadReceiverGroup(group, directory),
                      std::back_inserter(receivers));
  }
  if (receivers.empty()) Fail(json, "is empty");
  return receivers;
}

}  // namespace YSpec::InputSections
