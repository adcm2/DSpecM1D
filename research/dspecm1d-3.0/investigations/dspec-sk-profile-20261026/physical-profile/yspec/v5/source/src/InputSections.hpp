#pragma once

// The sections of the JSON input files, shared by the programs. Each
// program's reader (src/Input.cpp) checks the top level of its file, reads
// the sections it shares with the others, and then its own.

#include <YSpec/Calculation.hpp>
#include <YSpec/Receiver.hpp>
#include <YSpec/Source.hpp>
#include <filesystem>
#include <format>
#include <initializer_list>
#include <planetmodel/json.hpp>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace YSpec::InputSections {

using planetmodel::Json;

// Throws an error naming the file, the member and the fault.
[[noreturn]] void Fail(const Json& json, std::string_view message);

// Refuses any member of an object not among those named, so that a misspelt
// member is not silently ignored.
void Allow(const Json& json, std::span<const std::string_view> names);
void Allow(const Json& json, std::initializer_list<std::string_view> names);

// One of a set of words, each standing for a value.
template <class T>
T Choice(const Json& json,
         std::initializer_list<std::pair<std::string_view, T>> choices) {
  const auto& word = json.AsString();
  for (const auto& [name, value] : choices) {
    if (word == name) return value;
  }
  auto names = std::string{};
  for (const auto& [name, value] : choices) {
    names += (names.empty() ? "\"" : ", \"") + std::string(name) + "\"";
  }
  Fail(json, std::format("is \"{}\", which is none of {}", word, names));
}

// A number that must be there, and one that may not be.
double Number(const Json& object, std::string_view key);
double Number(const Json& object, std::string_view key, double fallback);

// A [latitude, longitude] pair, in degrees.
Location ReadLocation(const Json& json);

// Values given as one number, an array of them, or a range
// { "min": a, "max": b, "count": n } of n equally spaced, both ends included,
// each multiplied by a scale.
std::vector<double> ReadValues(const Json& json, double scale = 1);

// What every program reads: the model, how the equations are solved and
// summed, and the sources.
struct Shared {
  std::filesystem::path model;
  Calculation calculation;
  std::vector<Source> sources;
  std::vector<std::string> names;
};

// The members of the top level that every program reads.
std::vector<std::string_view> SharedMembers();

// Reads the shared members of the top level of a file in a directory, to
// which relative paths are relative.
Shared ReadShared(const Json& json, const std::filesystem::path& directory);

// The optional "quantity" of an output section.
Quantity ReadQuantity(const Json& output, Quantity fallback);

// A receiver of its own: a depth, a "location" or a "distance_deg" and
// "azimuth_deg" from a point, and an optional "instrument".
Receiver ReadReceiver(const Json& json, const Location& from,
                      const std::filesystem::path& directory);

// A "receivers" section.
std::vector<Receiver> ReadReceivers(const Json& json,
                                    const std::filesystem::path& directory);

}  // namespace YSpec::InputSections
