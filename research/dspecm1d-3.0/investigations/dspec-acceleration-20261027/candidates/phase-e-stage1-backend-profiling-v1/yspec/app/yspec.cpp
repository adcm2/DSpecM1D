// yspec: synthetic seismograms for a spherically symmetric earth model.
//
//   yspec input.json
//
// Reads the calculation from a JSON file (see YSpec/Input.hpp for its
// members), computes it and writes the seismograms (see YSpec/Run.hpp for the
// files), reporting progress on the standard error.

#include <YSpec/YSpec.hpp>
#include <filesystem>
#include <format>
#include <iostream>
#include <string>

#include "Program.hpp"

namespace {

const auto kUsage =
    "usage: yspec input.json\n"
    "\n"
    "Computes synthetic seismograms for the calculation described in\n"
    "input.json (see README.md) and writes one file per receiver,\n"
    "<prefix>.1, <prefix>.2, ... (<prefix>.<source>.1, ... for several\n"
    "sources), each line holding the time (s) and the vertical, north and\n"
    "east components in SI units, or in those of the receiver's instrument.\n"
    "The number of threads is set by OMP_NUM_THREADS.\n";

std::string Run(const std::filesystem::path& path) {
  using Program::Plural;
  const auto input = YSpec::ReadSeismogramInput(path);
  const auto& calculation = input.calculation;
  std::cerr << std::format(
      "yspec: model {}, degrees {} to {}, {} to {} mHz, {}, {}, {}\n",
      input.model.string(), calculation.degrees.min, calculation.degrees.max,
      calculation.band.min * 1e3, calculation.band.max * 1e3,
      Plural(input.sources.size(), "source"),
      Plural(input.receivers.size(), "receiver"),
      Plural(Program::Threads(), "thread"));
  const auto seismograms =
      YSpec::ComputeSeismograms(input, Program::Reporter("yspec"));
  YSpec::WriteSeismograms(seismograms, input);
  return std::format(
      "{} with prefix {}",
      Plural(seismograms.Sources() * seismograms.Receivers(), "file"),
      input.prefix.string());
}

}  // namespace

int main(int argc, char** argv) {
  return Program::Main(argc, argv, "yspec", kUsage, Run);
}
