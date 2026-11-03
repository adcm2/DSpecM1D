// yspec-snapshot: the wavefield at chosen times, on a grid around each source,
// for a spherically symmetric earth model.
//
//   yspec-snapshot input.json
//
// Reads the calculation from a JSON file (see YSpec/Input.hpp for its
// members), computes it and writes the snapshots (see YSpec/Run.hpp for the
// files), reporting progress on the standard error.

#include <YSpec/YSpec.hpp>
#include <filesystem>
#include <format>
#include <iostream>
#include <string>

#include "Program.hpp"

namespace {

const auto kUsage =
    "usage: yspec-snapshot input.json\n"
    "\n"
    "Computes snapshots of the wavefield for the calculation described in\n"
    "input.json (see README.md), at the times it lists and on a grid of\n"
    "depths, distances and azimuths around each source, and writes one file\n"
    "per time, <prefix>.1, <prefix>.2, ... (<prefix>.<source>.1, ... for\n"
    "several sources), each line holding a point of the grid, its latitude\n"
    "and longitude, and the up, away and transverse components in SI units.\n"
    "The number of threads is set by OMP_NUM_THREADS.\n";

std::string Run(const std::filesystem::path& path) {
  using Program::Plural;
  const auto input = YSpec::ReadSnapshotInput(path);
  const auto& calculation = input.calculation;
  const auto& grid = input.grid;
  std::cerr << std::format(
      "yspec-snapshot: model {}, degrees {} to {}, {} to {} mHz, {}, {}, "
      "{} x {} x {}, {}\n",
      input.model.string(), calculation.degrees.min, calculation.degrees.max,
      calculation.band.min * 1e3, calculation.band.max * 1e3,
      Plural(input.sources.size(), "source"),
      Plural(input.times.size(), "time"), Plural(grid.depths.size(), "depth"),
      Plural(grid.distances.size(), "distance"),
      Plural(grid.azimuths.size(), "azimuth"),
      Plural(Program::Threads(), "thread"));
  const auto snapshots =
      YSpec::ComputeSnapshots(input, Program::Reporter("yspec-snapshot"));
  YSpec::WriteSnapshots(snapshots, input);
  return std::format(
      "{} with prefix {}",
      Plural(snapshots.Sources() * snapshots.Times().size(), "file"),
      input.prefix.string());
}

}  // namespace

int main(int argc, char** argv) {
  return Program::Main(argc, argv, "yspec-snapshot", kUsage, Run);
}
