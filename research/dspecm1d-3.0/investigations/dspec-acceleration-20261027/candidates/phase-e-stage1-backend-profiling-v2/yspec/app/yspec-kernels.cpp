// yspec-kernels: the sensitivity of a traveltime measured on a seismogram to
// spherically symmetric perturbations of the earth model, by the adjoint
// method.
//
//   yspec-kernels input.json
//
// Reads the calculation from a JSON file (see YSpec/Input.hpp for its
// members), computes the kernels and writes them (see YSpec/Run.hpp for the
// files), reporting progress on the standard error.

#include <YSpec/YSpec.hpp>
#include <filesystem>
#include <format>
#include <iostream>
#include <string>

#include "Program.hpp"

namespace {

const auto kUsage =
    "usage: yspec-kernels input.json\n"
    "\n"
    "Computes the kernels of cross-correlation traveltimes measured on the\n"
    "seismogram of one source at one receiver, as input.json describes (see\n"
    "README.md), by the adjoint method: radial kernels, in <prefix>.kernels,\n"
    "and, if asked, three-dimensional kernels on sections along and across\n"
    "the path (<prefix>.sections), on maps around it (<prefix>.maps) and\n"
    "through the earth (<prefix>.nc, a planetmodel model file), with the\n"
    "seismogram, the window and the adjoint source in <prefix>.measurement.\n"
    "With several traveltimes each name has its own after the prefix. The\n"
    "number of threads is set by OMP_NUM_THREADS.\n";

std::string Run(const std::filesystem::path& path) {
  using Program::Plural;
  const auto input = YSpec::ReadKernelInput(path);
  const auto& calculation = input.calculation;
  auto names = std::string{};
  for (const auto& m : input.measurements) {
    names += (names.empty() ? "" : ", ") + m.name + " on " + m.component;
  }
  std::cerr << std::format(
      "yspec-kernels: model {}, degrees {} to {}, {} to {} mHz, {}: {}, {}\n",
      input.model.string(), calculation.degrees.min, calculation.degrees.max,
      calculation.band.min * 1e3, calculation.band.max * 1e3,
      Plural(input.measurements.size(), "traveltime"), names,
      Plural(Program::Threads(), "thread"));
  const auto result = YSpec::RunKernels(path, [](std::string_view stage) {
    return Program::Reporter(std::format("yspec-kernels: {}", stage));
  });
  return std::format("the kernels of {} with prefix {}",
                     Plural(result.measurements.size(), "traveltime"),
                     input.prefix.string());
}

}  // namespace

int main(int argc, char** argv) {
  return Program::Main(argc, argv, "yspec-kernels", kUsage, Run);
}
