#pragma once

// What the programs share: their command line, their reports on the standard
// error, and their handling of errors.

#include <YSpec/Calculation.hpp>
#include <chrono>
#include <cstddef>
#include <exception>
#include <format>
#include <iostream>
#include <memory>
#include <string>
#include <string_view>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace Program {

inline std::string Plural(std::size_t n, std::string_view noun) {
  return std::format("{} {}{}", n, noun, n == 1 ? "" : "s");
}

inline std::size_t Threads() {
#ifdef _OPENMP
  return static_cast<std::size_t>(omp_get_max_threads());
#else
  return 1;
#endif
}

// Reports progress at every tenth of the frequencies.
inline YSpec::Progress Reporter(std::string name) {
  auto reported = std::make_shared<int>(0);
  return
      [name = std::move(name), reported](std::size_t done, std::size_t total) {
        const auto percent = static_cast<int>(100 * done / total);
        if (percent >= *reported + 10 || done == total) {
          *reported = percent - percent % 10;
          std::cerr << std::format("{}: {:3}% of {} frequencies\n", name,
                                   percent, total);
        }
      };
}

// The program's main: handles --help, --version, a wrong number of arguments
// and errors, and reports the time taken. run(path) does the work, given the
// input file, and returns what was written.
template <class Run>
int Main(int argc, char** argv, std::string_view name, std::string_view usage,
         Run&& run) {
  const auto argument = argc == 2 ? std::string(argv[1]) : std::string{};
  if (argument == "--help" || argument == "-h") {
    std::cout << usage;
    return 0;
  }
  if (argument == "--version") {
    std::cout << name << " " << YSPEC_VERSION << "\n";
    return 0;
  }
  if (argc != 2) {
    std::cerr << usage;
    return 2;
  }
  try {
    const auto start = std::chrono::steady_clock::now();
    const auto written = run(argument);
    const auto seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
            .count();
    std::cerr << std::format("{}: wrote {} in {:.2f} s\n", name, written,
                             seconds);
  } catch (const std::exception& e) {
    std::cerr << name << ": " << e.what() << "\n";
    return 1;
  }
  return 0;
}

}  // namespace Program
