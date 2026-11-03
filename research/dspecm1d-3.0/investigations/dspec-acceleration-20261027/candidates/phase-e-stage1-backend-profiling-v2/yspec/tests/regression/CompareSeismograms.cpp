// Compares two seismogram files, columns t Z N E: the largest difference in
// each component, relative to the peak of that component in the reference,
// must be below a tolerance.
//
//   CompareSeismograms output reference tolerance

#include <algorithm>
#include <array>
#include <cmath>
#include <format>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace {

std::vector<std::array<double, 4>> Read(const std::string& path) {
  auto file = std::ifstream(path);
  if (!file) throw std::runtime_error("cannot read " + path);
  auto rows = std::vector<std::array<double, 4>>{};
  for (auto row = std::array<double, 4>{};
       file >> row[0] >> row[1] >> row[2] >> row[3];) {
    rows.push_back(row);
  }
  return rows;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 4) {
    std::cerr << "usage: CompareSeismograms output reference tolerance\n";
    return 2;
  }
  try {
    const auto output = Read(argv[1]);
    const auto reference = Read(argv[2]);
    const auto tolerance = std::stod(argv[3]);
    if (output.size() != reference.size()) {
      std::cerr << std::format("{}: {} samples, the reference has {}\n",
                               argv[1], output.size(), reference.size());
      return 1;
    }
    auto pass = true;
    auto line = std::format("{}:", argv[1]);
    for (auto c : {1, 2, 3}) {
      auto peak = 0.0;
      auto difference = 0.0;
      for (auto i = std::size_t{0}; i < output.size(); ++i) {
        peak = std::max(peak, std::abs(reference[i][c]));
        difference =
            std::max(difference, std::abs(output[i][c] - reference[i][c]));
      }
      const auto relative = peak > 0 ? difference / peak : difference;
      line += std::format(" {:.2e}", relative);
      pass = pass && relative <= tolerance;
    }
    std::cout << line << (pass ? "" : "  FAILED") << "\n";
    return pass ? 0 : 1;
  } catch (const std::exception& e) {
    std::cerr << e.what() << "\n";
    return 1;
  }
}
