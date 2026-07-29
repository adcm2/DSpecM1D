#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <DSpecM1D/LoveNumbers>

namespace {

using DSpecM1D::LoveNumbers::Config;
using DSpecM1D::LoveNumbers::DegreeResult;

const char *usage =
    "usage: dspecm1d-love MODEL_FILE LMAX OUTPUT_FILE "
    "[POLYNOMIAL_ORDER] [MAXIMUM_RADIAL_STEP]";

int parseInteger(const char *text, const std::string &name) {
  try {
    std::size_t end = 0;
    const int value = std::stoi(text, &end);
    if (end != std::string(text).size()) {
      throw std::invalid_argument("");
    }
    return value;
  } catch (const std::exception &) {
    throw std::invalid_argument(name + " must be an integer.");
  }
}

double parseDouble(const char *text, const std::string &name) {
  try {
    std::size_t end = 0;
    const double value = std::stod(text, &end);
    if (end != std::string(text).size()) {
      throw std::invalid_argument("");
    }
    return value;
  } catch (const std::exception &) {
    throw std::invalid_argument(name + " must be a number.");
  }
}

void writeResults(std::ostream &output, const std::string &model_path,
                  const Config &config,
                  const std::vector<DegreeResult> &results) {
  output << std::scientific << std::setprecision(16);
  output << "# DSpecM1D static Love numbers\n"
         << "# model: " << model_path << '\n'
         << "# maximum_degree: " << config.maximum_degree << '\n'
         << "# polynomial_order: " << config.polynomial_order << '\n'
         << "# maximum_radial_step: " << config.maximum_radial_step << '\n'
         << "# degree_one_convention: P(a) = 0\n"
         << "# columns: l h_u k_u h_phi k_phi h_load k_load h_t k_t\n"
         << "# units: h_u,h_phi,h_load=m^3 kg^-1; "
            "k_u,k_phi,k_load=m^4 kg^-1 s^-2; "
            "h_t=s^2 m^-1; k_t=dimensionless\n";

  for (const DegreeResult &result : results) {
    output << result.degree << ' ' << result.h_u << ' ' << result.k_u
           << ' ' << result.h_phi << ' ' << result.k_phi << ' '
           << result.h_load() << ' ' << result.k_load() << ' '
           << result.h_t << ' ' << result.k_t << '\n';
  }
}

}   // namespace

int main(int argc, char **argv) {
  if (argc < 4 || argc > 6) {
    std::cerr << usage << '\n';
    return 2;
  }

  try {
    const std::string model_path = argv[1];
    const int maximum_degree = parseInteger(argv[2], "LMAX");
    const Config config{
        .maximum_degree = maximum_degree,
        .polynomial_order =
            argc >= 5 ? parseInteger(argv[4], "POLYNOMIAL_ORDER") : 6,
        .maximum_radial_step =
            argc >= 6
                ? parseDouble(argv[5], "MAXIMUM_RADIAL_STEP")
                : 0.01,
    };

    const EarthModels::ModelInput<double> model(model_path);
    const std::vector<DegreeResult> results =
        DSpecM1D::LoveNumbers::calculate(model, config);

    if (std::string(argv[3]) == "-") {
      writeResults(std::cout, model_path, config, results);
      if (!std::cout) {
        throw std::runtime_error("Could not write Love numbers.");
      }
    } else {
      std::ofstream output(argv[3]);
      if (!output) {
        throw std::runtime_error(
            "Could not open output file: " + std::string(argv[3]));
      }
      writeResults(output, model_path, config, results);
      if (!output) {
        throw std::runtime_error(
            "Could not write output file: " + std::string(argv[3]));
      }
    }
  } catch (const std::exception &error) {
    std::cerr << "dspecm1d-love: " << error.what() << '\n';
    return 1;
  }

  return 0;
}
