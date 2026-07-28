#include <algorithm>
#include <cfenv>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/SparseLU>

// This keeps EarthMesh's header definitions in one translation unit while
// exposing the private diagnostics beside the public calculation.
#include "../src/LoveNumbers.cpp"

namespace {

using DSpecM1D::LoveNumbers::Config;
using DSpecM1D::LoveNumbers::DegreeResult;
using DSpecM1D::LoveNumbers::calculate;
using DSpecM1D::LoveNumbers::detail::LoveDofMap;
using DSpecM1D::LoveNumbers::detail::RadialState;
using DSpecM1D::LoveNumbers::detail::assembleStaticOperator;
using DSpecM1D::LoveNumbers::detail::assembleStaticRightHandSides;

void printDiagnostics(RadialState &state, int degree) {
  const LoveDofMap degrees_of_freedom(state.mesh(), degree);
  const Eigen::SparseMatrix<double> matrix =
      assembleStaticOperator(state, degrees_of_freedom, degree);
  const Eigen::MatrixXd right_hand_sides =
      assembleStaticRightHandSides(state, degrees_of_freedom, degree);
  for (int outer = 0; outer < matrix.outerSize(); ++outer) {
    for (Eigen::SparseMatrix<double>::InnerIterator entry(matrix, outer);
         entry; ++entry) {
      if (!std::isfinite(entry.value())) {
        throw std::runtime_error(
            "controlled matrix contains a non-finite entry");
      }
    }
  }
  if (!right_hand_sides.allFinite()) {
    throw std::runtime_error(
        "controlled right-hand side contains a non-finite entry");
  }

  Eigen::SparseLU<Eigen::SparseMatrix<double>> factorization;
  factorization.analyzePattern(matrix);
  factorization.factorize(matrix);
  if (factorization.info() != Eigen::Success) {
    throw std::runtime_error("controlled matrix factorization failed");
  }
  const Eigen::MatrixXd solutions =
      factorization.solve(right_hand_sides);
  if (factorization.info() != Eigen::Success) {
    throw std::runtime_error("controlled three-column solve failed");
  }
  if (!solutions.allFinite()) {
    throw std::runtime_error(
        "controlled solution contains a non-finite entry");
  }

  std::cerr << "diagnostic " << degree << ' '
            << degrees_of_freedom.size();
  for (int column = 0; column < right_hand_sides.cols(); ++column) {
    const double residual =
        (matrix * solutions.col(column) -
         right_hand_sides.col(column))
            .norm() /
        std::max(1.0, right_hand_sides.col(column).norm());
    std::cerr << ' ' << residual;
  }
  std::cerr << '\n';
}

}   // namespace

int main(int argc, char **argv) {
  if (argc < 4) {
    std::cerr
        << "usage: controlled_love_numbers model elements degree...\n";
    return 2;
  }

  const int requested_elements = std::stoi(argv[2]);
  if (requested_elements < 1) {
    throw std::invalid_argument("elements must be positive");
  }

  std::vector<int> degrees;
  int maximum_degree = 0;
  for (int argument = 3; argument < argc; ++argument) {
    const int degree = std::stoi(argv[argument]);
    if (degree < 0) {
      throw std::invalid_argument("degrees must be non-negative");
    }
    degrees.push_back(degree);
    maximum_degree = std::max(maximum_degree, degree);
  }

  std::feclearexcept(FE_INVALID | FE_DIVBYZERO);

  const EarthModels::ModelInput<double> model(argv[1]);
  const Config config{
      .maximum_degree = maximum_degree,
      .polynomial_order = 4,
      .maximum_radial_step = 1.0 / requested_elements,
  };
  const std::vector<DegreeResult> results = calculate(model, config);
  RadialState state(model, config);

  const int elements = state.mesh().NE();
  const int radial_nodes =
      elements * (state.mesh().NN() - 1) + 1;
  const double surface_radius =
      state.surfaceRadius() * model.LengthNorm();
  const double surface_gravity =
      state.surfaceGravity() * model.AccelerationNorm();
  const double legacy_surface_gravity =
      state.meshModel().Gravity(elements - 1, state.mesh().NN() - 1) *
      model.AccelerationNorm();

  std::cout << std::scientific << std::setprecision(17);
  std::cerr << std::scientific << std::setprecision(17);
  std::cerr << "metadata " << surface_radius << ' ' << surface_gravity
            << ' ' << elements << ' ' << radial_nodes << '\n';
  std::cerr << "legacy_surface_gravity " << legacy_surface_gravity
            << '\n';

  std::vector<int> layer_elements(model.NumberOfLayers(), 0);
  for (int element = 0; element < state.mesh().NE(); ++element) {
    ++layer_elements.at(state.mesh().LayerNumber(element));
  }
  std::cerr << "layers " << layer_elements.size();
  for (const int count : layer_elements) {
    std::cerr << ' ' << count;
  }
  std::cerr << '\n';

  for (int layer = 0; layer < model.NumberOfLayers(); ++layer) {
    if (!model.IsFluid(layer)) {
      continue;
    }

    int first_element = -1;
    int last_element = -1;
    double minimum_derivative =
        std::numeric_limits<double>::infinity();
    double maximum_derivative =
        -std::numeric_limits<double>::infinity();
    for (int element = 0; element < state.mesh().NE(); ++element) {
      if (state.mesh().LayerNumber(element) != layer) {
        continue;
      }
      if (first_element < 0) {
        first_element = element;
      }
      last_element = element;
      for (int node = 0; node < state.mesh().NN(); ++node) {
        const double derivative =
            state.densityDerivative(element, node) *
            model.DensityNorm() / model.LengthNorm();
        minimum_derivative = std::min(minimum_derivative, derivative);
        maximum_derivative = std::max(maximum_derivative, derivative);
      }
    }

    const double lower_density =
        state.meshModel().Density(first_element, 0) *
        model.DensityNorm();
    const double upper_density =
        state.meshModel().Density(
            last_element, state.mesh().NN() - 1) *
        model.DensityNorm();
    std::cerr << "fluid " << layer << ' ' << lower_density << ' '
              << upper_density << ' ' << minimum_derivative << ' '
              << maximum_derivative << '\n';
  }

  for (const int degree : degrees) {
    printDiagnostics(state, degree);
    const DegreeResult &result = results.at(degree);
    std::cout << degree << ' ' << result.h_u << ' ' << result.k_u
              << ' ' << result.h_phi << ' ' << result.k_phi << ' '
              << result.h_load() << ' ' << result.k_load() << ' '
              << result.h_t << ' ' << result.k_t << '\n';
  }

  const int exceptions = std::fetestexcept(FE_INVALID | FE_DIVBYZERO);
  const int invalid = (exceptions & FE_INVALID) != 0;
  const int divide_by_zero = (exceptions & FE_DIVBYZERO) != 0;
  std::cerr << "floating_point_flags " << invalid << ' '
            << divide_by_zero << '\n';
  if (invalid || divide_by_zero) {
    throw std::runtime_error(
        "controlled calculation raised an invalid or divide-by-zero "
        "floating-point exception");
  }
}
