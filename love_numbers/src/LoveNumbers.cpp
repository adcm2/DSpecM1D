#include <DSpecM1D/LoveNumbers>

#include "RadialState.h"
#include "StaticSolve.h"

namespace DSpecM1D::LoveNumbers {

std::vector<DegreeResult>
calculate(const EarthModels::ModelInput<double> &model, const Config &config) {
  detail::RadialState radial_state(model, config);
  std::vector<DegreeResult> results;
  results.reserve(config.maximum_degree + 1);

  for (int degree = 0; degree <= config.maximum_degree; ++degree) {
    results.push_back(
        detail::solveStaticDegree(model, radial_state, degree));
  }

  return results;
}

}   // namespace DSpecM1D::LoveNumbers
