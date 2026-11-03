// Synthetic seismograms from C++: the example in README.md, which the test
// suite builds and runs.

#include <YSpec/YSpec.hpp>
#include <cmath>
#include <iostream>
#include <vector>

int main() {
  const auto model = YSpec::Model::ReadDeck(YSPEC_EXAMPLES_DIR "/prem.200");

  auto calculation = YSpec::Calculation{};
  calculation.degrees = {.min = 0, .max = 200};
  calculation.band = {.min = 0.2e-3, .max = 15e-3};           // Hz
  calculation.sampling = {.duration = 3600, .timeStep = 20};  // s

  const auto sources = std::vector<YSpec::Source>{
      {.location = {0, 0}, .depth = 29.5e3, .moment = {.rr = 1e18}}};
  // At 3 km depth: the sea floor of this model.
  const auto receivers = std::vector<YSpec::Receiver>{
      {.location = {0, 30}, .depth = 3e3}, {.location = {0, 60}, .depth = 3e3}};

  const auto seismograms =
      YSpec::ComputeSeismograms(model, calculation, sources, receivers);
  const auto z = seismograms.Trace(0, 0, YSpec::Vertical);

  std::cout << "vertical displacement at 30 degrees, 2000 s: " << z[100]
            << " m\n";
  return std::isfinite(z[100]) && z[100] != 0 ? 0 : 1;
}
