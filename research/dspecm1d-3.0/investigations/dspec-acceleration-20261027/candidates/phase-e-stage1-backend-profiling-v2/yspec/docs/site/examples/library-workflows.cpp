#include <YSpec/YSpec.hpp>
#include <cmath>
#include <iostream>
#include <vector>

int main(int argc, char** argv) {
  if (argc != 2) return 2;
  const auto model = YSpec::Model::ReadDeck(argv[1]);
  auto calculation = YSpec::Calculation{};
  calculation.degrees = {.min = 0, .max = 40};
  calculation.band = {.min = 0.2e-3, .max = 5e-3};
  calculation.sampling = {.duration = 1200, .timeStep = 50};

  const auto sources = std::vector<YSpec::Source>{
      {.location = {0, 0}, .depth = 29.5e3, .moment = {.rr = 1e17}}};
  const auto receivers =
      std::vector<YSpec::Receiver>{{.location = {0, 30}, .depth = 3e3}};
  const auto strains =
      YSpec::ComputeStrains(model, calculation, sources, receivers);
  const auto radial_strain = strains.Trace(0, 0, YSpec::StrainComponent::rr);
  if (radial_strain.empty() || !std::isfinite(radial_strain[0])) return 1;

  const auto measurement =
      YSpec::Traveltime(model, calculation, sources[0], receivers[0],
                        YSpec::Along(YSpec::Vertical),
                        YSpec::Window{.start = 200, .end = 800, .taper = 30});
  const auto mesh = YSpec::RadialMesh{.order = 3, .maxElement = 100e3};
  const auto kernels = YSpec::ComputeRadialKernels(
      model, calculation, sources[0], measurement, mesh);
  const auto density = kernels[YSpec::RadialParameter::Rho];
  if (density.empty() || !std::isfinite(density.front())) return 1;
  std::cout << "radial strain at sample 0: " << radial_strain[0]
            << "\ndensity-kernel nodes: " << density.size() << '\n';
}
