#include <iostream>
#include <cmath>
#include <functional>
#include <filesystem>
#include <fstream>
// #include <new_coupling/All>
#include <PlanetaryModel/All>

int
main() {
  // input
  std::string pathtoprem =
      "/space/adcm2/mineos-1.0.2/DEMO/models/prem.200.no.noatten.txt";
  auto prem0 = EarthModels::ModelInput(pathtoprem);

  // output
  std::string pathtofile = "./work/toroidal/prem_pert.out";
  auto file = std::ofstream(pathtofile);
  std::size_t numnodes = prem0.LayerUpperIndex(prem0.NumberOfLayers() - 1) + 1;

  // make file
  file << "prem.pert.200" << std::endl;
  file << 0 << " " << std::setprecision(3) << 1.0 << " " << 1 << " " << 2
       << std::endl;
  file << numnodes << " " << prem0.LayerUpperIndex(0) + 1 << " "
       << prem0.LayerUpperIndex(1) + 1 << "\n";

  auto vec_radii = prem0.LayerRadii();
  for (int i = 0; i < vec_radii.size(); ++i) {
    double mf = 1.0;
    if (i == 3) {
      mf *= 1.02;
    }
    for (int j = 0; j < vec_radii[i].size(); ++j) {
      double cr = vec_radii[i][j];
      double rho = prem0.Density(i)(cr) * prem0.DensityNorm();
      double vpv = prem0.VPV(i)(cr) * prem0.VelocityNorm();
      double vsv = prem0.VSV(i)(cr) * prem0.VelocityNorm() * mf;
      double vph = prem0.VPH(i)(cr) * prem0.VelocityNorm();
      double vsh = prem0.VSH(i)(cr) * prem0.VelocityNorm() * mf;
      cr *= prem0.LengthNorm();
      file << std::setprecision(16) << cr << " " << rho << " " << vpv << " "
           << vsv << " " << 0.0 << " " << 0.0 << " " << vph << " " << vsh << " "
           << prem0.Eta(i)(cr) << std::endl;
    }
  }
  file.close();

  return 0;
};