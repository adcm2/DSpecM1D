#include <iostream>
#include <cmath>
#include <functional>
#include <filesystem>
#include <fstream>
#include <Eigen/Core>
#include <new_coupling/All>
#include <PlanetaryModel/All>
#include "toroidal_clean.h"
#include "toroidal_bench.h"

int
main() {

  // parameters
  int nfreq;
  int nsolid = 1;           // mantle mode
  double maxstep = 0.07;    // max step
  double sigshift = 0.01;   // shift for eigensolver
  double twopi = 2.0 * 3.1415926535897932;
  int pn = 8;   // order of spectral element method

  ///////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////
  // paths
  std::string pathtopremeig =
      "/space/adcm2/mineos-1.0.2/OUTPUT/eprem_noocean_S_IC";
  std::string pathtoprem =
      "/space/adcm2/mineos-1.0.2/DEMO/models/prem.200.no.noatten.txt";

  // models
  // auto prem = TestTools::EarthModel(pathtoprem);
  auto prem = EarthModels::ModelInput(pathtoprem);

  // our calculation
  std::vector<std::vector<double>> vec_diff;
  for (int lv = 1; lv < 11; ++lv) {

    // mineos result
    std::string path = pathtopremeig + "_" + std::to_string(lv);
    ModeCoupling::modecattoroidal modes_prem_lv(path, prem);

    // our spectral element method
    Toroidal::spectral_element_planet test_planet(prem, maxstep, pn, lv,
                                                  nsolid);

    // calculate eigenfrequencies
    int nfreq2 = modes_prem_lv.NumberOfModes() + (lv == 1);
    test_planet.CalculateEigenfrequencies(nfreq2, sigshift);

    Eigen::VectorXd vec_eig = test_planet.efrequencies_gen() / twopi;

    // storage vector:
    std::vector<double> vec_pb;
    for (int i = 0; i < modes_prem_lv.NumberOfModes(); ++i) {
      auto wf = modes_prem_lv.singlemodep(i).w();
      double tmp = std::abs(wf - vec_eig(i + (lv == 1))) / wf;
      vec_pb.push_back(tmp);
    }
    vec_diff.push_back(vec_pb);
  }

  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  // output
  for (int i = 0; i < vec_diff.size(); ++i) {
    std::string path1 = "./work/toroidal/mineos_bench.out";
    std::string pathtofile = path1 + "_" + std::to_string(i + 1);
    auto file = std::ofstream(pathtofile);
    for (int j = 0; j < vec_diff[i].size(); ++j) {
      file << std::setprecision(16) << j + (i == 0) << ";" << vec_diff[i][j]
           << std::endl;
    }
    file.close();
  }
  return 0;
}