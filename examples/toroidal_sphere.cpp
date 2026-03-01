#include <iostream>
#include <cmath>
#include <functional>
#include <filesystem>
#include <fstream>
#include <Eigen/Core>
#include <new_coupling/All>
#include "toroidal_clean.h"
#include "toroidal_bench.h"

int
main() {

  // parameters
  int nfreq = 5;
  int nsolid = 1;
  double maxstep;
  std::cout << "Enter maxstep: \n";
  std::cin >> maxstep;
  double sigshift;
  std::cout << "Enter shift:\n";
  std::cin >> sigshift;
  int lval;
  std::cout << "Enter l:\n";
  std::cin >> lval;
  std::cout << "Enter # of frequencies:\n";
  std::cin >> nfreq;
  std::cout << "Enter solid layer:\n";
  std::cin >> nsolid;

  // homosphere
  Toroidal_Bench::spherical_model homosphere(1.0, 1.0, 1.0);

  // exact solution
  Toroidal_Bench::sphere_bench testbench(lval, nfreq);
  auto wbench = testbench.w();

  // benchmark with different polynomial order and maxstep
  std::vector<std::vector<double>> vec_relerr(6, std::vector<double>(10, 0.0));
  std::vector<std::vector<double>> vec_step(6, std::vector<double>(10, 0.0));
  std::vector<int> vec_p{4, 5, 6, 7, 8, 9};

  // for averaging purposes
  double multfact = 1.0 / (9.0 + static_cast<double>((lval != 1)));
  int idxmin = (lval == 1);

  // do benchmarking
  for (int idxp = 0; idxp < 6; ++idxp) {
    int pval = vec_p[idxp];
    for (int idxi = 0; idxi < 10; ++idxi) {

      // values of maxstep
      double max_i = std::pow(
          10.0, -0.5 - 2.5 * std::log(pval) * static_cast<double>(idxi) /
                           static_cast<double>(pval * pval + 1));
      vec_step[idxp][idxi] = max_i;

      // planet for benchmark
      Toroidal::spectral_element_planet bench_planet(homosphere, max_i,
                                                     vec_p[idxp], lval);

      // calculate eigenfrequencies
      bench_planet.CalculateEigenfrequencies(nfreq, sigshift);
      Eigen::VectorXcd eig_bench = bench_planet.efrequencies_gen();

      // calculate relative error (average)
      double tmp = 0.0;

      // average
      for (int idxf = idxmin; idxf < nfreq; ++idxf) {
        tmp += std::abs(wbench[idxf] - eig_bench(idxf).real()) / wbench[idxf];
      }

      // save
      vec_relerr[idxp][idxi] = tmp * multfact;
    }
  }

  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  // output
  std::string pathtofile = "./work/toroidal/benchmark.out";
  auto file = std::ofstream(pathtofile);
  for (int idxi = 0; idxi < 10; ++idxi) {
    file << std::setprecision(16) << idxi;
    for (int idxp = 0; idxp < 6; ++idxp) {
      file << ";" << vec_step[idxp][idxi] << ";" << vec_relerr[idxp][idxi];
    }
    file << std::endl;
  }
  file.close();

  return 0;
}