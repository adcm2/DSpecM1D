#include <iostream>
#include <cmath>
#include <functional>
#include <GaussQuad/All>
#include <Interpolation/All>
#include <filesystem>
#include <fstream>
#include <Eigen/Core>
// #include <Gravitational_Field/Timer>
#include <PlanetaryModel/All>
#include <new_coupling/All>
#include "toroidal.h"
#include "toroidal_clean.h"
#include "toroidal_gal.h"
#include "toroidal_bench.h"

using namespace Toroidal;

class spherical_model {
public:
  // default
  spherical_model() {};

  // spherical_model(double, double, double, double, double);
  spherical_model(double density, double mu, double rad)
      : _length_norm{1.0}, _mass_norm{1.0}, _time_norm{1.0},
        _vec_layer_boundaries{{0.0, rad}}, _vec_layer_densities{{density}},
        _vec_layer_mu{{mu}} {};

  // return functions
  double LengthNorm() const { return _length_norm; };
  double MassNorm() const { return _mass_norm; };
  double TimeNorm() const { return _time_norm; };
  double DensityNorm() const {
    return _mass_norm / std::pow(_length_norm, 3.0);
  };
  double InertiaNorm() const {
    return _mass_norm * std::pow(_length_norm, 2.0);
  };
  double VelocityNorm() const { return _length_norm / _time_norm; };
  double AccelerationNorm() const {
    return _length_norm / std::pow(_time_norm, 2.0);
  };
  double ForceNorm() const {
    return _mass_norm * _length_norm / std::pow(_time_norm, 2.0);
  };
  double StressNorm() const {
    return _mass_norm / (std::pow(_time_norm, 2.0) * _length_norm);
  };
  int NumberOfLayers() const { return _number_of_layers; };
  auto LowerRadius(int i) const {
    assert(i < _number_of_layers && "Not in model");
    return _vec_layer_boundaries[i];
  };
  auto UpperRadius(int i) const {
    assert(i < _number_of_layers && "Not in model");
    return _vec_layer_boundaries[i + 1];
  };
  auto OuterRadius() const { return _vec_layer_boundaries.back(); };
  auto Density(int i) const {
    assert(i < _number_of_layers && "Not in model");
    auto densitylambda = [intdensity = this->_vec_layer_densities[i]](
                             double eval_rad) { return intdensity; };
    return densitylambda;
  };
  auto Mu(int i) const {
    assert(i < _number_of_layers && "Not in model");
    auto mulambda = [intmu = this->_vec_layer_mu[i]](double eval_rad) {
      return intmu;
    };
    return mulambda;
  };
  auto L(int i) const { return Mu(i); };
  auto N(int i) const { return Mu(i); };
  auto VS(int i) const {
    assert(i < _number_of_layers && "Not in model");
    auto mulambda = [intmu = this->_vec_layer_mu[i],
                     intrho = this->_vec_layer_densities[i]](double eval_rad) {
      return std::sqrt(intmu / intrho);
    };
    return mulambda;
  };

  bool IsSolid(int i) const { return true; };
  bool IsFluid(int i) const { return false; };

private:
  double _length_norm, _time_norm, _mass_norm;
  int _number_of_layers = 1;
  std::vector<double> _vec_layer_boundaries, _vec_layer_densities,
      _vec_layer_mu;

  // // general constructor
  // spherical_model(std::vector<double> &, std::vector<double> &, double,
  // double,
  //                 double);
};

int
main() {
  // testing model string
  std::vector<double> x_radii{0.0, 1.0};
  std::vector<double> mu_vec{1.0};
  std::vector<double> rho_vec{1.0};
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

  Timer timer;
  // timer.start();
  model_string def_str(x_radii, mu_vec, rho_vec);
  // timer.stop("Defining string");

  // timer.start();
  spectral_element def_spec(def_str, maxstep, 6, 1);
  // timer.stop("Defining spectral element scheme");

  // timer.start();
  def_spec.CalculateEigenfrequencies(nfreq);
  // timer.stop("Calculating eigenfrequncies");

  // testing newton method
  eigenfunction_catalogue eig_exact(10, 1.0, 1.0, 1.0);
  auto vec_eig = eig_exact.w();

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

  // new constructor with planetary model
  EarthModels::PREM<double, int> myprem;
  GaussQuad::Quadrature1D<double> q =
      GaussQuad::GaussLobattoLegendreQuadrature1D<double>(6);

  std::cout << "q.N(): " << q.N() << "\n\n";
  spherical_model homosphere(1.0, 1.0, 1.0);
  double d_nd = 13088.5 / prem.DensityNorm();
  double vs_nd = 3500.0 / prem.VelocityNorm();

  double mu_nd = vs_nd * vs_nd * d_nd;

  spherical_model innercore(d_nd, mu_nd, 0.191729);

  spectral_element test_planet(prem, maxstep, 6, lval);
  spectral_element_planet test_planet1(homosphere, 0.1, 9, lval);
  spectral_element_planet test_planet2(prem, maxstep, 8, lval, nsolid);

  test_planet1.CalculateEigenfrequencies(nfreq, sigshift);
  test_planet2.CalculateEigenfrequencies(nfreq, sigshift);

  Eigen::VectorXcd v_eig = test_planet1.efrequencies();
  Eigen::VectorXcd v_eiggen = test_planet1.efrequencies_gen();

  auto radmesh = test_planet2.mesh();

  int N = radmesh.NN();

  // modes and output containers
  ModeCoupling::modecattoroidal modes_prem(pathtopremeig, prem);
  auto vec_minval = modes_prem.AllModes(radmesh);
  auto vec_eigd = modes_prem.AllModesDeriv(radmesh);
  auto vec_eigval = test_planet2.evectors_std();
  auto vec_eigderiv = test_planet2.evector_deriv();
  auto vec_traction = test_planet2.traction_std(prem);

  // augment calc
  test_planet2.augment_basis_calculate();
  auto vec_augbasis = test_planet2.augment_basis();

  // test out gal
  gal_gen galtest(test_planet2, prem, prem, false);

  // eigenvalues and vectors
  // std::cout << "Homosphere: " << v_eiggen << "\n\n";
  // for (auto idx : vec_eig) {
  //   std::cout << idx << "\n";
  // }

  // bench
  Toroidal_Bench::sphere_bench testbench(lval, nfreq);
  auto wbench = testbench.w();
  // for (auto idx : wbench) {
  //   std::cout << idx << "\n";
  // }
  std::cout << "Homosphere, l = " << lval << ", up to n = " << nfreq << ":\n";
  for (int idx = 0; idx < nfreq; ++idx) {
    int idx2 = idx;
    if (lval == 1) {
      idx2 = idx;
    }
    std::cout << "Rel. diff: "
              << std::abs(wbench[idx] - v_eiggen(idx2).real()) / wbench[idx] *
                     100
              << "\n";
  }

  // proper benchmarks
  std::vector<std::vector<double>> vec_relerr(6, std::vector<double>(10, 0.0));
  for (int idxp = 3; idxp < 9; ++idxp) {
    for (int idxi = 0; idxi < 10; ++idxi) {
      double max_i = std::pow(
          10.0, -0.5 - 2.5 * std::log(idxp) * static_cast<double>(idxi) /
                           static_cast<double>(idxp * idxp + 1));
      spectral_element_planet bench_planet(homosphere, max_i, idxp + 1, lval);
      bench_planet.CalculateEigenfrequencies(nfreq, sigshift);
      Eigen::VectorXcd eig_bench = bench_planet.efrequencies_gen();

      // calculate relative error (average)
      double tmp = 0.0;
      double multfact = 0.1;
      int idxmin = 0;
      if (lval == 1) {
        idxmin = 1;
        multfact = 1.0 / 9.0;
      }
      for (int idxf = idxmin; idxf < nfreq; ++idxf) {
        tmp += std::abs(wbench[idxf] - eig_bench(idxf).real()) / wbench[idxf];
      }
      vec_relerr[idxp - 3][idxi] = tmp * multfact;
    }
  }
  // for (int idxp = 3; idxp < 9; ++idxp) {
  //   for (int idxi = 0; idxi < 10; ++idxi) {
  //     std::cout << "p: " << idxp << ", h: " << -1.0 - idxi / 5.0
  //               << ", err: " << vec_relerr[idxp - 3][idxi] << "\n";
  //   }
  // }
  // for (auto &idx : vec_relerr) {
  //   std::cout << idx << "\n";
  // }
  // double freqnorm = 1.0 / (prem.TimeNorm() * 2.0 * 3.1415926535);
  // std::cout << std::setprecision(10)
  //           << test_planet2.efrequencies() / (2.0 * 3.1415926535) << "\n\n";
  // std::cout << std::setprecision(10) << galtest.efrequencies() * freqnorm
  //           << "\n\n";
  // std::cout << modes_prem.singlemodep(0).w() << "\n\n";
  // std::cout << std::setprecision(3) << galtest.evectors() << "\n\n";
  // testing the derivative by approximate method
  // std::vector<std::vector<std::vector<double>>> vec_eigd2;
  // for (int idxe = 0; idxe < radmesh.NE(); ++idxe) {
  //   std::vector<std::vector<double>> vec_out;
  //   for (int idxq = 0; idxq < radmesh.NumNodesWithinElement() - 1; ++idxq) {
  //     std::vector<double> vec_x;
  //     for (int idx = 0; idx < vec_eigval.back().back().size(); ++idx) {
  //       double tmp =
  //           (vec_eigval[idxe][idxq + 1][idx] - vec_eigval[idxe][idxq][idx]) /
  //           (radmesh.NodeRadius(idxe, idxq + 1) -
  //            radmesh.NodeRadius(idxe, idxq));
  //       vec_x.push_back(tmp);
  //     }
  //     vec_out.push_back(vec_x);
  //   }
  //   auto vec_tmp = vec_out.back();
  //   vec_out.push_back(vec_tmp);
  //   vec_eigd2.push_back(vec_out);
  // }
  ///////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  // output
  // outputting result
  std::string pathtofile = "./work/toroidal/efunction.out";
  auto file = std::ofstream(pathtofile);

  for (int idxe = 0; idxe < radmesh.NE(); ++idxe) {
    for (int idxq = 0; idxq < N; ++idxq) {
      file << std::setprecision(16) << radmesh.NodeRadius(idxe, idxq);
      for (int idx = 0; idx < vec_eigval.back().back().size(); ++idx) {
        file << ";" << vec_eigval[idxe][idxq][idx];
      }
      file << std::endl;
    }
  }
  file.close();

  std::string pathtofile2 = "./work/toroidal/efunction_mineos.out";
  auto file2 = std::ofstream(pathtofile2);
  for (int idxe = 0; idxe < radmesh.NE(); ++idxe) {
    for (int idxq = 0; idxq < N; ++idxq) {
      file2 << std::setprecision(16) << radmesh.NodeRadius(idxe, idxq);
      for (int idx = 0; idx < modes_prem.NumberOfModes(); ++idx) {
        file2 << ";" << vec_minval[idxe][idxq][idx];
      }
      file2 << std::endl;
    }
  }
  file2.close();

  std::string pathtofile3 = "./work/toroidal/efunctiond.out";
  auto file3 = std::ofstream(pathtofile3);

  for (int idxe = 0; idxe < radmesh.NE(); ++idxe) {
    for (int idxq = 0; idxq < N; ++idxq) {
      file3 << std::setprecision(16) << radmesh.NodeRadius(idxe, idxq);
      for (int idx = 0; idx < vec_eigderiv.back().back().size(); ++idx) {
        file3 << ";" << vec_eigderiv[idxe][idxq][idx];
      }
      file3 << std::endl;
    }
  }
  file3.close();

  std::string pathtofile4 = "./work/toroidal/efunctiond2.out";
  auto file4 = std::ofstream(pathtofile4);

  for (int idxe = 0; idxe < radmesh.NE(); ++idxe) {
    for (int idxq = 0; idxq < N; ++idxq) {
      file4 << std::setprecision(16) << radmesh.NodeRadius(idxe, idxq);
      for (int idx = 0; idx < vec_eigd.back().back().size(); ++idx) {
        file4 << ";" << vec_eigd[idxe][idxq][idx];
      }
      file4 << std::endl;
    }
  }
  file4.close();

  std::string pathtofile5 = "./work/toroidal/traction.out";
  auto file5 = std::ofstream(pathtofile5);

  for (int idxe = 0; idxe < radmesh.NE(); ++idxe) {
    for (int idxq = 0; idxq < N; ++idxq) {
      file5 << std::setprecision(16) << radmesh.NodeRadius(idxe, idxq);
      for (int idx = 0; idx < vec_traction.back().back().size(); ++idx) {
        file5 << ";" << vec_traction[idxe][idxq][idx];
      }
      file5 << std::endl;
    }
  }
  file5.close();

  std::string pathtofile6 = "./work/toroidal/augbasis.out";
  auto file6 = std::ofstream(pathtofile6);

  for (int idxe = 0; idxe < radmesh.NE(); ++idxe) {
    for (int idxq = 0; idxq < N; ++idxq) {
      file6 << std::setprecision(16) << radmesh.NodeRadius(idxe, idxq);
      for (int idx = 0; idx < vec_augbasis.back().back().size(); ++idx) {
        file6 << ";" << vec_augbasis[idxe][idxq][idx];
      }
      file6 << std::endl;
    }
  }
  file6.close();

  std::string pathtofile7 = "./work/toroidal/benchmark.out";
  auto file7 = std::ofstream(pathtofile7);
  for (int idxi = 0; idxi < 10; ++idxi) {

    file7 << std::setprecision(16) << idxi;
    for (int idxp = 3; idxp < 9; ++idxp) {
      double max_i = std::pow(
          10.0, -0.5 - 2.5 * std::log(idxp) * static_cast<double>(idxi) /
                           static_cast<double>(idxp * idxp + 1));
      file7 << ";" << max_i << ";" << vec_relerr[idxp - 3][idxi];
    }
    file7 << std::endl;
  }
  file7.close();

  return 0;
}