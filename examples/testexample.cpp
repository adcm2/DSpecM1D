
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <FFTWpp/Ranges>
#include <GSHTrans/All>
#include <GaussQuad/All>
#include <new_coupling/All>
#include <Interpolation/All>
#include <PlanetaryModel/All>
#include <TomographyModels/All>
// #include <Gravitational_Field/All>
// #include <Gravitational_Field/Test>
// #include <Gravitational_Field/Timer>
#include <algorithm>
#include <chrono>
#include <complex>
#include <concepts>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <random>
#include <vector>

auto
modescalarproduct(ModeCoupling::mineos_eigenfunction_continuous &mode1,
                  ModeCoupling::mineos_eigenfunction_continuous &mode2,
                  EarthModels::ModelInput<double, int> &inp_model) {
  // Try to use integrate in gauss quad
  auto q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(4);

  double totint = 0.0;
  double densitynorm = 5515.0;
  double radiusnorm = inp_model.OuterRadius() * inp_model.LengthNorm();
  //    std::cout << radiusnorm << "\n\n";
  double pi_db = 3.14159265358979;
  double bigg_db = 6.6723 * std::pow(10.0, -11.0);
  double velocitynorm = radiusnorm / std::sqrt(pi_db * bigg_db * densitynorm);
  double frequencynorm = velocitynorm / radiusnorm;

  for (int idxlayer = 0; idxlayer < inp_model.NumberOfLayers(); ++idxlayer) {
    auto vec_radii = inp_model.LayerRadii(idxlayer);

    for (int idxpl = 0; idxpl < vec_radii.size() - 1; ++idxpl) {
      //  if (idxlayer == inp_model.NumberOfLayers() - 1) {
      //     std::cout << vec_radii[idxpl] << " " << vec_radii[idxpl + 1]
      //               << "\n";
      //  }

      auto fun = [&inp_model, &mode1, &mode2, &idxlayer, &vec_radii,
                  &idxpl](double x) {
        double x1 = vec_radii[idxpl];
        double x2 = vec_radii[idxpl + 1];

        double delta = (x2 - x1) / 2.0;
        double xpl = (x1 + x2) / 2.0;
        double xscale = delta * x + xpl;
        return delta * inp_model.Density(idxlayer)(xscale) * mode1.w() *
               mode2.w() *
               (mode1.u(idxlayer)(xscale) * mode2.u(idxlayer)(xscale) +
                mode1.v(idxlayer)(xscale) * mode2.v(idxlayer)(xscale)) *
               xscale * xscale;
      };
      totint += q.Integrate(fun);
    };
  };
  totint *= frequencynorm * frequencynorm;
  totint *= 4.0 * pi_db * pi_db;
  totint *= inp_model.DensityNorm() / densitynorm;
  return totint;
};

int
main() {
  std::cout << "Hello\n";
  Timer timer1;
  timer1.start();

  // paths
  std::string pathtopremeig =
      "/space/adcm2/mineos-1.0.2/OUTPUT/eprem_noocean_S_D";
  std::string pathtoprempert =
      "/space/adcm2/mineos-1.0.2/OUTPUT/eprem_noocean_S_D_perturb";
  std::string pathtoprem =
      "/space/adcm2/mineos-1.0.2/DEMO/models/prem_noocean_noattenuation.txt";
  std::string pathtopremperturb = "/space/adcm2/mineos-1.0.2/DEMO/models/"
                                  "prem_noocean_noattenuation_perturb.txt";

  // models
  auto prem = EarthModels::ModelInput(pathtoprem);
  auto premperturb = EarthModels::ModelInput(pathtopremperturb);

  // mode catalogues
  ModeCoupling::modecataloguecontinuous modes_prem(pathtopremeig, prem);
  ModeCoupling::modecataloguecontinuous modes_perturbed(pathtoprempert,
                                                        premperturb);

  timer1.stop("Time to read in modes");
  for (int idx = 0; idx < 5; ++idx) {
    auto testmode = modes_prem.singlemode(idx);
    std::cout << testmode.n() << " " << testmode.l() << " " << testmode.q()
              << " " << testmode.gv() << " " << testmode.w() << "\n";
  }
  std::cout << "Number of modes: " << modes_prem.NumberOfModes() << "\n";

  // layers
  int numlayers = prem.NumberOfLayers();

  std::cout << "Density at inner core: "
            << prem.Density(0)(0.0) * prem.DensityNorm() << "\n";
  std::cout << "Density norm: " << prem.DensityNorm() << "\n";
  std::cout << "Velocity prem at surface: "
            << prem.VPV(numlayers - 1)(1.0) * prem.VelocityNorm() << "\n";
  std::cout << "Velocity perturbed prem at surface: "
            << premperturb.VPV(numlayers - 1)(1.0) * premperturb.VelocityNorm()
            << "\n";

  // checking continuous
  auto mode1 = modes_prem.singlemode(0);
  auto mode2 = modes_prem.singlemode(2);
  auto mode3 = modes_perturbed.singlemode(0);

  auto myint1 = modescalarproduct(mode1, mode1, prem);
  auto myint2 = modescalarproduct(mode1, mode2, prem);
  auto myint21 = modescalarproduct(mode2, mode1, prem);
  auto myint3 = modescalarproduct(mode1, mode3, prem);
  auto myint4 = modescalarproduct(mode2, mode3, prem);
  auto myint5 = modescalarproduct(mode3, mode3, premperturb);
  std::cout << "Integral: " << myint1 << " " << myint2 << " " << myint21 << " "
            << myint3 << " " << myint4 << " " << myint5 << "\n";

  timer1.start();

  int maxmodes = 50;
  Eigen::MatrixXd mat_norms(maxmodes, maxmodes);
  for (int idx1 = 0; idx1 < maxmodes; ++idx1) {
    for (int idx2 = 0; idx2 < maxmodes; ++idx2) {
      auto mode1 = modes_prem.singlemode(idx1);
      auto mode2 = modes_prem.singlemode(idx2);
      mat_norms(idx1, idx2) = modescalarproduct(mode1, mode2, prem);
    }
  }
  Eigen::VectorXd vec_force(maxmodes);
  Eigen::VectorXd vec_force2(maxmodes);
  auto l21_pert = modes_perturbed.singlemode(0);
  auto l21_prem = modes_prem.singlemode(0);
  for (int idx1 = 0; idx1 < maxmodes; ++idx1) {
    //   for (int idx2 = 0; idx2 < modes_prem.NumberOfModes(); ++idx2) {
    auto mode1 = modes_prem.singlemode(idx1);

    vec_force(idx1) = modescalarproduct(mode1, l21_pert, prem);
    vec_force2(idx1) = modescalarproduct(mode1, l21_prem, prem);
    //   }
  }
  timer1.stop("Time for matrix assembly");

  timer1.start();
  Eigen::FullPivLU<Eigen::MatrixXd> solver(mat_norms);
  Eigen::VectorXd vec_sol = solver.solve(vec_force);
  Eigen::VectorXd vec_sol2 = solver.solve(vec_force2);
  timer1.stop("Time for LU");
  //    timer1.start();
  //    Eigen::LDLT<Eigen::MatrixXd> solver2(mat_norms);
  //    Eigen::VectorXd vec_sol2 = solver2.solve(vec_force);
  //    timer1.stop("Time for Cholesky");
  // for (int idx1 = 0; idx1 < modes_prem.NumberOfModes(); ++idx1) {
  //    std::cout << vec_sol(idx1) << " " << vec_force(idx1) << "\n";
  //    //   }
  // }
  // std::cout << "\n\n";
  // for (int idx1 = 0; idx1 < modes_prem.NumberOfModes(); ++idx1) {
  //    std::cout << vec_sol2(idx1) << " " << vec_force2(idx1) << "\n";
  //    //   }
  // }
  timer1.start();
  std::vector<double> xval, yval1, yval2, yval3;
  std::vector<int> xidx;
  double maxstep = 0.01;
  // xval.push_back(0.0);
  // xidx.push_back(0);
  double xtmp = 0.0;
  {
    for (int idx = 0; idx < prem.NumberOfLayers(); ++idx) {
      int numsteps =
          std::ceil((prem.UpperRadius(idx) - prem.LowerRadius(idx)) / maxstep);
      numsteps = std::max(numsteps, 3);
      double stepsize =
          (prem.UpperRadius(idx) - prem.LowerRadius(idx)) / numsteps;
      xval.push_back(xtmp);
      xidx.push_back(idx);
      for (int idx1 = 1; idx1 < numsteps + 1; ++idx1) {
        xval.push_back(xtmp + stepsize);
        xidx.push_back(idx);
        xtmp = xval.back();
      }
      // std::cout << "idx: " << idx << " " << stepsize << "\n";
    }
  }
  yval1.resize(xval.size());
  yval2.resize(xval.size());
  yval3.resize(xval.size());

  for (int idx = 0; idx < xval.size(); ++idx) {
    yval1[idx] = l21_pert.u(xidx[idx])(xval[idx]);
    yval3[idx] = l21_prem.u(xidx[idx])(xval[idx]);

    yval2[idx] = 0.0;
    for (int idx1 = 0; idx1 < maxmodes; ++idx1) {
      yval2[idx] += modes_prem.singlemode(idx1).u(xidx[idx])(xval[idx]) *
                    vec_sol(idx1) * modes_prem.singlemode(idx1).w() /
                    l21_pert.w();
      // std::cout << idx << " " << idx1 << "\n";
    }
    // std::cout << idx << " " << yval1[idx] << " " << yval2[idx] << "\n";
  }

  // output

  timer1.stop("Final part");

  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  // outputting result
  std::string pathtofile = "./work/prempert.out";
  auto file2 = std::ofstream(pathtofile);
  for (int i = 0; i < xval.size(); ++i) {
    file2 << std::setprecision(16) << xval[i] << ";" << yval1[i] << ";"
          << yval2[i] << ";" << yval3[i] << std::endl;
  };

  // for (int idx = 0; idx < xval.size(); ++idx) {
  //    std::cout << idx << " " << xval[idx] << " " << xidx[idx] << "\n";
  // }

  //    auto q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(4);
  //    auto fun = [](double x) { return 1.0; };
  //    std::cout << q.Integrate(fun) << "\n";
  //    std::cout << std::setprecision(10)
  //              << prem.Density(0)(190859 / prem.LengthNorm()) *
  //              prem.DensityNorm()
  //              << "\n";
  //    std::cout << q.Integrate(fun) << "\n";
  //    std::cout << "\n u at zero: " << mode1.u(0)(0.0) << " "
  //              << mode1.u(0)(prem.UpperRadius(0)) << "\n";
  //    auto uvec = mode1.u();
  //    for (int idx = 0; idx < prem.NumberOfLayers(); ++idx) {
  //       std::cout << "Difference: "
  //                 << uvec[prem.LayerLowerIndex(idx)] -
  //                        mode1.u(idx)(prem.LowerRadius(idx))
  //                 << " "
  //                 << uvec[prem.LayerUpperIndex(idx)] -
  //                        mode1.u(idx)(prem.UpperRadius(idx))
  //                 << "\n";
  //    }
  //    std::cout << "\nIndices: \n";
  //    for (int idx = 0; idx < 13; ++idx) {
  //       std::cout << prem.LayerLowerIndex(idx) << " " <<
  //       prem.LayerUpperIndex(idx)
  //                 << "\n";
  //    }
  return 0;
}