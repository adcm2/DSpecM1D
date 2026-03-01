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
#include "toroidal_gal.h"

int
main() {
  // parameters
  int nfreq = 5;
  int nsolid = 1;
  int pn = 8;
  double twopi = 2.0 * 3.1415926535897932;
  double maxstep;
  int lval;
  double sigshift;
  bool toaug;
  std::string pathpert;
  std::cout << "Enter maxstep: \n";
  std::cin >> maxstep;
  std::cout << "Enter shift:\n";
  std::cin >> sigshift;
  std::cout << "Enter l:\n";
  std::cin >> lval;
  std::cout << "Enter # of frequencies:\n";
  std::cin >> nfreq;
  std::cout << "Enter solid layer:\n";
  std::cin >> nsolid;
  std::cout << "To aug or not?\n";
  std::cin >> toaug;
  std::cout << "Path to perturbed model:\n";
  std::cin >> pathpert;

  ///////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////
  // paths
  std::string pathtopremeig =
      "/space/adcm2/mineos-1.0.2/OUTPUT/eprem_noocean_S_IC";
  std::string pathtoprem =
      "/space/adcm2/mineos-1.0.2/DEMO/models/prem.200.no.noatten.txt";
  std::string pathtoprempert =
      "/space/adcm2/mineos-1.0.2/DEMO/models/prem.200.no.noatten.perturb.txt";
  std::string pathtopremnew = "work/toroidal/prem_pert.out";

  std::cout << pathpert << std::endl;

  // models
  auto prem = EarthModels::ModelInput(pathtoprem);
  auto prempert = EarthModels::ModelInput(pathpert);
  // auto prempert = TestTools::EarthModel();

  // vs

  // our spectral element method
  Toroidal::spectral_element_planet test_planet(prem, maxstep, pn, lval,
                                                nsolid);
  test_planet.CalculateEigenfrequencies(nfreq, sigshift);
  test_planet.augment_basis_calculate();

  // perturbed
  Toroidal::spectral_element_planet test_planet_pert(prempert, maxstep, pn,
                                                     lval, nsolid);
  test_planet_pert.CalculateEigenfrequencies(nfreq, sigshift);
  test_planet_pert.augment_basis_calculate();

  // test out gal
  Toroidal::gal_gen galtest(test_planet, prem, prempert, 0);
  Toroidal::gal_gen galtest2(test_planet, prem, prempert, 1);

  // check values
  Eigen::VectorXd w1 = test_planet.efrequencies_gen() / twopi;
  Eigen::VectorXd w2 = galtest.efrequencies();
  Eigen::VectorXd w3 = test_planet_pert.efrequencies_gen() / twopi;
  Eigen::VectorXd w4 = galtest2.efrequencies();

  double freqnorm = 1.0 / (prem.TimeNorm() * twopi);
  // std::cout << std::setprecision(16) << w1 << "\n\n";
  // std::cout << w3 << "\n\n";
  // std::cout << w2 * freqnorm << "\n";

  Eigen::VectorXd vec_dp = Eigen::VectorXd::Zero(w1.size());
  Eigen::VectorXd vec_dg = Eigen::VectorXd::Zero(w1.size());
  Eigen::VectorXd vec_dg2 = Eigen::VectorXd::Zero(w1.size());
  for (int idx = 0; idx < w1.size(); ++idx) {
    vec_dp(idx) = std::abs((w1(idx) - w3(idx)) / w3(idx));
    vec_dg(idx) = std::abs((w2(idx) * freqnorm - w3(idx)) / w3(idx));
    vec_dg2(idx) = std::abs((w4(idx) * freqnorm - w3(idx)) / w3(idx));
  }

  // difference output
  std::cout << "Relative difference in frequencies: \n";
  std::cout << std::setprecision(16) << vec_dp.head(5) << "\n\n";
  std::cout << "Relative error in Galerkin expansion: \n";
  std::cout << std::setprecision(16) << vec_dg.head(5) << "\n\n";
  std::cout << "Relative error in Galerkin expansion 2: \n";
  std::cout << std::setprecision(16) << vec_dg2.head(5) << "\n\n";

  // traction:
  auto vec_traction = test_planet.traction_std(prem);
  auto vec_traction_p = test_planet.traction_std(prempert);
  auto vec_traction_pp = test_planet_pert.traction_std(prempert);
  auto vec_traction_gal_unaug = galtest.Traction(test_planet, prempert);
  auto vec_traction_gal_aug = galtest2.Traction(test_planet, prempert);

  // eigenvectors:
  auto eig0 = test_planet.evectors_std_ref();
  auto eig_pert = test_planet_pert.evectors_std_ref();
  auto eig_gal_unaug = galtest.Modes_Ref();
  auto eig_gal_aug = galtest2.Modes_Ref();

  // eigenvectors
  Eigen::MatrixXcd mateig = galtest.evectors();

  // augment basis
  auto vec_augbasis = test_planet.augment_basis();

  auto aug_val = test_planet.aug_std_ref();
  auto aug_deriv = test_planet.aug_deriv_ref();
  // std::vector<std::vector<std::vector<double>>> vec_traction_augment =
  // aug_deriv; for(int idx = 0; ) std::cout << std::setprecision(2) << mateig
  // << "\n\n";

  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////

  std::string pathtofile = "./work/toroidal/eigenvectors_pert.out";
  auto file = std::ofstream(pathtofile);
  auto radmesh = test_planet.mesh();
  for (int idxe = 0; idxe < radmesh.NE(); ++idxe) {
    for (int idxq = 0; idxq < radmesh.NN(); ++idxq) {
      file << std::setprecision(16) << radmesh.NodeRadius(idxe, idxq);
      for (int idx = 0; idx < vec_traction.back().back().size(); ++idx) {
        file << ";" << eig0[idxe][idxq][idx] << ";" << eig_pert[idxe][idxq][idx]
             << ";" << eig_gal_unaug[idxe][idxq][idx].real() << ";"
             << eig_gal_aug[idxe][idxq][idx].real();
      }
      file << std::endl;
    }
  }
  file.close();

  std::string pathtofile5 = "./work/toroidal/traction_pert.out";
  auto file5 = std::ofstream(pathtofile5);
  // auto radmesh = test_planet.mesh();
  for (int idxe = 0; idxe < radmesh.NE(); ++idxe) {
    for (int idxq = 0; idxq < radmesh.NN(); ++idxq) {
      file5 << std::setprecision(16) << radmesh.NodeRadius(idxe, idxq);
      for (int idx = 0; idx < vec_traction.back().back().size(); ++idx) {
        file5 << ";" << vec_traction[idxe][idxq][idx] << ";"
              << vec_traction_p[idxe][idxq][idx] << ";"
              << vec_traction_pp[idxe][idxq][idx] << ";"
              << vec_traction_gal_unaug[idxe][idxq][idx].real() << ";"
              << vec_traction_gal_aug[idxe][idxq][idx].real();
      }
      file5 << std::endl;
    }
  }
  file5.close();

  std::string pathtofile6 = "./work/toroidal/augbasis.out";
  auto file6 = std::ofstream(pathtofile6);

  for (int idxe = 0; idxe < radmesh.NE(); ++idxe) {
    for (int idxq = 0; idxq < radmesh.NN(); ++idxq) {
      double cr = radmesh.NodeRadius(idxe, idxq);
      file6 << std::setprecision(16) << cr;
      for (int idx = 0; idx < aug_val.back().back().size(); ++idx) {
        file6 << ";" << aug_val[idxe][idxq][idx] << ";"
              << aug_deriv[idxe][idxq][idx] << ";"
              << prem.L(radmesh.LayerNumber(idxe))(cr) *
                     (aug_deriv[idxe][idxq][idx] -
                      aug_val[idxe][idxq][idx] / cr);
      }
      file6 << std::endl;
    }
  }
  file6.close();
  return 0;
}