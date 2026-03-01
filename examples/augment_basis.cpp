#include <iostream>
#include <cmath>
#include <functional>
#include <GaussQuad/All>
#include <Interpolation/All>
#include <filesystem>
#include <fstream>
#include <Eigen/Core>
#include <Spectra/SymEigsSolver.h>
#include <Eigen/SparseCore>
#include <Eigen/Eigenvalues>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/SymGEigsSolver.h>
#include <Spectra/MatOp/DenseSymMatProd.h>
#include <Spectra/MatOp/SparseCholesky.h>
#include <gplspec/Timer>
#include "string_classdefinitions.h"
using namespace Spectra;

// main program
int
main() {

  // testing model string
  std::vector<double> x_radii{0.0, 0.2, 0.3, 1.0};
  std::vector<double> mu_vec{1.0, 1.0, 1.0};
  std::vector<double> rho_vec{1.0, 1.0, 1.0};

  // string and spectral element
  model_string test_string(x_radii, mu_vec, rho_vec);
  spectral_element test_spec(test_string, 0.05, 6);
  test_spec.CalculateEigenfrequencies(5);
  auto x = test_spec.xvalues();
  auto myres = test_spec.augment_basis(2);
  test_spec.augment_basis_calculate();
  auto myres2 = test_spec.augment_basis();
  gsab gsab_test(test_spec, test_string, test_string, true);
  gsab gsab_test2(test_spec, test_string, test_string, 1);

  std::cout << std::setprecision(16) << "Frequencies 1:\n"
            << gsab_test.efrequencies() << "\n\n";
  std::cout << std::setprecision(16) << "Frequencies 2:\n"
            << gsab_test2.efrequencies() << "\n\n";
  //   std::cout << myres << "\n\n";

  // ///////////////////////////////////////////////////////////////
  // ///////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  // outputting result
  std::string pathtofile = "./work/string/augment_f.out";
  auto file = std::ofstream(pathtofile);
  for (int i = 0; i < x.size(); ++i) {
    file << std::setprecision(16) << x[i] << ";" << myres(i) << std::endl;
  };
  file.close();

  // outputting result
  std::string pathtofile2 = "./work/string/augment_f2.out";
  auto file2 = std::ofstream(pathtofile2);
  for (int i = 0; i < x.size(); ++i) {
    file2 << std::setprecision(16) << x[i] << ";" << myres2(i, 0) << ";"
          << myres2(i, 1) << ";" << myres2(i, 2) << std::endl;
  };
  file2.close();
  return 0;
}