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
  using VECTOR = std::vector<double>;
  // declarations
  double mu0 = 1.0, mu1 = 1.3;
  double pi_db = 3.1415926535897932;
  int maxn = 100;

  // eigenfunction catalogues (analytical result)
  eigenfunction_catalogue cateig(maxn, mu0), cateig_p(maxn, mu1),
      cateig3(maxn, 2.0, 2.0);
  for (int i = 0; i < 5; ++i) {
    std::cout << std::setprecision(15) << i << " " << cateig.w(i) << "\n";
  }

  //////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////
  // augment function
  double aval =
      -(mu0 * cateig_p.fd(0)(1.0) + cateig_p.f(0)(1.0)) / (mu0 * pi_db);
  auto augmentfunc = [pi_db, aval](double x) {
    return aval * std::sin(pi_db * x);
  };

  // expansion coefficients
  auto vec_coeff = force(cateig_p.f(0), cateig);
  auto vec_coeff2 = force_augment(cateig_p.f(0), augmentfunc, cateig);

  //////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////
  Timer timer1;
  timer1.start();
  int num_elements = 120;
  int num_q = 8;
  // spectral element solution for eigenfunction
  auto spectest = spectral_element(mu0, num_elements, num_q);
  spectest.CalculateEigenfrequencies(5);

  auto xvals = spectest.xvalues();
  Eigen::VectorXcd eig_vec = spectest.evector_value(4);
  timer1.stop("Time for spectral element solution");

  std::cout << "\n xval first: " << xvals[0] << "\n";
  std::cout << "xval size: " << xvals.size() << " " << eig_vec.rows() << "\n\n";

  timer1.start();
  // output vectors
  int nvals = 1001;
  double xsep = 1.0 / (nvals - 1);
  // VECTOR x(nvals);   //, y1(nvals), y2(nvals), yapprox(nvals),
  // yapprox2(nvals); std::generate(x.begin(), x.end(),
  //               [n = 0, &xsep]() mutable { return n++ * xsep; });
  auto x = spectest.xvalues();
  auto y1 = cateig_p.f(0, x);
  auto y2 = cateig.f(0, x);
  auto yapprox = cateig.basis_expansion(vec_coeff, x);
  auto yapprox2 = cateig.basis_expansion_augment(vec_coeff2, x, augmentfunc);
  timer1.stop("Time for eigenexpansion");

  //////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////
  // galerkin solver
  timer1.start();
  int neig = 5;
  galerkin_solver gal_test(maxn, mu0, mu0, num_elements, num_q);
  galerkin_solver gal_mu(maxn, mu0, mu1, num_elements, num_q);

  galerkin_solver gal_mu2(maxn, mu0, mu1, num_elements, num_q, true);
  // std::cout << gal_mu2.efrequencies() << "\n\n";
  auto gal_eig = gal_test.evector_physical(0, xvals);
  auto gal_eig2 = gal_mu.evector_physical(0, x);
  auto gal_eig3 = gal_mu2.evector_physical(0, x);
  timer1.stop("Time for Galerkin normal mode");

  // outputting frequencies
  std::cout << "\nUnperturbed: " << cateig.w(0)
            << ", perturbed: " << cateig_p.w(0) << "\n";
  auto spectest2 = spectral_element(mu0, num_elements, num_q);
  spectest2.CalculateEigenfrequencies(5);
  Eigen::VectorXcd spec_pert = spectest2.efrequencies();
  std::cout << "Spec elem: " << spec_pert << "\n\n";
  /*

  Eigen::VectorXcd galpertfreq = gal_mu.efrequencies();
  Eigen::VectorXcd galpertfreq2 = gal_mu2.efrequencies();
  std::cout << "Galerkin perturbed: " << galpertfreq.segment(0, 5) << "\n\n";
  std::cout << "Galerkin augment perturbed: " << galpertfreq2.segment(0, 5)
            << "\n\n";

  std::cout << "Error in GalAug (rel %): "
            << std::abs((galpertfreq2(0) - cateig_p.w(0)) / cateig_p.w(0)) *
                   100.0
            << "\n";
  std::cout << "Error in Gal (rel %): "
            << std::abs((galpertfreq(0) - cateig_p.w(0)) / cateig_p.w(0)) *
                   100.0
            << "\n";
  std::cout << "Error in SpectralElement (rel %): "
            << std::abs((spec_pert(0) - cateig_p.w(0)) / cateig_p.w(0)) * 100.0
            << "\n";
*/

  // ///////////////////////////////////////////////////////////////

  // ///////////////////////////////////////////////////////////////

  // testing model string
  std::vector<double> x_radii{0.0, 0.2, 0.205, 0.3, 1.0};
  std::vector<double> mu_vec{1.0, 1.1, 1.0, 1.0}, mu_vec2{1.0, 1.0, 1.0, 1.0};

  std::vector<double> rho_vec{1.0, 1.0, 1.0, 1.0};
  // for (auto &idx : mu_vec) {
  //   idx *= 2.0;
  // }
  // for (auto &idx : rho_vec) {
  //   idx *= 2.0;
  // }
  // std::vector<double> x_radii{0.0, 1.0};
  // std::vector<double> mu_vec{1.0};
  // std::vector<double> rho_vec{1.0};

  model_string test_string(x_radii, mu_vec, rho_vec);
  model_string test_string2(x_radii, mu_vec2, rho_vec);
  // spectral_element test_spectral_string(test_string, 0.005, 6);
  // spectral_element test_spectral_string2(test_string2, 0.005, 6);
  spectral_element test_spectral_string(test_string, 0.01, 6);
  spectral_element test_spectral_string2(test_string2, 0.01, 6);
  // std::cout << "xelem: \n";
  // for (auto idx : test_spectral_string.xelem()) {
  //   std::cout << idx << "\n";
  // }
  auto veclayer = test_spectral_string.layers();
  auto x2 = test_spectral_string.xvalues();
  // for (auto idx : veclayer) {
  //    std::cout << idx << "\n";
  // }
  // std::cout << "Test 1\n";
  test_spectral_string.CalculateEigenfrequencies(30);
  test_spectral_string2.CalculateEigenfrequencies(5);
  Eigen::VectorXcd eig_vec_test = test_spectral_string.evector_value(4);
  auto eig_vec_deriv = test_spectral_string.evector_deriv(4);
  auto eig_vec_deriv2 = test_spectral_string2.evector_deriv(4);
  // std::cout << "Test 2\n";
  Eigen::VectorXcd spec_pert_test = test_spectral_string.efrequencies();
  // std::cout << "Spec elem new: " << spec_pert_test << "\n\n";

  // std::cout << "Error in old (rel %): "
  //           << std::abs((spec_pert(0) - cateig.w(0)) / cateig.w(0)) * 100.0
  //           << "\n";
  // std::cout << "Error in new (rel %): "
  //           << std::abs((spec_pert_test(0) - cateig3.w(0)) / cateig.w(0)) *
  //                  100.0
  //           << "\n";

  ///////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////
  // testing the gsab class
  timer1.start();
  test_spectral_string.augment_basis_calculate();
  gsab gsab_test0(test_spectral_string, test_string, test_string2);
  gsab gsab_test(test_spectral_string, test_string, test_string2, true);
  gsab gsab_test2(test_spectral_string, test_string, test_string2, 1);
  auto eig_val_test_gsab = gsab_test.evalues();
  auto eig_vec_test_gsab = gsab_test.evectors();
  // std::cout << std::setprecision(3) << eig_vec_test_gsab << "\n\n";
  // std::cout << std::setprecision(6) << eig_val_test_gsab << "\n\n";

  std::cout << "first result: " << test_spectral_string2.efrequencies()
            << "\n\n";

  std::cout << "no aug result: " << gsab_test0.efrequencies().head(5) << "\n\n";
  std::cout << "augment result 1: " << gsab_test.efrequencies().head(5)
            << "\n\n";
  std::cout << "augment result 2: " << gsab_test2.efrequencies().head(5)
            << "\n\n";

  timer1.stop("Time for gsab x 3");
  std::cout << "Unaug diff: "
            << test_spectral_string2.efrequencies() -
                   gsab_test0.efrequencies().head(5)
            << "\n\n";
  std::cout << "Aug diff: "
            << test_spectral_string2.efrequencies() -
                   gsab_test.efrequencies().head(5)
            << "\n\n";

  std::cout << "Aug diff 2: "
            << test_spectral_string2.efrequencies() -
                   gsab_test2.efrequencies().head(5)
            << "\n\n";

  // ///////////////////////////////////////////////////////////////
  // ///////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  // outputting result
  std::string pathtofile = "./work/string/eigenfunction.out";
  auto file = std::ofstream(pathtofile);
  for (int i = 0; i < x.size(); ++i) {
    file << std::setprecision(16) << x[i] << ";" << y1[i] << ";" << y2[i] << ";"
         << yapprox[i] << ";" << yapprox2[i] << ";" << std::sin(x[i] * pi_db)
         << ";" << gal_eig2[i].real() << ";" << gal_eig3[i].real() << std::endl;
  };
  file.close();

  // outputting result
  std::string pathtofile1 = "./work/string/eigenfunction_spec.out";
  auto file1 = std::ofstream(pathtofile1);
  for (int i = 0; i < xvals.size(); ++i) {
    file1 << std::setprecision(16) << xvals[i] << ";" << eig_vec(i).real()
          << ";" << eig_vec(i).imag() << ";" << gal_eig[i].real() << std::endl;
  };
  file1.close();

  // outputting result
  std::string pathtofile2 = "./work/string/eigenfunction_spec_2.out";
  auto file2 = std::ofstream(pathtofile2);
  for (int i = 0; i < x2.size(); ++i) {
    file2 << std::setprecision(16) << x2[i] << ";" << eig_vec_test(i).real()
          << ";" << eig_vec_test(i).imag() << std::endl;
  };
  file2.close();

  // outputting result
  auto xnodesout = test_spectral_string.xnodes();
  std::string pathtofile3 = "./work/string/eigenfunction_spec_3.out";
  auto file3 = std::ofstream(pathtofile3);
  for (int idx_i = 0; idx_i < xnodesout.size(); ++idx_i) {
    for (int idx_j = 0; idx_j < xnodesout[0].size(); ++idx_j) {
      file3 << std::setprecision(16) << xnodesout[idx_i][idx_j] << ";"
            << eig_vec_deriv[idx_i][idx_j].real() << ";"
            << eig_vec_deriv[idx_i][idx_j].imag() << ";"
            << eig_vec_deriv2[idx_i][idx_j].real() << ";"
            << eig_vec_deriv2[idx_i][idx_j].imag() << std::endl;
    }
  }
  file3.close();

  return 0;
}