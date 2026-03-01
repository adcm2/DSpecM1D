
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <FFTWpp/Ranges>
#include <GSHTrans/All>
#include <GaussQuad/All>
#include <new_coupling/All>
#include <Interpolation/All>
// #include <PlanetaryModel/All>
// #include <TomographyModels/All>
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

int
main() {
  using myvector = std::vector<double>;
  using myiter = myvector::iterator;
  using namespace Interpolation;
  using InterpA = Interpolation::CubicSpline<myiter, myiter>;
  auto newf = [](double x) { return std::sin(x); };

  int nsize = 21, nsize2 = 100, nsize3 = 3;
  std::vector<double> xval(nsize, 0.0), yval(nsize, 0.0);
  std::vector<double> xval2(nsize2, 0.0), yval2(nsize2, 0.0),
      yval3(nsize2, 0.0);
  std::vector<double> xval_test(nsize3, 0.0), yval_test(nsize3, 0.0);

  double pi_db = 3.1415926535897932;
  for (int idx = 0; idx < xval.size(); ++idx) {
    xval[idx] = pi_db * static_cast<double>(idx) /
                static_cast<double>((xval.size() - 1));
    yval[idx] = newf(xval[idx]);
  }

  for (int idx = 0; idx < xval_test.size(); ++idx) {
    xval_test[idx] = pi_db * static_cast<double>(idx) /
                     static_cast<double>((xval_test.size() - 1));
    // yval[idx] = newf(xval[idx]);
  }
  for (auto idx : xval_test) {
    std::cout << idx << "\n";
  }

  auto interpval =
      InterpA(xval.begin(), xval.end(), yval.begin(), CubicSplineBC::Clamped, 1,
              CubicSplineBC::Clamped, -1);
  auto q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(5);
  auto qx = q.Points();
  auto qw = q.Weights();
  auto q2 = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(8);
  auto qx2 = q2.Points();
  auto qw2 = q2.Weights();
  std::vector<std::vector<double>> vec_xq;

  for (int idx = 0; idx < xval2.size(); ++idx) {
    xval2[idx] = pi_db * static_cast<double>(idx) /
                 static_cast<double>((xval2.size() - 1));
    yval2[idx] = newf(xval2[idx]);
    yval3[idx] = interpval(xval2[idx]);
  }

  std::cout << "Hello\n";
  double totint = 0.0;
  for (int idxlayer = 0; idxlayer < xval.size() - 1; ++idxlayer) {
    double tmpint = 0.0;
    double x1 = xval[idxlayer];
    double x2 = xval[idxlayer + 1];
    double delta = (x2 - x1) / 2.0;
    double xpl = (x1 + x2) / 2.0;
    for (int idxq = 0; idxq < qx.size(); ++idxq) {
      double xscale = delta * qx[idxq] + xpl;
      tmpint += qw[idxq] * interpval(xscale);
    }
    tmpint *= delta;
    totint += tmpint;
  }

  double totint2 = 0.0, totint3 = 0.0;
  for (int idxlayer = 0; idxlayer < xval_test.size() - 1; ++idxlayer) {
    double tmpint = 0.0, tmpint1 = 0.0;
    double x1 = xval_test[idxlayer];
    double x2 = xval_test[idxlayer + 1];
    double delta = (x2 - x1) / 2.0;
    double xpl = (x1 + x2) / 2.0;
    std::vector<double> xtmp;
    for (int idxq = 0; idxq < qx2.size(); ++idxq) {
      double xscale = delta * qx2[idxq] + xpl;
      xtmp.push_back(xscale);
      tmpint1 += qw2[idxq] * newf(xscale);
      tmpint += qw2[idxq] * interpval(xscale);
    }
    tmpint *= delta;
    totint2 += tmpint;

    tmpint1 *= delta;
    totint3 += tmpint1;
    vec_xq.push_back(xtmp);
  }

  std::cout << std::setprecision(16) << "Exact: " << 2.0
            << ", interp: " << totint << ", 2 interp: " << totint2
            << ", 2 exact: " << totint3 << "\n\n";
  std::cout << std::setprecision(16)
            << "Exact 1: " << std::abs((2.0 - totint) / 2.0)
            << ", error 2: " << std::abs((2.0 - totint2) / 2.0)
            << ", error 3: " << std::abs((2.0 - totint3) / 2.0) << "\n\n";
  std::cout << std::setprecision(16) << "Error against interpolant: "
            << std::abs((totint2 - totint) / 2.0) << "\n\n";

  // outputting result
  std::string pathtofile1 = "./work/splinecheck.out";
  auto file1 = std::ofstream(pathtofile1);
  {
    for (int idxe = 0; idxe < xval2.size(); ++idxe) {
      file1 << std::setprecision(16) << xval2[idxe] << ";" << yval2[idxe] << ";"
            << yval3[idxe] << std::endl;
    };
  };

  file1.close();

  // outputting result
  std::string pathtofile2 = "./work/splinecheck2.out";
  auto file2 = std::ofstream(pathtofile2);
  {
    for (int idxe = 0; idxe < xval.size(); ++idxe) {
      file2 << std::setprecision(16) << xval[idxe] << ";" << yval[idxe]
            << std::endl;
    };
  };

  file2.close();

  // outputting result
  std::string pathtofile3 = "./work/splinecheck3.out";
  auto file3 = std::ofstream(pathtofile3);
  {
    for (int idx1 = 0; idx1 < vec_xq.size(); ++idx1) {
      for (int idxq = 0; idxq < qx2.size(); ++idxq) {
        file3 << std::setprecision(16) << vec_xq[idx1][idxq] << ";"
              << newf(vec_xq[idx1][idxq]) << std::endl;
      }
    };
  };

  file3.close();
  return 0;
}