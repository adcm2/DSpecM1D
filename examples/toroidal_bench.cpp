#include <iostream>
#include <cmath>
#include <functional>
#include <filesystem>
#include <fstream>
// #include "toroidal.h"

// auto
// sph_bessel_deriv(int n, double x, double h = std::pow(10.0, -6.0)) {
//   double xh = x + h;
//   double xhm = x - h;
//   double tmp1 = std::sph_bessel(n, xhm);
//   double tmp2 = std::sph_bessel(n, xh);
//   return (tmp2 - tmp1) / (2.0 * h);
// }

auto
sph_bessel_deriv(int n, double x) {
  return -std::sph_bessel(n + 1, x) + (n / x) * std::sph_bessel(n, x);
}

auto
tracfunc(int n, double x) {
  if (x <= 0) {
    std::cout << "Less than 0!\n";
  }
  return (n - 1) * std::sph_bessel(n, x) - x * std::sph_bessel(n + 1, x);
}
auto
funcderiv(int n, double x) {
  if (x <= 0) {
    std::cout << "Less than 0!\n";
  }
  double tmp = n * (n - 1) / x * std::sph_bessel(n, x);
  tmp -= (2.0 * n + 1) * std::sph_bessel(n + 1, x);
  tmp += x * std::sph_bessel(n + 2, x);
  return tmp;
}

auto
bisect(const std::function<double(double)> &f, double x0, double x1,
       int maxit = 50) {
  double f0 = f(x0);
  double f1 = f(x1);
  int fp0 = (f0 > 0), fp1 = (f1 > 0);
  //   std::cout << f0 << " " << f1 << " " << fp0 << " " << fp1 << "\n\n";
  // ((f1 > 0) && (f0 > 0)) || ((f1 < 0) && (f0 < 0))
  if (fp0 == fp1) {
    std::cout << "Both same sign\n";
    return x0;
  } else {
    double xh = x0;
    double xl = x0;
    double xu = x1;
    for (int idx = 0; idx < maxit; ++idx) {
      if (std::abs(f(xh)) < std::pow(10.0, -14.0)) {
        break;
      }
      xh = 0.5 * (xl + xu);
      //   std::cout << std::setprecision(3) << xl << " " << xh << " " << xu <<
      //   "\n"; std::cout << f(xl) << " " << f(xh) << " " << f(xu) << "\n";
      int fph = (f(xh) > 0);

      if (fph == fp0) {
        xl = xh;

      } else {
        xu = xh;
      }

      //   fp0 = (f(xl))
    }
    return xh;
  }
};

int
main() {
  /*
// spot check for n == 1
double x = 1.2345;
double func_sphb = std::sph_bessel(1, x);
double func_exact = std::sin(x) / (x * x) - std::cos(x) / x;
double deriv_sphb = sph_bessel_deriv(1, x);
double deriv_exact =
    std::sin(x) * (1 - 2.0 / (x * x)) / x + 2.0 * std::cos(x) / (x * x);

std::cout << "j_1(" << x << ") = " << std::setprecision(16) << func_sphb
          << '\n';

// exact solution for j_1
std::cout << "sin(x)/x² - cos(x)/x = " << func_exact << '\n';

// derivatives
std::cout << "\nj_1'(" << x << ") = " << deriv_sphb << '\n';

// exact solution for j_1

std::cout << "sin x(1 - 2/x^2)/x + 2 cos x/x^2= " << deriv_exact << '\n';

std::cout << "\nRel error function: "
          << std::abs(func_exact - func_sphb) / func_sphb * 100
          << ", derivative: "
          << std::abs(deriv_exact - deriv_sphb) / deriv_sphb * 100 << "\n\n";
*/
  //

  // test newton
  int l = 0;
  std::cout << "Enter l:\n";
  std::cin >> l;
  auto fv = [l](double x) { return tracfunc(l, x); };
  auto fd = [l](double x) { return funcderiv(l, x); };
  //   double xval = newton(fv, fd, 0.1, 20);
  //   double xval2 = newton(fv, fd, 5.0);
  //   double xval3 = newton(fv, fd, 5.0);
  double xl = 0.5, xu = 6.0, xval = l / 2.0;
  int numx = 10;
  std::vector<double> x_sol(numx);
  for (int idx = 0; idx < numx; ++idx) {
    // get values of lower and upper
    xl = xval + 0.05;
    xu = xl;
    int fpl = (fv(xl) > 0);
    while (((fv(xu) > 0) == fpl)) {
      xu += 0.3;
    }
    xval = bisect(fv, xl, xu);
    x_sol[idx] = xval;
  }
  std::vector<double> vec_x(5001, 0.0), vec_y(5001, 0.0);
  for (int idx = 1; idx < vec_x.size(); ++idx) {
    vec_x[idx] = idx * 0.01;
    vec_y[idx] = fv(vec_x[idx]);
  }

  for (auto &idx : x_sol) {
    std::cout << std::setprecision(15) << idx << "\n";
  }

  ///////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  // output
  // outputting result
  std::string pathtofile = "./work/toroidal/functiontozero.out";
  auto file = std::ofstream(pathtofile);

  for (int idx = 0; idx < vec_x.size(); ++idx) {

    file << std::setprecision(16) << vec_x[idx] << ";" << vec_y[idx]
         << std::endl;
  }
  file.close();
  //   std::cout << "\n\n" << xval2 << "\n";
  return 0;
}