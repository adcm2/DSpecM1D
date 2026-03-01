#ifndef TOROIDAL_BENCH_CLASS_GUARD_H
#define TOROIDAL_BENCH_CLASS_GUARD_H
#include <iostream>
#include <cmath>
#include <functional>
#include <filesystem>
#include <fstream>
#include <vector>
// #include "toroidal.h"

// auto
// sph_bessel_deriv(int n, double x, double h = std::pow(10.0, -6.0)) {
//   double xh = x + h;
//   double xhm = x - h;
//   double tmp1 = std::sph_bessel(n, xhm);
//   double tmp2 = std::sph_bessel(n, xh);
//   return (tmp2 - tmp1) / (2.0 * h);
// }

namespace Toroidal_Bench {

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
       int maxit = 100) {
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
      if (std::abs(f(xh)) < std::pow(10.0, -15.0)) {
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

class sphere_bench {
private:
  int _l, _nmax;
  std::vector<double> _xsol;

public:
  sphere_bench(int l, int nmax) : _l{l}, _nmax{nmax} {

    // function to find result
    auto fv = [lv = _l](double x) { return tracfunc(lv, x); };

    // variables to do bisection
    double xl, xu, xval = l / 2.0;
    // int numx = 10;
    int maxnum = _nmax;
    int cidx = 0;
    _xsol = std::vector<double>(_nmax, 0.0);
    if (_l == 1) {
      maxnum = _nmax - 1;
      cidx = 1;
    }

    for (int idx = 0; idx < maxnum; ++idx) {
      // get values of lower and upper
      xl = xval + 0.05;
      xu = xl;
      int fpl = (fv(xl) > 0);
      while (((fv(xu) > 0) == fpl)) {
        xu += 0.3;
      }
      xval = bisect(fv, xl, xu);
      _xsol[cidx++] = xval;
    }
  };

  auto w() const { return _xsol; };
};
/*
int
main() {

  // test newton
  int l = 1;
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
}*/
}   // namespace Toroidal_Bench

#endif