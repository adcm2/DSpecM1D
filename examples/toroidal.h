#ifndef TOROIDAL_CLASS_GUARD_H
#define TOROIDAL_CLASS_GUARD_H

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
#include <Eigen/SparseLU>
#include <Eigen/Eigenvalues>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/SymGEigsSolver.h>
#include <Spectra/MatOp/DenseSymMatProd.h>
#include <Spectra/MatOp/SparseCholesky.h>
#include <gplspec/Timer>
#include <PlanetaryModel/All>
#include <gplspec/All>

using namespace Spectra;

namespace Toroidal {
double
newton(const std::function<double(double)> &f,
       const std::function<double(double)> &fd, double x0, int maxit = 10) {
  double xi, xf;
  xi = x0;
  xf = x0;
  for (int idx = 0; idx < maxit; ++idx) {
    if (std::abs(f(xf)) < std::pow(10.0, -14.0)) {
      break;
    };
    xf = xi - f(xi) / fd(xi);
    xi = xf;
  }
  return xf;
};

class eigenfunction_catalogue {
public:
  eigenfunction_catalogue() {};
  eigenfunction_catalogue(int num_eig, double a, double mu, double rho = 1.0)
      : _num_eig{num_eig}, _a{a}, _mu{mu}, _rho{rho},
        _beta{std::sqrt(mu / rho)} {
    auto fval = [](double x) {
      auto tmp = (3.0 - x * x) * std::sin(x) - 3.0 * x * std::cos(x);
      return tmp;
    };

    auto fdval = [](double x) {
      auto tmp = x * std::sin(x) - x * x * std::cos(x);
      return tmp;
    };

    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////
    // solving for wn
    _vec_wn.resize(_num_eig);
    _vec_wn[0] = _beta / a * newton(fval, fdval, 0.0, 30);
    _vec_wn[1] = _beta / a * newton(fval, fdval, 5.0, 30);
    for (int i = 2; i < _num_eig; ++i) {
      // eigenfrequencies
      double xtmp = (i + 1) * pi_db;
      _vec_wn[i] = _beta / a * newton(fval, fdval, xtmp, 30);
    }
  };

  //////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////
  int NumberOfEigenfunctions() const { return _num_eig; };
  double Mu() const { return _mu; };
  double w(int i) const { return _vec_wn[i]; };
  auto w() const { return _vec_wn; };

private:
  std::vector<double> _vec_wn;
  double _mu, _rho, _a, _beta;
  int _num_eig;
  double pi_db = 3.1415926535897932;
};

class model_string {
public:
  model_string() {};
  model_string(std::vector<double> &x_radii, std::vector<double> &vec_mu,
               std::vector<double> &vec_rho)
      : _x{x_radii}, _mu{vec_mu}, _rho{vec_rho}, _numlayers{vec_mu.size()} {
    assert((vec_mu.size() == vec_rho.size()) && "Different sizes");
    assert(((x_radii.size() - 1) == _numlayers) && "Different sizes");
  };

  auto NumberOfLayers() const { return _numlayers; };

  auto mu(int i) const {
    auto rf = [&mu = _mu, i](double x) { return mu[i]; };
    return rf;
  };
  auto rho(int i) const {
    auto rf = [&rho = _rho, i](double x) { return rho[i]; };
    return rf;
  };

  auto LowerRadius(int i) const { return _x[i]; };
  auto UpperRadius(int i) const { return _x[i + 1]; };
  auto LayerWidth(int i) const { return _x[i + 1] - _x[i]; };

private:
  std::size_t _numlayers;
  std::vector<double> _mu, _rho, _x;
};

// spectral element solver
class spectral_element {
public:
  spectral_element(const model_string &, double, int, int);
  template <class model1d> spectral_element(const model1d &, double, int, int);

  void CalculateEigenfrequencies(int);

  // value of eigenvector
  auto evector_value(int i) const {
    return evectors.block(0, i, evectors.rows(), 1);
  };

  // derivative of eigenvector
  //  outputs in form of vector of vectors, corresponding to element and nodes
  //  within each element in the inner vector
  auto evector_deriv(int) const;

  // matrix of all eigenvectors
  Eigen::MatrixXcd evectors_all() const { return evectors; };

  // eigenvalues and eigenfrequencies
  auto evalues() const { return evalues_seig; };
  auto efrequencies() const;

  auto efrequency(int i) const {
    return std::sqrt(std::abs(evalues_seig(evalues_seig.rows() - 1 - i)));
  };

  // augmentation basis function calculation
  void augment_basis_calculate();
  auto augment_deriv(int) const;

  auto augment_basis() const { return _augbasis; };

  auto is_augmented() const { return _augcalc; };
  Eigen::VectorXd augment_basis(int) const;

  auto xvalues() const { return x_valnodes; };
  auto xnodes() const { return x_nodes; };
  auto xelem() const { return x_elem; };

  auto &layers() const { return vec_layer; };

  Eigen::SparseMatrix<double> MatEig() const { return mat_seig; };

private:
  int _NE, _l, _k2;
  double _mu;
  GaussQuad::Quadrature1D<double> _q;
  Eigen::SparseMatrix<double> mat_seig, mat_ke;
  std::vector<double> x_elem, x_valnodes;
  std::vector<std::vector<double>> x_nodes;
  std::vector<int> vec_layer;
  Eigen::VectorXcd evalues_seig;
  Eigen::MatrixXcd evectors;
  Eigen::MatrixXd _augbasis;
  bool _calc_eig = false;
  bool _augcalc = false;
  Radial_Tools::RadialMesh _mesh;

  std::size_t overallidx(int idxe, int idxq) const {
    return (_q.N() - 1) * idxe + idxq - 1;
  };

  std::size_t overallidxfinal(int idxe, int idxq) const {
    return (_q.N() - 1) * idxe + idxq;
  };
};

// constructor
spectral_element::spectral_element(const model_string &inp_model,
                                   double maxwidth, int NQ, int l)
    : _l{l}, _k2{l * (l + 1)} {   // quadrature and Lagrange polynomials
  _q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(NQ);
  auto pleg =
      Interpolation::LagrangePolynomial(_q.Points().begin(), _q.Points().end());

  int numlayers = inp_model.NumberOfLayers();
  // x_elem.push_back(0.0);
  {
    int totelem = 0;
    for (int idx = 0; idx < numlayers; ++idx) {
      totelem += std::ceil(inp_model.LayerWidth(idx) / maxwidth);
    }
    _NE = totelem;
    x_elem.resize(totelem + 1);
    vec_layer.resize(totelem);
    x_nodes.resize(totelem);

    x_elem[0] = 0.0;
    int totidx = 0;
    for (int idx = 0; idx < numlayers; ++idx) {
      int numelements = std::ceil(inp_model.LayerWidth(idx) / maxwidth);
      double elemwidth = inp_model.LayerWidth(idx) / numelements;
      for (int idx1 = 0; idx1 < numelements - 1; ++idx1) {
        vec_layer[totidx] = idx;
        x_elem[totidx + 1] = x_elem[totidx] + elemwidth;
        ++totidx;
      }
      x_elem[totidx + 1] = inp_model.UpperRadius(idx);
      vec_layer[totidx] = idx;
      ++totidx;
    }
    for (int idx = 0; idx < totelem; ++idx) {
      std::vector<double> xtmp(_q.N());
      double xl = x_elem[idx];
      double xu = x_elem[idx + 1];
      double xd = (xu - xl) / 2.0;
      double xp = (xl + xu) / 2.0;
      xtmp[0] = xl;
      for (int k = 1; k < _q.N() - 1; ++k) {
        xtmp[k] = xp + xd * _q.X(k);
      }
      xtmp[_q.N() - 1] = xu;
      x_nodes[idx] = xtmp;
    }
  }

  // vectors for filling out matrix
  mat_seig.resize((_q.N() - 1) * _NE, (_q.N() - 1) * _NE);
  mat_ke.resize((_q.N() - 1) * _NE, (_q.N() - 1) * _NE);
  std::vector<double> vec_nn((_q.N() - 1) * _NE, 0.0),
      vec_lm1((_q.N() - 1) * _NE, 0.0);
  {
    {
      x_valnodes.resize((_q.N() - 1) * _NE + 1);
      x_valnodes[0] = 0.0;
      int idxov = 1;
      for (int i = 0; i < _NE; ++i) {
        for (int k = 1; k < _q.N(); ++k) {
          x_valnodes[idxov++] = x_nodes[i][k];
        }
      }
    }

    // generating matrix of derivative values for Lagrange polynomials
    std::vector<std::vector<std::vector<double>>> mat_gauss_deriv;
    std::vector<std::vector<double>> vec_lag_deriv, vec_delta;
    mat_gauss_deriv.reserve(_q.N());
    for (int idxk = 0; idxk < _q.N(); ++idxk) {
      std::vector<std::vector<double>> mat_tmp(_q.N(),
                                               std::vector<double>(_q.N()));
      for (int idxi = 0; idxi < _q.N(); ++idxi) {
        std::vector<double> vec_tmp(_q.N());
        for (int idxj = 0; idxj < _q.N(); ++idxj) {
          vec_tmp[idxj] = pleg.Derivative(idxi, _q.X(idxk)) *
                          pleg.Derivative(idxj, _q.X(idxk));
        }
        mat_tmp[idxi] = vec_tmp;
      }
      mat_gauss_deriv.push_back(mat_tmp);
    };
    for (int idxk = 0; idxk < _q.N(); ++idxk) {
      std::vector<double> vec_tmp(_q.N(), 0.0), vec_tmp1(_q.N(), 0.0);
      for (int idxi = 0; idxi < _q.N(); ++idxi) {
        vec_tmp[idxi] = pleg.Derivative(idxi, _q.X(idxk));
        vec_tmp1[idxi] = pleg(idxi, _q.X(idxk));
      }
      vec_lag_deriv.push_back(vec_tmp);
      vec_delta.push_back(vec_tmp1);
    };

    //////////////////////////////////////////////////////////////
    // rhs matrix
    // we don't evaluate the matrix itself.
    //  THIS ASSUMES THAT THE MATRIX IS DIAGONAL. ALTHOUGH THIS SHOULD HOLD IN
    //  BOTH THE SPECTRAL ELEMENT FORMULATION AS WELL AS USING THE NORMAL
    //  MODES FROM THE UNPERTURBED PROBLEM form matrix, using tripletlist
    //  vector
    //   tripletlist vector
    using T = Eigen::Triplet<double>;
    std::vector<T> tripletlist3, tripletlist_ke;
    tripletlist3.reserve(_q.N() * _q.N() * _NE);
    tripletlist_ke.reserve(_q.N() * _q.N() * _NE);

    {
      int idxtl = 0;
      double elem_width = x_elem[1] - x_elem[0];

      // first element
      for (int i = 1; i < _q.N(); ++i) {

        vec_nn[idxtl + i - 1] += elem_width / 2.0 * _q.W(i) *
                                 inp_model.rho(vec_layer[0])(x_nodes[0][i]) *
                                 x_nodes[0][i] * x_nodes[0][i];
      }
      idxtl += _q.N() - 2;

      // middle and last elements
      for (int idxe = 1; idxe < _NE; ++idxe) {
        elem_width = x_elem[idxe + 1] - x_elem[idxe];
        for (int i = 0; i < _q.N(); ++i) {
          vec_nn[idxtl + i] +=
              elem_width / 2.0 * _q.W(i) *
              inp_model.rho(vec_layer[idxe])(x_nodes[idxe][i]) *
              x_nodes[idxe][i] * x_nodes[idxe][i];
        }
        idxtl += _q.N() - 1;
      }
    }

    // find lambda^{-1} (ie diagonal)
    for (int i = 0; i < vec_nn.size(); ++i) {
      vec_lm1[i] = 1.0 / vec_nn[i];
    }

    //////////////////////////////////////////////////////////////
    // get tripletlist for matrix
    {
      // std::cout << "k2-2: " << (_k2 - 2.0) << "\n";
      int idxtl = 0;
      // first element
      // loop over basis
      {
        double elem_width = x_elem[1] - x_elem[0];
        double d_val = 2.0 / elem_width;
        for (int i = 1; i < _q.N(); ++i) {
          for (int j = 1; j < _q.N(); ++j) {
            double tmp = 0.0;

            // loop over points
            for (int k = 0; k < _q.N(); ++k) {
              tmp += _q.W(k) *
                     ((x_nodes[0][k] * vec_lag_deriv[k][i] * d_val -
                       vec_delta[k][i]) *
                          (x_nodes[0][k] * vec_lag_deriv[k][j] * d_val -
                           vec_delta[k][j]) +
                      (_k2 - 2.0) * (vec_delta[k][i] * vec_delta[k][j])) *
                     inp_model.mu(vec_layer[0])(x_nodes[0][k]);
            }
            tmp *= elem_width / 2.0;

            long int idx_x = idxtl + i - 1;
            long int idx_y = idxtl + j - 1;
            tripletlist3.push_back(T(idx_x, idx_y, tmp * vec_lm1[idx_x]));
            tripletlist_ke.push_back(T(idx_x, idx_y, tmp));
          }
        }
      }
      idxtl += _q.N() - 2;

      // middle and last elements
      // over elements
      for (int idxe = 1; idxe < _NE; ++idxe) {
        double elem_width = x_elem[idxe + 1] - x_elem[idxe];
        double d_val = 2.0 / elem_width;

        // over basis
        for (int i = 0; i < _q.N(); ++i) {
          for (int j = 0; j < _q.N(); ++j) {
            double tmp = 0.0;

            // over nodes
            for (int k = 0; k < _q.N(); ++k) {
              tmp += _q.W(k) *
                     ((x_nodes[idxe][k] * vec_lag_deriv[k][i] * d_val -
                       vec_delta[k][i]) *
                          (x_nodes[idxe][k] * vec_lag_deriv[k][j] * d_val -
                           vec_delta[k][j]) +
                      (_k2 - 2.0) * (vec_delta[k][i] * vec_delta[k][j])) *
                     inp_model.mu(vec_layer[idxe])(x_nodes[idxe][k]);
            }
            tmp *= elem_width / 2.0;

            tripletlist3.push_back(
                T(idxtl + i, idxtl + j, tmp * vec_lm1[idxtl + i]));
            tripletlist_ke.push_back(T(idxtl + i, idxtl + j, tmp));
          }
        }
        idxtl += _q.N() - 1;
      }

      //////////////////////////////////////////////////////////////
      // boundary and continuity conditions
      // add to last
      // tripletlist3.push_back(T(idxtl, idxtl, vec_lm1[idxtl]));
      // tripletlist_ke.push_back(T(idxtl, idxtl, 1.0));

      int idxoverall = 0;

      int qmax = _q.N() - 1;
    }
    mat_seig.setFromTriplets(tripletlist3.begin(), tripletlist3.end());
    mat_ke.setFromTriplets(tripletlist_ke.begin(), tripletlist_ke.end());
  }
  mat_seig.makeCompressed();
  mat_ke.makeCompressed();
};

template <class model1d>
spectral_element::spectral_element(const model1d &inp_model, double maxstep,
                                   int NQ, int l) {
  _q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(NQ);
  bool incball = false;
  _mesh = Radial_Tools::RadialMesh(inp_model, _q, maxstep, 1.2, incball);

  // std::cout << "Num elements: " << _mesh.NumberOfElements() << "\n";
  // std::cout << "Outer radius: " << _mesh.OuterRadius() << "\n";
  // std::cout << "Planet radius: " << _mesh.PlanetRadius() << "\n";
  // std::cout << "HELLO!\n";
};

void
spectral_element::CalculateEigenfrequencies(int N) {
  // initiate matrix multiplicaiton wrapper
  SparseGenMatProd<double> op_seig(mat_seig);

  // get nc
  long int maxn;
  if ((5 * N) > mat_seig.rows()) {
    maxn = mat_seig.rows();
  } else {
    maxn = 5 * N;
  }

  // initiate eigensolver
  GenEigsSolver<SparseGenMatProd<double>> eigs(op_seig, N, maxn);
  eigs.init();
  int nconv_seig = eigs.compute(SortRule::SmallestMagn);

  if (eigs.info() == CompInfo::Successful) {
    std::cout << "Successful\n";
    evalues_seig = eigs.eigenvalues();
    evectors.resize(eigs.eigenvectors().rows() + 1, eigs.eigenvectors().cols());
    evectors.block(0, 0, 1, evectors.cols()) =
        Eigen::MatrixXcd::Zero(1, evectors.cols());
    evectors.block(1, 0, evectors.rows() - 1, evectors.cols()) =
        eigs.eigenvectors();
    _calc_eig = true;
    // evectors = eigs.eigenvectors();
  } else {
    std::cout << "Unsuccessful\n";
  }
};

auto
spectral_element::evector_deriv(int i) const {
  std::vector<std::vector<std::complex<double>>> ret_vec =
      std::vector<std::vector<std::complex<double>>>(
          _NE, std::vector<std::complex<double>>(_q.N(), 0.0));
  auto pleg =
      Interpolation::LagrangePolynomial(_q.Points().begin(), _q.Points().end());
  std::vector<std::vector<double>> vec_lagderiv(
      _q.N(), std::vector<double>(_q.N(), 0.0));
  for (int idx_i = 0; idx_i < _q.N(); ++idx_i) {
    for (int idx_j = 0; idx_j < _q.N(); ++idx_j) {
      vec_lagderiv[idx_i][idx_j] = pleg.Derivative(idx_i, _q.X(idx_j));
    }
  }
  for (int idx = 0; idx < _NE; ++idx) {
    double elem_width = x_nodes[idx][_q.N() - 1] - x_nodes[idx][0];
    for (int idxq = 0; idxq < _q.N(); ++idxq) {

      std::complex<double> tmp = 0.0;

      for (int idxq2 = 0; idxq2 < _q.N(); ++idxq2) {
        std::size_t ovidx = overallidxfinal(idx, idxq2);
        tmp += evectors(ovidx, i) * vec_lagderiv[idxq2][idxq];
      }
      tmp *= 2.0 / elem_width;
      ret_vec[idx][idxq] = tmp;
    }
  }
  return ret_vec;
};

auto
spectral_element::efrequencies() const {
  Eigen::VectorXd ret_eig(evalues_seig.rows());
  for (int i = 0; i < evalues_seig.rows(); ++i) {
    ret_eig(evalues_seig.rows() - 1 - i) = std::sqrt(std::abs(evalues_seig(i)));
  }
  return ret_eig;
};

void
spectral_element::augment_basis_calculate() {
  std::size_t nrows = mat_ke.cols() + 1;
  std::size_t ncols = vec_layer.back() + 1;
  _augbasis = Eigen::MatrixXd::Zero(nrows, ncols);
  Eigen::SparseLU<Eigen::SparseMatrix<double>> _solver;
  _solver.analyzePattern(mat_ke);
  _solver.factorize(mat_ke);

  for (int idx = 0; idx < vec_layer.back() + 1; ++idx) {

    // determine force vector
    Eigen::VectorXd vec_force = Eigen::VectorXd::Zero(mat_ke.cols());
    std::size_t ovidx;
    int idxe = 0;
    if (idx < vec_layer.back()) {
      for (int idx_e = 0; idx_e < vec_layer.size(); ++idx_e) {
        if (vec_layer[idx_e] == idx + 1) {
          idxe = idx_e;

          break;
        }
      }
      ovidx = overallidx(idxe, 0);
    } else {
      ovidx = mat_ke.cols() - 1;
    }
    vec_force(ovidx) = -1.0;
    _augbasis.block(1, idx, mat_ke.cols(), 1) = _solver.solve(vec_force);
  }

  _augcalc = true;
};

auto
spectral_element::augment_deriv(int i) const {
  std::vector<std::vector<std::complex<double>>> ret_vec =
      std::vector<std::vector<std::complex<double>>>(
          _NE, std::vector<std::complex<double>>(_q.N(), 0.0));
  auto pleg =
      Interpolation::LagrangePolynomial(_q.Points().begin(), _q.Points().end());
  std::vector<std::vector<double>> vec_lagderiv(
      _q.N(), std::vector<double>(_q.N(), 0.0));
  for (int idx_i = 0; idx_i < _q.N(); ++idx_i) {
    for (int idx_j = 0; idx_j < _q.N(); ++idx_j) {
      vec_lagderiv[idx_i][idx_j] = pleg.Derivative(idx_i, _q.X(idx_j));
    }
  }
  for (int idx = 0; idx < _NE; ++idx) {
    double elem_width = x_nodes[idx][_q.N() - 1] - x_nodes[idx][0];
    for (int idxq = 0; idxq < _q.N(); ++idxq) {

      std::complex<double> tmp = 0.0;

      for (int idxq2 = 0; idxq2 < _q.N(); ++idxq2) {
        std::size_t ovidx = overallidxfinal(idx, idxq2);
        tmp += _augbasis(ovidx, i) * vec_lagderiv[idxq2][idxq];
      }
      tmp *= 2.0 / elem_width;
      ret_vec[idx][idxq] = tmp;
    }
  }
  return ret_vec;
};

Eigen::VectorXd
spectral_element::augment_basis(int i) const {
  assert((i < vec_layer.back() + 1) && "Interface not in model");
  Eigen::VectorXd vec_force = Eigen::VectorXd::Zero(mat_ke.cols());
  int idxe = 0;
  std::size_t ovidx;
  if (i == (vec_layer.back())) {
    ovidx = mat_ke.cols() - 1;
  } else {
    for (int idx = 0; idx < vec_layer.size(); ++idx) {
      if (vec_layer[idx] == i + 1) {
        idxe = idx;
        break;
      }
    }
    ovidx = overallidx(idxe, 0);
  }
  vec_force(ovidx) = -1.0;

  Eigen::SparseLU<Eigen::SparseMatrix<double>> _solver;
  _solver.analyzePattern(mat_ke);
  _solver.factorize(mat_ke);
  Eigen::VectorXd vec_output = Eigen::VectorXd::Zero(mat_ke.cols() + 1);
  vec_output.tail(mat_ke.cols()) = _solver.solve(vec_force);
  return vec_output;
};
}   // namespace Toroidal

#endif