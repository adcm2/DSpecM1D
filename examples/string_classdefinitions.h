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

using namespace Spectra;

// functions to use
double
newton(const std::function<double(double)> &f,
       const std::function<double(double)> &fd, double x0, int maxit = 10) {
  double xi, xf;
  xi = x0;
  for (int idx = 0; idx < maxit; ++idx) {
    xf = xi - f(xi) / fd(xi);
    if (std::abs(f(xf)) < std::pow(10.0, -14.0)) {
      break;
    };
    xi = xf;
  }
  return xf;
};

double
orthoeigen(double mu, double wn, double x) {
  double x0 = std::sqrt(2.0 / (1.0 + mu * std::cos(wn) * std::cos(wn)));
  return x0 * std::sin(wn * x);
};

double
orthoeigenderiv(double mu, double wn, double x) {
  double x0 = std::sqrt(2.0 / (1.0 + mu * std::cos(wn) * std::cos(wn)));
  return x0 * wn * std::cos(wn * x);
};

double
innerproduct(const std::function<double(double)> &f1,
             const std::function<double(double)> &f2,
             const GaussQuad::Quadrature1D<double> &q) {
  int npoints = 101;
  double xw = 1.0 / (npoints - 1);
  double totint = 0.0;
  for (int idx = 0; idx < npoints - 1; ++idx) {
    double xl = idx * xw;
    double xu = (idx + 1) * xw;
    double inttmp = 0.0;
    for (int idxq = 0; idxq < q.N(); ++idxq) {
      double xval = 0.5 * (xl + xu) + 0.5 * xw * q.X(idxq);
      inttmp += q.W(idxq) * f1(xval) * f2(xval);
    }
    totint += inttmp * xw / 2.0;
  }
  return totint;
};

class eigenfunction_catalogue {
public:
  eigenfunction_catalogue() {};
  eigenfunction_catalogue(int num_eig, double mu, double rho = 1.0)
      : _num_eig{num_eig}, _mu{mu}, _rho{rho} {
    auto fval = [mu = _mu](double x) {
      return std::cos(x) / std::sin(x) + 1 / (mu * x);
    };

    auto fdval = [mu = _mu](double x) {
      return -1.0 / (std::sin(x) * std::sin(x)) - 1 / (mu * x * x);
    };

    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////
    // solving for wn
    _vec_wn.resize(_num_eig);
    _vec_a.resize(_num_eig);

    for (int i = 0; i < _num_eig; ++i) {
      // eigenfrequencies
      double xtmp = (2 * i + 1) * pi_db / 2.0;
      _vec_wn[i] = newton(fval, fdval, xtmp);
      _vec_a[i] = std::sqrt(
          2.0 / (1.0 + mu * std::cos(_vec_wn[i]) * std::cos(_vec_wn[i])));
    }
  };

  //////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////
  int NumberOfEigenfunctions() const { return _num_eig; };
  double Mu() const { return _mu; };
  double w(int i) const { return std::sqrt(_mu / _rho) * _vec_wn[i]; };
  auto f(int i) const {
    auto tmpf = [&wn = _vec_wn, &a = _vec_a, i](double x) {
      return a[i] * std::sin(wn[i] * x);
    };
    return tmpf;
  };
  auto fd(int i) const {
    auto tmpf = [&wn = _vec_wn, &a = _vec_a, i](double x) {
      return a[i] * wn[i] * std::cos(wn[i] * x);
    };
    return tmpf;
  };
  double basis_expansion(std::vector<double> &vec_a, double x) const {
    double tmp = 0.0;
    for (int i = 0; i < _num_eig; ++i) {
      tmp += vec_a[i] * f(i)(x);
    }
    return tmp;
  };

  // function, derivative and expansion at a point
  auto f(int i, std::vector<double> &x) const {
    std::vector<double> retvec(x.size());
    for (int idx = 0; idx < x.size(); ++idx) {
      retvec[idx] = f(i)(x[idx]);
    }
    return retvec;
  };

  auto fd(int i, std::vector<double> &x) const {
    std::vector<double> retvec(x.size());
    for (int idx = 0; idx < x.size(); ++idx) {
      retvec[idx] = fd(i)(x[idx]);
    }
    return retvec;
  };

  auto basis_expansion(std::vector<double> &vec_a,
                       std::vector<double> &x) const {
    std::vector<double> retvec(x.size());
    for (int i = 0; i < x.size(); ++i) {
      retvec[i] = basis_expansion(vec_a, x[i]);
    }
    return retvec;
  };

  auto basis_expansion_augment(std::vector<double> &vec_a,
                               std::vector<double> &x,
                               const std::function<double(double)> &f1) const {
    std::vector<double> retvec(x.size());
    for (int i = 0; i < x.size(); ++i) {
      retvec[i] = basis_expansion(vec_a, x[i]) + f1(x[i]);
    }
    return retvec;
  };

private:
  std::vector<double> _vec_wn, _vec_a;
  double _mu, _rho;
  int _num_eig;
  double pi_db = 3.1415926535897932;
};

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
auto
force(const std::function<double(double)> &f1,
      const eigenfunction_catalogue &cat_eig) {
  int N = cat_eig.NumberOfEigenfunctions();
  std::vector<double> vec_coeff(N);
  auto q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(10);
  for (int i = 0; i < N; ++i) {
    vec_coeff[i] = innerproduct(f1, cat_eig.f(i), q);
  }

  return vec_coeff;
};

auto
force_augment(const std::function<double(double)> &f1,
              const std::function<double(double)> &f2,
              const eigenfunction_catalogue &cat_eig) {
  auto f3 = [&f1, &f2](double x) { return f1(x) - f2(x); };

  return force(f3, cat_eig);
};

///////////////////////////////////////////////////////////////////////////////
// class defining the problem
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

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// spectral element solver
class spectral_element {
public:
  spectral_element(const model_string &, double, int);
  spectral_element(double mu, int NE, int NQ) : _mu{mu}, _NE{NE} {

    // quadrature and Lagrange polynomials
    _q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(NQ);
    auto pleg = Interpolation::LagrangePolynomial(_q.Points().begin(),
                                                  _q.Points().end());

    // vectors for filling out matrix
    mat_seig.resize((_q.N() - 1) * _NE, (_q.N() - 1) * _NE);
    mat_ke.resize((_q.N() - 1) * _NE, (_q.N() - 1) * _NE);
    std::vector<double> vec_nn((_q.N() - 1) * _NE, 0.0),
        vec_lm1((_q.N() - 1) * _NE, 0.0);
    {
      // generate elements and nodes
      x_elem.resize(_NE + 1);
      double elem_width = 1.0 / _NE;
      std::generate(
          x_elem.begin(), x_elem.end(),
          [n = 0, &elem_width]() mutable { return n++ * elem_width; });

      // fill out nodes
      x_nodes.resize(_NE);
      for (int i = 0; i < _NE; ++i) {
        std::vector<double> xtmp(_q.N());
        double xl = x_elem[i];
        double xu = x_elem[i + 1];
        double xd = (xu - xl) / 2.0;
        double xp = (xl + xu) / 2.0;
        for (int k = 0; k < _q.N(); ++k) {
          xtmp[k] = xp + xd * _q.X(k);
        }
        x_nodes[i] = xtmp;
      }
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

      // form matrix, using tripletlist vector
      //  tripletlist vector
      using T = Eigen::Triplet<double>;
      std::vector<T> tripletlist3, tripletlist_ke;
      tripletlist3.reserve(_q.N() * _q.N() * _NE);
      tripletlist_ke.reserve(_q.N() * _q.N() * _NE);

      {
        int idxtl = 0;
        // first element
        for (int i = 1; i < _q.N(); ++i) {
          vec_nn[idxtl + i - 1] += elem_width / 2.0 * _q.W(i);
        }
        idxtl += _q.N() - 2;

        // middle and last elements
        for (int idxe = 1; idxe < _NE; ++idxe) {
          for (int i = 0; i < _q.N(); ++i) {
            vec_nn[idxtl + i] += elem_width / 2.0 * _q.W(i);
          }
          idxtl += _q.N() - 1;
        }
      }

      // find lambda^{-1} (ie diagonal)
      for (int i = 0; i < vec_nn.size(); ++i) {
        vec_lm1[i] = 1.0 / vec_nn[i];
      }
      // get tripletlist for matrix
      {
        int idxtl = 0;
        // first element
        // loop over basis
        for (int i = 1; i < _q.N(); ++i) {
          for (int j = 1; j < _q.N(); ++j) {
            double tmp = 0.0;

            // loop over points
            for (int k = 0; k < _q.N(); ++k) {
              tmp += _q.W(k) * mat_gauss_deriv[k][i][j];
            }
            tmp *= 2.0 / elem_width;

            long int idx_x = idxtl + i - 1;
            long int idx_y = idxtl + j - 1;
            tripletlist3.push_back(T(idx_x, idx_y, _mu * tmp * vec_lm1[idx_x]));
            tripletlist_ke.push_back(T(idx_x, idx_y, _mu * tmp));
          }
        }
        idxtl += _q.N() - 2;

        // middle and last elements
        // over elements
        for (int idxe = 1; idxe < _NE; ++idxe) {

          // over basis
          for (int i = 0; i < _q.N(); ++i) {
            for (int j = 0; j < _q.N(); ++j) {
              double tmp = 0.0;

              // over nodes
              for (int k = 0; k < _q.N(); ++k) {
                tmp += _q.W(k) * mat_gauss_deriv[k][i][j];
              }
              tmp *= 2.0 / elem_width;

              tripletlist3.push_back(
                  T(idxtl + i, idxtl + j, _mu * tmp * vec_lm1[idxtl + i]));
              tripletlist_ke.push_back(T(idxtl + i, idxtl + j, _mu * tmp));
            }
          }
          idxtl += _q.N() - 1;
        }
        // add to last
        tripletlist3.push_back(T(idxtl, idxtl, vec_lm1[idxtl]));
        tripletlist_ke.push_back(T(idxtl, idxtl, 1.0));
      }
      mat_seig.setFromTriplets(tripletlist3.begin(), tripletlist3.end());
      mat_ke.setFromTriplets(tripletlist_ke.begin(), tripletlist_ke.end());
    }
    mat_seig.makeCompressed();
    mat_ke.makeCompressed();
  };

  void CalculateEigenfrequencies(int N) {
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
      evectors.resize(eigs.eigenvectors().rows() + 1,
                      eigs.eigenvectors().cols());
      evectors.block(0, 0, 1, evectors.cols()) =
          Eigen::MatrixXcd::Zero(1, evectors.cols());
      evectors.block(1, 0, evectors.rows() - 1, evectors.cols()) =
          eigs.eigenvectors();
      _calc_eig = true;
      // evectors = eigs.eigenvectors();
    } else {
      std::cout << "Unsuccessful\n";
    }

    // for (int i = 4; i > -1; i--) {
    //    std::cout << std::setprecision(15) << "Eigenfrequency:\n"
    //              << std::sqrt(evalues_seig(i)) << std::endl;
    // }

    // std::cout << "Rows: " << evectors.rows() << ", cols: " <<
    // evectors.cols()
    //           << "\n";
  };

  auto evector_value(int i) const {
    return evectors.block(0, i, evectors.rows(), 1);
  };
  auto evector_deriv(int i) const {
    std::vector<std::vector<std::complex<double>>> ret_vec =
        std::vector<std::vector<std::complex<double>>>(
            _NE, std::vector<std::complex<double>>(_q.N(), 0.0));
    auto pleg = Interpolation::LagrangePolynomial(_q.Points().begin(),
                                                  _q.Points().end());
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

  Eigen::MatrixXcd evectors_all() const { return evectors; };

  auto evalues() const { return evalues_seig; };
  auto efrequencies() const {
    Eigen::VectorXd ret_eig(evalues_seig.rows());
    for (int i = 0; i < evalues_seig.rows(); ++i) {
      ret_eig(evalues_seig.rows() - 1 - i) =
          std::sqrt(std::abs(evalues_seig(i)));
    }
    return ret_eig;
  };

  // augmentation basis function calculation
  void augment_basis_calculate() {
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
        // std::cout << "idx: " << idx << ", idxe: " << idxe << ", ovidx: " <<
        // ovidx << "\n";
      } else {
        ovidx = mat_ke.cols() - 1;
        // std::cout << "idx: " << idx << ", ovidx: " << ovidx << "\n";
        // std::cout << "length: " << vec_force.size() <<"\n";
      }
      vec_force(ovidx) = -1.0;
      _augbasis.block(1, idx, mat_ke.cols(), 1) = _solver.solve(vec_force);
    }

    _augcalc = true;
  };

  auto augment_deriv(int i) const {
    std::vector<std::vector<std::complex<double>>> ret_vec =
        std::vector<std::vector<std::complex<double>>>(
            _NE, std::vector<std::complex<double>>(_q.N(), 0.0));
    auto pleg = Interpolation::LagrangePolynomial(_q.Points().begin(),
                                                  _q.Points().end());
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

  auto augment_basis() const { return _augbasis; };

  auto is_augmented() const { return _augcalc; };

  Eigen::VectorXd augment_basis(int i) const {
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

  auto xvalues() const { return x_valnodes; };
  auto xnodes() const { return x_nodes; };
  auto xelem() const { return x_elem; };

  auto &layers() const { return vec_layer; };

private:
  int _NE;
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

  std::size_t overallidx(int idxe, int idxq) const {
    return (_q.N() - 1) * idxe + idxq - 1;
  };

  std::size_t overallidxfinal(int idxe, int idxq) const {
    return (_q.N() - 1) * idxe + idxq;
  };
};

spectral_element::spectral_element(
    const model_string &inp_model, double maxwidth,
    int NQ) {   // quadrature and Lagrange polynomials
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
  // for (auto idx : x_elem) {
  //    std::cout << idx << "\n";
  // }
  // for (auto idx : x_nodes) {
  //    for (auto idx1 : idx) {
  //       std::cout << idx1 << "\n";
  //    }
  // }
  // _NE = 50;
  // vectors for filling out matrix
  mat_seig.resize((_q.N() - 1) * _NE, (_q.N() - 1) * _NE);
  mat_ke.resize((_q.N() - 1) * _NE, (_q.N() - 1) * _NE);
  std::vector<double> vec_nn((_q.N() - 1) * _NE, 0.0),
      vec_lm1((_q.N() - 1) * _NE, 0.0);
  {
    // // generate elements and nodes
    // // x_elem.resize(_NE + 1);
    // double elem_width = 1.0 / _NE;
    // std::generate(
    //     x_elem.begin(), x_elem.end(),
    //     [n = 0, &elem_width]() mutable { return n++ * elem_width; });

    // // fill out nodes
    // x_nodes.resize(_NE);
    // for (int i = 0; i < _NE; ++i) {
    //    std::vector<double> xtmp(_q.N());
    //    double xl = x_elem[i];
    //    double xu = x_elem[i + 1];
    //    double xd = (xu - xl) / 2.0;
    //    double xp = (xl + xu) / 2.0;
    //    for (int k = 0; k < _q.N(); ++k) {
    //       xtmp[k] = xp + xd * _q.X(k);
    //    }
    //    x_nodes[i] = xtmp;
    // }
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
                                 inp_model.rho(vec_layer[0])(x_nodes[0][i]);
      }
      idxtl += _q.N() - 2;

      // middle and last elements
      for (int idxe = 1; idxe < _NE; ++idxe) {
        elem_width = x_elem[idxe + 1] - x_elem[idxe];
        for (int i = 0; i < _q.N(); ++i) {
          vec_nn[idxtl + i] += elem_width / 2.0 * _q.W(i) *
                               inp_model.rho(vec_layer[idxe])(x_nodes[idxe][i]);
          ;
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
      int idxtl = 0;
      // first element
      // loop over basis
      {
        double elem_width = x_elem[1] - x_elem[0];
        for (int i = 1; i < _q.N(); ++i) {
          for (int j = 1; j < _q.N(); ++j) {
            double tmp = 0.0;

            // loop over points
            for (int k = 0; k < _q.N(); ++k) {
              tmp += _q.W(k) * mat_gauss_deriv[k][i][j] *
                     inp_model.mu(vec_layer[0])(x_nodes[0][k]);
            }
            tmp *= 2.0 / elem_width;

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
        // over basis
        for (int i = 0; i < _q.N(); ++i) {
          for (int j = 0; j < _q.N(); ++j) {
            double tmp = 0.0;

            // over nodes
            for (int k = 0; k < _q.N(); ++k) {
              tmp += _q.W(k) * mat_gauss_deriv[k][i][j] *
                     inp_model.mu(vec_layer[idxe])(x_nodes[idxe][k]);
            }
            tmp *= 2.0 / elem_width;

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
      tripletlist3.push_back(T(idxtl, idxtl, vec_lm1[idxtl]));
      tripletlist_ke.push_back(T(idxtl, idxtl, 1.0));

      int idxoverall = 0;

      int qmax = _q.N() - 1;
    }
    mat_seig.setFromTriplets(tripletlist3.begin(), tripletlist3.end());
    mat_ke.setFromTriplets(tripletlist_ke.begin(), tripletlist_ke.end());
  }
  mat_seig.makeCompressed();
  mat_ke.makeCompressed();
};

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

class augment_function_cat {
private:
public:
  augment_function_cat() {};
  auto f(int i) {
    auto mylambda = [i](double x) {
      double pi_db = 3.1415926535897932;
      return std::sin((i + 1.0) * pi_db * x);
    };
    return mylambda;
  };
  auto fd(int i) {
    auto mylambda = [i](double x) {
      double pi_db = 3.1415926535897932;
      return pi_db * (i + 1.0) * std::cos((i + 1.0) * pi_db * x);
    };
    return mylambda;
  };
};

// spectral element solver
class galerkin_solver {

private:
  using VD = std::vector<double>;
  using VVD = std::vector<VD>;
  using VVVD = std::vector<VVD>;

  int _NE, _NQ, _numbasis, _numeig;
  double _mu;
  GaussQuad::Quadrature1D<double> _q;
  VD x_elem, x_width;
  VVD x_nodes;
  VVVD val_ef, val_efd;
  Eigen::VectorXcd vec_eig;
  Eigen::MatrixXcd mat_eig;
  Eigen::MatrixXd mat_lhs, mat_rhs;
  eigenfunction_catalogue _eigcat;
  augment_function_cat _augcat;
  bool _augmented;

public:
  galerkin_solver(int neig, double mu0, double mu, int NE, int NQ,
                  bool to_augment = false)
      : _NE{NE}, _NQ{NQ}, _mu{mu}, _numeig{neig}, _augmented{to_augment} {
    _eigcat = eigenfunction_catalogue(neig, mu0);

    // quadrature and Lagrange polynomials
    _q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(NQ);

    // number of eigenfunctions
    if (to_augment) {
      _numbasis = neig + 1;
    } else {
      _numbasis = neig;
    }
    {   // generate elements and nodes
      x_elem.resize(_NE + 1);

      double elem_width = 1.0 / _NE;
      x_width = VD(_NE, elem_width);
      std::generate(
          x_elem.begin(), x_elem.end(),
          [n = 0, &elem_width]() mutable { return n++ * elem_width; });

      // fill out nodes
      x_nodes.resize(_NE);
      for (int i = 0; i < _NE; ++i) {
        std::vector<double> xtmp(_q.N());
        double xl = x_elem[i];
        double xu = x_elem[i + 1];
        double xd = (xu - xl) / 2.0;
        double xp = (xl + xu) / 2.0;
        for (int k = 0; k < _q.N(); ++k) {
          xtmp[k] = xp + xd * _q.X(k);
        }
        x_nodes[i] = xtmp;
      }
    }

    //////////////////////////////////////////////////////////////
    // get values of the eigenfunctions at the nodes
    {
      val_ef = VVVD(_numbasis, VVD(NE, VD(NQ, 0.0)));
      val_efd = val_ef;
      for (int idx_f = 0; idx_f < neig; ++idx_f) {
        for (int idx_e = 0; idx_e < _NE; ++idx_e) {
          val_ef[idx_f][idx_e] = _eigcat.f(idx_f, x_nodes[idx_e]);
          val_efd[idx_f][idx_e] = _eigcat.fd(idx_f, x_nodes[idx_e]);
        }
      }
      if (to_augment) {
        for (int i = neig; i < _numbasis; ++i) {
          for (int idx_e = 0; idx_e < _NE; ++idx_e) {
            for (int idx_q = 0; idx_q < _NQ; ++idx_q) {
              val_ef[i][idx_e][idx_q] =
                  _augcat.f(i - neig)(x_nodes[idx_e][idx_q]);
              val_efd[i][idx_e][idx_q] =
                  _augcat.fd(i - neig)(x_nodes[idx_e][idx_q]);
            }
          }
        }
      }

      //////////////////////////////////////////////////////////////
      // integrate to get matrices
      mat_lhs.resize(_numbasis, _numbasis);
      mat_rhs.resize(_numbasis, _numbasis);
      {
        // mass matrix
        for (int idx_i = 0; idx_i < _numbasis; ++idx_i) {
          for (int idx_j = 0; idx_j < _numbasis; ++idx_j) {
            double totint = 0.0;
            for (int i = 0; i < _NE; ++i) {
              double tmp = 0.0;
              for (int k = 0; k < _q.N(); ++k) {
                tmp += _q.W(k) * val_ef[idx_i][i][k] * val_ef[idx_j][i][k];
              }
              tmp *= x_width[i] / 2.0;
              totint += tmp;
            }
            mat_rhs(idx_i, idx_j) = totint;
          }
        }
      }
      {
        // ke matrix
        for (int idx_i = 0; idx_i < _numbasis; ++idx_i) {
          for (int idx_j = 0; idx_j < _numbasis; ++idx_j) {
            double totint = 0.0;
            for (int i = 0; i < _NE; ++i) {
              double tmp = 0.0;
              for (int k = 0; k < _q.N(); ++k) {
                tmp += _q.W(k) * val_efd[idx_i][i][k] * val_efd[idx_j][i][k];
              }
              tmp *= x_width[i] / 2.0;
              totint += tmp;
            }
            mat_lhs(idx_i, idx_j) =
                totint * mu +
                val_ef[idx_i][_NE - 1][_q.N() - 1] *
                    val_ef[idx_j][_NE - 1][_q.N() - 1] +
                mu * val_ef[idx_i][0][0] * val_efd[idx_j][0][0];
          }
        }
      }

      // eigenvalues
      Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> ges;
      ges.compute(mat_lhs, mat_rhs);
      Eigen::VectorXcd tmp_eig = ges.eigenvalues();
      Eigen::MatrixXcd tmp_vec = ges.eigenvectors();
      vec_eig.resize(tmp_eig.rows());
      mat_eig.resize(tmp_vec.rows(), tmp_vec.cols());

      // find index permutation for ascending order
      std::vector<int> vec_idx(tmp_eig.rows());
      for (int i = 0; i < vec_idx.size(); ++i) {
        vec_idx[i] = i;
      }
      auto myfun = [&vec_comp = tmp_eig](int a, int b) {
        return std::abs(vec_comp(a)) < std::abs(vec_comp(b));
      };
      std::sort(vec_idx.begin(), vec_idx.end(), myfun);

      // get final result in order of increasing absolute magnitude
      for (int i = 0; i < vec_eig.rows(); ++i) {
        vec_eig(i) = tmp_eig(vec_idx[i]);
      }
      for (int i = 0; i < tmp_vec.cols(); ++i) {
        mat_eig.block(0, i, tmp_vec.rows(), 1) =
            tmp_vec.block(0, vec_idx[i], tmp_vec.rows(), 1);
      }
    };
  }
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  Eigen::VectorXd efrequencies() {
    Eigen::VectorXd ret_eig(vec_eig.rows());
    for (int i = 0; i < vec_eig.rows(); ++i) {
      ret_eig(i) = std::sqrt(std::abs(vec_eig(i)));
    }
    return ret_eig;
  };
  Eigen::VectorXcd evalues() { return vec_eig; };
  Eigen::MatrixXcd evectors() { return mat_eig; };

  // physical outputs

  auto evector_physical(int i) {
    // Eigen::VectorXcd ret_vec = Eigen::VectorXcd::Zero(mat_eig.rows());
    std::vector<std::complex<double>> ret_vec;
    for (int idx_e = 0; idx_e < _NE; ++idx_e) {
      for (int idx_q = 0; idx_q < _NQ - 1; ++idx_q) {
        std::complex<double> tmp = 0.0;
        for (int idx = 0; idx < _numbasis; ++idx) {
          tmp += val_ef[idx][idx_e][idx_q] * mat_eig(idx, i);
        }
        ret_vec.push_back(tmp);
      }
    }
    {
      std::complex<double> tmp = 0.0;
      for (int idx = 0; idx < mat_eig.rows(); ++idx) {
        tmp += val_ef[idx][_NE - 1][_NQ - 1] * mat_eig(idx, i);
      }
      ret_vec.push_back(tmp);
    }
    return ret_vec;
  };

  auto evector_physical(int i, std::vector<double> &x) {
    // Eigen::VectorXcd ret_vec = Eigen::VectorXcd::Zero(mat_eig.rows());
    std::vector<std::complex<double>> ret_vec(x.size(), 0.0);
    for (int idx = 0; idx < x.size(); ++idx) {
      for (int idx_e = 0; idx_e < _numeig; ++idx_e) {
        ret_vec[idx] += mat_eig(idx_e, i) * _eigcat.f(idx_e)(x[idx]);
      }
      if (_augmented) {
        for (int idx_b = 0; idx_b < _numbasis - _numeig; ++idx_b) {
          ret_vec[idx] +=
              mat_eig(idx_b + _numeig, i) * _augcat.f(idx_b)(x[idx]);
        }
      }
    }
    return ret_vec;
  };
};

class gsab {
private:
  Eigen::VectorXcd _evalues;
  Eigen::MatrixXcd _evectors, mat_inertia, mat_ke;
  std::vector<std::vector<double>> _nodes;
  GaussQuad::Quadrature1D<double> _q;
  std::size_t _numbasis, _numeig;
  std::size_t index_evec(std::size_t i, std::size_t k) {
    return (_q.N() - 1) * i + k;
  };
  Eigen::VectorXcd vec_eig;
  Eigen::MatrixXcd mat_eig;

  using VD = std::vector<std::complex<double>>;
  using VVD = std::vector<VD>;
  using VVVD = std::vector<VVD>;

  VVVD val_ef, val_efd;

public:
  gsab() {};
  gsab(const spectral_element &, const model_string &, const model_string &,
       int);
  gsab(const spectral_element &, const model_string &, const model_string &,
       bool);

  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  Eigen::VectorXd efrequencies() {
    Eigen::VectorXd ret_eig(vec_eig.rows());
    for (int i = 0; i < vec_eig.rows(); ++i) {
      ret_eig(i) = std::sqrt(std::abs(vec_eig(i)));
    }
    return ret_eig;
  };
  Eigen::VectorXcd evalues() { return vec_eig; };
  Eigen::MatrixXcd evectors() { return mat_eig; };
};

gsab::gsab(const spectral_element &inp_basis, const model_string &mod_init,
           const model_string &mod_new, int idxrand)
    : _evalues{inp_basis.evalues()}, _evectors{inp_basis.evectors_all()} {
  // check that the new model has the same discontinuities
  if (mod_init.NumberOfLayers() == mod_new.NumberOfLayers()) {
    for (int idx = 0; idx < mod_init.NumberOfLayers(); ++idx) {
      if ((mod_init.UpperRadius(idx) - mod_new.UpperRadius(idx)) >
          std::pow(10.0, -12.0)) {
        std::cout << "Models have different discontinuities!\n";
      }
    }
  } else {
    std::cout << "Models have different number of layers!\n";
  }

  // checked that the models have the same discontinuities
  // we now generate the matrices
  _nodes = inp_basis.xnodes();
  _q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(_nodes[0].size());
  auto pleg =
      Interpolation::LagrangePolynomial(_q.Points().begin(), _q.Points().end());

  // generating matrix of derivative values for Lagrange polynomials
  std::vector<std::vector<std::vector<double>>> mat_gauss_deriv;
  std::vector<std::vector<double>> mat_lag_deriv;
  //  std::cout << "HELLO 1\n";
  mat_gauss_deriv.reserve(_q.N());
  for (int idxk = 0; idxk < _q.N(); ++idxk) {
    std::vector<std::vector<double>> mat_tmp(_q.N(),
                                             std::vector<double>(_q.N()));
    std::vector<double> tmp_lagderiv(_q.N(), 0.0);
    for (int idxi = 0; idxi < _q.N(); ++idxi) {
      std::vector<double> vec_tmp(_q.N());

      // derivative of i-th function at k-th point:
      tmp_lagderiv[idxi] = pleg.Derivative(idxi, _q.X(idxk));
      for (int idxj = 0; idxj < _q.N(); ++idxj) {
        vec_tmp[idxj] = pleg.Derivative(idxi, _q.X(idxk)) *
                        pleg.Derivative(idxj, _q.X(idxk));
      }
      mat_tmp[idxi] = vec_tmp;
    }
    mat_gauss_deriv.push_back(mat_tmp);
    mat_lag_deriv.push_back(tmp_lagderiv);
  };

  //   lag deriv has:
  // vec[i][k] is derivative of k-th function at i-th point

  //  std::cout << "HELLO 2\n";
  //////////////////////////////////////////////////////////////
  // generate vectors containing the values and derivatives of the basis
  // functions
  int num_elements = _nodes.size();
  _numeig = _evectors.cols();

  if (inp_basis.is_augmented()) {
    // _numbasis = _numeig;
    _numbasis = _numeig + mod_new.NumberOfLayers();
    std::cout << "HELLO AGAIN!!\n";
  } else {
    _numbasis = _numeig;
  }
  val_ef = VVVD(_numbasis, VVD(num_elements, VD(_q.N(), 0.0)));
  val_efd = val_ef;
  //  std::cout << "HELLO 2.1\n";
  for (int idx_f = 0; idx_f < _numeig; ++idx_f) {
    auto tmp_deriv = inp_basis.evector_deriv(idx_f);
    for (int idx_e = 0; idx_e < num_elements; ++idx_e) {
      for (int idx_q = 0; idx_q < _q.N(); ++idx_q) {
        std::size_t ovidx = index_evec(idx_e, idx_q);
        val_ef[idx_f][idx_e][idx_q] = _evectors(ovidx, idx_f);

        // // evaluating derivative of function
        // std::complex<double> tmp = 0.0;
        // for (int idx = 0; idx < _q.N(); ++idx) {
        //   std::size_t ovidx1 = index_evec(idx_e, idx);
        //   tmp += _evectors(ovidx1, idx_f) * mat_lag_deriv[idx_q][idx];
        // }
        // tmp *= 2.0 / (_nodes[idx_e][_q.N() - 1] - _nodes[idx_e][0]);
        val_efd[idx_f][idx_e][idx_q] = tmp_deriv[idx_e][idx_q];
        // val_efd[idx_f][idx_e][idx_q] = tmp;
        //  val_efd[idx_f][idx_e]{idx_q} = _eigcat.fd(idx_f, x_nodes[idx_e]);
      }
    }
  }
  if (inp_basis.is_augmented()) {
    std::cout << "HELLO WHAT YOU SAYING\n";
    // if () {
    //   std::cout << "HELLO AGAIN!!\n";
    // }
    auto augbasis = inp_basis.augment_basis();
    int ovidx = _numeig;

    auto eig_vectors = inp_basis.augment_basis();
    // std::cout << "Rows: " << eig_vectors.rows()
    //           << ", cols: " << eig_vectors.cols() << "\n";
    for (int idx_l = 0; idx_l < mod_new.NumberOfLayers(); ++idx_l) {
      double layer_width =
          mod_new.UpperRadius(idx_l) - mod_new.LowerRadius(idx_l);
      // for (int iaug = 1; iaug < 3; ++iaug) {
      // double sc = 3.1415926535 * iaug;
      for (int idx_e = 0; idx_e < num_elements; ++idx_e) {
        // if (inp_basis.layers()[idx_e] == idx_l) {
        for (int idx_q = 0; idx_q < _q.N(); ++idx_q) {
          // double xtmp =
          //     (_nodes[idx_e][idx_q] - mod_new.LowerRadius(idx_l)) /
          //     layer_width;
          std::size_t idxfull = idx_e * (_q.N() - 1) + idx_q;
          // std::cout << "idxfull: " << idxfull << "\n";
          val_ef[ovidx][idx_e][idx_q] = eig_vectors(idxfull, idx_l);
          // val_efd[ovidx][idx_e][idx_q] =
          //     sc / layer_width * std::cos(sc * xtmp);
        }
        // }
      }
      // std::cout << "Test 1\n";
      val_efd[ovidx] = inp_basis.augment_deriv(idx_l);
      ++ovidx;
      // }
    }
  }

  //////////////////////////////////////////////////////////////
  // integrate to get matrices
  // now we need to generate the matrices
  //  _numbasis = _evectors.cols();
  std::cout << "Numbasis: " << _numbasis << "\n";
  std::cout << "evectors: " << _evalues.size() << "\n";
  mat_inertia.resize(_numbasis, _numbasis);
  mat_ke.resize(_numbasis, _numbasis);
  {
    // mass matrix
    // outer two for loops are for the different basis functions
    for (int idx_i = 0; idx_i < _numbasis; ++idx_i) {
      for (int idx_j = 0; idx_j < _numbasis; ++idx_j) {

        // integrate for each set of basis functions
        std::complex<double> totint = 0.0;
        for (int i = 0; i < _nodes.size(); ++i) {
          std::complex<double> tmp = 0.0;
          double elem_width = _nodes[i][_q.N() - 1] - _nodes[i][0];
          for (int k = 0; k < _q.N(); ++k) {
            std::size_t ovidx = index_evec(i, k);
            tmp += _q.W(k) * val_ef[idx_i][i][k] * val_ef[idx_j][i][k] *
                   mod_new.rho(inp_basis.layers()[i])(_nodes[i][k]);
          }
          tmp *= elem_width / 2.0;
          totint += tmp;
        }
        mat_inertia(idx_i, idx_j) = totint;
      }
    }
  }

  {
    // ke matrix
    for (int idx_i = 0; idx_i < _numbasis; ++idx_i) {
      for (int idx_j = 0; idx_j < _numbasis; ++idx_j) {
        std::complex<double> totint = 0.0;
        for (int i = 0; i < num_elements; ++i) {
          std::complex<double> tmp = 0.0;
          for (int k = 0; k < _q.N(); ++k) {
            tmp += _q.W(k) * val_efd[idx_i][i][k] * val_efd[idx_j][i][k] *
                   mod_new.mu(inp_basis.layers()[i])(_nodes[i][k]);
          }
          tmp *= (_nodes[i][_q.N() - 1] - _nodes[i][0]) / 2.0;
          totint += tmp;
        }
        mat_ke(idx_i, idx_j) =
            totint + val_ef[idx_i][num_elements - 1][_q.N() - 1] *
                         val_ef[idx_j][num_elements - 1][_q.N() - 1];
        //      +
        // mod_new.mu(0)(0.0) * val_ef[idx_i][0][0] * val_efd[idx_j][0][0]
      }
    }
  }

  // std::cout << std::setprecision(5) << mat_ke << "\n\n" << mat_inertia <<
  // "\n";
  //
  // setting up the problem, we need to find the inverse of the inertia matrix
  // and pre-multiply it through against the kinetic energy matrix
  Eigen::FullPivLU<Eigen::MatrixXcd> lu(mat_inertia);
  Eigen::MatrixXcd mat_decomp = lu.solve(mat_ke);
  // std::cout << std::setprecision(3) << mat_decomp << "\n\n";

  // eigendecompose
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(mat_decomp);
  vec_eig = es.eigenvalues();
  mat_eig = es.eigenvectors();

  // std::cout << "Size of vec_eig: " << vec_eig.rows() << "\n";
};

gsab::gsab(const spectral_element &inp_basis, const model_string &mod_init,
           const model_string &mod_new, bool aug_inc = false)
    : _evalues{inp_basis.evalues()}, _evectors{inp_basis.evectors_all()} {
  // check that the new model has the same discontinuities
  if (mod_init.NumberOfLayers() == mod_new.NumberOfLayers()) {
    for (int idx = 0; idx < mod_init.NumberOfLayers(); ++idx) {
      if ((mod_init.UpperRadius(idx) - mod_new.UpperRadius(idx)) >
          std::pow(10.0, -12.0)) {
        std::cout << "Models have different discontinuities!\n";
      }
    }
  } else {
    std::cout << "Models have different number of layers!\n";
  }

  // checked that the models have the same discontinuities
  // we now generate the matrices
  _nodes = inp_basis.xnodes();
  _q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(_nodes[0].size());
  auto pleg =
      Interpolation::LagrangePolynomial(_q.Points().begin(), _q.Points().end());

  // generating matrix of derivative values for Lagrange polynomials
  std::vector<std::vector<std::vector<double>>> mat_gauss_deriv;
  std::vector<std::vector<double>> mat_lag_deriv;
  //  std::cout << "HELLO 1\n";
  mat_gauss_deriv.reserve(_q.N());
  for (int idxk = 0; idxk < _q.N(); ++idxk) {
    std::vector<std::vector<double>> mat_tmp(_q.N(),
                                             std::vector<double>(_q.N()));
    std::vector<double> tmp_lagderiv(_q.N(), 0.0);
    for (int idxi = 0; idxi < _q.N(); ++idxi) {
      std::vector<double> vec_tmp(_q.N());

      // derivative of i-th function at k-th point:
      tmp_lagderiv[idxi] = pleg.Derivative(idxi, _q.X(idxk));
      for (int idxj = 0; idxj < _q.N(); ++idxj) {
        vec_tmp[idxj] = pleg.Derivative(idxi, _q.X(idxk)) *
                        pleg.Derivative(idxj, _q.X(idxk));
      }
      mat_tmp[idxi] = vec_tmp;
    }
    mat_gauss_deriv.push_back(mat_tmp);
    mat_lag_deriv.push_back(tmp_lagderiv);
  };

  //   lag deriv has:
  // vec[i][k] is derivative of k-th function at i-th point

  //  std::cout << "HELLO 2\n";
  //////////////////////////////////////////////////////////////
  // generate vectors containing the values and derivatives of the basis
  // functions
  int num_elements = _nodes.size();
  _numeig = _evectors.cols();
  if (aug_inc) {
    _numbasis = _numeig + 2 * mod_new.NumberOfLayers();
  } else {
    _numbasis = _numeig;
  }
  val_ef = VVVD(_numbasis, VVD(num_elements, VD(_q.N(), 0.0)));
  val_efd = val_ef;
  //  std::cout << "HELLO 2.1\n";
  for (int idx_f = 0; idx_f < _numeig; ++idx_f) {
    auto tmp_deriv = inp_basis.evector_deriv(idx_f);
    for (int idx_e = 0; idx_e < num_elements; ++idx_e) {
      for (int idx_q = 0; idx_q < _q.N(); ++idx_q) {
        std::size_t ovidx = index_evec(idx_e, idx_q);
        val_ef[idx_f][idx_e][idx_q] = _evectors(ovidx, idx_f);

        // // evaluating derivative of function
        // std::complex<double> tmp = 0.0;
        // for (int idx = 0; idx < _q.N(); ++idx) {
        //   std::size_t ovidx1 = index_evec(idx_e, idx);
        //   tmp += _evectors(ovidx1, idx_f) * mat_lag_deriv[idx_q][idx];
        // }
        // tmp *= 2.0 / (_nodes[idx_e][_q.N() - 1] - _nodes[idx_e][0]);
        val_efd[idx_f][idx_e][idx_q] = tmp_deriv[idx_e][idx_q];
        // val_efd[idx_f][idx_e][idx_q] = tmp;
        //  val_efd[idx_f][idx_e]{idx_q} = _eigcat.fd(idx_f, x_nodes[idx_e]);
      }
    }
  }
  if (aug_inc) {
    // if (inp_basis.is_augmented()) {
    //   std::cout << "HELLO AGAIN!!\n";
    // }
    int ovidx = _numeig;
    for (int idx_l = 0; idx_l < mod_new.NumberOfLayers(); ++idx_l) {
      double layer_width =
          mod_new.UpperRadius(idx_l) - mod_new.LowerRadius(idx_l);
      for (int iaug = 1; iaug < 3; ++iaug) {
        double sc = 3.1415926535 * iaug;
        for (int idx_e = 0; idx_e < num_elements; ++idx_e) {
          if (inp_basis.layers()[idx_e] == idx_l) {
            for (int idx_q = 0; idx_q < _q.N(); ++idx_q) {
              double xtmp =
                  (_nodes[idx_e][idx_q] - mod_new.LowerRadius(idx_l)) /
                  layer_width;
              val_ef[ovidx][idx_e][idx_q] = std::sin(sc * xtmp);
              val_efd[ovidx][idx_e][idx_q] =
                  sc / layer_width * std::cos(sc * xtmp);
            }
          }
        }

        ++ovidx;
      }
    }
  }
  //  std::cout << "HELLO 3\n";

  //////////////////////////////////////////////////////////////
  // integrate to get matrices
  // now we need to generate the matrices
  //  _numbasis = _evectors.cols();
  mat_inertia.resize(_numbasis, _numbasis);
  mat_ke.resize(_numbasis, _numbasis);
  {
    // mass matrix
    // outer two for loops are for the different basis functions
    for (int idx_i = 0; idx_i < _numbasis; ++idx_i) {
      for (int idx_j = 0; idx_j < _numbasis; ++idx_j) {

        // integrate for each set of basis functions
        std::complex<double> totint = 0.0;
        for (int i = 0; i < _nodes.size(); ++i) {
          std::complex<double> tmp = 0.0;
          double elem_width = _nodes[i][_q.N() - 1] - _nodes[i][0];
          for (int k = 0; k < _q.N(); ++k) {
            std::size_t ovidx = index_evec(i, k);
            tmp += _q.W(k) * val_ef[idx_i][i][k] * val_ef[idx_j][i][k] *
                   mod_new.rho(inp_basis.layers()[i])(_nodes[i][k]);
          }
          tmp *= elem_width / 2.0;
          totint += tmp;
        }
        mat_inertia(idx_i, idx_j) = totint;
      }
    }
  }

  {
    // ke matrix
    for (int idx_i = 0; idx_i < _numbasis; ++idx_i) {
      for (int idx_j = 0; idx_j < _numbasis; ++idx_j) {
        std::complex<double> totint = 0.0;
        for (int i = 0; i < num_elements; ++i) {
          std::complex<double> tmp = 0.0;
          for (int k = 0; k < _q.N(); ++k) {
            tmp += _q.W(k) * val_efd[idx_i][i][k] * val_efd[idx_j][i][k] *
                   mod_new.mu(inp_basis.layers()[i])(_nodes[i][k]);
          }
          tmp *= (_nodes[i][_q.N() - 1] - _nodes[i][0]) / 2.0;
          totint += tmp;
        }
        mat_ke(idx_i, idx_j) =
            totint + val_ef[idx_i][num_elements - 1][_q.N() - 1] *
                         val_ef[idx_j][num_elements - 1][_q.N() - 1];
        //      +
        // mod_new.mu(0)(0.0) * val_ef[idx_i][0][0] * val_efd[idx_j][0][0]
      }
    }
  }

  // setting up the problem, we need to find the inverse of the inertia matrix
  // and pre-multiply it through against the kinetic energy matrix
  Eigen::FullPivLU<Eigen::MatrixXcd> lu(mat_inertia);
  Eigen::MatrixXcd mat_decomp = lu.solve(mat_ke);
  // std::cout << std::setprecision(3) << mat_decomp << "\n\n";

  // eigendecompose
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(mat_decomp);
  vec_eig = es.eigenvalues();
  mat_eig = es.eigenvectors();
  // std::cout << "Size of vec_eig: " << vec_eig.rows() << "\n";

  // std::cout << vec_eig << "\n\n";
  /*
    int qmax = _q.N() - 1;
    for (int idxe = 0; idxe < num_elements - 1; ++idxe) {
      std::cout << "Layer: " << inp_basis.layers()[idxe] << "\n";
      if ((inp_basis.layers()[idxe + 1] - inp_basis.layers()[idxe]) == 1) {
        {
          int laynum = inp_basis.layers()[idxe];
          // double elem_width = x_elem[idxe + 1] - x_elem[idxe];
          double tmpval = mod_new.mu(laynum)(_nodes[laynum][qmax]);
          for (int idxi = 0; idxi < _numbasis; ++idxi) {
            for (int idxj = 0; idxj < _numbasis; ++idxj) {
              std::complex<double> tmp2 = tmpval;
              tmp2 *= val_ef[idxi][laynum][qmax];
              tmp2 *= val_efd[idxj][laynum][qmax];
              mat_ke(idxi, idxj) += tmp2;
            }
          }
        }

        // do minus from next layer as well
        int idxup = idxe + 1;
        {
          int laynum = inp_basis.layers()[idxup];
          // double elem_width = x_elem[idxe + 1] - x_elem[idxe];
          double tmpval = mod_new.mu(laynum)(_nodes[laynum][0]);
          for (int idxi = 0; idxi < _numbasis; ++idxi) {
            for (int idxj = 0; idxj < _numbasis; ++idxj) {
              std::complex<double> tmp2 = -tmpval;
              tmp2 *= val_ef[idxi][laynum][0];
              tmp2 *= val_efd[idxj][laynum][0];
              mat_ke(idxi, idxj) += tmp2;
            }
          }
        }
      }
    }
    */

  // output
  // std::cout << std::setprecision(3) << mat_ke << "\n\n";
  // std::cout << std::setprecision(3) << mat_inertia << "\n\n";

  //   // eigenvalues
  //   Eigen::GeneralizedEigenSolver<Eigen::MatrixXcd> ges;
  //   ges.compute(mat_ke, mat_inertia);
  //   Eigen::VectorXcd tmp_eig = ges.eigenvalues();
  //   Eigen::MatrixXcd tmp_vec = ges.eigenvectors();
  //   vec_eig.resize(tmp_eig.rows());
  //   mat_eig.resize(tmp_vec.rows(), tmp_vec.cols());

  //   // find index permutation for ascending order
  //   std::vector<int> vec_idx(tmp_eig.rows());
  //   for (int i = 0; i < vec_idx.size(); ++i) {
  //     vec_idx[i] = i;
  //   }
  //   auto myfun = [&vec_comp = tmp_eig](int a, int b) {
  //     return std::abs(vec_comp(a)) < std::abs(vec_comp(b));
  //   };
  //   std::sort(vec_idx.begin(), vec_idx.end(), myfun);

  //   // get final result in order of increasing absolute magnitude
  //   for (int i = 0; i < vec_eig.rows(); ++i) {
  //     vec_eig(i) = tmp_eig(vec_idx[i]);
  //   }
  //   for (int i = 0; i < tmp_vec.cols(); ++i) {
  //     mat_eig.block(0, i, tmp_vec.rows(), 1) =
  //         tmp_vec.block(0, vec_idx[i], tmp_vec.rows(), 1);
  //   }
};

/*
      if ((vec_layer[1] - vec_layer[0]) == 1) {
        double tmpval = inp_model.mu(0)(x_nodes[0][qmax]);
        for (int idxq = 1; idxq < _q.N(); ++idxq) {
          double tmp2 = tmpval * pleg.Derivative(idxq, _q.X(qmax));
          int idxy = overallidx(0, idxq);
          int idxx = overallidx(0, qmax);
          // std::cout << "x: " << idxx << " " << idxy << "\n";
          tripletlist3.push_back(T(idxx, idxy, tmp2));
        }

        // do next layer as well
        int idxup = 1;
        {
          int laynum = vec_layer[idxup];
          // double elem_width = x_elem[idxe + 1] - x_elem[idxe];
          double tmpval = inp_model.mu(laynum)(x_nodes[idxup][0]);
          for (int idxq = 0; idxq < _q.N(); ++idxq) {
            double tmp2 = -tmpval * pleg.Derivative(idxq, -1.0);
            int idxy = overallidx(1, idxq);
            int idxx = overallidx(1, 0);
            // std::cout << "x: " << idxx << " " << idxy << "\n";
            tripletlist3.push_back(T(idxx, idxy, tmp2));
          }
        }
      }
      for (int idxe = 1; idxe < _NE - 1; ++idxe) {
        if ((vec_layer[idxe + 1] - vec_layer[idxe]) == 1) {
          {
            int laynum = vec_layer[idxe];
            // double elem_width = x_elem[idxe + 1] - x_elem[idxe];
            double tmpval = inp_model.mu(laynum)(x_nodes[idxe][qmax]);
            for (int idxq = 0; idxq < _q.N(); ++idxq) {
              double tmp2 = tmpval * pleg.Derivative(idxq, 1.0);
              int idxy = overallidx(idxe, idxq);
              int idxx = overallidx(idxe, qmax);

              // std::cout << "x: " << idxx << " " << idxy << "\n";
              tripletlist3.push_back(T(idxx, idxy, tmp2));
            }
          }
          // do next layer as well
          int idxup = idxe + 1;
          {
            int laynum = vec_layer[idxup];
            // double elem_width = x_elem[idxe + 1] - x_elem[idxe];
            double tmpval = inp_model.mu(laynum)(x_nodes[idxup][0]);
            for (int idxq = 0; idxq < _q.N(); ++idxq) {
              double tmp2 = -tmpval * pleg.Derivative(idxq, -1.0);
              int idxy = overallidx(idxup, idxq);
              int idxx = overallidx(idxup, 0);
              // std::cout << "x: " << idxx << " " << idxy << "\n";
              tripletlist3.push_back(T(idxx, idxy, tmp2));
            }
          }
        }
      }
      */