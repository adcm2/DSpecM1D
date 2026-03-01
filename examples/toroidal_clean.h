#ifndef TOROIDAL_CLEAN_CLASS_GUARD_H
#define TOROIDAL_CLEAN_CLASS_GUARD_H
#include <iostream>
#include <cmath>
#include <functional>
#include <GaussQuad/All>
#include <Interpolation/All>
#include <filesystem>
#include <fstream>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/SymGEigsSolver.h>
#include <PlanetaryModel/All>
#include <GSHTrans/All>
#include <Spectra/SymGEigsShiftSolver.h>
#include "SourceInfo.h"

using namespace Spectra;

namespace Toroidal {

// spectral element solver
class spectral_element_planet {
private:
  using stdvec = std::vector<double>;
  using stdvvec = std::vector<stdvec>;
  using stdvvvec = std::vector<stdvvec>;
  int _NE, _l, _k2, _el, _eu, _en, _solint;
  int nummodes, numaug = 0;
  double _mu;
  GaussQuad::Quadrature1D<double> _q;
  Eigen::SparseMatrix<double> mat_seig, mat_ke, mat_inertia;
  std::vector<double> x_elem, x_valnodes;
  stdvvec x_nodes;
  stdvvvec vec_eigval, vec_eigderiv;
  stdvvvec vec_augval, vec_augderiv;
  // stdvvvec vec_fullval, vec_fullderiv;

  std::vector<int> vec_layer;
  Eigen::VectorXcd evalues_seig, evalues_gen;
  Eigen::MatrixXcd evectors, evectors_gen;
  Eigen::MatrixXd _augbasis;
  bool _calc_eig = false, _calc_gen = false;
  bool _augcalc = false;
  bool _has_fluid = false;
  EarthMesh::RadialMesh _mesh;
  std::vector<int> _vec_fluid;
  std::vector<bool> _vec_dof;
  double _freq_norm;
  std::size_t mlen;
  stdvvvec mat_gauss_deriv;
  stdvvec vec_lag_deriv, vec_delta;
  std::size_t overallidx(int idxe, int idxq) const {
    return (_q.N() - 1) * idxe + idxq - 1;
  };

  std::size_t overallidxfinal(int idxe, int idxq) const {
    return (_q.N() - 1) * idxe + idxq;
  };

  std::size_t idx_fullmesh(int idxe, int idxq) const {
    return (_q.N() - 1) * idxe + idxq;
  };

  double densitynorm = 5515.0;
  double pi_db = 3.14159265358979;
  double bigg_db = 6.6723 * std::pow(10.0, -11.0);
  double frequencynorm = std::sqrt(pi_db * bigg_db * densitynorm);
  double normint;

  std::size_t numlen;

  // forced problem
  Eigen::VectorXcd vec_force;

public:
  template <class model1d>
  spectral_element_planet(const model1d &, double, int, int, int = 0);

  void CalculateEigenfrequencies(int, double = 0.0);
  template <class model1d>
  void CalculateForce(const model1d &inp_model, SourceInfo::EarthquakeCMT &,
                      int);

  // value of eigenvector
  auto evector_value(int i) const {
    return evectors.block(0, i, evectors.rows(), 1);
  };

  // matrix of all eigenvectors
  Eigen::MatrixXcd evectors_all() const { return evectors; };

  // eigenvalues and eigenfrequencies
  auto evalues() const { return evalues_seig; };
  auto efrequencies() const;
  auto efrequencies_gen() const;

  auto efrequency(int i) const {
    return std::sqrt(std::abs(evalues_seig(evalues_seig.rows() - 1 - i)));
  };

  auto efunction(int i) const {
    return evectors_gen.block(0, evalues_seig.rows() - 1 - i, evectors.rows(),
                              1);
  };
  auto efunction_all() const { return evectors_gen; };
  auto evectors_std() const;
  stdvvvec &evectors_std_ref();

  template <class model1d> auto traction_std(const model1d &inp_model) const;
  // derivative of eigenvector
  //  outputs in form of vector of vectors, corresponding to element and nodes
  //  within each element in the inner vector
  auto evector_deriv() const;
  stdvvvec &evectors_deriv_ref();

  // number of modes
  auto NumberOfModes() const { return nummodes; };
  auto NumberOfAugment() const { return numaug; };
  // Eigen::VectorXd efunction(int) const;

  // augmentation basis function calculation
  void augment_basis_calculate();
  auto augment_deriv() const;
  auto augment_basis() const { return vec_augval; };
  stdvvvec &aug_std_ref() { return vec_augval; };
  stdvvvec &aug_deriv_ref() { return vec_augderiv; };
  auto is_augmented() const { return _augcalc; };

  // full basis
  // auto full_basis() const { return vec_fullbasis; };
  // auto full_deriv() const { return vec_fullderiv; };
  // Eigen::VectorXd augment_basis(int) const;

  auto xvalues() const { return x_valnodes; };
  auto xnodes() const { return x_nodes; };
  auto xelem() const { return x_elem; };

  auto &layers() const { return vec_layer; };

  auto &mesh() const { return _mesh; };
  auto el() const { return _el; };
  auto eu() const { return _eu; };
  Eigen::SparseMatrix<double> MatEig() const { return mat_seig; };
  Eigen::SparseMatrix<double> MatInertia() const { return mat_inertia; };

  auto EigenCalc() { return _calc_eig; };
  // auto normint() const { return normint; };
  auto q() const { return _q; };

  // k2
  auto k2() const { return _k2; };

  std::size_t idx_submesh(int idxe, int idxq) const {
    if ((_has_fluid && (_solint == 0)) || (!_has_fluid)) {
      return (_q.N() - 1) * (idxe - _el) + idxq - 1;
    } else {
      return (_q.N() - 1) * (idxe - _el) + idxq;
    }
  };
};

template <class model1d>
spectral_element_planet::spectral_element_planet(const model1d &inp_model,
                                                 double maxstep, int NQ, int l,
                                                 int idx_l)
    : _l{l}, _k2{l * (l + 1)}, _freq_norm{1.0 / inp_model.TimeNorm()},
      _solint{idx_l}, normint{1.0 / (inp_model.TimeNorm() * frequencynorm) *
                              sqrt(inp_model.DensityNorm() / densitynorm)} {

  // std::cout << "\n\n\n" << normint << "\n\n\n";
  // double normint = frequencynorm * sqrt(prem.DensityNorm() / densitynorm);
  // get node points for integral and declare radial mesh
  _q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(NQ);
  bool incball = false;
  _mesh = EarthMesh::RadialMesh(inp_model, NQ, maxstep, 1.2, incball);
  _NE = _mesh.NE();
  _vec_fluid = std::vector<int>(_NE, 0);
  _vec_dof = std::vector<bool>(_NE, false);
  for (int idxe = 0; idxe < _NE; ++idxe) {
    if (inp_model.IsFluid(_mesh.LayerNumber(idxe))) {
      _vec_fluid[idxe] = 1;
      _has_fluid = true;
    }

    if (idxe > 0) {
      if ((std::abs(_vec_fluid[idxe] - _vec_fluid[idxe - 1]) == 1)) {
        _vec_dof[idxe - 1] = true;
      }
    }
  }
  // {
  //   int idx1 = 0;
  //   for (auto idx : _vec_dof) {
  //     std::cout << idx1++ << " " << idx << "\n";
  //   }
  // }
  auto pleg =
      Interpolation::LagrangePolynomial(_q.Points().begin(), _q.Points().end());

  // vectors for filling out matrix
  // size of matrix:
  mlen = (_q.N() - 1) * _NE;

  if (_has_fluid) {
    int idx_solid = (_vec_fluid[0] == 0);
    int idx_fluid = (!(_vec_fluid[0] == 0));
    // std::cout << idx_solid << " " << idx_fluid << "\n";
    bool setlow = false;
    // if (_vec_fluid[0] == 0) {
    //   idx_solid = 1;
    // } else {
    //   idx_fluid = 1;
    // }

    for (int idxe = 0; idxe < _NE; ++idxe) {
      if ((!setlow) && ((idx_solid - 1) == idx_l)) {
        // if ((idx_solid - 1) == idx_l) {
        _el = idxe;
        setlow = true;
        // }
      } else if (setlow && (_vec_dof[idxe] || (idxe == (_NE - 1)))) {
        // if (_vec_dof[idxe] || (idxe == (_NE - 1))) {
        // std::cout << "Current index: " << idxe << ", dof: " << _vec_dof[idxe]
        //           << "\n";
        _eu = idxe + 1;
        break;
        // }
      }
      if (_vec_dof[idxe]) {
        idx_solid += _vec_fluid[idxe];
        idx_fluid += (1 - _vec_fluid[idxe]);
        // if (_vec_fluid[idxe] == 1) {
        //   ++idx_solid;
        // } else {
        //   ++idx_fluid;
        // }
      }
    }
  } else {
    _el = 0;
    _eu = _NE;
  }
  // _el = 0;
  // _eu = _NE;
  _en = _eu - _el;
  numlen = _en * (_q.N() - 1) + 1;
  if (_el == 0) {
    numlen -= 1;
  }
  /*
  std::cout << "\n#elements: " << _mesh.NE() << ", _el: " << _el
            << ", _eu: " << _eu << "\n";
  std::cout << "mlen: " << mlen << ", numlen: " << numlen << "\n";
  std::cout << "Minimum radius: " << _mesh.NodeRadius(_el, 0)
            << ", max radius: " << _mesh.NodeRadius(_eu - 1, _q.N() - 1)
            << "\n";
            */
  // std::cout << "\nNE: " << _NE << ", el: " << _el << ", eu: " << _eu <<
  // "\n\n"; std::cout << "\\\\\\\\\\\\\\\\\\\\\\\\\n";
  // {
  //   int idxout = 0;
  //   for (auto idx : _vec_fluid) {
  //     std::cout << idxout++ << " " << idx << "\n";
  //   }
  // }
  // std::cout << "\\\\\\\\\\\\\\\\\\\\\\\\\n";

  // for (auto idx : _vec_dof) {
  //   // if (idx) {
  //   //   ++mlen;
  //   // }
  // }

  // std::cout << "Check 1\n";
  {

    // generating matrix of derivative values for Lagrange polynomials

    mat_gauss_deriv.reserve(_q.N());
    for (int idxk = 0; idxk < _q.N(); ++idxk) {
      stdvvec mat_tmp(_q.N(), std::vector<double>(_q.N()));
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
    // std::cout << "Check 2\n";
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

    // declaring
    mat_seig.resize(numlen, numlen);
    mat_ke.resize(numlen, numlen);
    mat_inertia.resize(numlen, numlen);
    std::vector<double> vec_nn(numlen, 0.0), vec_lm1(numlen, 0.0);
    using T = Eigen::Triplet<double>;
    std::vector<T> tpl_se, tpl_ke, tpl_in;
    tpl_se.reserve(_q.N() * _q.N() * _NE);
    tpl_ke.reserve(_q.N() * _q.N() * _NE);
    tpl_in.reserve(_q.N() * _q.N() * _NE);
    // std::cout << "Check 3\n";
    {

      // new:
      {
        std::size_t idxtl = 0;
        for (int idxe = _el; idxe < _eu; ++idxe) {
          int imin = (idxe == 0);
          double elem_width = _mesh.EW(idxe);
          int laynum = _mesh.LayerNumber(idxe);
          std::size_t idxu = idxtl;
          for (int i = imin; i < _q.N(); ++i) {
            // std::cout << idxu << "\n";
            double xrad = _mesh.NodeRadius(idxe, i);
            double tmp = elem_width / 2.0 * _q.W(i) *
                         inp_model.Density(laynum)(xrad) * xrad * xrad;
            vec_nn[idxu] += tmp;
            tpl_in.push_back(T(idxu, idxu, tmp));
            ++idxu;
          }
          idxtl += _q.N() - imin - 1;
        }
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

      // new:
      {
        std::size_t idxtl = 0;
        for (int idxe = _el; idxe < _eu; ++idxe) {
          int imin = (idxe == 0);
          double elem_width = _mesh.EW(idxe);

          int laynum = _mesh.LayerNumber(idxe);
          double d_val = 2.0 / elem_width;
          std::size_t idxti = idxtl;
          for (int i = imin; i < _q.N(); ++i) {
            std::size_t idxtj = idxtl;

            for (int j = imin; j < _q.N(); ++j) {
              double tmp = 0.0;

              // loop over points
              for (int k = 0; k < _q.N(); ++k) {
                double xrad = _mesh.NodeRadius(idxe, k);   // radius
                double vdi = (i == k);
                double vdj = (j == k);
                double tmp1 = xrad * vec_lag_deriv[k][i] * d_val - vdi;
                tmp1 *= xrad * vec_lag_deriv[k][j] * d_val - vdj;
                tmp1 *= inp_model.L(laynum)(xrad);
                double tmp2 = (_k2 - 2.0) * vdi * vdj;
                tmp2 *= inp_model.N(laynum)(xrad);
                tmp += (tmp1 + tmp2) * _q.W(k);
              }
              tmp *= elem_width / 2.0;

              tpl_se.push_back(T(idxti, idxtj, tmp * vec_lm1[idxti]));
              tpl_ke.push_back(T(idxti, idxtj, tmp));

              ++idxtj;
            }
            ++idxti;
          }
          idxtl += _q.N() - imin - 1;
        }
      }
    }
    mat_seig.setFromTriplets(tpl_se.begin(), tpl_se.end());
    mat_ke.setFromTriplets(tpl_ke.begin(), tpl_ke.end());
    mat_inertia.setFromTriplets(tpl_in.begin(), tpl_in.end());
    // std::cout << "Check before setting inertia matrix\n";
    // for (int k = 0; k < mat_inertia.outerSize(); ++k) {
    //   for (SparseMatrix<double>::InnerIterator it(mat_inertia, k); it; ++it)
    //   {
    //     std::cout << it.row() << " " << it.col() << " "
    //               << std::abs(it.value() - vec_nn[k]) << "\n";
    //     // it.row();   // row index
    //     // it.col();   // col index (here it is equal to k)
    //     // it.index(); // inner index, here it is equal to it.row()
    //   }
    // }
  }
  mat_seig.makeCompressed();
  mat_ke.makeCompressed();
  mat_inertia.makeCompressed();
};

void
spectral_element_planet::CalculateEigenfrequencies(int N, double sigshift) {
  // set number of modes
  nummodes = N;
  // initiate matrix multiplicaiton wrapper
  // std::cout << "Eig Check 1\n";
  SparseGenMatProd<double> op_seig(mat_seig);
  // std::cout << "Eig Check 2\n";
  using OpType = SymShiftInvert<double, Eigen::Sparse, Eigen::Sparse>;
  using BOpType = SparseSymMatProd<double>;
  OpType op_mult(mat_ke, mat_inertia);
  BOpType op_ke(mat_ke), op_inertia(mat_inertia);
  // std::cout << "Eig Check 3\n";

  // std::cout << "Eig Check 4\n";

  // get nc
  long int maxn, maxn2;
  if ((5 * N) > mat_seig.rows()) {
    maxn = mat_seig.rows();
  } else {
    maxn = 5 * N;
  }
  if ((3 * N) > mat_seig.rows()) {
    maxn2 = mat_seig.rows();
  } else {
    maxn2 = 3 * N;
  }

  // initiate eigensolver
  // std::cout << "Eig Check 5\n";
  GenEigsSolver<SparseGenMatProd<double>> eigs(op_seig, N, maxn);
  // std::cout << "Eig Check 6\n";
  // std::cout << "Enter shift:\n";
  // double sigshift;
  // std::cin >> sigshift;
  SymGEigsShiftSolver<OpType, BOpType, GEigsMode::ShiftInvert> eig_gen(
      op_mult, op_inertia, N, maxn2, sigshift);
  eigs.init();
  eig_gen.init();
  int nconv_seig = eigs.compute(SortRule::SmallestMagn);
  int nconv_gen = eig_gen.compute(SortRule::LargestMagn);

  if (eigs.info() == CompInfo::Successful) {
    std::cout << "Successful\n";
    evalues_seig = eigs.eigenvalues();

    // check if we are computing for the IC
    if ((_has_fluid && (_solint == 0)) || (!_has_fluid)) {
      std::size_t nrow = eigs.eigenvectors().rows() + 1;
      std::size_t ncol = eigs.eigenvectors().cols();
      evectors_gen.resize(nrow, ncol);
      evectors_gen.block(0, 0, 1, ncol) = Eigen::MatrixXcd::Zero(1, ncol);
      evectors_gen.block(1, 0, nrow - 1, ncol) = eigs.eigenvectors();
    } else {
      std::size_t nrow = eigs.eigenvectors().rows();
      std::size_t ncol = eigs.eigenvectors().cols();
      evectors_gen.resize(nrow, ncol);
      evectors = eigs.eigenvectors();
    }

    _calc_eig = true;

    // std::cout << "Eigenvalues from non-generalised: \n"
    //           << eigs.eigenvalues() << "\n\n";

    // evectors = eigs.eigenvectors();
  } else {
    std::cout << "Unsuccessful non-generalised\n";
  }

  if (eig_gen.info() == CompInfo::Successful) {
    // std::cout << "Eigenvalues from generalised: \n"
    //           << eig_gen.eigenvalues() * _freq_norm * _freq_norm << "\n\n";

    evalues_gen = eig_gen.eigenvalues();
    Eigen::MatrixXcd normval = eig_gen.eigenvectors().transpose() *
                               mat_inertia * eig_gen.eigenvectors();
    // std::cout << normval << "\n\n";

    // check if we are computing for the IC
    if ((_has_fluid && (_solint == 0)) || (!_has_fluid)) {
      std::size_t nrow = eig_gen.eigenvectors().rows() + 1;
      std::size_t ncol = eig_gen.eigenvectors().cols();
      evectors_gen.resize(nrow, ncol);
      evectors_gen.block(0, 0, 1, ncol) = Eigen::MatrixXcd::Zero(1, ncol);
      evectors_gen.block(1, 0, nrow - 1, ncol) = eig_gen.eigenvectors();

    } else {
      std::size_t nrow = eig_gen.eigenvectors().rows();
      std::size_t ncol = eig_gen.eigenvectors().cols();
      evectors_gen.resize(nrow, ncol);
      evectors_gen = eig_gen.eigenvectors();
    }

    // normalise:
    for (int i = 0; i < evectors_gen.rows(); ++i) {
      for (int j = 0; j < evectors_gen.cols(); ++j) {
        evectors_gen(i, j) *= 1.0 / (normint * sqrt(evalues_gen(j).real()));
      }
    }
    _calc_gen = true;
    // std::cout << "POST\n";
  } else {
    std::cout << "Unsuccessful generalised\n";
  }

  if (_calc_gen) {
    // save the eigenvectors in the "standard format":
    int NQ = _mesh.NN();
    for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
      stdvvec vec_out;
      if ((idxe < _el) || (idxe > _eu - 1)) {
        vec_eigval.push_back(stdvvec(NQ, stdvec(evectors_gen.cols(), 0.0)));
      } else {
        for (int idxq = 0; idxq < NQ; ++idxq) {
          std::vector<double> vec_x;
          std::size_t ovidx = (idxe - _el) * (NQ - 1) + idxq;
          for (int idx = 0; idx < evectors_gen.cols(); ++idx) {
            std::size_t colidx = evectors_gen.cols() - 1 - idx;
            std::complex<double> eigint = evectors_gen(ovidx, colidx);
            vec_x.push_back(eigint.real());
          }
          vec_out.push_back(vec_x);
        }
        vec_eigval.push_back(vec_out);
      }
    }

    // calculate the derivatives of the eigenvectors:
    vec_eigderiv =
        stdvvvec(_mesh.NE(), stdvvec(NQ, stdvec(evectors_gen.cols(), 0.0)));
    for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
      stdvvec vec_out;
      double elem_width = _mesh.EW(idxe);
      for (int idxq = 0; idxq < NQ; ++idxq) {
        for (int idx = 0; idx < evectors_gen.cols(); ++idx) {
          double tmp = 0.0;
          for (int idxq2 = 0; idxq2 < NQ; ++idxq2) {
            // std::size_t ovidx = overallidxfinal(idx, idxq2);
            // tmp += evectors(ovidx, i) * vec_lagderiv[idxq2][idxq];
            tmp += vec_eigval[idxe][idxq2][idx] * vec_lag_deriv[idxq][idxq2];
            // tmp += vec_lag_deriv[idxq][idxq2];
          }
          tmp *= 2.0 / elem_width;
          vec_eigderiv[idxe][idxq][idx] = tmp;
        }
      }
    }
  }

  // // calculate augmented
  // this->augment_basis_calculate();
  // vec_fullbasis = stdvvvec(_mesh.NE(),
  //                          stdvvec(NQ, stdvec(nummodes + numaug, 0.0)));
  // vec_fullderiv = stdvvvec(_mesh.NE(),
  //                          stdvvec(NQ, stdvec(nummodes + numaug, 0.0)));
  // for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
  //   for (int idxq = 0; idxq < NQ; ++idxq) {
  //     for (int idx = 0; idx < evectors_gen.cols(); ++idx) {
  //     }
  //   }
  // }
};

auto
spectral_element_planet::evectors_std() const {
  return vec_eigval;
};
spectral_element_planet::stdvvvec &
spectral_element_planet::evectors_std_ref() {
  return vec_eigval;
};

auto
spectral_element_planet::evector_deriv() const {

  return vec_eigderiv;
};

spectral_element_planet::stdvvvec &
spectral_element_planet::evectors_deriv_ref() {
  return vec_eigderiv;
};

template <class model1d>
auto
spectral_element_planet::traction_std(const model1d &inp_model) const {
  std::size_t nummodes0 = vec_eigderiv.back().back().size();
  int NQ = _mesh.NN();
  stdvvvec vec_traction(_mesh.NE(), stdvvec(NQ, stdvec(nummodes0, 0.0)));
  for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
    // stdvvec vec_out;
    int laynum = _mesh.LayerNumber(idxe);
    double elem_width = _mesh.EW(idxe);
    int qmin = 0;
    if (idxe == 0) {
      qmin = 1;
    }
    for (int idxq = 0; idxq < NQ; ++idxq) {
      // stdvec vec_x(nummodes,0.0);
      double xrad = _mesh.NodeRadius(idxe, idxq);
      for (int idx = 0; idx < nummodes0; ++idx) {
        double tmp =
            inp_model.L(laynum)(xrad) *
            (vec_eigderiv[idxe][idxq][idx] -
             vec_eigval[idxe][idxq][idx] / _mesh.NodeRadius(idxe, idxq));
        vec_traction[idxe][idxq][idx] = tmp;
      };
    };
  };
  return vec_traction;
};

auto
spectral_element_planet::efrequencies() const {
  Eigen::VectorXd ret_eig(evalues_seig.rows());
  for (int i = 0; i < evalues_seig.rows(); ++i) {
    ret_eig(evalues_seig.rows() - 1 - i) =
        std::sqrt(std::abs(evalues_seig(i))) * _freq_norm;
  }
  return ret_eig;
};

auto
spectral_element_planet::efrequencies_gen() const {
  Eigen::VectorXd ret_eig(evalues_gen.rows());
  for (int i = 0; i < evalues_gen.rows(); ++i) {
    ret_eig(evalues_gen.rows() - 1 - i) =
        std::sqrt(std::abs(evalues_gen(i))) * _freq_norm;
  }
  return ret_eig;
};

void
spectral_element_planet::augment_basis_calculate() {
  std::size_t nrows = evectors_gen.rows();
  std::size_t ncols = _mesh.NL() + 1;

  Eigen::SparseLU<Eigen::SparseMatrix<double>> _solver;
  _solver.analyzePattern(mat_ke);
  _solver.factorize(mat_ke);

  // number of layers within layer
  int num_aug = _mesh.LayerNumber(_eu - 1) - _mesh.LayerNumber(_el);
  // std::cout << "\n_el: " << _el << ", _eu:" << _eu - 1
  //           << ", # of layers: " << num_aug << "\n";
  numaug = num_aug;
  _augbasis = Eigen::MatrixXd::Zero(nrows, num_aug);
  // need to account for fact that in solid we only find the result in one solid
  // section
  // whether or not we're looking at IC or purely solid planet
  int lowidx = 0;
  if ((_has_fluid && (_solint == 0)) || (!_has_fluid)) {
    lowidx = 1;
  } else {
    lowidx = 0;
  }

  // force matrix, initialise and set (0,0) element (for first boundary)
  Eigen::MatrixXd mat_force = Eigen::MatrixXd::Zero(mat_ke.cols(), num_aug);
  // mat_force(0, 0) = -1;
  int colnum = 0;
  for (int idxe = _el + 1; idxe < _eu; ++idxe) {
    if ((_mesh.LayerNumber(idxe) - _mesh.LayerNumber(idxe - 1)) == 1) {
      std::size_t ovidx;
      ovidx = idx_submesh(idxe, 0);
      mat_force(ovidx, colnum++) = -1.0;
    }
  }
  // mat_force(mat_ke.cols() - 1, colnum) = -1;

  // calculate basis
  _augbasis.block(lowidx, 0, mat_ke.cols(), num_aug) = _solver.solve(mat_force);
  _augcalc = true;
  // Eigen::MatrixXcd normval =
  //     _augbasis.block(lowidx, 0, mat_ke.cols(), num_aug).transpose() *
  //     mat_inertia * _augbasis.block(lowidx, 0, mat_ke.cols(), num_aug);

  // std::cout << std::setprecision(4) << normval << "\n";
  // fill out augval and augderiv
  //  save the augmented basis in the "standard format":
  int NQ = _mesh.NN();
  vec_augval = stdvvvec(_mesh.NE(), stdvvec(NQ, stdvec(num_aug, 0.0)));
  for (int idxe = _el; idxe < _eu; ++idxe) {
    for (int idxq = 0; idxq < NQ; ++idxq) {
      std::size_t ovidx = (idxe - _el) * (NQ - 1) + idxq;
      for (int idx = 0; idx < num_aug; ++idx) {
        vec_augval[idxe][idxq][idx] = _augbasis(ovidx, idx);
      }
    }
  }
  // calculate the derivatives of the augmented basis:
  vec_augderiv = stdvvvec(_mesh.NE(), stdvvec(NQ, stdvec(num_aug, 0.0)));
  for (int idxe = _el; idxe < _eu; ++idxe) {
    double elem_width = _mesh.EW(idxe);
    for (int idxq = 0; idxq < NQ; ++idxq) {
      for (int idx = 0; idx < num_aug; ++idx) {
        double tmp = 0.0;
        for (int idxq2 = 0; idxq2 < NQ; ++idxq2) {
          tmp += vec_augval[idxe][idxq2][idx] * vec_lag_deriv[idxq][idxq2];
        }
        tmp *= 2.0 / elem_width;
        vec_augderiv[idxe][idxq][idx] = tmp;
      }
    }
  }
};

auto
spectral_element_planet::augment_deriv() const {
  return vec_augderiv;
};

template <class model1d>
void
spectral_element_planet::CalculateForce(const model1d &inp_model,
                                        SourceInfo::EarthquakeCMT &cmt, int m) {

  // resize force
  vec_force.resize(numlen);
  vec_force = Eigen::VectorXcd::Zero(numlen);

  // find element within which the source sits
  double depth = cmt.Depth();
  double rad_source = _mesh.PR() - 1000.0 * depth / inp_model.LengthNorm();
  int idxsource = 0;
  for (int idx = _el; idx < _eu; ++idx) {
    if ((_mesh.ELR(idx) < rad_source) && (_mesh.EUR(idx) > rad_source)) {
      idxsource = idx;
      std::cout << idx << " " << rad_source << " " << _mesh.ELR(idx) << " "
                << _mesh.EUR(idx) << "\n";
    };
  };

  int NQ = _mesh.NN();
  std::vector<double> vec_nodes(NQ, 0.0);
  for (int idx = 0; idx < NQ; ++idx) {
    vec_nodes[idx] = _mesh.NodeRadius(idxsource, idx);
  }
  auto pleg =
      Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());

  // get the y0-, y0+ values etc at the source location
  double theta_s = (90.0 - cmt.Latitude()) * EIGEN_PI / (180.0);
  double phi_s = cmt.Longitude() * EIGEN_PI / (180.0);

  auto wigdmat =
      GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                       GSHTrans::Single, GSHTrans::ColumnMajor>(_l, _l, 2,
                                                                theta_s);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
    return ylm;
  };

  std::complex<double> ymc = std::conj(ylmn(_l, m, -1, phi_s));
  std::complex<double> ypc = std::conj(ylmn(_l, m, 1, phi_s));
  std::complex<double> ymmc = std::conj(ylmn(_l, m, -2, phi_s));
  std::complex<double> yppc = std::conj(ylmn(_l, m, 2, phi_s));

  // need indices for the force, so loop through the nodes in the element
  double invsqrt2 = 1.0 / std::sqrt(2.0);
  std::complex<double> isq2 = std::complex<double>(0.0, invsqrt2);
  double omegal2 = (_l + 2) * (_l - 1) / 2.0;

  // size
  std::cout << "size: " << vec_force.rows() << "\n";
  for (int idxq = 0; idxq < NQ; ++idxq) {
    auto w_val = pleg(idxq, rad_source) / rad_source;
    auto w_prefactor =
        pleg.Derivative(idxq, rad_source) - w_val;   // prefactor for first term
    auto w_valc = w_val * isq2 * omegal2;
    auto w_prec = w_prefactor * isq2;

    // index
    std::size_t ovidx = idx_submesh(idxsource, idxq);
    std::cout << "ovidx: " << ovidx << "\n";
    vec_force(ovidx) = w_prec * (cmt.MC0m() * ymc - cmt.MC0p() * ypc) +
                       w_valc * (cmt.MCmm() * ymmc - cmt.MCpp() * yppc);
  }

  // fill out force vector
};

}   // namespace Toroidal

#endif