#ifndef SEM_TOROIDAL_GUARD_H
#define SEM_TOROIDAL_GUARD_H
#include <iostream>
#include <cmath>
#include <Eigen/Core>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/SymGEigsShiftSolver.h>
#include <GSHTrans/Core>
// #include <GaussQuad/All>
#include <Interpolation/Lagrange>
#include "SourceInfo.h"
#include <EarthMesh/All>
#include "MatrixIndices.h"

using namespace Spectra;

namespace Toroidal {

class TStore {
private:
  // private members
  int l, n;
  double f;

public:
  TStore(int l, int n, double f) : l(l), n(n), f(f) {};
  double GetFreq() const { return f; };
  int GetL() const { return l; };
  int GetN() const { return n; };
};

// spectral element solver
class sem {
private:
  // abbreviations
  using stdvec = std::vector<double>;
  using stdvvec = std::vector<stdvec>;
  using stdvvvec = std::vector<stdvvec>;
  using Complex = std::complex<double>;

  // mesh
  EarthMesh::RadialMesh _mesh;

  // integers for calc
  int _lmax, _il;

  // dealing with fluid regions
  bool _has_fluid = false;
  std::vector<int> _vec_fluid;
  std::vector<bool> _vec_dof;

  // matrix parameters
  std::size_t mlen, totlen;
  int _el, _eu, _en, numlen, _k2, _solint;

  // derivatives etc of basis functions
  stdvvvec mat_gauss_deriv;
  stdvvec vec_lag_deriv, vec_delta;

  // values and derivatives of eigenfunctions
  stdvvvec vec_eigval, vec_eigderiv;
  std::vector<stdvvvec> vec_eigval_l, vec_eigderiv_l;
  std::vector<stdvvvec> vec_augval_l, vec_augderiv_l;
  std::vector<stdvvvec> vec_all_l, vec_allderiv_l;
  std::vector<Eigen::MatrixXd> _augbasis_l;

  // matrices
  Eigen::SparseMatrix<double> mat_seig, mat_ke, mat_inertia;
  std::vector<Eigen::SparseMatrix<double>> mat_l_ke, mat_l_inertia;

  // eigenfunction calculation parameters
  int _num_modes, numaug, _num_modes_l;
  bool _calc_eig = false, _calc_gen = false;
  std::vector<bool> _calc_l;

  // normalisation factors
  double densitynorm = 5515.0;
  double pi_db = 3.14159265358979;
  double bigg_db = 6.6723 * std::pow(10.0, -11.0);
  double frequencynorm = std::sqrt(pi_db * bigg_db * densitynorm);
  double normint;   // normalisation factor for inertia matrix
  double _freq_norm, _length_norm;

  // eigenfunction storage
  Eigen::VectorXcd evalues_seig, evalues_gen;
  Eigen::MatrixXcd evectors, evectors_gen;
  std::vector<Eigen::VectorXcd> evalues_gen_l;
  std::vector<Eigen::MatrixXcd> evectors_gen_l;
  // stdvvvec vec_augval, vec_augderiv;

  // indexing class
  MatrixIndices mat_indices;

  // forced problem
  Eigen::VectorXcd vec_force;

  // store for all modes up to f
  int nmodes = 0;
  std::vector<TStore> vec_store;
  std::vector<std::size_t> vec_indices;

public:
  sem() {};
  template <class model1d> sem(const model1d &, double, int, int, int = 0);
  void CalculateEigenfrequencies(int, double);
  void CalculateEigenfrequenciesSeparate(int, double);
  auto efrequencies_gen() const;
  auto efrequencies_gen_separate(int) const;
  auto efunctions(int) const;
  auto efunctions_ref(int idxl) { return evectors_gen_l[idxl - 1]; };
  auto mesh() const { return _mesh; };
  auto el() const { return _el; };
  auto eu() const { return _eu; };

  auto modes_coupled() const { return vec_store; };

  Eigen::VectorXcd CalculateForce(SourceInfo::EarthquakeCMT &);
  Eigen::VectorXcd ReceiverVectorTheta(double, double, double);
  Eigen::VectorXcd ReceiverVectorThetaSurface(double, double, double);
  Eigen::VectorXcd ReceiverVectorPhi(double, double, double);
  Eigen::VectorXcd ReceiverVectorPhiSurface(double, double, double);
  Eigen::SparseMatrix<double> GetStiffnessMatrix() const { return mat_ke; };
  Eigen::SparseMatrix<double> GetInertiaMatrix() const { return mat_inertia; };

  Eigen::SparseMatrix<double> GetStiffnessMatrixL(int idxl) const {
    return mat_l_ke[idxl - 1];
  };
  Eigen::SparseMatrix<double> GetInertiaMatrixL(int idxl) const {
    return mat_l_inertia[idxl - 1];
  };
  Eigen::MatrixXcd CalculateForce(SourceInfo::EarthquakeCMT &, int);
  // Eigen::MatrixXcd ReceiverVectorTheta(double, double, double);
  Eigen::MatrixXcd ReceiverVectorThetaSurfaceL(double, double, double, int);
  // Eigen::VectorXcd ReceiverVectorPhi(double, double, double);
  Eigen::MatrixXcd ReceiverVectorPhiSurfaceL(double, double, double, int);

  auto PrintModesUpToFreq(double) const;

  // function that gets the matrices for modes up to a certain frequency
  auto FindModesForCoupling(double);
  template <class model1d>
  Eigen::MatrixXcd NMC_INERTIA(const model1d &, bool = false);
  template <class model1d>
  Eigen::MatrixXcd NMC_KE(const model1d &, bool = false);
  auto NMC_FORCE(SourceInfo::EarthquakeCMT &, bool = false);
  Eigen::VectorXcd ReceiverVectorThetaSurfaceCoupling(double, double, double,
                                                      bool = false);
  Eigen::VectorXcd ReceiverVectorPhiSurfaceCoupling(double, double, double,
                                                    bool = false);

  // augmentation
  void augment_basis_calculate();

  // nmc for particular l
  template <class model1d>
  Eigen::MatrixXcd NMC_INERTIA(const model1d &, int, bool = false);
  template <class model1d>
  Eigen::MatrixXcd NMC_KE(const model1d &, int, bool = false);
  Eigen::MatrixXcd ReceiverVectorThetaSurface_NMCL(double, double, double, int,
                                                   bool = false);
  Eigen::MatrixXcd ReceiverVectorPhiSurface_NMCL(double, double, double, int,
                                                 bool = false);
  Eigen::MatrixXcd CalculateForceNMC(SourceInfo::EarthquakeCMT &, int,
                                     bool = false);
  auto NumModesL(int idxl) const { return _num_modes_l; };
};

template <class model1d>
sem::sem(const model1d &inp_model, double maxstep, int NQ, int lmax, int idx_l)
    : _mesh(inp_model, NQ, 1.0, maxstep, false),
      _freq_norm{1.0 / inp_model.TimeNorm()}, _lmax{lmax}, _il{idx_l},
      _k2{lmax * (lmax + 1)},
      normint{1.0 / (inp_model.TimeNorm() * frequencynorm) *
              sqrt(inp_model.DensityNorm() / densitynorm)},
      _length_norm{inp_model.LengthNorm()} {
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  // vector to store whether an element is fluid or not
  _vec_fluid = std::vector<int>(_mesh.NE(), 0);
  _vec_dof = std::vector<bool>(_mesh.NE(), false);
  for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
    if (inp_model.IsFluid(_mesh.LayerNumber(idxe))) {
      _vec_fluid[idxe] = 1;
      _has_fluid = true;
    }

    // if there is a change from fluid to solid or vice versa, we set the dof to
    // true at the lower element
    if (idxe > 0) {
      if ((std::abs(_vec_fluid[idxe] - _vec_fluid[idxe - 1]) == 1)) {
        _vec_dof[idxe - 1] = true;
      }
    }
  }

  // determining if there is a fluid layer and then finding the indices of the
  // elements that correspond to the bottom and top of the appropriate solid
  // layer (as given by _il)
  if (_has_fluid) {
    int idx_solid = (_vec_fluid[0] == 0);
    int idx_fluid = (!(_vec_fluid[0] == 0));
    bool setlow = false;

    // looping through the elements of the mesh
    //  to find the indices of the solid layer
    for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {

      if ((!setlow) && ((idx_solid - 1) == idx_l)) {
        // if we haven't set the lower index yet and the current solid layer is
        // the same as idx_l
        _el = idxe;
        setlow = true;
      } else if (setlow && (_vec_dof[idxe] || (idxe == (_mesh.NE() - 1)))) {
        // if we have set the lower index and either the current element is a
        // change from fluid to solid or we are at the last element
        _eu = idxe + 1;
        break;
      }
      if (_vec_dof[idxe]) {
        idx_solid += _vec_fluid[idxe];
        idx_fluid += (1 - _vec_fluid[idxe]);
      }
    }
  } else {
    _el = 0;
    _eu = _mesh.NE();
  }
  _en = _eu - _el;
  numlen = _en * (_mesh.NN() - 1) + 1;
  if (_el == 0) {
    numlen -= 1;
  }

  // numaug
  if (_el == 0) {
    numaug = _mesh.LayerNumber(_eu - 1) - _mesh.LayerNumber(_el) + 1;
  } else {
    numaug = _mesh.LayerNumber(_eu - 1) - _mesh.LayerNumber(_el) + 2;
  }

  // matrix indices class
  mat_indices = MatrixIndices(_el, _eu, 1, _lmax, _mesh.NN());
  // toroidal modes have lmin = 1

  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////

  // getting matrices of derivative values for Lagrange polynomials
  auto q = _mesh.GLL();
  auto pleg =
      Interpolation::LagrangePolynomial(q.Points().begin(), q.Points().end());
  {
    // storing h_i'(x_k) * h_j'(x_k)
    mat_gauss_deriv.reserve(q.N());
    for (int idxk = 0; idxk < q.N(); ++idxk) {
      stdvvec mat_tmp(q.N(), std::vector<double>(q.N()));
      for (int idxi = 0; idxi < q.N(); ++idxi) {
        std::vector<double> vec_tmp(q.N());
        for (int idxj = 0; idxj < q.N(); ++idxj) {
          vec_tmp[idxj] = pleg.Derivative(idxi, q.X(idxk)) *
                          pleg.Derivative(idxj, q.X(idxk));
        }
        mat_tmp[idxi] = vec_tmp;
      }
      mat_gauss_deriv.push_back(mat_tmp);
    };

    // storing h_i'(x_k)
    for (int idxk = 0; idxk < q.N(); ++idxk) {
      std::vector<double> vec_tmp(q.N(), 0.0), vec_tmp1(q.N(), 0.0);
      for (int idxi = 0; idxi < q.N(); ++idxi) {
        vec_tmp[idxi] = pleg.Derivative(idxi, q.X(idxk));
        vec_tmp1[idxi] = pleg(idxi, q.X(idxk));
      }
      vec_lag_deriv.push_back(vec_tmp);
      vec_delta.push_back(vec_tmp1);
    };
  }

  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  {

    // size of matrix for full problem
    totlen = mat_indices.mat_size();
    std::cout << "Total size of matrix: " << totlen << "\n";
    std::cout << "Number of elements in the solid layer: " << _en << "\n";
    std::cout << "Number of nodes per element: " << _mesh.NN() << "\n";
    std::cout << "lmax: " << _lmax << "\n";

    // matrices
    mat_seig.resize(totlen, totlen);
    mat_ke.resize(totlen, totlen);
    mat_inertia.resize(totlen, totlen);

    // sizes:
    std::cout << "Size of seig: " << mat_seig.rows() << " x " << mat_seig.cols()
              << "\n";
    std::cout << "Size of ke: " << totlen << "\n";
    std::cout << "Return value: "
              << mat_indices.mat_index(_eu - 1, q.N() - 1, _lmax, _lmax)
              << "\n";

    // diagonal elements for inertia matrix
    std::vector<double> vec_nn(numlen, 0.0), vec_lm1(numlen, 0.0);

    // triplet list for matrices
    using T = Eigen::Triplet<double>;
    std::vector<T> tpl_se, tpl_ke, tpl_in;
    tpl_se.reserve(mat_indices.mat_nonzero());
    tpl_ke.reserve(mat_indices.mat_nonzero());
    tpl_in.reserve(mat_indices.mat_nonzero());
    {

      // new:
      {
        std::size_t idxtl = 0;
        for (int idxe = _el; idxe < _eu; ++idxe) {
          int imin = (idxe == 0);
          double elem_width = _mesh.EW(idxe);
          int laynum = _mesh.LayerNumber(idxe);
          std::size_t idxu = idxtl;
          for (int i = imin; i < q.N(); ++i) {
            // std::cout << idxu << "\n";
            double xrad = _mesh.NodeRadius(idxe, i);
            double tmp = elem_width / 2.0 * q.W(i) *
                         inp_model.Density(laynum)(xrad) * xrad * xrad;
            vec_nn[idxu] += tmp;

            // through lm
            for (int idxl = 1; idxl < _lmax + 1; ++idxl) {
              for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
                tpl_in.push_back(T(mat_indices.mat_index(idxe, i, idxl, idxm),
                                   mat_indices.mat_index(idxe, i, idxl, idxm),
                                   tmp));
              }
            }
            ++idxu;
          }
          idxtl += q.N() - imin - 1;
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
          for (int i = imin; i < q.N(); ++i) {
            // std::size_t idxtj = idxtl;

            for (int j = imin; j < q.N(); ++j) {

              for (int idxl = 1; idxl < _lmax + 1; ++idxl) {
                double kval = (idxl + 2) * (idxl - 1);
                double tmp = 0.0;
                // loop over points
                for (int k = 0; k < q.N(); ++k) {
                  double xrad = _mesh.NodeRadius(idxe, k);   // radius
                  double vdi = (i == k);
                  double vdj = (j == k);
                  double tmp1 = xrad * vec_lag_deriv[k][i] * d_val - vdi;
                  tmp1 *= xrad * vec_lag_deriv[k][j] * d_val - vdj;
                  tmp1 *= inp_model.L(laynum)(xrad);
                  double tmp2 = kval * vdi * vdj;
                  tmp2 *= inp_model.N(laynum)(xrad);
                  tmp += (tmp1 + tmp2) * q.W(k);
                }
                tmp *= elem_width / 2.0;

                for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
                  // get the index of the spherical harmonic
                  std::size_t idxtim =
                      mat_indices.mat_index(idxe, i, idxl, idxm);
                  std::size_t idxtjm =
                      mat_indices.mat_index(idxe, j, idxl, idxm);

                  tpl_ke.push_back(T(idxtim, idxtjm, tmp));
                  tpl_se.push_back(T(idxtim, idxtjm, tmp * vec_lm1[idxti]));
                }
              }

              // ++idxtj;
            }
            ++idxti;
          }
          idxtl += q.N() - imin - 1;
        }
      }
    }
    mat_seig.setFromTriplets(tpl_se.begin(), tpl_se.end());
    mat_ke.setFromTriplets(tpl_ke.begin(), tpl_ke.end());
    mat_inertia.setFromTriplets(tpl_in.begin(), tpl_in.end());
  }
  mat_seig.makeCompressed();
  mat_ke.makeCompressed();
  mat_inertia.makeCompressed();

  // make the vector of matrices for the different ls
  {
    using T = Eigen::Triplet<double>;
    std::size_t matlen = mat_indices.mat_size_l();
    mat_l_ke.reserve(_lmax);
    mat_l_inertia.reserve(_lmax);
    std::cout << "Number of matrices for ls: " << mat_l_ke.size() << "\n";
    for (int idxl = 1; idxl < _lmax + 1; ++idxl) {

      std::vector<T> tpl_ls, tpl_in, tpl_ke;
      tpl_ls.reserve(mat_indices.mat_nonzero());
      tpl_in.reserve(mat_indices.mat_nonzero_l());
      tpl_ke.reserve(mat_indices.mat_nonzero_l());
      {
        for (int idxe = _el; idxe < _eu; ++idxe) {
          int imin = (idxe == 0);
          double elem_width = _mesh.EW(idxe);
          int laynum = _mesh.LayerNumber(idxe);
          for (int i = imin; i < q.N(); ++i) {
            double xrad = _mesh.NodeRadius(idxe, i);
            double tmp = elem_width / 2.0 * q.W(i) *
                         inp_model.Density(laynum)(xrad) * xrad * xrad;
            // through lm
            tpl_in.push_back(T(mat_indices.mat_index_l(idxe, i),
                               mat_indices.mat_index_l(idxe, i), tmp));
          }
        }
        Eigen::SparseMatrix<double> mat_tmp(mat_indices.mat_size_l(),
                                            mat_indices.mat_size_l());
        mat_tmp.setFromTriplets(tpl_in.begin(), tpl_in.end());
        mat_tmp.makeCompressed();
        mat_l_inertia.push_back(mat_tmp);
        mat_l_inertia.back().makeCompressed();
      }

      {

        // new:
        {
          for (int idxe = _el; idxe < _eu; ++idxe) {
            int imin = (idxe == 0);
            double elem_width = _mesh.EW(idxe);

            int laynum = _mesh.LayerNumber(idxe);
            double d_val = 2.0 / elem_width;
            for (int i = imin; i < q.N(); ++i) {

              for (int j = imin; j < q.N(); ++j) {

                double kval = (idxl + 2) * (idxl - 1);
                double tmp = 0.0;
                // loop over points
                for (int k = 0; k < q.N(); ++k) {
                  double xrad = _mesh.NodeRadius(idxe, k);   // radius
                  double vdi = (i == k);
                  double vdj = (j == k);
                  double tmp1 = xrad * vec_lag_deriv[k][i] * d_val - vdi;
                  tmp1 *= xrad * vec_lag_deriv[k][j] * d_val - vdj;
                  tmp1 *= inp_model.L(laynum)(xrad);
                  double tmp2 = kval * vdi * vdj;
                  tmp2 *= inp_model.N(laynum)(xrad);
                  tmp += (tmp1 + tmp2) * q.W(k);
                }
                tmp *= elem_width / 2.0;

                // get the index of the spherical harmonic
                std::size_t idxtim = mat_indices.mat_index_l(idxe, i);
                std::size_t idxtjm = mat_indices.mat_index_l(idxe, j);

                tpl_ke.push_back(T(idxtim, idxtjm, tmp));
              }
            }
          }
        }
        {
          Eigen::SparseMatrix<double> mat_tmp(mat_indices.mat_size_l(),
                                              mat_indices.mat_size_l());
          mat_tmp.setFromTriplets(tpl_ke.begin(), tpl_ke.end());
          mat_tmp.makeCompressed();
          mat_l_ke.push_back(mat_tmp);
          mat_l_ke.back().makeCompressed();
        }
      }
    }
  }
};

void
sem::CalculateEigenfrequencies(int N, double sigshift) {
  // set number of modes
  _num_modes = N;
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
  // GenEigsSolver<SparseGenMatProd<double>> eigs(op_seig, N, maxn);
  // std::cout << "Eig Check 6\n";
  // std::cout << "Enter shift:\n";
  // double sigshift;
  // std::cin >> sigshift;
  SymGEigsShiftSolver<OpType, BOpType, GEigsMode::ShiftInvert> eig_gen(
      op_mult, op_inertia, N, maxn2, sigshift);
  // eigs.init();
  eig_gen.init();
  // int nconv_seig = eigs.compute(SortRule::SmallestMagn);
  int nconv_gen = eig_gen.compute(SortRule::LargestMagn);

  /*
  if (eigs.info() == CompInfo::Successful) {
    std::cout << "Successful\n";
    evalues_seig = eigs.eigenvalues();

    // check if we are computing for the IC
    if ((_has_fluid && (_il == 0)) || (!_has_fluid)) {
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
  */

  if (eig_gen.info() == CompInfo::Successful) {
    // std::cout << "Eigenvalues from generalised: \n"
    //           << eig_gen.eigenvalues() * _freq_norm * _freq_norm <<
    // "\n\n";

    evalues_gen = eig_gen.eigenvalues();
    Eigen::MatrixXcd normval = eig_gen.eigenvectors().transpose() *
                               mat_inertia * eig_gen.eigenvectors();
    // std::cout << normval << "\n\n";

    // check if we are computing for the IC
    if ((_has_fluid && (_il == 0)) || (!_has_fluid)) {
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
  // vec_fullbasis = stdvvvec(_mesh.NumberOfElements(),
  //                          stdvvec(NQ, stdvec(nummodes + numaug,
  //  0.0)));
  // vec_fullderiv = stdvvvec(_mesh.NumberOfElements(),
  //                          stdvvec(NQ, stdvec(nummodes + numaug,
  //  0.0)));
  // for (int idxe = 0; idxe < _mesh.NumberOfElements(); ++idxe) {
  //   for (int idxq = 0; idxq < NQ; ++idxq) {
  //     for (int idx = 0; idx < evectors_gen.cols(); ++idx) {
  //     }
  //   }
  // }
};

void
sem::CalculateEigenfrequenciesSeparate(int N, double sigshift) {
  // set number of modes
  _num_modes_l = N;
  evalues_gen_l.resize(_lmax);
  evectors_gen_l.resize(_lmax);
  vec_eigval_l.resize(_lmax);
  vec_eigderiv_l.resize(_lmax);
  vec_all_l = std::vector<stdvvvec>(
      _lmax,
      stdvvvec(_mesh.NE(), stdvvec(_mesh.NN(), stdvec(N + numaug, 0.0))));
  vec_allderiv_l = std::vector<stdvvvec>(
      _lmax,
      stdvvvec(_mesh.NE(), stdvvec(_mesh.NN(), stdvec(N + numaug, 0.0))));
  // initiate matrix multiplicaiton wrapper
  // std::cout << "Eig Check 1\n";
  // SparseGenMatProd<double> op_seig(mat_seig);
  // std::cout << "Eig Check 2\n";
  using OpType = SymShiftInvert<double, Eigen::Sparse, Eigen::Sparse>;
  using BOpType = SparseSymMatProd<double>;
  for (int idxl = 0; idxl < _lmax; ++idxl) {
    OpType op_mult(mat_l_ke[idxl], mat_l_inertia[idxl]);
    BOpType op_ke(mat_l_ke[idxl]), op_inertia(mat_l_inertia[idxl]);
    // std::cout << "Eig Check 3\n";

    // std::cout << "Eig Check 4\n";

    // get nc
    long int maxn, maxn2;
    if ((5 * N) > mat_l_ke[idxl].rows()) {
      maxn = mat_l_ke[idxl].rows();
    } else {
      maxn = 5 * N;
    }
    if ((3 * N) > mat_l_ke[idxl].rows()) {
      maxn2 = mat_l_ke[idxl].rows();
    } else {
      maxn2 = 3 * N;
    }

    // initiate eigensolver
    SymGEigsShiftSolver<OpType, BOpType, GEigsMode::ShiftInvert> eig_gen(
        op_mult, op_inertia, N, maxn2, sigshift);
    eig_gen.init();
    int nconv_gen = eig_gen.compute(SortRule::LargestMagn);

    if (eig_gen.info() == CompInfo::Successful) {
      // std::cout << "Eigenvalues from generalised: \n"
      //           << eig_gen.eigenvalues() * _freq_norm * _freq_norm <<
      // "\n\n";

      evalues_gen_l[idxl] = eig_gen.eigenvalues();
      // std::cout << "\n\n" << evalues_gen_l[idxl] << "\n\n";
      Eigen::MatrixXcd normval = eig_gen.eigenvectors().transpose() *
                                 mat_l_inertia[idxl] * eig_gen.eigenvectors();
      // std::cout << normval << "\n\n";

      // check if we are computing for the IC
      if ((_has_fluid && (_il == 0)) || (!_has_fluid)) {
        std::size_t nrow = eig_gen.eigenvectors().rows() + 1;
        std::size_t ncol = eig_gen.eigenvectors().cols();
        evectors_gen_l[idxl].resize(nrow, ncol);
        evectors_gen_l[idxl].block(0, 0, 1, ncol) =
            Eigen::MatrixXcd::Zero(1, ncol);
        evectors_gen_l[idxl].block(1, 0, nrow - 1, ncol) =
            eig_gen.eigenvectors();

      } else {
        std::size_t nrow = eig_gen.eigenvectors().rows();
        std::size_t ncol = eig_gen.eigenvectors().cols();
        evectors_gen_l[idxl].resize(nrow, ncol);
        evectors_gen_l[idxl] = eig_gen.eigenvectors();
      }

      // normalise:
      for (int i = 0; i < evectors_gen_l[idxl].rows(); ++i) {
        for (int j = 0; j < evectors_gen_l[idxl].cols(); ++j) {
          evectors_gen_l[idxl](i, j) *=
              1.0 / (normint * sqrt(evalues_gen_l[idxl](j).real()));
        }
      }
      _calc_gen = true;
      // std::cout << "POST\n";
    } else {
      std::cout << "Unsuccessful generalised\n";
    }

    if (_calc_gen) {
      // std::cout << "Saving eigenvectors for l = " << idxl + 1 << "\n";
      // save the eigenvectors in the "standard format":
      int NQ = _mesh.NN();
      for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
        stdvvec vec_out;
        if ((idxe < _el) || (idxe > _eu - 1)) {
          vec_eigval_l[idxl].push_back(
              stdvvec(NQ, stdvec(evectors_gen_l[idxl].cols(), 0.0)));
        } else {
          for (int idxq = 0; idxq < NQ; ++idxq) {
            std::vector<double> vec_x;
            std::size_t ovidx = (idxe - _el) * (NQ - 1) + idxq;
            for (int idx = 0; idx < evectors_gen_l[idxl].cols(); ++idx) {
              std::size_t colidx = evectors_gen_l[idxl].cols() - 1 - idx;
              std::complex<double> eigint = evectors_gen_l[idxl](ovidx, colidx);
              vec_x.push_back(eigint.real());
              vec_all_l[idxl][idxe][idxq][idx] = eigint.real();
            }
            vec_out.push_back(vec_x);
          }
          vec_eigval_l[idxl].push_back(vec_out);
        }
      }
      // std::cout << "Check 2\n";
      // calculate the derivatives of the eigenvectors:
      vec_eigderiv_l[idxl] = stdvvvec(
          _mesh.NE(), stdvvec(NQ, stdvec(evectors_gen_l[idxl].cols(), 0.0)));
      for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
        stdvvec vec_out;
        double elem_width = _mesh.EW(idxe);
        for (int idxq = 0; idxq < NQ; ++idxq) {
          for (int idx = 0; idx < evectors_gen_l[idxl].cols(); ++idx) {
            double tmp = 0.0;
            for (int idxq2 = 0; idxq2 < NQ; ++idxq2) {
              // std::size_t ovidx = overallidxfinal(idx, idxq2);
              // tmp += evectors(ovidx, i) * vec_lagderiv[idxq2][idxq];
              tmp += vec_eigval_l[idxl][idxe][idxq2][idx] *
                     vec_lag_deriv[idxq][idxq2];
              // tmp += vec_lag_deriv[idxq][idxq2];
            }
            tmp *= 2.0 / elem_width;
            vec_eigderiv_l[idxl][idxe][idxq][idx] = tmp;
            vec_allderiv_l[idxl][idxe][idxq][idx] = tmp;
          }
        }
      }
    }
  }

  // add to all
  // vec_all_l = vec_eigval_l;
  // vec_allderiv_l = vec_eigderiv_l;
};

auto
sem::efrequencies_gen() const {
  Eigen::VectorXd ret_eig(evalues_gen.rows());
  for (int i = 0; i < evalues_gen.rows(); ++i) {
    ret_eig(evalues_gen.rows() - 1 - i) =
        std::sqrt(std::abs(evalues_gen(i))) * _freq_norm;
  }
  return ret_eig;
};

auto
sem::efrequencies_gen_separate(int idxl) const {
  Eigen::VectorXd ret_eig(evalues_gen_l[idxl - 1].rows());
  for (int i = 0; i < evalues_gen_l[idxl - 1].rows(); ++i) {
    ret_eig(evalues_gen_l[idxl - 1].rows() - 1 - i) =
        std::sqrt(std::abs(evalues_gen_l[idxl - 1](i))) * _freq_norm;
  }
  return ret_eig;
};

auto
sem::efunctions(int idxl) const {
  return vec_eigval_l[idxl - 1];
};

Eigen::VectorXcd
sem::CalculateForce(SourceInfo::EarthquakeCMT &cmt) {

  // resize force
  vec_force.resize(totlen);
  vec_force = Eigen::VectorXcd::Zero(totlen);

  // find element within which the source sits
  double depth = cmt.Depth();
  double rad_source = _mesh.PR() - 1000.0 * depth / _length_norm;
  // std::cout << "Source radius: " << rad_source << "\n";
  int idxsource = 0;

  int NQ = _mesh.NN();

  // get the y0-, y0+ values etc at the source location
  double theta_s = (90.0 - cmt.Latitude()) * EIGEN_PI / (180.0);
  double phi_s = cmt.Longitude() * EIGEN_PI / (180.0);

  // wigner d matrix
  auto wigdmat =
      GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                       GSHTrans::Single, GSHTrans::ColumnMajor>(_lmax, _lmax, 2,
                                                                theta_s);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
    return ylm;
  };

  // std::cout << "Source location: " << theta_s << " " << phi_s << "\n";
  // int m = 0, _l = 2;
  // parameters
  double invsqrt2 = 1.0 / std::sqrt(2.0);
  std::complex<double> isq2 = std::complex<double>(0.0, invsqrt2);

  // loop through the elements
  //   to find the element that contains the source
  for (int idx = _el; idx < _eu; ++idx) {
    // std::cout << _mesh.ELR(idx) << " " << _mesh.EUR(idx) << " " <<
    // rad_source
    //           << "\n";
    if ((_mesh.ELR(idx) <= rad_source) && (_mesh.EUR(idx) >= rad_source)) {
      // std::cout << idx << " " << rad_source << " " << _mesh.ELR(idx) << " "
      //           << _mesh.EUR(idx) << "\n";
      std::vector<double> vec_nodes(NQ, 0.0);
      for (int idxn = 0; idxn < NQ; ++idxn) {
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      }
      // std::cout << "Check 1\n";
      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());
      // std::cout << "Check 2\n";
      for (int idxq = 0; idxq < NQ; ++idxq) {
        auto w_val = pleg(idxq, rad_source) / rad_source;
        auto w_prefactor = pleg.Derivative(idxq, rad_source) -
                           w_val;   // prefactor for first term
        for (int idxl = 1; idxl < _lmax + 1; ++idxl) {
          double omegal2 = (idxl + 2) * (idxl - 1) / 2.0;
          double lprefac =
              std::exp(-2.0 * 3.141592653589793 * (idxl + 1) / (1 + 0.5));
          // lprefac = 1.0;
          for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
            // std::cout << "Check 3\n";
            // spherical harmonic
            Complex ymc = std::conj(ylmn(idxl, idxm, -1, phi_s));
            // std::cout << "Check 4\n";
            Complex ypc = std::conj(ylmn(idxl, idxm, 1, phi_s));
            // std::cout << "Check 5\n";
            Complex ymmc = 0.0, yppc = 0.0;
            if (idxl > 1) {
              ymmc = std::conj(ylmn(idxl, idxm, -2, phi_s));
              // std::cout << "Check 6\n";
              yppc = std::conj(ylmn(idxl, idxm, 2, phi_s));
              //  std::cout << "Check 7\n";
            }

            // get the index of the spherical harmonic
            std::size_t ovidx = mat_indices.mat_index(idx, idxq, idxl, idxm);

            // std::cout << "Size: " << vec_force.rows() << ", ovidx: " <<
            // ovidx
            //           << "\n";
            Complex tmp = w_prefactor * (cmt.MC0m() * ymc - cmt.MC0p() * ypc);
            tmp += w_val * omegal2 * (cmt.MCmm() * ymmc - cmt.MCpp() * yppc);
            tmp *= isq2;

            vec_force(ovidx) = tmp * lprefac;
          };
        };
      };
    };
  };
  return vec_force;
};

Eigen::MatrixXcd
sem::CalculateForce(SourceInfo::EarthquakeCMT &cmt, int idxl) {

  std::size_t flen = mat_indices.mat_size_l();
  std::size_t fcols = 2 * idxl + 1;
  Eigen::MatrixXcd vec_lforce = Eigen::MatrixXcd::Zero(flen, fcols);

  // resize force
  // vec_force.resize(totlen);
  // vec_force = Eigen::VectorXcd::Zero(totlen);
  // auto ylmn = [](int l, int m, int N, double theta, double phi) {
  //   auto wigtemp = GSHTrans::Wigner(l, l, N, theta);   // Wigner D-matrix
  //   auto dl = wigtemp(l);
  //   auto tmp = dl(m);
  //   auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
  //   return ylm;
  // };

  // find element within which the source sits
  double depth = cmt.Depth();
  double rad_source = _mesh.PR() - 1000.0 * depth / _length_norm;
  // std::cout << "\nSource radius: " << rad_source << "\n";
  int idxsource = 0;

  int NQ = _mesh.NN();

  // get the y0-, y0+ values etc at the source location
  double theta_s = (90.0 - cmt.Latitude()) * EIGEN_PI / (180.0);
  double phi_s = cmt.Longitude() * EIGEN_PI / (180.0);

  auto wigdmat =
      GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                       GSHTrans::Single, GSHTrans::ColumnMajor>(_lmax, _lmax, 2,
                                                                theta_s);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
    return ylm;
  };
  // std::cout << "\nSource location: " << theta_s << " " << phi_s << "\n";
  // int m = 0, _l = 2;
  // parameters
  double invsqrt2 = 1.0 / std::sqrt(2.0);
  std::complex<double> isq2 = std::complex<double>(0.0, invsqrt2);
  double omegal2 = std::sqrt((idxl + 2) * (idxl - 1) / 2.0);
  double lprefac = std::exp(-2.0 * 3.141592653589793 * (idxl + 1) / (1 + 0.5));
  // lprefac = 1.0;
  // loop through the elements
  //   to find the element that contains the source
  for (int idx = _el; idx < _eu; ++idx) {
    // std::cout << _mesh.ELR(idx) << " " << _mesh.EUR(idx) << " " <<
    // rad_source
    //           << "\n";
    if ((_mesh.ELR(idx) <= rad_source) && (_mesh.EUR(idx) >= rad_source)) {
      // std::cout << idx << " " << rad_source << " " << _mesh.ELR(idx) << " "
      //           << _mesh.EUR(idx) << "\n";
      std::vector<double> vec_nodes(NQ, 0.0);
      for (int idxn = 0; idxn < NQ; ++idxn) {
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      }

      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());

      for (int idxq = 0; idxq < NQ; ++idxq) {
        auto w_val = pleg(idxq, rad_source) / rad_source;
        auto w_prefactor = pleg.Derivative(idxq, rad_source) -
                           w_val;   // prefactor for first term
        // for (int idxl = 1; idxl < _lmax + 1; ++idxl) {

        // get the index of the spherical harmonic
        std::size_t ridx = mat_indices.mat_index_l(idx, idxq);
        for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
          // spherical harmonic
          Complex ymc = std::conj(ylmn(idxl, idxm, -1, phi_s));
          Complex ypc = std::conj(ylmn(idxl, idxm, 1, phi_s));
          Complex ymmc = 0.0, yppc = 0.0;

          if (idxl > 1) {
            ymmc = std::conj(ylmn(idxl, idxm, -2, phi_s));
            yppc = std::conj(ylmn(idxl, idxm, 2, phi_s));
          }

          // std::cout << "Size: " << vec_force.rows() << ", ovidx: " << ovidx
          //           << "\n";
          Complex tmp = w_prefactor * (cmt.MC0m() * ymc - cmt.MC0p() * ypc);
          tmp += w_val * omegal2 * (cmt.MCmm() * ymmc - cmt.MCpp() * yppc);
          tmp *= isq2;

          vec_lforce(ridx, idxm + idxl) = tmp * lprefac;
        };
        // };
      };
    };
  };
  return vec_lforce;
};

Eigen::VectorXcd
sem::ReceiverVectorTheta(double rad_r, double theta_r, double phi_r) {
  // create the receiver vector
  Eigen::VectorXcd vec_receiver(totlen);
  vec_receiver = Eigen::VectorXcd::Zero(totlen);
  auto wigdmat =
      GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                       GSHTrans::Single, GSHTrans::ColumnMajor>(_lmax, _lmax, 2,
                                                                theta_r);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
    return ylm;
  };
  Complex i1 = std::complex<double>(0.0, 1.0);
  // loop through the elements
  for (int idx = _el; idx < _eu; ++idx) {
    if ((_mesh.ELR(idx) <= rad_r) && (_mesh.EUR(idx) >= rad_r)) {
      std::vector<double> vec_nodes(_mesh.NN(), 0.0);
      for (int idxn = 0; idxn < _mesh.NN(); ++idxn) {
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      }
      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());
      for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
        for (int idxl = 1; idxl < _lmax + 1; ++idxl) {
          for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
            Complex ylm = ylmn(idxl, idxm, -1, phi_r);
            Complex ylp = ylmn(idxl, idxm, 1, phi_r);
            std::size_t ovidx = mat_indices.mat_index(idx, idxq, idxl, idxm);
            vec_receiver(ovidx) = -i1 * pleg(idxq, rad_r) * (ylm + ylp) / 2.0;
          }
        }
      }
    }
  }
  return vec_receiver;
};

Eigen::VectorXcd
sem::ReceiverVectorThetaSurface(double rad_r, double theta_r, double phi_r) {
  // create the receiver vector
  Eigen::VectorXcd vec_receiver(totlen);
  vec_receiver = Eigen::VectorXcd::Zero(totlen);
  auto wigdmat =
      GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                       GSHTrans::Single, GSHTrans::ColumnMajor>(_lmax, _lmax, 2,
                                                                theta_r);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
    return ylm;
  };
  Complex i1 = std::complex<double>(0.0, 1.0);

  int idx = _eu - 1;
  int idxq = _mesh.NN() - 1;
  for (int idxl = 1; idxl < _lmax + 1; ++idxl) {
    for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
      Complex ylm = ylmn(idxl, idxm, -1, phi_r);
      Complex ylp = ylmn(idxl, idxm, 1, phi_r);
      std::size_t ovidx = mat_indices.mat_index(idx, idxq, idxl, idxm);
      vec_receiver(ovidx) = -i1 * (ylm + ylp) / 2.0;
    }
  }
  return vec_receiver;
};

Eigen::MatrixXcd
sem::ReceiverVectorThetaSurfaceL(double rad_r, double theta_r, double phi_r,
                                 int idxl) {

  std::size_t flen = mat_indices.mat_size_l();
  std::size_t fcols = 2 * idxl + 1;
  // create the receiver vector
  Eigen::MatrixXcd vec_receiver(flen, fcols);
  vec_receiver = Eigen::MatrixXcd::Zero(flen, fcols);
  auto wigdmat =
      GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                       GSHTrans::Single, GSHTrans::ColumnMajor>(_lmax, _lmax, 2,
                                                                theta_r);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
    return ylm;
  };
  Complex i1 = std::complex<double>(0.0, 1.0);

  int idx = _eu - 1;
  int idxq = _mesh.NN() - 1;
  std::size_t ovidx = mat_indices.mat_index_l(idx, idxq);
  for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
    Complex ylm = ylmn(idxl, idxm, -1, phi_r);
    Complex ylp = ylmn(idxl, idxm, 1, phi_r);
    vec_receiver(ovidx, idxm + idxl) = -i1 * (ylm + ylp) / 2.0;
  }
  return vec_receiver;
};

Eigen::VectorXcd
sem::ReceiverVectorPhi(double rad_r, double theta_r, double phi_r) {
  // create the receiver vector
  Eigen::VectorXcd vec_receiver(totlen);
  vec_receiver = Eigen::VectorXcd::Zero(totlen);
  auto wigdmat =
      GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                       GSHTrans::Single, GSHTrans::ColumnMajor>(_lmax, _lmax, 2,
                                                                theta_r);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
    return ylm;
  };
  Complex i1 = std::complex<double>(0.0, 1.0);
  // loop through the elements
  for (int idx = _el; idx < _eu; ++idx) {
    if ((_mesh.ELR(idx) <= rad_r) && (_mesh.EUR(idx) >= rad_r)) {
      std::vector<double> vec_nodes(_mesh.NN(), 0.0);
      for (int idxn = 0; idxn < _mesh.NN(); ++idxn) {
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      }
      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());
      for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
        for (int idxl = 1; idxl < _lmax + 1; ++idxl) {
          for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
            Complex ylm = ylmn(idxl, idxm, -1, phi_r);
            Complex ylp = ylmn(idxl, idxm, 1, phi_r);
            std::size_t ovidx = mat_indices.mat_index(idx, idxq, idxl, idxm);
            vec_receiver(ovidx) = pleg(idxq, rad_r) * (ylp - ylm) / 2.0;
          }
        }
      }
    }
  }
  return vec_receiver;
};

Eigen::VectorXcd
sem::ReceiverVectorPhiSurface(double rad_r, double theta_r, double phi_r) {
  // create the receiver vector
  Eigen::VectorXcd vec_receiver(totlen);
  vec_receiver = Eigen::VectorXcd::Zero(totlen);
  auto wigdmat =
      GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                       GSHTrans::Single, GSHTrans::ColumnMajor>(_lmax, _lmax, 2,
                                                                theta_r);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
    return ylm;
  };
  Complex i1 = std::complex<double>(0.0, 1.0);

  int idx = _eu - 1;
  int idxq = _mesh.NN() - 1;
  for (int idxl = 1; idxl < _lmax + 1; ++idxl) {
    for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
      Complex ylm = ylmn(idxl, idxm, -1, phi_r);
      Complex ylp = ylmn(idxl, idxm, 1, phi_r);
      std::size_t ovidx = mat_indices.mat_index(idx, idxq, idxl, idxm);
      vec_receiver(ovidx) = (ylp - ylm) / 2.0;
    }
  }
  return vec_receiver;
};

Eigen::MatrixXcd
sem::ReceiverVectorPhiSurfaceL(double rad_r, double theta_r, double phi_r,
                               int idxl) {
  std::size_t flen = mat_indices.mat_size_l();
  std::size_t fcols = 2 * idxl + 1;
  // create the receiver vector
  Eigen::MatrixXcd vec_receiver(flen, fcols);
  vec_receiver = Eigen::MatrixXcd::Zero(flen, fcols);
  auto wigdmat =
      GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                       GSHTrans::Single, GSHTrans::ColumnMajor>(_lmax, _lmax, 2,
                                                                theta_r);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
    return ylm;
  };
  Complex i1 = std::complex<double>(0.0, 1.0);

  int idx = _eu - 1;
  int idxq = _mesh.NN() - 1;
  std::size_t ovidx = mat_indices.mat_index_l(idx, idxq);
  // for (int idxl = 1; idxl < _lmax + 1; ++idxl) {
  for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
    Complex ylm = ylmn(idxl, idxm, -1, phi_r);
    Complex ylp = ylmn(idxl, idxm, 1, phi_r);

    vec_receiver(ovidx, idxm + idxl) = (ylp - ylm) / 2.0;
  }
  // }
  return vec_receiver;
};

auto
sem::PrintModesUpToFreq(double fmax) const {
  for (int idxl = 0; idxl < _lmax; ++idxl) {
    Eigen::VectorXcd eigvals = evalues_gen_l[idxl];
    for (int i = 0; i < eigvals.size(); ++i) {
      double freq = std::sqrt(std::abs(eigvals(eigvals.size() - 1 - i))) *
                    _freq_norm / (2.0 * EIGEN_PI);
      if (freq <= fmax) {
        std::cout << "l = " << (idxl + 1) << " , overtone = " << i
                  << ", freq = " << freq << "\n";
      }
    }
  }
};
// ...existing code...

auto
sem::FindModesForCoupling(double fmax) {
  for (int idxl = 0; idxl < _lmax; ++idxl) {
    Eigen::VectorXcd eigvals = evalues_gen_l[idxl];
    for (int i = 0; i < eigvals.size(); ++i) {
      double freq = std::sqrt(std::abs(eigvals(eigvals.size() - 1 - i))) *
                    _freq_norm / (2.0 * EIGEN_PI);
      if (freq <= fmax) {
        if (!((idxl == 0) && (i == 0))) {
          nmodes += 2 * (idxl + 1) + 1;   // number of modes for this l
          vec_store.push_back(TStore(idxl + 1, i, freq));
        }
      }
    }
  }
  std::sort(vec_store.begin(), vec_store.end(),
            [](const TStore &a, const TStore &b) {
              return a.GetFreq() < b.GetFreq();
            });
  for (auto idx : vec_store) {
    // std::cout << "l = " << idx.GetL() << ", n = " << idx.GetN()
    //           << ", freq = " << idx.GetFreq() << "\n";
  }
  vec_indices.push_back(0);   // first index is 0
  for (int idx = 0; idx < vec_store.size() - 1; ++idx) {
    int idxl = vec_store[idx].GetL();
    int tmp = vec_indices[idx] + (2 * idxl + 1);
    vec_indices.push_back(tmp);   // number of modes for this l
  }
};

// function that gets the matrices for modes up to a certain frequency
template <class model1d>
Eigen::MatrixXcd
sem::NMC_INERTIA(const model1d &mod_new, bool aug_inertia) {

  Eigen::MatrixXcd retmat;

  // total size
  auto totnum = nmodes;
  if (aug_inertia) {
    totnum += ((_lmax + 1) * (_lmax + 1) - 1) * numaug;
  }

  // resize the matrix
  retmat.resize(totnum, totnum);
  retmat.setZero();

  // fill the matrix
  // inertia matrix
  std::size_t ovidx = 0;
  auto q = _mesh.GLL();
  // couple normal modes together
  for (int is = 0; is < vec_store.size(); ++is) {
    int idxl = vec_store[is].GetL();
    int idxn = vec_store[is].GetN();
    for (int js = 0; js < vec_store.size(); ++js) {
      int idxl2 = vec_store[js].GetL();
      int idxn2 = vec_store[js].GetN();

      // if the l values are the same, we can fill the matrix
      if (idxl == idxl2) {
        double tmp = 0.0;
        for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
          double elem_width = _mesh.EW(idxe);
          double d_val = 2.0 / elem_width;
          int laynum = _mesh.LayerNumber(idxe);
          double tmp1 = 0.0;
          for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
            double cr = _mesh.NodeRadius(idxe, idxq);
            double tmp2 = vec_eigval_l[idxl - 1][idxe][idxq][idxn] *
                          vec_eigval_l[idxl - 1][idxe][idxq][idxn2];
            tmp2 *= mod_new.Density(laynum)(cr) * cr * cr;
            tmp1 += q.W(idxq) * tmp2;
          }
          tmp += tmp1 * elem_width * 0.5;
        }

        // go through m values
        for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
          int idx1 = vec_indices[is] + idxm + idxl;
          int idx2 = vec_indices[js] + idxm + idxl;
          if (idx1 != idx2) {
            retmat(idx1, idx2) = tmp;
            retmat(idx2, idx1) = tmp;
          } else {
            retmat(idx1, idx2) = tmp;
          }
        }
      }
    }
  }

  auto idxouter = vec_indices.back() + (2 * vec_store.back().GetL() + 1);
  if (aug_inertia) {   // cross augmentation with modes
    for (int idxl = 0; idxl < _lmax; ++idxl) {
      for (int idxn = 0; idxn < numaug; ++idxn) {
        auto lval = idxl + 1;
        for (int is = 0; is < vec_store.size(); ++is) {
          if (lval == vec_store[is].GetL()) {
            double tmp = 0.0;
            for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
              double elem_width = _mesh.EW(idxe);
              double d_val = 2.0 / elem_width;
              int laynum = _mesh.LayerNumber(idxe);
              double tmp1 = 0.0;
              for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
                double cr = _mesh.NodeRadius(idxe, idxq);
                double tmp2 =
                    vec_augval_l[idxl][idxe][idxq][idxn] *
                    vec_eigval_l[idxl][idxe][idxq][vec_store[is].GetN()];
                tmp2 *= mod_new.Density(laynum)(cr) * cr * cr;
                tmp1 += q.W(idxq) * tmp2;
              }
              tmp += tmp1 * elem_width * 0.5;
            }

            // go through m values
            for (int idxm = -lval; idxm < lval + 1; ++idxm) {
              int idx1 = vec_indices[is] + idxm + lval;
              int idx2 = idxouter;
              idx2 += (lval * lval - 1) * numaug +
                      (2 * lval + 1) * idxn;   // offset for previous l values
              idx2 += idxm + lval;
              if (idx1 != idx2) {
                retmat(idx1, idx2) = tmp;
                retmat(idx2, idx1) = tmp;
              } else {
                retmat(idx1, idx2) = tmp;
              }
            }
          }
        }
      }
    }

    // cross augmentation functions against each other
    for (int idxl = 0; idxl < _lmax; ++idxl) {
      for (int idxn = 0; idxn < numaug; ++idxn) {
        for (int idxn2 = 0; idxn2 < numaug; ++idxn2) {
          auto lval = idxl + 1;

          double tmp = 0.0;
          for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
            double elem_width = _mesh.EW(idxe);
            double d_val = 2.0 / elem_width;
            int laynum = _mesh.LayerNumber(idxe);
            double tmp1 = 0.0;
            for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
              double cr = _mesh.NodeRadius(idxe, idxq);
              double tmp2 = vec_augval_l[idxl][idxe][idxq][idxn] *
                            vec_augval_l[idxl][idxe][idxq][idxn2];
              tmp2 *= mod_new.Density(laynum)(cr) * cr * cr;
              tmp1 += q.W(idxq) * tmp2;
            }
            tmp += tmp1 * elem_width * 0.5;
          }

          // go through m values
          for (int idxm = -lval; idxm < lval + 1; ++idxm) {
            int idx1 = idxouter;
            idx1 += (lval * lval - 1) * numaug +
                    (2 * lval + 1) * idxn;   // offset for previous l values
            idx1 += idxm + lval;
            int idx2 = idxouter;
            idx2 += (lval * lval - 1) * numaug +
                    (2 * lval + 1) * idxn2;   // offset for previous l values
            idx2 += idxm + lval;
            if (idx1 != idx2) {
              retmat(idx1, idx2) = tmp;
              retmat(idx2, idx1) = tmp;
            } else {
              retmat(idx1, idx2) = tmp;
            }
          }
        }
      }
    }
  }

  return retmat;
};

// function that gets the matrices for modes up to a certain frequency
template <class model1d>
Eigen::MatrixXcd
sem::NMC_KE(const model1d &mod_new, bool aug_ke) {

  Eigen::MatrixXcd retmat;

  // total size
  auto totnum = nmodes;
  if (aug_ke) {
    totnum += ((_lmax + 1) * (_lmax + 1) - 1) * numaug;
  }

  // resize the matrix
  retmat.resize(totnum, totnum);
  retmat.setZero();

  // fill the matrix
  // ke matrix
  std::size_t ovidx = 0;
  auto q = _mesh.GLL();
  for (int is = 0; is < vec_store.size(); ++is) {
    int idxl = vec_store[is].GetL();
    int idxn = vec_store[is].GetN();
    for (int js = 0; js < vec_store.size(); ++js) {
      int idxl2 = vec_store[js].GetL();
      int idxn2 = vec_store[js].GetN();

      // if the l values are the same, we can fill the matrix
      if (idxl == idxl2) {
        double tmp = 0.0;
        auto k2 = (idxl * (idxl + 1));
        for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
          double elem_width = _mesh.EW(idxe);
          double d_val = 2.0 / elem_width;
          int laynum = _mesh.LayerNumber(idxe);
          double tmp1 = 0.0;
          for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
            double cr = _mesh.NodeRadius(idxe, idxq);

            // first integral
            double tmp2 = cr * vec_eigderiv_l[idxl - 1][idxe][idxq][idxn] -
                          vec_eigval_l[idxl - 1][idxe][idxq][idxn];
            double tmp3 = cr * vec_eigderiv_l[idxl - 1][idxe][idxq][idxn2] -
                          vec_eigval_l[idxl - 1][idxe][idxq][idxn2];
            double tmp4 = tmp2 * tmp3 * mod_new.L(laynum)(cr);

            // second integral
            double tmp5 = (k2 - 2.0) *
                          vec_eigval_l[idxl - 1][idxe][idxq][idxn] *
                          vec_eigval_l[idxl - 1][idxe][idxq][idxn2];
            tmp5 *= mod_new.N(laynum)(cr);

            // sum
            tmp1 += q.W(idxq) * (tmp4 + tmp5);
          }
          tmp += tmp1 * elem_width * 0.5;
        }

        // go through m values
        for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
          int idx1 = vec_indices[is] + idxm + idxl;
          int idx2 = vec_indices[js] + idxm + idxl;
          if (idx1 != idx2) {
            retmat(idx1, idx2) = tmp;
            retmat(idx2, idx1) = tmp;
          } else {
            retmat(idx1, idx2) = tmp;
          }
        }
      }
    }
  }

  auto idxouter = vec_indices.back() + (2 * vec_store.back().GetL() + 1);
  if (aug_ke) {   // cross augmentation with modes
    for (int idxl = 0; idxl < _lmax; ++idxl) {
      for (int idxn = 0; idxn < numaug; ++idxn) {
        auto lval = idxl + 1;
        for (int is = 0; is < vec_store.size(); ++is) {
          int idxn2 = vec_store[is].GetN();
          if (lval == vec_store[is].GetL()) {

            double tmp = 0.0;
            auto k2 = (lval * (lval + 1));
            for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
              double elem_width = _mesh.EW(idxe);
              double d_val = 2.0 / elem_width;
              int laynum = _mesh.LayerNumber(idxe);
              double tmp1 = 0.0;
              for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
                double cr = _mesh.NodeRadius(idxe, idxq);

                // first integral
                double tmp2 = cr * vec_augderiv_l[idxl][idxe][idxq][idxn] -
                              vec_augval_l[idxl][idxe][idxq][idxn];
                double tmp3 = cr * vec_eigderiv_l[idxl][idxe][idxq][idxn2] -
                              vec_eigval_l[idxl][idxe][idxq][idxn2];
                double tmp4 = tmp2 * tmp3 * mod_new.L(laynum)(cr);

                // second integral
                double tmp5 = (k2 - 2.0) *
                              vec_augval_l[idxl][idxe][idxq][idxn] *
                              vec_eigval_l[idxl][idxe][idxq][idxn2];
                tmp5 *= mod_new.N(laynum)(cr);

                // sum
                tmp1 += q.W(idxq) * (tmp4 + tmp5);
              }
              tmp += tmp1 * elem_width * 0.5;
            }

            // go through m values
            for (int idxm = -lval; idxm < lval + 1; ++idxm) {
              int idx1 = vec_indices[is] + idxm + lval;
              int idx2 = idxouter;
              idx2 += (lval * lval - 1) * numaug +
                      (2 * lval + 1) * idxn;   // offset for previous l values
              idx2 += idxm + lval;
              if (idx1 != idx2) {
                retmat(idx1, idx2) = tmp;
                retmat(idx2, idx1) = tmp;
              } else {
                retmat(idx1, idx2) = tmp;
              }
            }
          }
        }
      }
    }

    // cross augmentation functions against each other
    for (int idxl = 0; idxl < _lmax; ++idxl) {
      for (int idxn = 0; idxn < numaug; ++idxn) {
        for (int idxn2 = 0; idxn2 < numaug; ++idxn2) {
          auto lval = idxl + 1;

          double tmp = 0.0;
          auto k2 = (lval * (lval + 1));
          for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
            double elem_width = _mesh.EW(idxe);
            double d_val = 2.0 / elem_width;
            int laynum = _mesh.LayerNumber(idxe);
            double tmp1 = 0.0;
            for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
              double cr = _mesh.NodeRadius(idxe, idxq);

              // first integral
              double tmp2 = cr * vec_augderiv_l[idxl][idxe][idxq][idxn] -
                            vec_augval_l[idxl][idxe][idxq][idxn];
              double tmp3 = cr * vec_augderiv_l[idxl][idxe][idxq][idxn2] -
                            vec_augval_l[idxl][idxe][idxq][idxn2];
              double tmp4 = tmp2 * tmp3 * mod_new.L(laynum)(cr);

              // second integral
              double tmp5 = (k2 - 2.0) * vec_augval_l[idxl][idxe][idxq][idxn] *
                            vec_augval_l[idxl][idxe][idxq][idxn2];
              tmp5 *= mod_new.N(laynum)(cr);

              // sum
              tmp1 += q.W(idxq) * (tmp4 + tmp5);
            }
            tmp += tmp1 * elem_width * 0.5;
          }

          // go through m values
          for (int idxm = -lval; idxm < lval + 1; ++idxm) {
            int idx1 = idxouter;
            idx1 += (lval * lval - 1) * numaug +
                    (2 * lval + 1) * idxn;   // offset for previous l values
            idx1 += idxm + lval;
            int idx2 = idxouter;
            idx2 += (lval * lval - 1) * numaug +
                    (2 * lval + 1) * idxn2;   // offset for previous l values
            idx2 += idxm + lval;
            if (idx1 != idx2) {
              retmat(idx1, idx2) = tmp;
              retmat(idx2, idx1) = tmp;
            } else {
              retmat(idx1, idx2) = tmp;
            }
          }
        }
      }
    }
  }

  return retmat;
};

auto
sem::NMC_FORCE(SourceInfo::EarthquakeCMT &cmt, bool aug_force) {

  auto totnum = nmodes;
  if (aug_force) {
    totnum += ((_lmax + 1) * (_lmax + 1) - 1) * numaug;
  }

  Eigen::VectorXcd vec_force = Eigen::VectorXcd::Zero(totnum);

  // resize force
  // vec_force.resize(totlen);
  // vec_force = Eigen::VectorXcd::Zero(totlen);

  // find element within which the source sits
  double depth = cmt.Depth();
  double rad_source = _mesh.PR() - 1000.0 * depth / _length_norm;
  // std::cout << "Source radius: " << rad_source << "\n";
  int idxsource = 0;

  int NQ = _mesh.NN();

  // get the y0-, y0+ values etc at the source location
  double theta_s = (90.0 - cmt.Latitude()) * EIGEN_PI / (180.0);
  double phi_s = cmt.Longitude() * EIGEN_PI / (180.0);

  auto wigdmat =
      GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                       GSHTrans::Single, GSHTrans::ColumnMajor>(_lmax, _lmax, 2,
                                                                theta_s);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
    return ylm;
  };
  // std::cout << "Source location: " << theta_s << " " << phi_s << "\n";
  // int m = 0, _l = 2;
  // parameters
  double invsqrt2 = 1.0 / std::sqrt(2.0);
  std::complex<double> isq2 = std::complex<double>(0.0, invsqrt2);

  // loop through the elements
  //   to find the element that contains the source
  for (int is = 0; is < vec_store.size(); ++is) {
    int idxl = vec_store[is].GetL();
    int idxn = vec_store[is].GetN();
    for (int idx = _el; idx < _eu; ++idx) {
      // std::cout << _mesh.ELR(idx) << " " << _mesh.EUR(idx) << " " <<
      // rad_source
      //           << "\n";
      if ((_mesh.ELR(idx) <= rad_source) && (_mesh.EUR(idx) >= rad_source)) {
        // std::cout << idx << " " << rad_source << " " << _mesh.ELR(idx) << "
        // "
        //           << _mesh.EUR(idx) << "\n";
        std::vector<double> vec_nodes(NQ, 0.0);
        for (int idxv = 0; idxv < NQ; ++idxv) {
          vec_nodes[idxv] = _mesh.NodeRadius(idx, idxv);
        }

        auto pleg = Interpolation::LagrangePolynomial(vec_nodes.begin(),
                                                      vec_nodes.end());

        auto w_val = 0.0;
        auto w_deriv = 0.0;
        for (int idxq = 0; idxq < NQ; ++idxq) {
          w_val +=
              vec_eigval_l[idxl - 1][idx][idxq][idxn] * pleg(idxq, rad_source);
          w_deriv += vec_eigval_l[idxl - 1][idx][idxq][idxn] *
                     pleg.Derivative(idxq, rad_source);
        };
        auto w_vald = w_val / rad_source;
        auto w_prefactor = w_deriv - w_vald;
        double omegal2 = (idxl + 2) * (idxl - 1) / 2.0;
        double lprefac =
            std::exp(-2.0 * 3.141592653589793 * (idxl + 1) / (1 + 0.5));
        // lprefac = 1.0;
        for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
          // spherical harmonic
          Complex ymc = std::conj(ylmn(idxl, idxm, -1, phi_s));
          Complex ypc = std::conj(ylmn(idxl, idxm, 1, phi_s));
          Complex ymmc = 0.0, yppc = 0.0;

          if (idxl > 1) {
            ymmc = std::conj(ylmn(idxl, idxm, -2, phi_s));
            yppc = std::conj(ylmn(idxl, idxm, 2, phi_s));
          }

          Complex tmp = w_prefactor * (cmt.MC0m() * ymc - cmt.MC0p() * ypc);
          tmp += w_vald * omegal2 * (cmt.MCmm() * ymmc - cmt.MCpp() * yppc);
          tmp *= isq2 * lprefac;

          vec_force(vec_indices[is] + idxm + idxl) = tmp;
        };
      };
    };
  };

  auto idxouter = vec_indices.back() + (2 * vec_store.back().GetL() + 1);
  if (aug_force) {   // cross augmentation with modes

    for (int idxl = 0; idxl < _lmax; ++idxl) {
      for (int idxn = 0; idxn < numaug; ++idxn) {
        auto lval = idxl + 1;
        for (int idx = _el; idx < _eu; ++idx) {
          // std::cout << _mesh.ELR(idx) << " " << _mesh.EUR(idx) << " " <<
          // rad_source
          //           << "\n";
          if ((_mesh.ELR(idx) <= rad_source) &&
              (_mesh.EUR(idx) >= rad_source)) {
            // std::cout << idx << " " << rad_source << " " << _mesh.ELR(idx)
            // << "
            // "
            //           << _mesh.EUR(idx) << "\n";
            std::vector<double> vec_nodes(NQ, 0.0);
            for (int idxv = 0; idxv < NQ; ++idxv) {
              vec_nodes[idxv] = _mesh.NodeRadius(idx, idxv);
            }

            auto pleg = Interpolation::LagrangePolynomial(vec_nodes.begin(),
                                                          vec_nodes.end());

            auto w_val = 0.0;
            auto w_deriv = 0.0;
            for (int idxq = 0; idxq < NQ; ++idxq) {
              w_val +=
                  vec_augval_l[idxl][idx][idxq][idxn] * pleg(idxq, rad_source);
              w_deriv += vec_augval_l[idxl][idx][idxq][idxn] *
                         pleg.Derivative(idxq, rad_source);
            };
            auto w_vald = w_val / rad_source;
            auto w_prefactor = w_deriv - w_vald;
            double omegal2 = (lval + 2) * (lval - 1) / 2.0;
            double lprefac =
                std::exp(-2.0 * 3.141592653589793 * (idxl + 1) / (1 + 0.5));
            // lprefac = 1.0;
            for (int idxm = -lval; idxm < lval + 1; ++idxm) {
              // std::cout << "l = " << lval << ", n = " << idxn
              //           << ", m = " << idxm << "\n";
              // spherical harmonic
              Complex ymc = std::conj(ylmn(lval, idxm, -1, phi_s));
              Complex ypc = std::conj(ylmn(lval, idxm, 1, phi_s));
              Complex ymmc = 0.0, yppc = 0.0;

              if (lval > 1) {
                ymmc = std::conj(ylmn(lval, idxm, -2, phi_s));
                yppc = std::conj(ylmn(lval, idxm, 2, phi_s));
              }

              Complex tmp = w_prefactor * (cmt.MC0m() * ymc - cmt.MC0p() * ypc);
              tmp += w_vald * omegal2 * (cmt.MCmm() * ymmc - cmt.MCpp() * yppc);
              tmp *= isq2;

              int idx1 = idxouter;
              idx1 += (lval * lval - 1) * numaug +
                      (2 * lval + 1) * idxn;   // offset for previous l values
              idx1 += idxm + lval;
              vec_force(idx1) = tmp * lprefac;
            };
          };
        }
      };
    };
  }

  return vec_force;
};

Eigen::VectorXcd
sem::ReceiverVectorThetaSurfaceCoupling(double rad_r, double theta_r,
                                        double phi_r, bool toaug) {
  // std::size_t flen = mat_indices.mat_size_l();
  // std::size_t fcols = 2 * idxl + 1;
  // create the receiver vector
  auto totnum = nmodes;
  if (toaug) {
    totnum += ((_lmax + 1) * (_lmax + 1) - 1) * numaug;
  }
  Eigen::VectorXcd vec_receiver(totnum);
  vec_receiver = Eigen::VectorXcd::Zero(totnum);
  auto wigdmat =
      GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                       GSHTrans::Single, GSHTrans::ColumnMajor>(_lmax, _lmax, 2,
                                                                theta_r);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
    return ylm;
  };
  Complex i1 = std::complex<double>(0.0, 1.0);

  int idx = _eu - 1;
  int idxq = _mesh.NN() - 1;
  // std::size_t ovidx = mat_indices.mat_index_l(idx, idxq);
  for (int is = 0; is < vec_store.size(); ++is) {
    int idxl = vec_store[is].GetL();
    int idxn = vec_store[is].GetN();
    // for (int idxl = 1; idxl < _lmax + 1; ++idxl) {
    for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
      Complex ylm = ylmn(idxl, idxm, -1, phi_r);
      Complex ylp = ylmn(idxl, idxm, 1, phi_r);

      vec_receiver(vec_indices[is] + idxm + idxl) =
          vec_eigval_l[idxl - 1][idx][idxq][idxn] * -i1 * (ylm + ylp) / 2.0;
    }
    // }
  }
  auto idxouter = vec_indices.back() + (2 * vec_store.back().GetL() + 1);
  if (toaug) {   // cross augmentation with modes

    for (int idxl = 0; idxl < _lmax; ++idxl) {
      for (int idxn = 0; idxn < numaug; ++idxn) {
        auto lval = idxl + 1;
        for (int idxm = -lval; idxm < lval + 1; ++idxm) {
          Complex ylm = ylmn(lval, idxm, -1, phi_r);
          Complex ylp = ylmn(lval, idxm, 1, phi_r);

          int idx1 = idxouter;
          idx1 += (lval * lval - 1) * numaug +
                  (2 * lval + 1) * idxn;   // offset for previous l values
          idx1 += idxm + lval;
          vec_receiver(idx1) =
              vec_augval_l[idxl][idx][idxq][idxn] * -i1 * (ylm + ylp) / 2.0;
        }
      }
    }
  }
  return vec_receiver;
};

Eigen::VectorXcd
sem::ReceiverVectorPhiSurfaceCoupling(double rad_r, double theta_r,
                                      double phi_r, bool toaug) {
  // std::size_t flen = mat_indices.mat_size_l();
  // std::size_t fcols = 2 * idxl + 1;
  // create the receiver vector
  auto totnum = nmodes;
  if (toaug) {
    totnum += ((_lmax + 1) * (_lmax + 1) - 1) * numaug;
  }

  Eigen::VectorXcd vec_receiver(totnum);
  vec_receiver = Eigen::VectorXcd::Zero(totnum);
  auto wigdmat =
      GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                       GSHTrans::Single, GSHTrans::ColumnMajor>(_lmax, _lmax, 2,
                                                                theta_r);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
    return ylm;
  };
  Complex i1 = std::complex<double>(0.0, 1.0);

  int idx = _eu - 1;
  int idxq = _mesh.NN() - 1;
  // std::size_t ovidx = mat_indices.mat_index_l(idx, idxq);
  for (int is = 0; is < vec_store.size(); ++is) {
    int idxl = vec_store[is].GetL();
    int idxn = vec_store[is].GetN();
    // for (int idxl = 1; idxl < _lmax + 1; ++idxl) {
    for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
      Complex ylm = ylmn(idxl, idxm, -1, phi_r);
      Complex ylp = ylmn(idxl, idxm, 1, phi_r);

      vec_receiver(vec_indices[is] + idxm + idxl) =
          vec_eigval_l[idxl - 1][idx][idxq][idxn] * (ylp - ylm) / 2.0;
    }
    // }
  }

  auto idxouter = vec_indices.back() + (2 * vec_store.back().GetL() + 1);
  if (toaug) {   // cross augmentation with modes
    for (int idxl = 0; idxl < _lmax; ++idxl) {
      for (int idxn = 0; idxn < numaug; ++idxn) {
        auto lval = idxl + 1;
        for (int idxm = -lval; idxm < lval + 1; ++idxm) {
          Complex ylm = ylmn(lval, idxm, -1, phi_r);
          Complex ylp = ylmn(lval, idxm, 1, phi_r);

          int idx1 = idxouter;
          idx1 += (lval * lval - 1) * numaug +
                  (2 * lval + 1) * idxn;   // offset for previous l values
          idx1 += idxm + lval;
          vec_receiver(idx1) =
              vec_augval_l[idxl][idx][idxq][idxn] * (ylp - ylm) / 2.0;
        }
      }
    }
  }
  return vec_receiver;
};

void
sem::augment_basis_calculate() {
  std::cout << "Calculating augmented basis\n";

  std::size_t ncols = _mesh.NL() + 1;
  _augbasis_l.resize(_lmax);
  vec_augval_l = std::vector<stdvvvec>(_lmax);
  vec_augderiv_l = std::vector<stdvvvec>(_lmax);
  _calc_l = std::vector<bool>(_lmax, false);
  int lowidx = 0;
  if ((_has_fluid && (_il == 0)) || (!_has_fluid)) {
    lowidx = 1;
  } else {
    lowidx = 0;
  }
  // std::cout << "Check 1\n";
  // number of layers within layer
  // int num_aug = _mesh.LayerNumber(_eu - 1) - _mesh.LayerNumber(_el) + 2;
  // numaug = num_aug;
  // std::cout << "\nNumber of augmentation functions: " << num_aug << "\n\n";
  // std::cout << "Check 2\n";
  for (int idx = 0; idx < _lmax; ++idx) {
    // std::cout << "Check 3, l = " << (idx + 1) << "\n";
    auto mat_keuse = mat_l_ke[idx];
    Eigen::SparseLU<Eigen::SparseMatrix<double>> _solver;
    _solver.analyzePattern(mat_keuse);
    _solver.factorize(mat_keuse);
    std::size_t nrows = lowidx + mat_keuse.rows();
    _augbasis_l[idx] = Eigen::MatrixXd::Zero(nrows, numaug);
    // need to account for fact that in solid we only find the result in one
    // solid
    // section
    // whether or not we're looking at IC or purely solid planet

    // std::cout << "Check 4\n";
    // force matrix, initialise and set (0,0) element (for first boundary)
    Eigen::MatrixXd mat_force = Eigen::MatrixXd::Zero(mat_keuse.cols(), numaug);
    // std::cout << "mat_force size: " << mat_force.rows() << " "
    //           << mat_force.cols() << "\n\n";
    // mat_force(0, 0) = -1;
    int colnum = 0;
    auto idxlow = _el;
    if (idxlow == 0) {
      idxlow += 1;
    }
    for (int idxe = idxlow; idxe < _eu; ++idxe) {
      if ((_mesh.LayerNumber(idxe) - _mesh.LayerNumber(idxe - 1)) == 1) {
        std::size_t ovidx;
        ovidx = mat_indices.mat_index_l(idxe, 0);
        // std::cout << "Setting force at element " << idxe << " row " << ovidx
        //           << " col " << colnum << "\n";
        mat_force(ovidx, colnum) = -1.0;
        colnum += 1;
      }
    }
    // mat_force(mat_ke.cols() - 1, colnum) = -1;
    // std::cout << "Check 5\n";
    // std::cout << "augbasis size: " << _augbasis_l[idx].rows() << " "
    //           << _augbasis_l[idx].cols() << "\n";
    // std::cout << "Block trying to make: " << lowidx << " " <<
    // mat_keuse.cols()
    //           << " " << num_aug << "\n";
    // calculate basis
    _augbasis_l[idx].block(lowidx, 0, mat_keuse.cols(), numaug) =
        _solver.solve(mat_force);
    _calc_l[idx] = true;
    // Eigen::MatrixXcd normval =
    //     _augbasis.block(lowidx, 0, mat_ke.cols(), num_aug).transpose() *
    //     mat_inertia * _augbasis.block(lowidx, 0, mat_ke.cols(), num_aug);

    // std::cout << std::setprecision(4) << normval << "\n";
    // fill out augval and augderiv
    //  save the augmented basis in the "standard format":
    // std::cout << "Check 6\n";
    int NQ = _mesh.NN();
    // std::cout << "num_aug: " << num_aug << "\n";
    // std::cout << "NQ: " << NQ << "\n";
    // std::cout << "_mesh.NE(): " << _mesh.NE() << "\n";
    vec_augval_l[idx] = stdvvvec(_mesh.NE(), stdvvec(NQ, stdvec(numaug, 0.0)));
    // std::cout << "Check 6.1\n";
    for (int idxe = _el; idxe < _eu; ++idxe) {
      for (int idxq = 0; idxq < NQ; ++idxq) {
        std::size_t ovidx = (idxe - _el) * (NQ - 1) + idxq;
        for (int idxa = 0; idxa < numaug; ++idxa) {
          // std::cout << "idxe, idxq, idxa: " << idxe << " " << idxq << " "
          //           << idxa << "\n";
          vec_augval_l[idx][idxe][idxq][idxa] = _augbasis_l[idx](ovidx, idxa);
        }
      }
    }
    // calculate the derivatives of the augmented basis:
    vec_augderiv_l[idx] =
        stdvvvec(_mesh.NE(), stdvvec(NQ, stdvec(numaug, 0.0)));
    // std::cout << "Check 7\n";
    for (int idxe = _el; idxe < _eu; ++idxe) {
      double elem_width = _mesh.EW(idxe);
      for (int idxq = 0; idxq < NQ; ++idxq) {
        for (int idxa = 0; idxa < numaug; ++idxa) {
          double tmp = 0.0;
          for (int idxq2 = 0; idxq2 < NQ; ++idxq2) {
            tmp += vec_augval_l[idx][idxe][idxq2][idxa] *
                   vec_lag_deriv[idxq][idxq2];
          }
          tmp *= 2.0 / elem_width;
          vec_augderiv_l[idx][idxe][idxq][idxa] = tmp;
        }
      }
    }
  }

  // we now need to add the augmentation to the eigenvector matrix
  for (int idxl = 0; idxl < _lmax; ++idxl) {
    for (int idxe = _el; idxe < _eu; ++idxe) {
      int NQ = _mesh.NN();
      for (int idxq = 0; idxq < NQ; ++idxq) {
        for (int idxa = 0; idxa < numaug; ++idxa) {
          vec_all_l[idxl][idxe][idxq][idxa + _num_modes_l] =
              vec_augval_l[idxl][idxe][idxq][idxa];
          vec_allderiv_l[idxl][idxe][idxq][idxa + _num_modes_l] =
              vec_augderiv_l[idxl][idxe][idxq][idxa];
        }
      }
    }
  }
};

template <class model1d>
Eigen::MatrixXcd
sem::NMC_INERTIA(const model1d &mod_new, int idxl, bool toaug) {
  Eigen::MatrixXcd retmat;

  // total size
  auto totnum = _num_modes_l;
  if (toaug) {
    totnum += numaug;
  }

  // resize the matrix
  retmat.resize(totnum, totnum);
  retmat.setZero();

  // fill the matrix
  // inertia matrix
  std::size_t ovidx = 0;
  auto q = _mesh.GLL();
  // couple normal modes together

  if (!toaug) {
    for (int idxn = 0; idxn < _num_modes_l; ++idxn) {

      for (int idxn2 = 0; idxn2 < idxn + 1; ++idxn2) {

        double tmp = 0.0;
        for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
          double elem_width = _mesh.EW(idxe);
          double d_val = 2.0 / elem_width;
          int laynum = _mesh.LayerNumber(idxe);
          double tmp1 = 0.0;
          for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
            double cr = _mesh.NodeRadius(idxe, idxq);
            double tmp2 = vec_eigval_l[idxl - 1][idxe][idxq][idxn] *
                          vec_eigval_l[idxl - 1][idxe][idxq][idxn2];
            tmp2 *= mod_new.Density(laynum)(cr) * cr * cr;
            tmp1 += q.W(idxq) * tmp2;
          }
          tmp += tmp1 * elem_width * 0.5;
        }

        // go through m values

        if (idxn != idxn2) {
          retmat(idxn, idxn2) = tmp;
          retmat(idxn2, idxn) = tmp;
        } else {
          retmat(idxn, idxn2) = tmp;
        }
      }
    }
  } else if (toaug) {
    for (int idxn = 0; idxn < _num_modes_l + numaug; ++idxn) {

      for (int idxn2 = 0; idxn2 < idxn + 1; ++idxn2) {

        double tmp = 0.0;
        for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
          double elem_width = _mesh.EW(idxe);
          double d_val = 2.0 / elem_width;
          int laynum = _mesh.LayerNumber(idxe);
          double tmp1 = 0.0;
          for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
            double cr = _mesh.NodeRadius(idxe, idxq);
            double tmp2 = vec_all_l[idxl - 1][idxe][idxq][idxn] *
                          vec_all_l[idxl - 1][idxe][idxq][idxn2];
            tmp2 *= mod_new.Density(laynum)(cr) * cr * cr;
            tmp1 += q.W(idxq) * tmp2;
          }
          tmp += tmp1 * elem_width * 0.5;
        }

        // go through m values

        if (idxn != idxn2) {
          retmat(idxn, idxn2) = tmp;
          retmat(idxn2, idxn) = tmp;
        } else {
          retmat(idxn, idxn2) = tmp;
        }
      }
    }
  }
  return retmat;
};

template <class model1d>
Eigen::MatrixXcd
sem::NMC_KE(const model1d &mod_new, int idxl, bool toaug) {

  // total size
  auto totnum = _num_modes_l;
  if (toaug) {
    totnum += numaug;
  }

  // resize the matrix
  Eigen::MatrixXcd retmat = Eigen::MatrixXcd::Zero(totnum, totnum);
  // retmat.setZero();

  // fill the matrix
  // inertia matrix
  std::size_t ovidx = 0;
  auto q = _mesh.GLL();

  // couple normal modes together
  if (!toaug) {
    std::cout << "\n Not augmenting\n\n";
    for (int idxn = 0; idxn < _num_modes_l; ++idxn) {
      for (int idxn2 = 0; idxn2 < idxn + 1; ++idxn2) {

        double tmp = 0.0;
        auto k2 = (idxl * (idxl + 1));
        for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
          double elem_width = _mesh.EW(idxe);
          double d_val = 2.0 / elem_width;
          int laynum = _mesh.LayerNumber(idxe);
          double tmp1 = 0.0;
          for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
            double cr = _mesh.NodeRadius(idxe, idxq);

            // first integral
            double tmp2 = cr * vec_eigderiv_l[idxl - 1][idxe][idxq][idxn] -
                          vec_eigval_l[idxl - 1][idxe][idxq][idxn];
            double tmp3 = cr * vec_eigderiv_l[idxl - 1][idxe][idxq][idxn2] -
                          vec_eigval_l[idxl - 1][idxe][idxq][idxn2];
            double tmp4 = tmp2 * tmp3 * mod_new.L(laynum)(cr);

            // second integral
            double tmp5 = (k2 - 2.0) *
                          vec_eigval_l[idxl - 1][idxe][idxq][idxn] *
                          vec_eigval_l[idxl - 1][idxe][idxq][idxn2];
            tmp5 *= mod_new.N(laynum)(cr);

            // sum
            tmp1 += q.W(idxq) * (tmp4 + tmp5);
          }
          tmp += tmp1 * elem_width * 0.5;
        }

        // go through m values

        if (idxn != idxn2) {
          retmat(idxn, idxn2) = tmp;
          retmat(idxn2, idxn) = tmp;
        } else {
          retmat(idxn, idxn2) = tmp;
        }

        // }
        // }
      }
    }
  } else if (toaug) {
    std::cout << "\n Augmenting\n\n";
    for (int idxn = 0; idxn < totnum; ++idxn) {
      for (int idxn2 = 0; idxn2 < idxn + 1; ++idxn2) {

        double tmp = 0.0;
        auto k2 = (idxl * (idxl + 1));
        for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
          double elem_width = _mesh.EW(idxe);
          double d_val = 2.0 / elem_width;
          int laynum = _mesh.LayerNumber(idxe);
          double tmp1 = 0.0;
          for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
            double cr = _mesh.NodeRadius(idxe, idxq);

            // first integral
            double tmp2 = cr * vec_allderiv_l[idxl - 1][idxe][idxq][idxn] -
                          vec_all_l[idxl - 1][idxe][idxq][idxn];
            double tmp3 = cr * vec_allderiv_l[idxl - 1][idxe][idxq][idxn2] -
                          vec_all_l[idxl - 1][idxe][idxq][idxn2];
            double tmp4 = tmp2 * tmp3 * mod_new.L(laynum)(cr);

            // second integral
            double tmp5 = (k2 - 2.0) * vec_all_l[idxl - 1][idxe][idxq][idxn] *
                          vec_all_l[idxl - 1][idxe][idxq][idxn2];
            tmp5 *= mod_new.N(laynum)(cr);

            // sum
            tmp1 += q.W(idxq) * (tmp4 + tmp5);
          }
          tmp += tmp1 * elem_width * 0.5;
        }

        // go through m values

        if (idxn != idxn2) {
          retmat(idxn, idxn2) = tmp;
          retmat(idxn2, idxn) = tmp;
        } else {
          retmat(idxn, idxn2) = tmp;
        }

        // }
        // }
      }
    }
  }
  return retmat;
};

Eigen::MatrixXcd
sem::ReceiverVectorThetaSurface_NMCL(double rad_r, double theta_r, double phi_r,
                                     int idxl, bool toaug) {

  auto totnum = _num_modes_l;
  if (toaug) {
    totnum += numaug;
  }
  std::size_t flen = _num_modes_l;
  std::size_t fcols = 2 * idxl + 1;
  // create the receiver vector
  Eigen::MatrixXcd vec_receiver(totnum, fcols);
  vec_receiver = Eigen::MatrixXcd::Zero(totnum, fcols);
  auto wigdmat =
      GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                       GSHTrans::Single, GSHTrans::ColumnMajor>(_lmax, _lmax, 2,
                                                                theta_r);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
    return ylm;
  };
  Complex i1 = std::complex<double>(0.0, 1.0);

  int idx = _eu - 1;
  int idxq = _mesh.NN() - 1;
  if (!toaug) {
    for (int idxn = 0; idxn < _num_modes_l; ++idxn) {
      for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
        Complex ylm = ylmn(idxl, idxm, -1, phi_r);
        Complex ylp = ylmn(idxl, idxm, 1, phi_r);
        vec_receiver(idxn, idxm + idxl) =
            vec_eigval_l[idxl - 1][idx][idxq][idxn] * -i1 * (ylm + ylp) / 2.0;
      }
    }
  } else if (toaug) {
    for (int idxn = 0; idxn < totnum; ++idxn) {
      for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
        Complex ylm = ylmn(idxl, idxm, -1, phi_r);
        Complex ylp = ylmn(idxl, idxm, 1, phi_r);
        vec_receiver(idxn, idxm + idxl) =
            vec_all_l[idxl - 1][idx][idxq][idxn] * -i1 * (ylm + ylp) / 2.0;
      }
    }
  }
  return vec_receiver;
};

Eigen::MatrixXcd
sem::ReceiverVectorPhiSurface_NMCL(double rad_r, double theta_r, double phi_r,
                                   int idxl, bool toaug) {
  auto totnum = _num_modes_l;
  if (toaug) {
    totnum += numaug;
  }
  // std::size_t flen = _num_modes_l;
  std::size_t fcols = 2 * idxl + 1;
  // create the receiver vector
  Eigen::MatrixXcd vec_receiver(totnum, fcols);
  vec_receiver = Eigen::MatrixXcd::Zero(totnum, fcols);
  auto wigdmat =
      GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                       GSHTrans::Single, GSHTrans::ColumnMajor>(_lmax, _lmax, 2,
                                                                theta_r);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
    return ylm;
  };
  Complex i1 = std::complex<double>(0.0, 1.0);

  int idx = _eu - 1;
  int idxq = _mesh.NN() - 1;
  if (!toaug) {
    for (int idxn = 0; idxn < _num_modes_l; ++idxn) {
      for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
        Complex ylm = ylmn(idxl, idxm, -1, phi_r);
        Complex ylp = ylmn(idxl, idxm, 1, phi_r);
        vec_receiver(idxn, idxm + idxl) =
            vec_eigval_l[idxl - 1][idx][idxq][idxn] * (ylp - ylm) / 2.0;
      }
    }
  } else if (toaug) {
    for (int idxn = 0; idxn < totnum; ++idxn) {
      for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
        Complex ylm = ylmn(idxl, idxm, -1, phi_r);
        Complex ylp = ylmn(idxl, idxm, 1, phi_r);
        vec_receiver(idxn, idxm + idxl) =
            vec_all_l[idxl - 1][idx][idxq][idxn] * (ylp - ylm) / 2.0;
      }
    }
  }
  return vec_receiver;
};

Eigen::MatrixXcd
sem::CalculateForceNMC(SourceInfo::EarthquakeCMT &cmt, int idxl, bool toaug) {
  auto totnum = _num_modes_l;
  if (toaug) {
    totnum += numaug;
  }
  // std::size_t flen = _num_modes_l;
  std::size_t fcols = 2 * idxl + 1;
  Eigen::MatrixXcd vec_lforce = Eigen::MatrixXcd::Zero(totnum, fcols);

  // resize force
  // vec_force.resize(totlen);
  // vec_force = Eigen::VectorXcd::Zero(totlen);

  // find element within which the source sits
  double depth = cmt.Depth();
  double rad_source = _mesh.PR() - 1000.0 * depth / _length_norm;
  // std::cout << "\nSource radius: " << rad_source << "\n";
  int idxsource = 0;

  int NQ = _mesh.NN();

  // get the y0-, y0+ values etc at the source location
  double theta_s = (90.0 - cmt.Latitude()) * EIGEN_PI / (180.0);
  double phi_s = cmt.Longitude() * EIGEN_PI / (180.0);
  auto wigdmat =
      GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                       GSHTrans::Single, GSHTrans::ColumnMajor>(_lmax, _lmax, 2,
                                                                theta_s);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
    return ylm;
  };
  // std::cout << "\nSource location: " << theta_s << " " << phi_s << "\n";
  // int m = 0, _l = 2;
  // parameters
  double invsqrt2 = 1.0 / std::sqrt(2.0);
  std::complex<double> isq2 = std::complex<double>(0.0, invsqrt2);
  double omegal2 = std::sqrt((idxl + 2) * (idxl - 1) / 2.0);
  double lprefac = std::exp(-2.0 * 3.141592653589793 * (idxl + 1) / (1 + 0.5));
  // loop through the elements
  //   to find the element that contains the source
  for (int idx = _el; idx < _eu; ++idx) {
    // std::cout << _mesh.ELR(idx) << " " << _mesh.EUR(idx) << " " <<
    // rad_source
    //           << "\n";
    if ((_mesh.ELR(idx) <= rad_source) && (_mesh.EUR(idx) >= rad_source)) {
      // std::cout << idx << " " << rad_source << " " << _mesh.ELR(idx) << " "
      //           << _mesh.EUR(idx) << "\n";
      std::vector<double> vec_nodes(NQ, 0.0);
      for (int idxi = 0; idxi < NQ; ++idxi) {
        vec_nodes[idxi] = _mesh.NodeRadius(idx, idxi);
      }
      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());
      for (int idxn = 0; idxn < totnum; ++idxn) {

        auto w_val = 0.0;
        auto w_deriv = 0.0;
        for (int idxq = 0; idxq < NQ; ++idxq) {
          w_val +=
              vec_all_l[idxl - 1][idx][idxq][idxn] * pleg(idxq, rad_source);
          w_deriv += vec_all_l[idxl - 1][idx][idxq][idxn] *
                     pleg.Derivative(idxq, rad_source);
          if (idxn == 0) {
            // std::cout << "vec_eigval_l: " << "idxl, idx, idxq, idxn: " <<
            // idxl
            //           << " " << idx << " " << idxq << " " << idxn << " "
            //           << vec_eigval_l[idxl - 1][idx][idxq][idxn] << "\n";
          }
        };
        auto w_vald = w_val / rad_source;
        auto w_prefactor = w_deriv - w_vald;
        // std::cout << "w_vald, w_deriv, w_prefactor: " << w_vald << " "
        //           << w_deriv << " " << w_prefactor << "\n";
        // double omegal2 = (idxl + 2) * (idxl - 1) / 2.0;

        // get the index of the spherical harmonic
        for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
          // spherical harmonic
          Complex ymc = std::conj(ylmn(idxl, idxm, -1, phi_s));
          Complex ypc = std::conj(ylmn(idxl, idxm, 1, phi_s));
          Complex ymmc = 0.0, yppc = 0.0;

          if (idxl > 1) {
            ymmc = std::conj(ylmn(idxl, idxm, -2, phi_s));
            yppc = std::conj(ylmn(idxl, idxm, 2, phi_s));
          }

          Complex tmp = w_prefactor * (cmt.MC0m() * ymc - cmt.MC0p() * ypc);
          tmp += w_vald * omegal2 * (cmt.MCmm() * ymmc - cmt.MCpp() * yppc);
          tmp *= isq2 * lprefac;

          vec_lforce(idxn, idxm + idxl) = tmp;
          // };
          // };
        };
      };
    };
  }
  return vec_lforce;
};

}   // namespace Toroidal

#endif