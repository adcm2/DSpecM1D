#ifndef SEM_FULL_GUARD_H
#define SEM_FULL_GUARD_H
#include <iostream>
#include <iomanip>
#include <cmath>
#include <Eigen/Core>
// #include <Spectra/MatOp/SparseGenMatProd.h>
// #include <Spectra/SymGEigsShiftSolver.h>
#include <GSHTrans/Core>
// #include <GaussQuad/All>
#include <Interpolation/Lagrange>
#include "SourceInfo.h"
#include <EarthMesh/All>
#include "input_parser.h"
#include "mesh_model.h"
// #include "MatrixIndices.h"

// using namespace Spectra;

namespace Full1D {

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
  using SMAT = Eigen::SparseMatrix<double>;

  // mesh
  EarthMesh::RadialMesh _mesh;
  MeshModel _mesh_model;

  // integers for calc
  int _lmax, _il;

  // dealing with fluid regions
  bool _has_fluid = false;
  std::vector<int> _vec_fluid;
  std::vector<bool> _vec_dof;

  // matrix parameters
  std::size_t mlen, totlen;
  int _el = 0, _eu = 0, _en = 0, numlen = 0, _k2 = 0, _solint = 0;

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
  Eigen::SparseMatrix<double> mat_seig, mat_ke, mat_inertia, mat_inertia_0,
      mat_ke_0;
  std::vector<Eigen::SparseMatrix<double>> mat_l_ke, mat_l_inertia, mat_l_ke_t,
      mat_l_inertia_t;

  // eigenfunction calculation parameters
  int _num_modes, numaug, _num_modes_l;
  bool _calc_eig = false, _calc_gen = false;
  std::vector<bool> _calc_l;

  // normalisation factors
  double densitynorm = 5515.0;
  double pi_db = 3.141592653589793238462643383279502884197;
  double bigg_nd = 6.6723 * std::pow(10.0, -11.0);
  double bigg_db;
  double frequencynorm = std::sqrt(pi_db * bigg_nd * densitynorm);
  double normint;   // normalisation factor for inertia matrix
  double _freq_norm, _length_norm, _moment_norm;

  // eigenfunction storage
  Eigen::VectorXcd evalues_seig, evalues_gen;
  Eigen::MatrixXcd evectors, evectors_gen;
  std::vector<Eigen::VectorXcd> evalues_gen_l;
  std::vector<Eigen::MatrixXcd> evectors_gen_l;
  std::vector<std::vector<Eigen::Triplet<double>>> tp_ke_s, tp_in_s, tp_ke_t,
      tp_in_t;
  std::vector<std::vector<std::vector<std::size_t>>> vec_ke_s, vec_in_s,
      vec_ke_t, vec_in_t;
  std::vector<std::vector<std::complex<double>>> vec_ke_s_val, vec_in_s_val,
      vec_ke_t_val, vec_in_t_val;
  // stdvvvec vec_augval, vec_augderiv;

  // indexing class
  // MatrixIndices mat_indices;

  // forced problem
  Eigen::VectorXcd vec_force;

  // store for all modes up to f
  int nmodes = 0;
  std::vector<TStore> vec_store;
  std::vector<std::size_t> vec_indices;

  // fsb
  std::vector<int> fsb, vec_offset{0};

public:
  sem() {};
  template <class model1d> sem(const model1d &, double, int, int);
  auto LtG_S(int, int, int) const;
  auto LtG_T(int, int) const;
  auto LtG_R(int, int, int) const;
  auto EL() const { return _el; };
  auto EU() const { return _eu; };

  Eigen::MatrixXcd CalculateForce(SourceInfo::EarthquakeCMT &, int);
  Eigen::MatrixXcd CalculateForce_T(SourceInfo::EarthquakeCMT &, int);
  Eigen::MatrixXcd CalculateForce_R(SourceInfo::EarthquakeCMT &);

  Eigen::MatrixXcd CalculateForce_SPH_U(SourceInfo::EarthquakeCMT &, int);
  Eigen::MatrixXcd CalculateForce_SPH_UP(SourceInfo::EarthquakeCMT &, int);
  Eigen::MatrixXcd CalculateForce_SPH_V(SourceInfo::EarthquakeCMT &, int);
  Eigen::MatrixXcd CalculateForce_SPH_VP(SourceInfo::EarthquakeCMT &, int);

  // cleaner method
  Eigen::MatrixXcd CalculateForce_All(SourceInfo::EarthquakeCMT &, int);
  Eigen::MatrixXcd CalculateForce_Coefficients(SourceInfo::EarthquakeCMT &,
                                               int);
  Eigen::MatrixXcd CalculateForce_All_T(SourceInfo::EarthquakeCMT &, int);
  Eigen::MatrixXcd CalculateForce_Coefficients_T(SourceInfo::EarthquakeCMT &,
                                                 int);

  const EarthMesh::RadialMesh &mesh() const { return _mesh; };
  const MeshModel &mesh_model() const { return _mesh_model; };

  auto Receiver_Elements(InputParameters &) const;

  // Matrices
  SMAT MAT_KE(int idxl) const { return mat_l_ke[idxl - 1]; };
  SMAT MAT_IN(int idxl) const { return mat_l_inertia[idxl - 1]; };

  SMAT MAT_KE_T(int idxl) const { return mat_l_ke_t[idxl - 1]; };
  SMAT MAT_IN_T(int idxl) const { return mat_l_inertia_t[idxl - 1]; };

  SMAT MAT_KE_R() const { return mat_ke_0; };
  SMAT MAT_IN_R() const { return mat_inertia_0; };

  // tripletlists:
  auto tripletlist_ke_s(int idxl) const { return tp_ke_s[idxl - 1]; };
  auto tripletlist_in_s(int idxl) const { return tp_in_s[idxl - 1]; };
  auto tripletlist_ke_t(int idxl) const { return tp_ke_t[idxl - 1]; };
  auto tripletlist_in_t(int idxl) const { return tp_in_t[idxl - 1]; };

  // vectors of values
  auto tripletlist_ke_s_idx(int idxl) const { return vec_ke_s[idxl - 1]; };
  auto tripletlist_in_s_idx(int idxl) const { return vec_in_s[idxl - 1]; };
  auto tripletlist_ke_t_idx(int idxl) const { return vec_ke_t[idxl - 1]; };
  auto tripletlist_in_t_idx(int idxl) const { return vec_in_t[idxl - 1]; };
  auto tripletlist_ke_s_val(int idxl) const { return vec_ke_s_val[idxl - 1]; };
  auto tripletlist_in_s_val(int idxl) const { return vec_in_s_val[idxl - 1]; };
  auto tripletlist_ke_t_val(int idxl) const { return vec_ke_t_val[idxl - 1]; };
  auto tripletlist_in_t_val(int idxl) const { return vec_in_t_val[idxl - 1]; };

  // Receiver vectors
  Eigen::MatrixXcd RV_Z(double, double, int);
  Eigen::MatrixXcd RV_THETA(double, double, int);
  Eigen::MatrixXcd RV_PHI(double, double, int);

  // spheroidal receiver vectors
  Eigen::MatrixXcd RV_Z(InputParameters &, int, int);
  Eigen::MatrixXcd RV_BASE_Z(InputParameters &, int, int);
  Eigen::MatrixXcd RV_VAL_Z(InputParameters &, int, int);
  // Eigen::MatrixXcd RV_BASE_Z_T(InputParameters &, int, int);
  // Eigen::MatrixXcd RV_VAL_Z_T(InputParameters &, int, int);
  Eigen::MatrixXcd RV_THETA(InputParameters &, int, int);
  Eigen::MatrixXcd RV_BASE_THETA(InputParameters &, int, int);
  Eigen::MatrixXcd RV_VAL_THETA(InputParameters &, int, int);
  Eigen::MatrixXcd RV_BASE_THETA_T(InputParameters &, int, int);
  Eigen::MatrixXcd RV_VAL_THETA_T(InputParameters &, int, int);
  Eigen::MatrixXcd RV_PHI(InputParameters &, int, int);
  // Eigen::MatrixXcd RV_BASE_PHI(InputParameters &, int, int);
  Eigen::MatrixXcd RV_VAL_PHI(InputParameters &, int, int);
  Eigen::MatrixXcd RV_BASE_PHI_T(InputParameters &, int, int);
  Eigen::MatrixXcd RV_VAL_PHI_T(InputParameters &, int, int);

  // toroidal receiver vectors
  Eigen::MatrixXcd RV_THETA_T(InputParameters &, int, int);
  Eigen::MatrixXcd RV_PHI_T(InputParameters &, int, int);

  // radial receiver vector
  Eigen::MatrixXcd RV_Z_R(InputParameters &, int);

  // radial receiver vector
  // Eigen::MatrixXcd RV_Z_R(InputParameters &, int, int);

  // void CalculateEigenfrequencies(int, double);
  // void CalculateEigenfrequenciesSeparate(int, double);
  // auto efrequencies_gen() const;
  // auto efrequencies_gen_separate(int) const;
  // auto efunctions(int) const;
  // auto efunctions_ref(int idxl) { return evectors_gen_l[idxl - 1]; };
  // auto el() const { return _el; };
  // auto eu() const { return _eu; };

  // auto modes_coupled() const { return vec_store; };

  // Eigen::VectorXcd CalculateForce(SourceInfo::EarthquakeCMT &);
  // Eigen::VectorXcd ReceiverVectorTheta(double, double, double);
  // Eigen::VectorXcd ReceiverVectorThetaSurface(double, double, double);
  // Eigen::VectorXcd ReceiverVectorPhi(double, double, double);
  // Eigen::VectorXcd ReceiverVectorPhiSurface(double, double, double);
  // Eigen::SparseMatrix<double> GetStiffnessMatrix() const { return mat_ke; };
  // Eigen::SparseMatrix<double> GetInertiaMatrix() const { return mat_inertia;
  // };

  // Eigen::SparseMatrix<double> GetStiffnessMatrixL(int idxl) const {
  //   return mat_l_ke[idxl - 1];
  // };
  // Eigen::SparseMatrix<double> GetInertiaMatrixL(int idxl) const {
  //   return mat_l_inertia[idxl - 1];
  // };
  // Eigen::MatrixXcd CalculateForce(SourceInfo::EarthquakeCMT &, int);
  // // Eigen::MatrixXcd ReceiverVectorTheta(double, double, double);
  // Eigen::MatrixXcd ReceiverVectorThetaSurfaceL(double, double, double, int);
  // // Eigen::VectorXcd ReceiverVectorPhi(double, double, double);
  // Eigen::MatrixXcd ReceiverVectorPhiSurfaceL(double, double, double, int);

  // auto PrintModesUpToFreq(double) const;

  // // function that gets the matrices for modes up to a certain frequency
  // auto FindModesForCoupling(double);
  // template <class model1d>
  // Eigen::MatrixXcd NMC_INERTIA(const model1d &, bool = false);
  // template <class model1d>
  // Eigen::MatrixXcd NMC_KE(const model1d &, bool = false);
  // auto NMC_FORCE(SourceInfo::EarthquakeCMT &, bool = false);
  // Eigen::VectorXcd ReceiverVectorThetaSurfaceCoupling(double, double, double,
  //                                                     bool = false);
  // Eigen::VectorXcd ReceiverVectorPhiSurfaceCoupling(double, double, double,
  //                                                   bool = false);

  // // augmentation
  // void augment_basis_calculate();

  // // nmc for particular l
  // template <class model1d>
  // Eigen::MatrixXcd NMC_INERTIA(const model1d &, int, bool = false);
  // template <class model1d>
  // Eigen::MatrixXcd NMC_KE(const model1d &, int, bool = false);
  // Eigen::MatrixXcd ReceiverVectorThetaSurface_NMCL(double, double, double,
  // int,
  //                                                  bool = false);
  // Eigen::MatrixXcd ReceiverVectorPhiSurface_NMCL(double, double, double, int,
  //                                                bool = false);
  // Eigen::MatrixXcd CalculateForceNMC(SourceInfo::EarthquakeCMT &, int,
  //                                    bool = false);
  // auto NumModesL(int idxl) const { return _num_modes_l; };
};

template <class model1d>
sem::sem(const model1d &inp_model, double maxstep, int NQ, int lmax)
    : _mesh(inp_model, NQ, 1.0, maxstep, false),
      _freq_norm{1.0 / inp_model.TimeNorm()}, _lmax{lmax},
      _k2{lmax * (lmax + 1)},
      normint{1.0 / (inp_model.TimeNorm() * frequencynorm) *
              sqrt(inp_model.DensityNorm() / densitynorm)},
      _length_norm{inp_model.LengthNorm()},
      bigg_db{6.67230 * std::pow(10.0, -11.0) /
              inp_model.GravitationalConstant()},
      _moment_norm{
          inp_model.MassNorm() *
          std::pow(inp_model.LengthNorm() / inp_model.TimeNorm(), 2.0)} {
  std::cout << "Moment norm: " << _moment_norm << "\n";
  _mesh_model = MeshModel(_mesh, inp_model);
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  // getting boundary and offset information for the local to global map
  fsb = _mesh.FS_Boundaries();
  {
    std::size_t totnum = 0;
    for (int idx = 0; idx < _mesh.NE(); ++idx) {
      auto tmp = vec_offset[idx];
      if (idx == fsb[totnum]) {
        tmp += 1;
        totnum += 1;
      }
      vec_offset.push_back(tmp);
    }
  }

  // numaug
  // if (_el == 0) {
  //   numaug = _mesh.LayerNumber(_eu - 1) - _mesh.LayerNumber(_el) + 1;
  // } else {
  numaug = _mesh.NL();
  // }

  // matrix indices class
  // mat_indices = MatrixIndices(_el, _eu, 1, _lmax, _mesh.NN());
  // Spheroidal modes have lmin = 1

  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  // std::cout << "Calculating basis function derivatives...\n";
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

  // std::cout << "Got derivative matrices.\n";
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////

  // get gravitational acceleration in 1D model
  // this will be updated to use gplspec but for moment using simple integral
  // method
  stdvvec vec_allradii_g(_mesh.NE(), stdvec(NQ, 0.0));
  // double bigg_db = 6.6743015 * std::pow(10.0, -11.0);
  // bigg_db *= 1.0 / inp_model.GravitationalConstant();

  // std::cout << "Calculating gravitational acceleration profile...\n";
  // integrate
  for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
    int laynum = _mesh.LayerNumber(idxe);
    if (idxe != 0) {
      vec_allradii_g[idxe][0] = vec_allradii_g[idxe - 1][NQ - 1];
    }
    int idxlow = 1;
    for (int idxn = idxlow; idxn < NQ; ++idxn) {
      double urad = _mesh.NodeRadius(idxe, idxn);
      double lrad = _mesh.NodeRadius(idxe, idxn - 1);
      double tmp = 0.0;
      double ewidth2 = 0.5 * (urad - lrad);
      double crad2 = 0.5 * (urad + lrad);
      for (int idxi = 0; idxi < NQ; ++idxi) {
        double cradi = ewidth2 * q.X(idxi) + crad2;
        tmp += q.W(idxi) * _mesh_model.Density(idxe, idxi) * cradi * cradi;
      }
      tmp *= ewidth2;
      vec_allradii_g[idxe][idxn] = vec_allradii_g[idxe][idxn - 1] + tmp;
    }
  }

  // divide by r^2:
  for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
    int idxlow = (idxe == 0);
    for (int idxn = idxlow; idxn < NQ; ++idxn) {
      double crad = _mesh.NodeRadius(idxe, idxn);
      vec_allradii_g[idxe][idxn] *= 4.0 * pi_db * bigg_db / (crad * crad);
    }
  }

  // std::cout << "Got gravitational acceleration profile.\n";
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////

  {

    // size of matrix for full problem
    totlen = this->LtG_S(2, _mesh.NE() - 1, NQ - 1) + 1;

    // std::cout << "Total matrix length: " << totlen << "\n";

    // matrices
    mat_seig.resize(totlen, totlen);
    mat_ke.resize(totlen, totlen);
    mat_inertia.resize(totlen, totlen);

    // triplet list for matrices
    using T = Eigen::Triplet<double>;

    // make the vector of matrices for the different ls
    {
      using T = Eigen::Triplet<double>;
      std::size_t matlen = totlen;

      mat_l_ke_t.reserve(_lmax);
      mat_l_inertia_t.reserve(_lmax);
      ///////////////////////////////////////////////////////////////////////

      ///////////////////////////////////////////////////////////////////////
      // toroidals
      // vector to store whether an element is fluid or not
      _vec_fluid = std::vector<int>(_mesh.NE(), 0);
      _vec_dof = std::vector<bool>(_mesh.NE(), false);
      bool _has_fluid = false;
      for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
        if (inp_model.IsFluid(_mesh.LayerNumber(idxe))) {
          _vec_fluid[idxe] = 1;
          _has_fluid = true;
        }

        // if there is a change from fluid to solid or vice versa, we set the
        // dof to true at the lower element
        if (idxe > 0) {
          if ((std::abs(_vec_fluid[idxe] - _vec_fluid[idxe - 1]) == 1)) {
            _vec_dof[idxe - 1] = true;
          }
        }
      }
      // std::cout << "Has fluid: " << _has_fluid << "\n";

      // int _el = 0;
      {
        bool found = false;
        bool prev_fluid = false;
        while (!found) {
          ++_el;
          // std::cout << "Checking element: " << _el << "\n";
          auto laynum = _mesh.LayerNumber(_el);
          // std::cout << "Checking element: " << _el << " layer: " << laynum
          //           << "\n";
          if (prev_fluid) {
            if (inp_model.IsSolid(laynum)) {
              found = true;
              break;
            }
          }
          if (inp_model.IsFluid(laynum)) {
            prev_fluid = true;
          }
        }
      }

      // std::cout << "First solid element index: " << _el << "\n";
      // check if there is an ocean layer
      _eu = _el;

      {
        bool found = false;
        while (!found) {
          auto laynum = _mesh.LayerNumber(_eu);
          if (inp_model.IsFluid(laynum) || (++_eu == _mesh.NE())) {
            found = true;
            break;
          }
        }
      }
      // std::cout << "First solid element index: " << _el << "\n";
      // std::cout << "Upper element index: " << _eu << "\n";
      // std::cout << "Number of elements: " << _mesh.NE() << "\n";

      auto toroidal_l_g = [NQ, _el = this->_el](int idxe, int idxq) {
        return (idxe - _el) * (NQ - 1) + idxq;
      };
      auto num_toroidal_dof = toroidal_l_g(_eu - 1, NQ - 1) + 1;

      // std::cout << "Before looping through l\n";
      for (int idxl = 1; idxl < _lmax + 1; ++idxl) {
        auto k2 = idxl * (idxl + 1.0);
        std::vector<T> tpl_in, tpl_ke;
        std::vector<std::vector<std::size_t>> vec_in_idx, vec_ke_idx;
        std::vector<Complex> vec_in_val, vec_ke_val;
        {
          for (int idxe = _el; idxe < _eu; ++idxe) {
            // int imin = (idxe == 0);
            double elem_width = _mesh.EW(idxe);
            int laynum = _mesh.LayerNumber(idxe);
            for (int i = 0; i < q.N(); ++i) {
              // pre-factor
              double xrad = _mesh.NodeRadius(idxe, i);
              double tmp = elem_width / 2.0 * q.W(i) *
                           _mesh_model.Density(idxe, i) * xrad * xrad;

              // indices
              auto idx_ww = this->LtG_T(idxe, i);

              // U'U term
              tpl_in.push_back(T(idx_ww, idx_ww, tmp));
              vec_in_idx.push_back({idx_ww, idx_ww});
              vec_in_val.push_back(tmp);
            }
          }
          vec_in_t.push_back(vec_in_idx);
          vec_in_t_val.push_back(vec_in_val);
          tp_in_t.push_back(tpl_in);
          Eigen::SparseMatrix<double> mat_tmp(num_toroidal_dof,
                                              num_toroidal_dof);
          mat_tmp.setFromTriplets(tpl_in.begin(), tpl_in.end());
          mat_tmp.makeCompressed();
          mat_l_inertia_t.push_back(mat_tmp);
          mat_l_inertia_t.back().makeCompressed();
        }
        // std::cout << "idxl: " << idxl << ", post inertia matrix setting.\n";
        {

          for (int idxe = _el; idxe < _eu; ++idxe) {
            // int imin = (idxe == 0);
            stdvvec mat_d(q.N(), stdvec(q.N(), 0.0));
            // stdvvec mat_n = mat_m;
            double elem_width = _mesh.EW(idxe);
            int laynum = _mesh.LayerNumber(idxe);
            double d_val = 2.0 / elem_width;

            for (int i = 0; i < q.N(); ++i) {
              for (int j = 0; j < q.N(); ++j) {
                mat_d[i][j] = d_val * vec_lag_deriv[j][i];
              }
            }

            for (int i = 0; i < q.N(); ++i) {
              // auto ri = inp_model.NodeRadius(idxe, i);
              double ri = _mesh.NodeRadius(idxe, i);
              auto Li = _mesh_model.L(idxe, i);
              auto Ni = _mesh_model.N(idxe, i);
              double tmp0 = elem_width / 2.0 * q.W(i);
              auto tmp_ww = tmp0 * (Li + (k2 - 2.0) * Ni);

              // indices
              auto idx_ww = this->LtG_T(idxe, i);
              // auto idx_vv = this->LtG_S(1, idxe, i);

              // U'U term
              tpl_ke.push_back(T(idx_ww, idx_ww, tmp_ww));
              vec_ke_idx.push_back({idx_ww, idx_ww});
              vec_ke_val.push_back(tmp_ww);
            }

            for (int i = 0; i < q.N(); ++i) {
              auto ri = _mesh.NodeRadius(idxe, i);
              auto Li = _mesh_model.L(idxe, i);
              for (int j = 0; j < q.N(); ++j) {
                auto rj = _mesh.NodeRadius(idxe, j);
                auto Lj = _mesh_model.L(idxe, j);

                auto idx_i = this->LtG_T(idxe, i);
                auto idx_j = this->LtG_T(idxe, j);

                // finding values to push back
                auto tmp_wdw =
                    -0.5 * elem_width * q.W(j) * Lj * rj * mat_d[i][j];

                // push back
                tpl_ke.push_back(T(idx_i, idx_j, tmp_wdw));
                tpl_ke.push_back(T(idx_j, idx_i, tmp_wdw));
                vec_ke_idx.push_back({idx_i, idx_j});
                vec_ke_val.push_back(tmp_wdw);
                vec_ke_idx.push_back({idx_j, idx_i});
                vec_ke_val.push_back(tmp_wdw);
              }
            }

            for (int i = 0; i < q.N(); ++i) {
              for (int j = 0; j < q.N(); ++j) {
                auto tmp_ww = 0.0;

                // go over element
                for (int k = 0; k < q.N(); ++k) {
                  auto rk = _mesh.NodeRadius(idxe, k);
                  auto tmp_ddrr = rk * rk * mat_d[i][k] * mat_d[j][k];
                  tmp_ww += q.W(k) * _mesh_model.L(idxe, k) * tmp_ddrr;
                }

                // multiply
                tmp_ww *= elem_width / 2.0;

                auto idx_i = this->LtG_T(idxe, i);
                auto idx_j = this->LtG_T(idxe, j);
                tpl_ke.push_back(T(idx_i, idx_j, tmp_ww));
                vec_ke_idx.push_back({idx_i, idx_j});
                vec_ke_val.push_back(tmp_ww);
              }
            }
          }
          tp_ke_t.push_back(tpl_ke);
          vec_ke_t.push_back(vec_ke_idx);
          vec_ke_t_val.push_back(vec_ke_val);
          Eigen::SparseMatrix<double> mat_tmp(num_toroidal_dof,
                                              num_toroidal_dof);
          mat_tmp.setFromTriplets(tpl_ke.begin(), tpl_ke.end());
          mat_tmp.makeCompressed();
          mat_l_ke_t.push_back(mat_tmp);
          mat_l_ke_t.back().makeCompressed();
        }
      }
      // for (auto idx : _vec_fluid) {
      //   std::cout << idx << "\n";
      // }

      ///////////////////////////////////////////////////////////////////////
      // spheroidals
      mat_l_ke.reserve(_lmax);
      mat_l_inertia.reserve(_lmax);
      for (int idxl = 1; idxl < _lmax + 1; ++idxl) {
        // auto k = std::sqrt(static_cast<double>(idxl) *
        //                    (static_cast<double>(idxl) + 1.0));
        std::vector<std::vector<std::size_t>> vec_in_idx, vec_ke_idx;
        std::vector<Complex> vec_in_val, vec_ke_val;
        auto k2 = idxl * (idxl + 1.0);
        std::vector<T> tpl_in, tpl_ke;
        {
          for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
            // int imin = (idxe == 0);
            double elem_width = _mesh.EW(idxe);
            int laynum = _mesh.LayerNumber(idxe);
            for (int i = 0; i < q.N(); ++i) {
              // pre-factor
              double xrad = _mesh.NodeRadius(idxe, i);
              double tmp = elem_width / 2.0 * q.W(i) *
                           _mesh_model.Density(idxe, i) * xrad * xrad;

              // indices
              auto idx_uu = this->LtG_S(0, idxe, i);
              auto idx_vv = this->LtG_S(1, idxe, i);

              // U'U term
              tpl_in.push_back(T(idx_uu, idx_uu, tmp));
              vec_in_idx.push_back({idx_uu, idx_uu});
              vec_in_val.push_back(tmp);

              // V'V term
              tpl_in.push_back(T(idx_vv, idx_vv, tmp * k2));
              vec_in_idx.push_back({idx_vv, idx_vv});
              vec_in_val.push_back(tmp * k2);
            }
          }
          vec_in_s.push_back(vec_in_idx);
          vec_in_s_val.push_back(vec_in_val);
          tp_in_s.push_back(tpl_in);
          Eigen::SparseMatrix<double> mat_tmp(totlen, totlen);
          mat_tmp.setFromTriplets(tpl_in.begin(), tpl_in.end());
          mat_tmp.makeCompressed();
          mat_l_inertia.push_back(mat_tmp);
          mat_l_inertia.back().makeCompressed();
        }

        //  kinetic energy matrix
        {

          // diagonal terms
          {
            // do PP at 00
            // {
            //   int idxpb = this->LtG_S(2, 0, 0);
            //   double tmp0 = _mesh.EW(0) / 2.0 * q.W(0);
            //   tpl_ke.push_back(
            //       T(idxpb, idxpb, k2 / (4.0 * pi_db * bigg_db) * tmp0));
            // }
            for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
              // int imin = (idxe == 0);
              double elem_width = _mesh.EW(idxe);

              int laynum = _mesh.LayerNumber(idxe);
              double d_val = 2.0 / elem_width;
              for (int i = 0; i < q.N(); ++i) {
                auto crad = _mesh.NodeRadius(idxe, i);
                auto crho = _mesh_model.Density(idxe, i);
                auto gi = vec_allradii_g[idxe][i];
                // universal factors
                double tmp0 = elem_width / 2.0 * q.W(i);
                double tmpp = tmp0 * k2 * crho * crad;

                // material parameters
                auto Li = _mesh_model.L(idxe, i);
                auto Ai = _mesh_model.A(idxe, i);
                auto Ni = _mesh_model.N(idxe, i);
                auto Fi = _mesh_model.F(idxe, i);

                //////////////////////////////
                // U'U term
                double tmp_u =
                    tmp0 *
                    (4.0 * crho * (pi_db * bigg_db * crho * crad - gi) * crad +
                     k2 * Li + 4 * (Ai - Ni));

                // V'V term
                double tmp_v = tmp0 * k2 * (k2 * Ai - 2 * Ni + Li);

                // P'P term
                double tmp_pp = k2 / (4.0 * pi_db * bigg_db) * tmp0;

                // P'V term
                double tmp_pv = tmp0 * k2 * crho * crad;

                // U'V term
                double tmp_uv =
                    tmp0 * k2 * (crho * gi * crad - Li - 2 * (Ai - Ni));

                //////////////////////////////

                // indices
                std::size_t idxtiu = this->LtG_S(0, idxe, i);
                std::size_t idxtiv = this->LtG_S(1, idxe, i);
                std::size_t idxtip = this->LtG_S(2, idxe, i);

                //////////////////////////////
                // pushback
                tpl_ke.push_back(T(idxtiu, idxtiu, tmp_u));    // U'U term
                tpl_ke.push_back(T(idxtiv, idxtiv, tmp_v));    // V'V term
                tpl_ke.push_back(T(idxtip, idxtip, tmp_pp));   // P'P term
                tpl_ke.push_back(T(idxtip, idxtiv, tmp_pv));   // P'V term
                tpl_ke.push_back(T(idxtiv, idxtip, tmp_pv));   // V'P term
                tpl_ke.push_back(T(idxtiu, idxtiv, tmp_uv));   // U'V term
                tpl_ke.push_back(T(idxtiv, idxtiu, tmp_uv));   // V'U term
                //////////////////////////////
                vec_ke_idx.push_back({idxtiu, idxtiu});
                vec_ke_val.push_back(tmp_u);
                vec_ke_idx.push_back({idxtiv, idxtiv});
                vec_ke_val.push_back(tmp_v);
                vec_ke_idx.push_back({idxtip, idxtip});
                vec_ke_val.push_back(tmp_pp);
                vec_ke_idx.push_back({idxtip, idxtiv});
                vec_ke_val.push_back(tmp_pv);
                vec_ke_idx.push_back({idxtiv, idxtip});
                vec_ke_val.push_back(tmp_pv);
                vec_ke_idx.push_back({idxtiu, idxtiv});
                vec_ke_val.push_back(tmp_uv);
                vec_ke_idx.push_back({idxtiv, idxtiu});
                vec_ke_val.push_back(tmp_uv);
              }
            }

            // add in boundary term for gravity
            std::size_t idxpb = this->LtG_S(2, _mesh.NE() - 1, NQ - 1);
            double rpb = _mesh.NodeRadius(_mesh.NE() - 1, NQ - 1);
            tpl_ke.push_back(
                T(idxpb, idxpb, rpb * (idxl + 1) / (4.0 * pi_db * bigg_db)));
            vec_ke_idx.push_back({idxpb, idxpb});
            vec_ke_val.push_back(rpb * (idxl + 1) / (4.0 * pi_db * bigg_db));
          }

          // coupling terms
          {
            for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
              // int imin = (idxe == 0);
              double elem_width = _mesh.EW(idxe);
              auto e2 = 0.5 * _mesh.EW(idxe);
              int laynum = _mesh.LayerNumber(idxe);
              double d_val = 2.0 / _mesh.EW(idxe);

              // calculate M and N
              stdvvec mat_d(q.N(), stdvec(q.N(), 0.0));
              // stdvvec mat_n = mat_m;

              for (int i = 0; i < q.N(); ++i) {
                for (int j = 0; j < q.N(); ++j) {
                  mat_d[i][j] = d_val * vec_lag_deriv[j][i];
                }
              }

              for (int i = 0; i < q.N(); ++i) {
                auto ri = _mesh.NodeRadius(idxe, i);
                auto rhoi = _mesh_model.Density(idxe, i);
                auto Ai = _mesh_model.A(idxe, i);
                auto Ni = _mesh_model.N(idxe, i);
                auto Fi = _mesh_model.F(idxe, i);
                auto Li = _mesh_model.L(idxe, i);
                for (int j = 0; j < q.N(); ++j) {
                  auto rj = _mesh.NodeRadius(idxe, j);
                  auto rhoj = _mesh_model.Density(idxe, j);
                  auto Aj = _mesh_model.A(idxe, j);
                  auto Nj = _mesh_model.N(idxe, j);
                  auto Fj = _mesh_model.F(idxe, j);
                  auto Lj = _mesh_model.L(idxe, j);

                  // finding values to push back
                  auto tmp_uu = elem_width * q.W(j) * Fj * rj * mat_d[i][j];

                  auto tmp_vv = -e2 * k2 * q.W(j) * Lj * rj * mat_d[i][j];

                  auto tmp_uvd = e2 * k2 * q.W(i) * Li * ri * mat_d[j][i];

                  auto tmp_udv = -e2 * k2 * q.W(j) * Fj * rj * mat_d[i][j];

                  auto tmp_pdu = q.W(j) * rhoj * rj * rj * vec_lag_deriv[j][i];

                  // indices
                  auto idx_u_i = this->LtG_S(0, idxe, i);
                  auto idx_u_j = this->LtG_S(0, idxe, j);
                  auto idx_v_i = this->LtG_S(1, idxe, i);
                  auto idx_v_j = this->LtG_S(1, idxe, j);
                  auto idx_p_i = this->LtG_S(2, idxe, i);
                  auto idx_p_j = this->LtG_S(2, idxe, j);

                  // if (((idxe == 0) && ((i == 4) && (j == 4))) ||
                  //     ((idxe == 1) && ((i == 0) && (j == 0)))) {
                  //   std::cout << "i: " << i << " j: " << j
                  //             << " idx_p_i: " << idx_p_i
                  //             << " idx_u_j: " << idx_u_j << "\n";
                  //   std::cout << "rhoj: " << rhoj << " rj: " << rj
                  //             << " vec_lag_deriv[j][i]: " <<
                  //             vec_lag_deriv[j][i]
                  //             << "\n";
                  //   std::cout << "pdu: " << std::setprecision(16) <<
                  //   tmp_pdu
                  //             << "\n\n";
                  // }
                  // pushbacks
                  tpl_ke.push_back(T(idx_u_i, idx_u_j, tmp_uu));    // U'U
                  tpl_ke.push_back(T(idx_u_j, idx_u_i, tmp_uu));    // UU'
                  tpl_ke.push_back(T(idx_v_i, idx_v_j, tmp_vv));    // V'V
                  tpl_ke.push_back(T(idx_v_j, idx_v_i, tmp_vv));    // VV'
                  tpl_ke.push_back(T(idx_u_i, idx_v_j, tmp_uvd));   // U'V
                  tpl_ke.push_back(T(idx_v_j, idx_u_i, tmp_uvd));   // VU'
                  tpl_ke.push_back(T(idx_u_i, idx_v_j, tmp_udv));   // UV'
                  tpl_ke.push_back(T(idx_v_j, idx_u_i, tmp_udv));   // V'U
                  tpl_ke.push_back(T(idx_p_i, idx_u_j, tmp_pdu));   // P'U
                  tpl_ke.push_back(T(idx_u_j, idx_p_i, tmp_pdu));   // UP'

                  vec_ke_idx.push_back({idx_u_j, idx_u_j});
                  vec_ke_val.push_back(tmp_uu);
                  vec_ke_idx.push_back({idx_u_j, idx_u_i});
                  vec_ke_val.push_back(tmp_uu);
                  vec_ke_idx.push_back({idx_v_i, idx_v_j});
                  vec_ke_val.push_back(tmp_vv);
                  vec_ke_idx.push_back({idx_v_j, idx_v_i});
                  vec_ke_val.push_back(tmp_vv);
                  vec_ke_idx.push_back({idx_u_i, idx_v_j});
                  vec_ke_val.push_back(tmp_uvd);
                  vec_ke_idx.push_back({idx_v_j, idx_u_i});
                  vec_ke_val.push_back(tmp_uvd);
                  vec_ke_idx.push_back({idx_u_i, idx_v_j});
                  vec_ke_val.push_back(tmp_udv);
                  vec_ke_idx.push_back({idx_v_j, idx_u_i});
                  vec_ke_val.push_back(tmp_udv);
                  vec_ke_idx.push_back({idx_p_i, idx_u_j});
                  vec_ke_val.push_back(tmp_pdu);
                  vec_ke_idx.push_back({idx_u_j, idx_p_i});
                  vec_ke_val.push_back(tmp_pdu);
                }
              }

              for (int i = 0; i < q.N(); ++i) {
                for (int j = 0; j < q.N(); ++j) {
                  auto tmp_uu = 0.0;
                  auto tmp_vv = 0.0;
                  auto tmp_pp = 0.0;

                  // go over element
                  for (int k = 0; k < q.N(); ++k) {
                    auto rk = _mesh.NodeRadius(idxe, k);
                    auto tmp_ddrr = rk * rk * mat_d[i][k] * mat_d[j][k];
                    tmp_uu += q.W(k) * _mesh_model.C(idxe, k) * tmp_ddrr;

                    tmp_vv += q.W(k) * _mesh_model.L(idxe, k) * tmp_ddrr;

                    tmp_pp += q.W(k) * tmp_ddrr;
                  }

                  // multiply
                  tmp_uu *= e2;
                  tmp_vv *= e2 * k2;
                  tmp_pp *= e2 / (4.0 * pi_db * bigg_db);

                  // indices
                  auto idx_u_i = this->LtG_S(0, idxe, i);
                  auto idx_u_j = this->LtG_S(0, idxe, j);
                  auto idx_v_i = this->LtG_S(1, idxe, i);
                  auto idx_v_j = this->LtG_S(1, idxe, j);
                  auto idx_p_i = this->LtG_S(2, idxe, i);
                  auto idx_p_j = this->LtG_S(2, idxe, j);

                  // pushbacks
                  tpl_ke.push_back(T(idx_u_i, idx_u_j, tmp_uu));   // U'U
                  tpl_ke.push_back(T(idx_v_i, idx_v_j, tmp_vv));   // V'V
                  tpl_ke.push_back(T(idx_p_i, idx_p_j, tmp_pp));   // P'P
                  vec_ke_idx.push_back({idx_u_i, idx_u_j});
                  vec_ke_val.push_back(tmp_uu);
                  vec_ke_idx.push_back({idx_v_i, idx_v_j});
                  vec_ke_val.push_back(tmp_vv);
                  vec_ke_idx.push_back({idx_p_i, idx_p_j});
                  vec_ke_val.push_back(tmp_pp);
                }
              }
            }
          }
          vec_ke_s.push_back(vec_ke_idx);
          vec_ke_s_val.push_back(vec_ke_val);
          tp_ke_s.push_back(tpl_ke);
          {
            Eigen::SparseMatrix<double> mat_tmp(totlen, totlen);
            mat_tmp.setFromTriplets(tpl_ke.begin(), tpl_ke.end());
            mat_tmp.makeCompressed();
            mat_l_ke.push_back(mat_tmp);
            mat_l_ke.back().makeCompressed();
          }

          // std::cout << "Completed L = " << idxl << " matrices\n";
        }
      }
      ///////////////////////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////////////////
      // radial matrices

      // for (int idxl = 0; idxl < _lmax + 1; ++idxl) {
      //   // auto k = std::sqrt(static_cast<double>(idxl) *
      //   //                    (static_cast<double>(idxl) + 1.0));
      //   auto k2 = idxl * (idxl + 1.0);
      // std::cout << "Calculating radial matrices...\n";
      {
        auto num_radial_dof = this->LtG_R(1, _mesh.NE() - 1, NQ - 1) + 1;
        std::vector<T> tpl_in, tpl_ke;
        {
          for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
            // int imin = (idxe == 0);
            double elem_width = _mesh.EW(idxe);
            int laynum = _mesh.LayerNumber(idxe);
            for (int i = 0; i < q.N(); ++i) {
              // pre-factor
              double xrad = _mesh.NodeRadius(idxe, i);
              double tmp = elem_width / 2.0 * q.W(i) *
                           _mesh_model.Density(idxe, i) * xrad * xrad;

              // indices
              auto idx_uu = this->LtG_R(0, idxe, i);

              // U'U term
              tpl_in.push_back(T(idx_uu, idx_uu, tmp));
            }
          }
          // Eigen::SparseMatrix<double> mat_tmp(totlen, totlen);
          mat_inertia_0.resize(num_radial_dof, num_radial_dof);
          mat_inertia_0.setFromTriplets(tpl_in.begin(), tpl_in.end());
          mat_inertia_0.makeCompressed();
        }
        // std::cout << "Post radial inertia matrix setting.\n";

        //  kinetic energy matrix
        {

          // diagonal terms
          {
            // do PP at 00
            // {
            //   int idxpb = this->LtG_S(2, 0, 0);
            //   double tmp0 = _mesh.EW(0) / 2.0 * q.W(0);
            //   tpl_ke.push_back(
            //       T(idxpb, idxpb, k2 / (4.0 * pi_db * bigg_db) * tmp0));
            // }
            for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
              // int imin = (idxe == 0);
              double elem_width = _mesh.EW(idxe);

              int laynum = _mesh.LayerNumber(idxe);
              double d_val = 2.0 / elem_width;
              for (int i = 0; i < q.N(); ++i) {
                auto crad = _mesh.NodeRadius(idxe, i);
                auto crho = _mesh_model.Density(idxe, i);
                auto gi = vec_allradii_g[idxe][i];
                // universal factors
                double tmp0 = elem_width / 2.0 * q.W(i);

                // material parameters
                auto Li = _mesh_model.L(idxe, i);
                auto Fi = _mesh_model.F(idxe, i);
                auto Ai = _mesh_model.A(idxe, i);
                auto Ni = _mesh_model.N(idxe, i);

                //////////////////////////////
                // U'U term
                double tmp_u =
                    tmp0 *
                    (4.0 * crho * (pi_db * bigg_db * crho * crad - gi) * crad +
                     4 * (Ai - Ni));

                //////////////////////////////

                // indices
                std::size_t idxtiu = this->LtG_R(0, idxe, i);

                //////////////////////////////
                // pushback
                tpl_ke.push_back(T(idxtiu, idxtiu, tmp_u));   // U'U term
                //////////////////////////////
              }
            }

            // add in boundary term for gravity
            int idxpb = this->LtG_R(1, _mesh.NE() - 1, NQ - 1);
            double rpb = _mesh.NodeRadius(_mesh.NE() - 1, NQ - 1);
            tpl_ke.push_back(T(idxpb, idxpb, rpb / (4.0 * pi_db * bigg_db)));
          }
          // std::cout << "Post radial KE diagonal terms.\n";
          // coupling terms
          {
            for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
              // int imin = (idxe == 0);
              double elem_width = _mesh.EW(idxe);
              auto e2 = 0.5 * _mesh.EW(idxe);
              int laynum = _mesh.LayerNumber(idxe);
              double d_val = 2.0 / _mesh.EW(idxe);

              // calculate M and N
              stdvvec mat_d(q.N(), stdvec(q.N(), 0.0));
              // stdvvec mat_n = mat_m;

              for (int i = 0; i < q.N(); ++i) {
                for (int j = 0; j < q.N(); ++j) {
                  mat_d[i][j] = d_val * vec_lag_deriv[j][i];
                }
              }

              for (int i = 0; i < q.N(); ++i) {
                auto ri = _mesh.NodeRadius(idxe, i);
                auto rhoi = _mesh_model.Density(idxe, i);
                auto Ai = _mesh_model.A(idxe, i);
                auto Ni = _mesh_model.N(idxe, i);
                auto Fi = _mesh_model.F(idxe, i);
                auto Li = _mesh_model.L(idxe, i);
                for (int j = 0; j < q.N(); ++j) {
                  auto rj = _mesh.NodeRadius(idxe, j);
                  auto rhoj = _mesh_model.Density(idxe, j);
                  auto Aj = _mesh_model.A(idxe, j);
                  auto Nj = _mesh_model.N(idxe, j);
                  auto Fj = _mesh_model.F(idxe, j);
                  auto Lj = _mesh_model.L(idxe, j);

                  // finding values to push back
                  auto tmp_uu = elem_width * q.W(j) * Fj * rj * mat_d[i][j];

                  auto tmp_pdu = q.W(j) * rhoj * rj * rj * vec_lag_deriv[j][i];

                  // indices
                  auto idx_u_i = this->LtG_R(0, idxe, i);
                  auto idx_u_j = this->LtG_R(0, idxe, j);
                  auto idx_p_i = this->LtG_R(1, idxe, i);
                  auto idx_p_j = this->LtG_R(1, idxe, j);

                  // pushbacks
                  tpl_ke.push_back(T(idx_u_i, idx_u_j, tmp_uu));    // U'U
                  tpl_ke.push_back(T(idx_u_j, idx_u_i, tmp_uu));    // UU'
                  tpl_ke.push_back(T(idx_p_i, idx_u_j, tmp_pdu));   // P'U
                  tpl_ke.push_back(T(idx_u_j, idx_p_i, tmp_pdu));   // UP'
                }
              }

              for (int i = 0; i < q.N(); ++i) {
                for (int j = 0; j < q.N(); ++j) {
                  auto tmp_uu = 0.0;
                  auto tmp_pp = 0.0;

                  // go over element
                  for (int k = 0; k < q.N(); ++k) {
                    auto rk = _mesh.NodeRadius(idxe, k);
                    auto tmp_ddrr = rk * rk * mat_d[i][k] * mat_d[j][k];
                    tmp_uu += q.W(k) * _mesh_model.C(idxe, k) * tmp_ddrr;

                    tmp_pp += q.W(k) * tmp_ddrr;
                  }

                  // multiply
                  tmp_uu *= e2;
                  tmp_pp *= e2 / (4.0 * pi_db * bigg_db);

                  // indices
                  auto idx_u_i = this->LtG_R(0, idxe, i);
                  auto idx_u_j = this->LtG_R(0, idxe, j);
                  auto idx_p_i = this->LtG_R(1, idxe, i);
                  auto idx_p_j = this->LtG_R(1, idxe, j);

                  // pushbacks
                  tpl_ke.push_back(T(idx_u_i, idx_u_j, tmp_uu));   // U'U
                  tpl_ke.push_back(T(idx_p_i, idx_p_j, tmp_pp));   // P'P
                }
              }
            }
          }

          {
            mat_ke_0.resize(num_radial_dof, num_radial_dof);
            mat_ke_0.setFromTriplets(tpl_ke.begin(), tpl_ke.end());
            mat_ke_0.makeCompressed();
          }
          // std::cout << "Post radial KE FINAL.\n";
          // std::cout << "Completed L = " << idxl << " matrices\n";
          // }
        }
      }
    }
  }
};

auto
sem::LtG_S(int neig, int idx_e, int idx_n) const {
  assert((neig >= 0) && (neig <= 2) &&
         "Error: neig must be 0, 1 or 2 in LtG_S");
  assert((idx_e >= 0) && (idx_e < _mesh.NE()) &&
         "Error: idx_e out of range in LtG_S");
  assert((idx_n >= 0) && (idx_n < _mesh.NN()) &&
         "Error: idx_n out of range in LtG_S");
  std::size_t retval = 0;

  int offset_val = vec_offset[idx_e];
  if (idx_n == 0) {
    if ((std::find(fsb.begin(), fsb.end(), idx_e - 1) != fsb.end())) {
      if (neig == 1) {
        offset_val += 1;
      } else {
        offset_val -= 1;
      }
    }
  }
  retval = (3 * idx_e * (_mesh.NN() - 1) + idx_n * 3 + neig) + offset_val;
  return retval;
};

auto
sem::LtG_R(int neig, int idx_e, int idx_n) const {
  assert((neig >= 0) && (neig < 2) && "Error: neig must be 0 or 1 in LtG_R");
  assert((idx_e >= 0) && (idx_e < _mesh.NE()) &&
         "Error: idx_e out of range in LtG_R");
  assert((idx_n >= 0) && (idx_n < _mesh.NN()) &&
         "Error: idx_n out of range in LtG_R");

  // first index in element
  std::size_t retval = 2 * idx_e * (_mesh.NN() - 1);

  // then node
  retval += idx_n * 2;

  // then function
  retval += neig;

  // return
  return retval;
};

auto
sem::LtG_T(int idx_e, int idx_n) const {
  assert((idx_e >= _el) && (idx_e < _eu) &&
         "Error: idx_e out of range in LtG_T");
  assert((idx_n >= 0) && (idx_n < _mesh.NN()) &&
         "Error: idx_n out of range in LtG_T");
  std::size_t retval = (idx_e - _el) * (_mesh.NN() - 1) + idx_n;
  return retval;
};

Eigen::MatrixXcd
sem::CalculateForce(SourceInfo::EarthquakeCMT &cmt, int idxl) {

  //////////////////////////////
  int NQ = _mesh.NN();
  totlen = this->LtG_S(2, _mesh.NE() - 1, NQ - 1) + 1;
  Eigen::MatrixXcd vec_lforce = Eigen::MatrixXcd::Zero(totlen, 2 * idxl + 1);
  double kval =
      std::sqrt(static_cast<double>(idxl) * (static_cast<double>(idxl) + 1.0));
  double kd2 = kval / std::sqrt(2.0);

  // find element within which the source sits
  double rad_source = _mesh.PR() - 1000.0 * cmt.Depth() / _length_norm;

  // source location in spherical coordinates
  double theta_s = (90.0 - cmt.Latitude()) * EIGEN_PI / (180.0);
  double phi_s = cmt.Longitude() * EIGEN_PI / (180.0);

  //////////////////////////////
  // wigner matrix for evaluation
  int maxn = 2;
  if (maxn > _lmax) {
    maxn = _lmax;
  }
  auto wigdmat =
      GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                       GSHTrans::Single, GSHTrans::ColumnMajor>(_lmax, _lmax,
                                                                maxn, theta_s);
  // ylmn lambda
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(Complex(0.0, m * phi));
    return ylm;
  };

  //////////////////////////////
  // parameters
  double invsq2 = 1.0 / std::sqrt(2.0);
  Complex isq2 = Complex(0.0, 1.0 / std::sqrt(2.0));
  double omegal2 = std::sqrt((idxl + 2) * (idxl - 1) / 2.0);
  // double lprefac = std::exp(-2.0 * this->pi_db * (idxl + 1) / (1 + 0.5));
  double lprefac = 1.0;

  // loop through the elements
  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_source) && (_mesh.EUR(idx) >= rad_source)) {
      // std::cout << "\nSource located in element " << idx << "\n\n";
      //////////////////////////////
      stdvec vec_nodes(NQ, 0.0);
      for (int idxn = 0; idxn < NQ; ++idxn) {
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      }

      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());
      //////////////////////////////

      for (int idxq = 0; idxq < NQ; ++idxq) {
        auto w_val = pleg(idxq, rad_source) / rad_source;
        auto w_deriv = pleg.Derivative(idxq, rad_source);

        //////////////////////////////
        // indices
        auto idx_u = this->LtG_S(0, idx, idxq);
        auto idx_v = this->LtG_S(1, idx, idxq);

        // loop over m
        for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
          // spherical harmonic values
          Complex y0c = std::conj(ylmn(idxl, idxm, 0, phi_s));
          Complex ymc = std::conj(ylmn(idxl, idxm, -1, phi_s));
          Complex ypc = std::conj(ylmn(idxl, idxm, 1, phi_s));
          Complex ymmc = 0.0, yppc = 0.0;
          if (idxl > 1) {
            ymmc = std::conj(ylmn(idxl, idxm, -2, phi_s));
            yppc = std::conj(ylmn(idxl, idxm, 2, phi_s));
          }

          //////////////////////////////
          // common values
          Complex tmp_0pm = cmt.MC0p() * ypc + cmt.MC0m() * ymc;
          Complex tmp_pm = 2.0 * cmt.MCmp() * y0c;
          Complex tmp_ppmm = cmt.MCpp() * yppc + cmt.MCmm() * ymmc;

          // u component
          Complex tmp_u = w_val * (kd2 * tmp_0pm - tmp_pm);
          tmp_u += cmt.MC00() * w_deriv * y0c;

          // v component
          Complex tmp_v = w_deriv * tmp_0pm;
          tmp_v += w_val * (kd2 * tmp_pm - tmp_0pm + omegal2 * tmp_ppmm);
          tmp_v *= kd2;

          //////////////////////////////
          // put in force
          vec_lforce(idx_u, idxm + idxl) = tmp_u * lprefac;
          vec_lforce(idx_v, idxm + idxl) = tmp_v * lprefac;

          //////////////////////////////
        };
      };
    };
  };

  vec_lforce *= (1.0 / _moment_norm);
  return vec_lforce;
};

Eigen::MatrixXcd
sem::CalculateForce_All(SourceInfo::EarthquakeCMT &cmt, int idxl) {

  //////////////////////////////
  int NQ = _mesh.NN();
  totlen = this->LtG_S(2, _mesh.NE() - 1, NQ - 1) + 1;
  Eigen::MatrixXcd vec_lforce = Eigen::MatrixXcd::Zero(totlen, 4);
  double kval =
      std::sqrt(static_cast<double>(idxl) * (static_cast<double>(idxl) + 1.0));
  double kd2 = kval / std::sqrt(2.0);

  // find element within which the source sits
  double rad_source = _mesh.PR() - 1000.0 * cmt.Depth() / _length_norm;

  // loop through the elements
  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_source) && (_mesh.EUR(idx) >= rad_source)) {
      // std::cout << "\nSource located in element " << idx << "\n\n";
      //////////////////////////////
      stdvec vec_nodes(NQ, 0.0);
      for (int idxn = 0; idxn < NQ; ++idxn) {
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      }

      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());
      //////////////////////////////

      for (int idxq = 0; idxq < NQ; ++idxq) {
        auto w_val = pleg(idxq, rad_source) / rad_source;
        auto w_deriv = pleg.Derivative(idxq, rad_source);

        //////////////////////////////
        // indices
        auto idx_u = this->LtG_S(0, idx, idxq);
        auto idx_v = this->LtG_S(1, idx, idxq);

        // put in force
        vec_lforce(idx_u, 0) = w_deriv;
        vec_lforce(idx_u, 1) = w_val;
        vec_lforce(idx_v, 2) = kd2 * w_deriv;
        vec_lforce(idx_v, 3) = kd2 * w_val;
      };
    };
  };

  return vec_lforce;
};

Eigen::MatrixXcd
sem::CalculateForce_All_T(SourceInfo::EarthquakeCMT &cmt, int idxl) {

  //////////////////////////////
  int NQ = _mesh.NN();
  totlen = this->LtG_T(_mesh.NE() - 1, NQ - 1) + 1;
  Eigen::MatrixXcd vec_lforce = Eigen::MatrixXcd::Zero(totlen, 2);
  double kval =
      std::sqrt(static_cast<double>(idxl) * (static_cast<double>(idxl) + 1.0));
  double kd2 = kval / std::sqrt(2.0);

  // find element within which the source sits
  double rad_source = _mesh.PR() - 1000.0 * cmt.Depth() / _length_norm;

  // loop through the elements
  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_source) && (_mesh.EUR(idx) >= rad_source)) {
      // std::cout << "\nSource located in element " << idx << "\n\n";
      //////////////////////////////
      stdvec vec_nodes(NQ, 0.0);
      for (int idxn = 0; idxn < NQ; ++idxn) {
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      }

      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());
      //////////////////////////////

      for (int idxq = 0; idxq < NQ; ++idxq) {
        auto w_val = pleg(idxq, rad_source) / rad_source;
        auto w_deriv = pleg.Derivative(idxq, rad_source);

        //////////////////////////////
        // indices
        auto idx_v = this->LtG_T(idx, idxq);

        // put in force
        {
          vec_lforce(idx_v, 0) = w_deriv;
          vec_lforce(idx_v, 1) = w_val;
        };
      };
    };
  };

  return vec_lforce;
};

Eigen::MatrixXcd
sem::CalculateForce_Coefficients(SourceInfo::EarthquakeCMT &cmt, int idxl) {

  //////////////////////////////
  int NQ = _mesh.NN();
  totlen = this->LtG_S(2, _mesh.NE() - 1, NQ - 1) + 1;
  Eigen::MatrixXcd vec_lforce = Eigen::MatrixXcd::Zero(2 * idxl + 1, 4);
  double kval =
      std::sqrt(static_cast<double>(idxl) * (static_cast<double>(idxl) + 1.0));
  double kd2 = kval / std::sqrt(2.0);

  // source location in spherical coordinates
  double theta_s = (90.0 - cmt.Latitude()) * EIGEN_PI / (180.0);
  double phi_s = cmt.Longitude() * EIGEN_PI / (180.0);

  //////////////////////////////
  // wigner matrix for evaluation
  int maxn = 2;
  if (maxn > _lmax) {
    maxn = _lmax;
  }
  auto wigdmat =
      GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                       GSHTrans::Single, GSHTrans::ColumnMajor>(_lmax, _lmax,
                                                                maxn, theta_s);
  // ylmn lambda
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(Complex(0.0, m * phi));
    return ylm;
  };

  //////////////////////////////
  // parameters
  // double invsq2 = 1.0 / std::sqrt(2.0);
  // Complex isq2 = Complex(0.0, 1.0 / std::sqrt(2.0));
  double omegal2 = std::sqrt((idxl + 2) * (idxl - 1) / 2.0);
  // double lprefac = std::exp(-2.0 * this->pi_db * (idxl + 1) / (1 + 0.5));
  // double lprefac = 1.0;

  // loop through the elements

  // loop over m
  for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
    // spherical harmonic values
    Complex y0c = std::conj(ylmn(idxl, idxm, 0, phi_s));
    Complex ymc = std::conj(ylmn(idxl, idxm, -1, phi_s));
    Complex ypc = std::conj(ylmn(idxl, idxm, 1, phi_s));
    Complex ymmc = 0.0, yppc = 0.0;
    if (idxl > 1) {
      ymmc = std::conj(ylmn(idxl, idxm, -2, phi_s));
      yppc = std::conj(ylmn(idxl, idxm, 2, phi_s));
    }

    //////////////////////////////
    // common values
    Complex tmp_0pm = cmt.MC0p() * ypc + cmt.MC0m() * ymc;
    Complex tmp_pm = 2.0 * cmt.MCmp() * y0c;
    Complex tmp_ppmm = cmt.MCpp() * yppc + cmt.MCmm() * ymmc;

    //////////////////////////////
    // put in force
    vec_lforce(idxm + idxl, 0) = cmt.MC00() * y0c;
    vec_lforce(idxm + idxl, 1) = kd2 * tmp_0pm - tmp_pm;
    vec_lforce(idxm + idxl, 2) = tmp_0pm;
    vec_lforce(idxm + idxl, 3) = kd2 * tmp_pm - tmp_0pm + omegal2 * tmp_ppmm;

    //////////////////////////////
  };

  vec_lforce *= (1.0 / _moment_norm);
  return vec_lforce;
};

Eigen::MatrixXcd
sem::CalculateForce_Coefficients_T(SourceInfo::EarthquakeCMT &cmt, int idxl) {

  //////////////////////////////
  int NQ = _mesh.NN();
  totlen = this->LtG_T(_mesh.NE() - 1, NQ - 1) + 1;
  Eigen::MatrixXcd vec_lforce = Eigen::MatrixXcd::Zero(2 * idxl + 1, 2);
  double kval =
      std::sqrt(static_cast<double>(idxl) * (static_cast<double>(idxl) + 1.0));
  double kd2 = kval / std::sqrt(2.0);

  // source location in spherical coordinates
  double theta_s = (90.0 - cmt.Latitude()) * EIGEN_PI / (180.0);
  double phi_s = cmt.Longitude() * EIGEN_PI / (180.0);

  //////////////////////////////
  // wigner matrix for evaluation
  int maxn = 2;
  if (maxn > _lmax) {
    maxn = _lmax;
  }
  auto wigdmat =
      GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                       GSHTrans::Single, GSHTrans::ColumnMajor>(_lmax, _lmax,
                                                                maxn, theta_s);
  // ylmn lambda
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(Complex(0.0, m * phi));
    return ylm;
  };

  //////////////////////////////
  // parameters
  Complex isq2 = Complex(0.0, 1.0 / std::sqrt(2.0));
  double omegal2 = std::sqrt((idxl + 2) * (idxl - 1) / 2.0);

  // loop over m
  for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
    // spherical harmonic values
    Complex ymc = std::conj(ylmn(idxl, idxm, -1, phi_s));
    Complex ypc = std::conj(ylmn(idxl, idxm, 1, phi_s));
    Complex ymmc = 0.0, yppc = 0.0;
    if (idxl > 1) {
      ymmc = std::conj(ylmn(idxl, idxm, -2, phi_s));
      yppc = std::conj(ylmn(idxl, idxm, 2, phi_s));
    }

    auto tmp1 = cmt.MC0m() * ymc - cmt.MC0p() * ypc;
    auto tmp2 = cmt.MCmm() * ymmc - cmt.MCpp() * yppc;
    //////////////////////////////
    // put in force
    vec_lforce(idxm + idxl, 0) = isq2 * tmp1;
    vec_lforce(idxm + idxl, 1) = isq2 * (omegal2 * tmp2 - tmp1);

    //////////////////////////////
  };

  vec_lforce *= (1.0 / _moment_norm);
  return vec_lforce;
};

Eigen::MatrixXcd
sem::CalculateForce_SPH_U(SourceInfo::EarthquakeCMT &cmt, int idxl) {

  //////////////////////////////
  int NQ = _mesh.NN();
  totlen = this->LtG_S(2, _mesh.NE() - 1, NQ - 1) + 1;
  Eigen::MatrixXcd vec_lforce = Eigen::MatrixXcd::Zero(totlen, 1);

  // find element within which the source sits
  double rad_source = _mesh.PR() - 1000.0 * cmt.Depth() / _length_norm;

  // loop through the elements
  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_source) && (_mesh.EUR(idx) >= rad_source)) {

      //////////////////////////////
      stdvec vec_nodes(NQ, 0.0);
      for (int idxn = 0; idxn < NQ; ++idxn) {
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      }

      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());
      //////////////////////////////

      for (int idxq = 0; idxq < NQ; ++idxq) {

        //////////////////////////////
        // indices
        auto idx_u = this->LtG_S(0, idx, idxq);

        //////////////////////////////
        // put in force
        auto w_val = pleg(idxq, rad_source) / rad_source;
        vec_lforce(idx_u, 0) = w_val;
      };
    };
  };

  vec_lforce *= (1.0 / _moment_norm);
  return vec_lforce;
};

Eigen::MatrixXcd
sem::CalculateForce_SPH_UP(SourceInfo::EarthquakeCMT &cmt, int idxl) {

  //////////////////////////////
  int NQ = _mesh.NN();
  totlen = this->LtG_S(2, _mesh.NE() - 1, NQ - 1) + 1;
  Eigen::MatrixXcd vec_lforce = Eigen::MatrixXcd::Zero(totlen, 1);

  // find element within which the source sits
  double rad_source = _mesh.PR() - 1000.0 * cmt.Depth() / _length_norm;

  // loop through the elements
  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_source) && (_mesh.EUR(idx) >= rad_source)) {

      //////////////////////////////
      stdvec vec_nodes(NQ, 0.0);
      for (int idxn = 0; idxn < NQ; ++idxn) {
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      }

      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());
      //////////////////////////////

      for (int idxq = 0; idxq < NQ; ++idxq) {

        //////////////////////////////
        // indices
        auto idx_u = this->LtG_S(0, idx, idxq);

        //////////////////////////////
        // put in force
        auto w_deriv = pleg.Derivative(idxq, rad_source);
        vec_lforce(idx_u, 0) = w_deriv;
      };
    };
  };

  vec_lforce *= (1.0 / _moment_norm);
  return vec_lforce;
};

Eigen::MatrixXcd
sem::CalculateForce_SPH_V(SourceInfo::EarthquakeCMT &cmt, int idxl) {

  //////////////////////////////
  int NQ = _mesh.NN();
  totlen = this->LtG_S(2, _mesh.NE() - 1, NQ - 1) + 1;
  Eigen::MatrixXcd vec_lforce = Eigen::MatrixXcd::Zero(totlen, 1);
  double kd2 = std::sqrt(static_cast<double>(idxl) *
                         (static_cast<double>(idxl) + 1.0) / 2.0);
  // find element within which the source sits
  double rad_source = _mesh.PR() - 1000.0 * cmt.Depth() / _length_norm;

  // loop through the elements
  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_source) && (_mesh.EUR(idx) >= rad_source)) {

      //////////////////////////////
      stdvec vec_nodes(NQ, 0.0);
      for (int idxn = 0; idxn < NQ; ++idxn) {
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      }

      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());
      //////////////////////////////

      for (int idxq = 0; idxq < NQ; ++idxq) {

        //////////////////////////////
        // indices
        auto idx_v = this->LtG_S(1, idx, idxq);

        //////////////////////////////
        // put in force
        auto w_val = kd2 * pleg(idxq, rad_source) / rad_source;
        vec_lforce(idx_v, 0) = w_val;
      };
    };
  };

  vec_lforce *= (1.0 / _moment_norm);
  return vec_lforce;
};

Eigen::MatrixXcd
sem::CalculateForce_SPH_VP(SourceInfo::EarthquakeCMT &cmt, int idxl) {

  //////////////////////////////
  int NQ = _mesh.NN();
  totlen = this->LtG_S(2, _mesh.NE() - 1, NQ - 1) + 1;
  Eigen::MatrixXcd vec_lforce = Eigen::MatrixXcd::Zero(totlen, 1);
  double kd2 = std::sqrt(static_cast<double>(idxl) *
                         (static_cast<double>(idxl) + 1.0) / 2.0);
  // find element within which the source sits
  double rad_source = _mesh.PR() - 1000.0 * cmt.Depth() / _length_norm;

  // loop through the elements
  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_source) && (_mesh.EUR(idx) >= rad_source)) {

      //////////////////////////////
      stdvec vec_nodes(NQ, 0.0);
      for (int idxn = 0; idxn < NQ; ++idxn) {
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      }

      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());
      //////////////////////////////

      for (int idxq = 0; idxq < NQ; ++idxq) {

        //////////////////////////////
        // indices
        auto idx_v = this->LtG_S(1, idx, idxq);

        //////////////////////////////
        // put in force
        auto w_deriv = kd2 * pleg.Derivative(idxq, rad_source);
        vec_lforce(idx_v, 0) = w_deriv;
      };
    };
  };

  vec_lforce *= (1.0 / _moment_norm);
  return vec_lforce;
};

Eigen::MatrixXcd
sem::CalculateForce_T(SourceInfo::EarthquakeCMT &cmt, int idxl) {

  //////////////////////////////
  int NQ = _mesh.NN();
  totlen = this->LtG_T(_eu - 1, NQ - 1) + 1;
  Eigen::MatrixXcd vec_lforce = Eigen::MatrixXcd::Zero(totlen, 2 * idxl + 1);
  double kval =
      std::sqrt(static_cast<double>(idxl) * (static_cast<double>(idxl) + 1.0));
  double kd2 = kval / std::sqrt(2.0);

  // find element within which the source sits
  double rad_source = _mesh.PR() - 1000.0 * cmt.Depth() / _length_norm;

  // source location in spherical coordinates
  double theta_s = (90.0 - cmt.Latitude()) * EIGEN_PI / (180.0);
  double phi_s = cmt.Longitude() * EIGEN_PI / (180.0);

  //////////////////////////////
  // wigner matrix for evaluation
  int maxn = 2;
  if (maxn > _lmax) {
    maxn = _lmax;
  }
  auto wigdmat =
      GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                       GSHTrans::Single, GSHTrans::ColumnMajor>(_lmax, _lmax,
                                                                maxn, theta_s);
  // ylmn lambda
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(Complex(0.0, m * phi));
    return ylm;
  };

  //////////////////////////////
  // parameters
  double invsq2 = 1.0 / std::sqrt(2.0);
  Complex isq2 = Complex(0.0, 1.0 / std::sqrt(2.0));
  double omegal2 = std::sqrt((idxl + 2) * (idxl - 1) / 2.0);
  // double lprefac = std::exp(-2.0 * this->pi_db * (idxl + 1) / (1 + 0.5));
  double lprefac = 1.0;

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
        std::size_t ridx = this->LtG_T(idx, idxq);
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

          vec_lforce(ridx, idxm + idxl) = tmp;
        };
        // };
      };
    };
  };
  vec_lforce *= (1.0 / _moment_norm);
  return vec_lforce;
};

Eigen::MatrixXcd
sem::CalculateForce_R(SourceInfo::EarthquakeCMT &cmt) {

  //////////////////////////////
  int NQ = _mesh.NN();
  totlen = this->LtG_R(1, _mesh.NE() - 1, NQ - 1) + 1;
  Eigen::MatrixXcd vec_lforce = Eigen::MatrixXcd::Zero(totlen, 1);
  // double kval =
  //     std::sqrt(static_cast<double>(idxl) * ( 1.0));
  // double kd2 = kval / std::sqrt(2.0);

  // find element within which the source sits
  double rad_source = _mesh.PR() - 1000.0 * cmt.Depth() / _length_norm;

  // source location in spherical coordinates
  double theta_s = (90.0 - cmt.Latitude()) * EIGEN_PI / (180.0);
  double phi_s = cmt.Longitude() * EIGEN_PI / (180.0);

  //////////////////////////////
  // wigner matrix for evaluation
  int maxn = 2;
  if (maxn > _lmax) {
    maxn = _lmax;
  }
  auto wigdmat =
      GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                       GSHTrans::Single, GSHTrans::ColumnMajor>(_lmax, _lmax,
                                                                maxn, theta_s);
  // ylmn lambda
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(Complex(0.0, m * phi));
    return ylm;
  };

  //////////////////////////////
  // parameters
  double invsq2 = 1.0 / std::sqrt(2.0);
  Complex isq2 = Complex(0.0, 1.0 / std::sqrt(2.0));
  // double omegal2 = std::sqrt((idxl + 2) * (idxl - 1) / 2.0);
  // double lprefac = std::exp(-2.0 * this->pi_db * (idxl + 1) / (1 + 0.5));
  double lprefac = 1.0;

  // loop through the elements
  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_source) && (_mesh.EUR(idx) >= rad_source)) {
      //////////////////////////////
      stdvec vec_nodes(NQ, 0.0);
      for (int idxn = 0; idxn < NQ; ++idxn) {
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      }

      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());
      //////////////////////////////

      for (int idxq = 0; idxq < NQ; ++idxq) {
        auto w_val = pleg(idxq, rad_source) / rad_source;
        auto w_deriv = pleg.Derivative(idxq, rad_source);

        //////////////////////////////
        // indices
        auto idx_u = this->LtG_R(0, idx, idxq);

        // values to use
        Complex y0c = std::conj(ylmn(0, 0, 0, phi_s));
        Complex tmp_pm = 2.0 * cmt.MCmp() * y0c;

        // u component
        Complex tmp_u = -w_val * tmp_pm;
        tmp_u += cmt.MC00() * w_deriv * y0c;

        //////////////////////////////
        // put in force
        vec_lforce(idx_u, 0) = tmp_u * lprefac;

        //////////////////////////////
      };
    };
  };

  vec_lforce *= (1.0 / _moment_norm);
  return vec_lforce;
};

Eigen::MatrixXcd
sem::RV_Z(double theta_r, double phi_r, int idxl) {
  // create the receiver vector
  std::size_t flen = this->LtG_S(2, _mesh.NE() - 1, _mesh.NN() - 1) + 1;
  std::size_t fcols = 2 * idxl + 1;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(flen, fcols);

  //////////////////////////////
  // wigner matrix for evaluation
  using namespace GSHTrans;
  int maxn = 2;
  if (maxn > _lmax) {
    maxn = _lmax;
  }
  auto wigdmat = Wigner<double, Ortho, All, All, Single, ColumnMajor>(
      _lmax, _lmax, maxn, theta_r);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
    return ylm;
  };

  //////////////////////////////
  // indices of outermost nodes
  auto idx_u = this->LtG_S(0, _mesh.NE() - 1, _mesh.NN() - 1);

  //////////////////////////////
  // fill out receiver vector for U
  for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
    // spherical harmonic values
    Complex yl0 = ylmn(idxl, idxm, 0, phi_r);

    // receiver vector
    vec_receiver(idx_u, idxm + idxl) = yl0;
  }

  //////////////////////////////
  //  return receiver vector
  return vec_receiver;
};

Eigen::MatrixXcd
sem::RV_Z(InputParameters &param, int idxl, int idxr) {
  auto rec = param.receivers()[idxr];
  double theta_r = (90.0 - rec.first) * EIGEN_PI / 180.0;
  double phi_r = rec.second * EIGEN_PI / 180.0;
  double rad_r = _mesh.PR() - param.receiver_depth() * 1000.0 / _length_norm;

  // create the receiver vector
  std::size_t flen = this->LtG_S(2, _mesh.NE() - 1, _mesh.NN() - 1) + 1;
  std::size_t fcols = 2 * idxl + 1;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(flen, fcols);

  //////////////////////////////
  // wigner matrix for evaluation
  using namespace GSHTrans;
  int maxn = 2;
  if (maxn > _lmax) {
    maxn = _lmax;
  }
  auto wigdmat = Wigner<double, Ortho, All, All, Single, ColumnMajor>(
      _lmax, _lmax, maxn, theta_r);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
    return ylm;
  };

  // int elem_idx = 0;
  bool evaluated = false;
  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_r) && (_mesh.EUR(idx) >= rad_r) &&
        (!evaluated)) {
      evaluated = true;

      // get the Lagrange polynomial for this element
      std::vector<double> vec_nodes(_mesh.NN(), 0.0);
      for (int idxn = 0; idxn < _mesh.NN(); ++idxn) {
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      }
      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());

      // loop through the nodes and evaluate
      for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
        auto idx_u = this->LtG_S(0, idx, idxq);
        for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
          // spherical harmonic values
          Complex yl0 = ylmn(idxl, idxm, 0, phi_r);

          // receiver vector
          vec_receiver(idx_u, idxm + idxl) = yl0 * pleg(idxq, rad_r);
        }
      }
    }
  }

  return vec_receiver;
};

Eigen::MatrixXcd
sem::RV_BASE_Z(InputParameters &param, int idxl, int idxr) {
  double rad_r = _mesh.PR() - param.receiver_depth() * 1000.0 / _length_norm;

  // create the receiver vector
  std::size_t flen = this->LtG_S(2, _mesh.NE() - 1, _mesh.NN() - 1) + 1;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(1, flen);

  // int elem_idx = 0;
  bool evaluated = false;
  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_r) && (_mesh.EUR(idx) >= rad_r) &&
        (!evaluated)) {
      evaluated = true;

      // get the Lagrange polynomial for this element
      std::vector<double> vec_nodes(_mesh.NN(), 0.0);
      for (int idxn = 0; idxn < _mesh.NN(); ++idxn) {
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      }
      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());

      // loop through the nodes and evaluate
      for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
        auto idx_u = this->LtG_S(0, idx, idxq);
        // receiver vector
        vec_receiver(0, idx_u) = pleg(idxq, rad_r);
      }
    }
  }

  return vec_receiver;
};

Eigen::MatrixXcd
sem::RV_VAL_Z(InputParameters &param, int idxl, int idxr) {
  auto rec = param.receivers()[idxr];
  double theta_r = (90.0 - rec.first) * EIGEN_PI / 180.0;
  double phi_r = rec.second * EIGEN_PI / 180.0;

  // create the receiver vector
  std::size_t fcols = 2 * idxl + 1;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(fcols, 1);

  //////////////////////////////
  // wigner matrix for evaluation
  using namespace GSHTrans;
  int maxn = 2;
  if (maxn > _lmax) {
    maxn = _lmax;
  }
  auto wigdmat = Wigner<double, Ortho, All, All, Single, ColumnMajor>(
      _lmax, _lmax, maxn, theta_r);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
    return ylm;
  };

  for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
    // spherical harmonic values
    Complex yl0 = ylmn(idxl, idxm, 0, phi_r);
    vec_receiver(idxm + idxl, 0) = yl0;
  }
  // Complex yl0 = ylmn(idxl, 0, 0, phi_r);
  // vec_receiver(idxl, 0) = yl0;

  return vec_receiver;
};

Eigen::MatrixXcd
sem::RV_Z_R(InputParameters &param, int idxr) {
  auto rec = param.receivers()[idxr];
  double theta_r = (90.0 - rec.first) * EIGEN_PI / 180.0;
  double phi_r = rec.second * EIGEN_PI / 180.0;
  double rad_r = _mesh.PR() - param.receiver_depth() * 1000.0 / _length_norm;

  // create the receiver vector
  std::size_t flen = this->LtG_R(1, _mesh.NE() - 1, _mesh.NN() - 1) + 1;
  // std::size_t fcols = 2 * idxl + 1;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(flen, 1);

  //////////////////////////////
  // wigner matrix for evaluation
  using namespace GSHTrans;
  int maxn = 2;
  if (maxn > _lmax) {
    maxn = _lmax;
  }
  auto wigdmat = Wigner<double, Ortho, All, All, Single, ColumnMajor>(
      _lmax, _lmax, maxn, theta_r);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
    return ylm;
  };

  // int elem_idx = 0;
  bool evaluated = false;
  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_r) && (_mesh.EUR(idx) >= rad_r) &&
        (!evaluated)) {
      evaluated = true;

      // get the Lagrange polynomial for this element
      std::vector<double> vec_nodes(_mesh.NN(), 0.0);
      for (int idxn = 0; idxn < _mesh.NN(); ++idxn) {
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      }
      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());

      // loop through the nodes and evaluate
      for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
        auto idx_u = this->LtG_R(0, idx, idxq);
        // for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
        // spherical harmonic values
        Complex yl0 = ylmn(0, 0, 0, phi_r);

        // receiver vector
        vec_receiver(idx_u, 0) = yl0 * pleg(idxq, rad_r);
      }
    }
  }

  return vec_receiver;
};

Eigen::MatrixXcd
sem::RV_THETA(double theta_r, double phi_r, int idxl) {
  // create the receiver vector
  std::size_t flen = this->LtG_S(2, _mesh.NE() - 1, _mesh.NN() - 1) + 1;
  std::size_t fcols = 2 * idxl + 1;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(flen, fcols);

  auto k =
      std::sqrt(static_cast<double>(idxl) * (static_cast<double>(idxl) + 1.0));

  //////////////////////////////
  // wigner matrix for evaluation
  using namespace GSHTrans;
  int maxn = 2;
  if (maxn > _lmax) {
    maxn = _lmax;
  }
  auto wigdmat = Wigner<double, Ortho, All, All, Single, ColumnMajor>(
      _lmax, _lmax, maxn, theta_r);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
    return ylm;
  };

  //////////////////////////////
  // indices of outermost nodes
  auto idx_v = this->LtG_S(1, _mesh.NE() - 1, _mesh.NN() - 1);

  //////////////////////////////
  // receiver vector for U and v
  for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
    // spherical harmonic values
    Complex ylm = ylmn(idxl, idxm, -1, phi_r);
    Complex ylp = ylmn(idxl, idxm, 1, phi_r);

    // receiver vector
    vec_receiver(idx_v, idxm + idxl) = k * (ylm - ylp) / 2.0;
  }
  //////////////////////////////

  // return receiver vector
  return vec_receiver;
};

Eigen::MatrixXcd
sem::RV_THETA(InputParameters &param, int idxl, int idxr) {
  auto rec = param.receivers()[idxr];
  double theta_r = (90.0 - rec.first) * EIGEN_PI / 180.0;
  double phi_r = rec.second * EIGEN_PI / 180.0;
  double rad_r = _mesh.PR() - param.receiver_depth() * 1000.0 / _length_norm;

  auto k =
      std::sqrt(static_cast<double>(idxl) * (static_cast<double>(idxl) + 1.0));
  // create the receiver vector
  std::size_t flen = this->LtG_S(2, _mesh.NE() - 1, _mesh.NN() - 1) + 1;
  std::size_t fcols = 2 * idxl + 1;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(flen, fcols);

  //////////////////////////////
  // wigner matrix for evaluation
  using namespace GSHTrans;
  int maxn = 2;
  if (maxn > _lmax) {
    maxn = _lmax;
  }
  auto wigdmat = Wigner<double, Ortho, All, All, Single, ColumnMajor>(
      _lmax, _lmax, maxn, theta_r);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
    return ylm;
  };

  // int elem_idx = 0;
  bool evaluated = false;
  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_r) && (_mesh.EUR(idx) >= rad_r) &&
        (!evaluated)) {
      evaluated = true;

      // get the Lagrange polynomial for this element
      std::vector<double> vec_nodes(_mesh.NN(), 0.0);
      for (int idxn = 0; idxn < _mesh.NN(); ++idxn) {
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      }
      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());

      // loop through the nodes and evaluate
      for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
        auto idx_v = this->LtG_S(1, idx, idxq);
        for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
          // spherical harmonic values
          Complex ylm = ylmn(idxl, idxm, -1, phi_r);
          Complex ylp = ylmn(idxl, idxm, 1, phi_r);

          // receiver vector
          vec_receiver(idx_v, idxm + idxl) =
              k * (ylm - ylp) / 2.0 * pleg(idxq, rad_r);
        }
      }
    }
  }

  return vec_receiver;
};

Eigen::MatrixXcd
sem::RV_BASE_THETA(InputParameters &param, int idxl, int idxr) {

  double rad_r = _mesh.PR() - param.receiver_depth() * 1000.0 / _length_norm;

  auto k =
      std::sqrt(static_cast<double>(idxl) * (static_cast<double>(idxl) + 1.0));
  // create the receiver vector
  std::size_t flen = this->LtG_S(2, _mesh.NE() - 1, _mesh.NN() - 1) + 1;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(1, flen);

  // int elem_idx = 0;
  bool evaluated = false;
  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_r) && (_mesh.EUR(idx) >= rad_r) &&
        (!evaluated)) {
      evaluated = true;

      // get the Lagrange polynomial for this element
      std::vector<double> vec_nodes(_mesh.NN(), 0.0);
      for (int idxn = 0; idxn < _mesh.NN(); ++idxn) {
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      }
      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());

      // loop through the nodes and evaluate
      for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
        auto idx_v = this->LtG_S(1, idx, idxq);

        // receiver vector
        vec_receiver(0, idx_v) = k / 2.0 * pleg(idxq, rad_r);
      }
    }
  }

  return vec_receiver;
};

Eigen::MatrixXcd
sem::RV_BASE_THETA_T(InputParameters &param, int idxl, int idxr) {

  double rad_r = _mesh.PR() - param.receiver_depth() * 1000.0 / _length_norm;

  auto k =
      std::sqrt(static_cast<double>(idxl) * (static_cast<double>(idxl) + 1.0));
  // create the receiver vector
  std::size_t flen = this->LtG_T(_mesh.NE() - 1, _mesh.NN() - 1) + 1;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(1, flen);
  // std::complex<double> i1 = std::complex<double>(0.0, 1.0);
  // int elem_idx = 0;
  bool evaluated = false;
  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_r) && (_mesh.EUR(idx) >= rad_r) &&
        (!evaluated)) {
      evaluated = true;

      // get the Lagrange polynomial for this element
      std::vector<double> vec_nodes(_mesh.NN(), 0.0);
      for (int idxn = 0; idxn < _mesh.NN(); ++idxn) {
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      }
      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());

      // loop through the nodes and evaluate
      for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
        auto idx_v = this->LtG_T(idx, idxq);

        // receiver vector
        vec_receiver(0, idx_v) = 0.5 * pleg(idxq, rad_r);
      }
    }
  }

  return vec_receiver;
};

Eigen::MatrixXcd
sem::RV_VAL_THETA(InputParameters &param, int idxl, int idxr) {
  auto rec = param.receivers()[idxr];
  double theta_r = (90.0 - rec.first) * EIGEN_PI / 180.0;
  double phi_r = rec.second * EIGEN_PI / 180.0;
  std::size_t fcols = 2 * idxl + 1;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(fcols, 1);

  //////////////////////////////
  // wigner matrix for evaluation
  using namespace GSHTrans;
  int maxn = 2;
  if (maxn > _lmax) {
    maxn = _lmax;
  }
  auto wigdmat = Wigner<double, Ortho, All, All, Single, ColumnMajor>(
      _lmax, _lmax, maxn, theta_r);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
    return ylm;
  };

  for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
    // spherical harmonic values
    Complex ylm = ylmn(idxl, idxm, -1, phi_r);
    Complex ylp = ylmn(idxl, idxm, 1, phi_r);
    vec_receiver(idxm + idxl, 0) = ylm - ylp;
  }

  return vec_receiver;
};

Eigen::MatrixXcd
sem::RV_VAL_THETA_T(InputParameters &param, int idxl, int idxr) {
  auto rec = param.receivers()[idxr];
  double theta_r = (90.0 - rec.first) * EIGEN_PI / 180.0;
  double phi_r = rec.second * EIGEN_PI / 180.0;
  std::size_t fcols = 2 * idxl + 1;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(fcols, 1);

  //////////////////////////////
  // wigner matrix for evaluation
  using namespace GSHTrans;
  int maxn = 2;
  if (maxn > _lmax) {
    maxn = _lmax;
  }
  auto wigdmat = Wigner<double, Ortho, All, All, Single, ColumnMajor>(
      _lmax, _lmax, maxn, theta_r);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
    return ylm;
  };
  auto i1 = std::complex<double>(0.0, 1.0);

  for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
    // spherical harmonic values
    Complex ylm = ylmn(idxl, idxm, -1, phi_r);
    Complex ylp = ylmn(idxl, idxm, 1, phi_r);
    vec_receiver(idxm + idxl, 0) = -i1 * (ylm + ylp);
  }

  return vec_receiver;
};

Eigen::MatrixXcd
sem::RV_THETA_T(InputParameters &param, int idxl, int idxr) {
  auto rec = param.receivers()[idxr];
  double theta_r = (90.0 - rec.first) * EIGEN_PI / 180.0;
  double phi_r = rec.second * EIGEN_PI / 180.0;
  double rad_r = _mesh.PR() - param.receiver_depth() * 1000.0 / _length_norm;

  auto k =
      std::sqrt(static_cast<double>(idxl) * (static_cast<double>(idxl) + 1.0));
  // create the receiver vector
  std::size_t flen = this->LtG_T(_eu - 1, _mesh.NN() - 1) + 1;
  std::size_t fcols = 2 * idxl + 1;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(flen, fcols);

  //////////////////////////////
  // wigner matrix for evaluation
  using namespace GSHTrans;
  int maxn = 2;
  if (maxn > _lmax) {
    maxn = _lmax;
  }
  auto wigdmat = Wigner<double, Ortho, All, All, Single, ColumnMajor>(
      _lmax, _lmax, maxn, theta_r);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    std::complex<double> ylm =
        tmp * std::exp(std::complex<double>(0.0, m * phi));
    return ylm;
  };
  Complex i1 = std::complex<double>(0.0, 1.0);
  // int elem_idx = 0;
  bool evaluated = false;
  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_r) && (_mesh.EUR(idx) >= rad_r) &&
        (!evaluated)) {
      evaluated = true;

      // get the Lagrange polynomial for this element
      std::vector<double> vec_nodes(_mesh.NN(), 0.0);
      for (int idxn = 0; idxn < _mesh.NN(); ++idxn) {
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      }
      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());

      // loop through the nodes and evaluate
      for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
        auto idx_v = this->LtG_T(idx, idxq);
        for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
          // spherical harmonic values
          Complex ylm = ylmn(idxl, idxm, -1, phi_r);
          Complex ylp = ylmn(idxl, idxm, 1, phi_r);

          // receiver vector
          vec_receiver(idx_v, idxm + idxl) =
              -i1 * pleg(idxq, rad_r) * (ylp + ylm) / 2.0;
        }
      }
    }
  }

  return vec_receiver;
};

Eigen::MatrixXcd
sem::RV_PHI(double theta_r, double phi_r, int idxl) {
  // create the receiver vector
  std::size_t flen = this->LtG_S(2, _mesh.NE() - 1, _mesh.NN() - 1) + 1;
  std::size_t fcols = 2 * idxl + 1;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(flen, fcols);

  auto k =
      std::sqrt(static_cast<double>(idxl) * (static_cast<double>(idxl) + 1.0));

  //////////////////////////////
  // wigner matrix for evaluation
  using namespace GSHTrans;
  int maxn = 2;
  if (maxn > _lmax) {
    maxn = _lmax;
  }
  auto wigdmat = Wigner<double, Ortho, All, All, Single, ColumnMajor>(
      _lmax, _lmax, maxn, theta_r);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
    return ylm;
  };
  Complex i1 = std::complex<double>(0.0, 1.0);

  //////////////////////////////
  // indices of outermost nodes
  auto idx_v = this->LtG_S(1, _mesh.NE() - 1, _mesh.NN() - 1);

  //////////////////////////////
  // receiver vector for U and v
  for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
    // spherical harmonic values
    Complex ylm = ylmn(idxl, idxm, -1, phi_r);
    Complex ylp = ylmn(idxl, idxm, 1, phi_r);

    // receiver vector
    vec_receiver(idx_v, idxm + idxl) = -i1 * k * (ylm + ylp) / 2.0;
  }

  //////////////////////////////
  // return receiver vector
  return vec_receiver;
};

Eigen::MatrixXcd
sem::RV_PHI(InputParameters &param, int idxl, int idxr) {
  auto rec = param.receivers()[idxr];
  double theta_r = (90.0 - rec.first) * EIGEN_PI / 180.0;
  double phi_r = rec.second * EIGEN_PI / 180.0;
  double rad_r = _mesh.PR() - param.receiver_depth() * 1000.0 / _length_norm;

  auto k =
      std::sqrt(static_cast<double>(idxl) * (static_cast<double>(idxl) + 1.0));

  // create the receiver vector
  std::size_t flen = this->LtG_S(2, _mesh.NE() - 1, _mesh.NN() - 1) + 1;
  std::size_t fcols = 2 * idxl + 1;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(flen, fcols);

  //////////////////////////////
  // wigner matrix for evaluation
  using namespace GSHTrans;
  int maxn = 2;
  if (maxn > _lmax) {
    maxn = _lmax;
  }
  auto wigdmat = Wigner<double, Ortho, All, All, Single, ColumnMajor>(
      _lmax, _lmax, maxn, theta_r);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
    return ylm;
  };
  Complex i1 = std::complex<double>(0.0, 1.0);

  // int elem_idx = 0;
  bool evaluated = false;
  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_r) && (_mesh.EUR(idx) >= rad_r) &&
        (!evaluated)) {
      evaluated = true;

      // get the Lagrange polynomial for this element
      std::vector<double> vec_nodes(_mesh.NN(), 0.0);
      for (int idxn = 0; idxn < _mesh.NN(); ++idxn) {
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      }
      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());

      // loop through the nodes and evaluate
      for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
        auto idx_v = this->LtG_S(1, idx, idxq);
        for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
          // spherical harmonic values
          Complex ylm = ylmn(idxl, idxm, -1, phi_r);
          Complex ylp = ylmn(idxl, idxm, 1, phi_r);

          // receiver vector
          vec_receiver(idx_v, idxm + idxl) =
              -k * i1 * (ylm + ylp) / 2.0 * pleg(idxq, rad_r);
        }
      }
    }
  }

  return vec_receiver;
};

Eigen::MatrixXcd
sem::RV_BASE_PHI_T(InputParameters &param, int idxl, int idxr) {
  return this->RV_BASE_THETA_T(param, idxl, idxr);
};

Eigen::MatrixXcd
sem::RV_VAL_PHI(InputParameters &param, int idxl, int idxr) {
  auto rec = param.receivers()[idxr];
  double theta_r = (90.0 - rec.first) * EIGEN_PI / 180.0;
  double phi_r = rec.second * EIGEN_PI / 180.0;
  std::size_t fcols = 2 * idxl + 1;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(fcols, 1);

  //////////////////////////////
  // wigner matrix for evaluation
  using namespace GSHTrans;
  int maxn = 2;
  if (maxn > _lmax) {
    maxn = _lmax;
  }
  auto wigdmat = Wigner<double, Ortho, All, All, Single, ColumnMajor>(
      _lmax, _lmax, maxn, theta_r);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
    return ylm;
  };
  Complex i1 = std::complex<double>(0.0, 1.0);
  for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
    // spherical harmonic values
    Complex ylm = ylmn(idxl, idxm, -1, phi_r);
    Complex ylp = ylmn(idxl, idxm, 1, phi_r);
    vec_receiver(idxm + idxl, 0) = -i1 * (ylm + ylp);
  }

  return vec_receiver;
};

Eigen::MatrixXcd
sem::RV_VAL_PHI_T(InputParameters &param, int idxl, int idxr) {
  auto rec = param.receivers()[idxr];
  double theta_r = (90.0 - rec.first) * EIGEN_PI / 180.0;
  double phi_r = rec.second * EIGEN_PI / 180.0;
  std::size_t fcols = 2 * idxl + 1;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(fcols, 1);

  //////////////////////////////
  // wigner matrix for evaluation
  using namespace GSHTrans;
  int maxn = 2;
  if (maxn > _lmax) {
    maxn = _lmax;
  }
  auto wigdmat = Wigner<double, Ortho, All, All, Single, ColumnMajor>(
      _lmax, _lmax, maxn, theta_r);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
    return ylm;
  };

  for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
    // spherical harmonic values
    Complex ylm = ylmn(idxl, idxm, -1, phi_r);
    Complex ylp = ylmn(idxl, idxm, 1, phi_r);
    vec_receiver(idxm + idxl, 0) = (ylp - ylm);
  }

  return vec_receiver;
};

Eigen::MatrixXcd
sem::RV_PHI_T(InputParameters &param, int idxl, int idxr) {
  auto rec = param.receivers()[idxr];
  double theta_r = (90.0 - rec.first) * EIGEN_PI / 180.0;
  double phi_r = rec.second * EIGEN_PI / 180.0;
  double rad_r = _mesh.PR() - param.receiver_depth() * 1000.0 / _length_norm;

  auto k =
      std::sqrt(static_cast<double>(idxl) * (static_cast<double>(idxl) + 1.0));
  // create the receiver vector
  std::size_t flen = this->LtG_T(_eu - 1, _mesh.NN() - 1) + 1;
  std::size_t fcols = 2 * idxl + 1;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(flen, fcols);

  //////////////////////////////
  // wigner matrix for evaluation
  using namespace GSHTrans;
  int maxn = 2;
  if (maxn > _lmax) {
    maxn = _lmax;
  }
  auto wigdmat = Wigner<double, Ortho, All, All, Single, ColumnMajor>(
      _lmax, _lmax, maxn, theta_r);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    std::complex<double> ylm =
        tmp * std::exp(std::complex<double>(0.0, m * phi));
    return ylm;
  };
  Complex i1 = std::complex<double>(0.0, 1.0);
  // int elem_idx = 0;
  bool evaluated = false;
  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_r) && (_mesh.EUR(idx) >= rad_r) &&
        (!evaluated)) {
      evaluated = true;

      // get the Lagrange polynomial for this element
      std::vector<double> vec_nodes(_mesh.NN(), 0.0);
      for (int idxn = 0; idxn < _mesh.NN(); ++idxn) {
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      }
      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());

      // loop through the nodes and evaluate
      for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
        auto idx_v = this->LtG_T(idx, idxq);
        for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
          // spherical harmonic values
          Complex ylm = ylmn(idxl, idxm, -1, phi_r);
          Complex ylp = ylmn(idxl, idxm, 1, phi_r);

          // receiver vector
          vec_receiver(idx_v, idxm + idxl) =
              pleg(idxq, rad_r) * (ylp - ylm) / 2.0;
        }
      }
    }
  }

  return vec_receiver;
};

auto
sem::Receiver_Elements(InputParameters &param) const {
  std::vector<int> receiver_elems;
  double rad_receiver =
      _mesh.PR() - 1000.0 * param.receiver_depth() / _length_norm;

  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_receiver) && (_mesh.EUR(idx) >= rad_receiver)) {
      receiver_elems.push_back(idx);
    }
  }

  return receiver_elems;
};

}   // namespace Full1D

// void
// sem::CalculateEigenfrequencies(int N, double sigshift) {
//   // set number of modes
//   _num_modes = N;
//   // initiate matrix multiplicaiton wrapper
//   // std::cout << "Eig Check 1\n";
//   SparseGenMatProd<double> op_seig(mat_seig);
//   // std::cout << "Eig Check 2\n";
//   using OpType = SymShiftInvert<double, Eigen::Sparse, Eigen::Sparse>;
//   using BOpType = SparseSymMatProd<double>;
//   OpType op_mult(mat_ke, mat_inertia);
//   BOpType op_ke(mat_ke), op_inertia(mat_inertia);
//   // std::cout << "Eig Check 3\n";

//   // std::cout << "Eig Check 4\n";

//   // get nc
//   long int maxn, maxn2;
//   if ((5 * N) > mat_seig.rows()) {
//     maxn = mat_seig.rows();
//   } else {
//     maxn = 5 * N;
//   }
//   if ((3 * N) > mat_seig.rows()) {
//     maxn2 = mat_seig.rows();
//   } else {
//     maxn2 = 3 * N;
//   }

//   // initiate eigensolver
//   // std::cout << "Eig Check 5\n";
//   // GenEigsSolver<SparseGenMatProd<double>> eigs(op_seig, N, maxn);
//   // std::cout << "Eig Check 6\n";
//   // std::cout << "Enter shift:\n";
//   // double sigshift;
//   // std::cin >> sigshift;
//   SymGEigsShiftSolver<OpType, BOpType, GEigsMode::ShiftInvert> eig_gen(
//       op_mult, op_inertia, N, maxn2, sigshift);
//   // eigs.init();
//   eig_gen.init();
//   // int nconv_seig = eigs.compute(SortRule::SmallestMagn);
//   int nconv_gen = eig_gen.compute(SortRule::LargestMagn);

//   /*
//   if (eigs.info() == CompInfo::Successful) {
//     std::cout << "Successful\n";
//     evalues_seig = eigs.eigenvalues();

//     // check if we are computing for the IC
//     if ((_has_fluid && (_il == 0)) || (!_has_fluid)) {
//       std::size_t nrow = eigs.eigenvectors().rows() + 1;
//       std::size_t ncol = eigs.eigenvectors().cols();
//       evectors_gen.resize(nrow, ncol);
//       evectors_gen.block(0, 0, 1, ncol) = Eigen::MatrixXcd::Zero(1, ncol);
//       evectors_gen.block(1, 0, nrow - 1, ncol) = eigs.eigenvectors();
//     } else {
//       std::size_t nrow = eigs.eigenvectors().rows();
//       std::size_t ncol = eigs.eigenvectors().cols();
//       evectors_gen.resize(nrow, ncol);
//       evectors = eigs.eigenvectors();
//     }

//     _calc_eig = true;

//     // std::cout << "Eigenvalues from non-generalised: \n"
//     //           << eigs.eigenvalues() << "\n\n";

//     // evectors = eigs.eigenvectors();
//   } else {
//     std::cout << "Unsuccessful non-generalised\n";
//   }
//   */

//   if (eig_gen.info() == CompInfo::Successful) {
//     // std::cout << "Eigenvalues from generalised: \n"
//     //           << eig_gen.eigenvalues() * _freq_norm * _freq_norm <<
//     // "\n\n";

//     evalues_gen = eig_gen.eigenvalues();
//     Eigen::MatrixXcd normval = eig_gen.eigenvectors().transpose() *
//                                mat_inertia * eig_gen.eigenvectors();
//     // std::cout << normval << "\n\n";

//     // check if we are computing for the IC
//     if ((_has_fluid && (_il == 0)) || (!_has_fluid)) {
//       std::size_t nrow = eig_gen.eigenvectors().rows() + 1;
//       std::size_t ncol = eig_gen.eigenvectors().cols();
//       evectors_gen.resize(nrow, ncol);
//       evectors_gen.block(0, 0, 1, ncol) = Eigen::MatrixXcd::Zero(1, ncol);
//       evectors_gen.block(1, 0, nrow - 1, ncol) = eig_gen.eigenvectors();

//     } else {
//       std::size_t nrow = eig_gen.eigenvectors().rows();
//       std::size_t ncol = eig_gen.eigenvectors().cols();
//       evectors_gen.resize(nrow, ncol);
//       evectors_gen = eig_gen.eigenvectors();
//     }

//     // normalise:
//     for (int i = 0; i < evectors_gen.rows(); ++i) {
//       for (int j = 0; j < evectors_gen.cols(); ++j) {
//         evectors_gen(i, j) *= 1.0 / (normint * sqrt(evalues_gen(j).real()));
//       }
//     }
//     _calc_gen = true;
//     // std::cout << "POST\n";
//   } else {
//     std::cout << "Unsuccessful generalised\n";
//   }

//   if (_calc_gen) {
//     // save the eigenvectors in the "standard format":
//     int NQ = _mesh.NN();
//     for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
//       stdvvec vec_out;
//       if ((idxe < _el) || (idxe > _eu - 1)) {
//         vec_eigval.push_back(stdvvec(NQ, stdvec(evectors_gen.cols(), 0.0)));
//       } else {
//         for (int idxq = 0; idxq < NQ; ++idxq) {
//           std::vector<double> vec_x;
//           std::size_t ovidx = (idxe - _el) * (NQ - 1) + idxq;
//           for (int idx = 0; idx < evectors_gen.cols(); ++idx) {
//             std::size_t colidx = evectors_gen.cols() - 1 - idx;
//             std::complex<double> eigint = evectors_gen(ovidx, colidx);
//             vec_x.push_back(eigint.real());
//           }
//           vec_out.push_back(vec_x);
//         }
//         vec_eigval.push_back(vec_out);
//       }
//     }

//     // calculate the derivatives of the eigenvectors:
//     vec_eigderiv =
//         stdvvvec(_mesh.NE(), stdvvec(NQ, stdvec(evectors_gen.cols(), 0.0)));
//     for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
//       stdvvec vec_out;
//       double elem_width = _mesh.EW(idxe);
//       for (int idxq = 0; idxq < NQ; ++idxq) {
//         for (int idx = 0; idx < evectors_gen.cols(); ++idx) {
//           double tmp = 0.0;
//           for (int idxq2 = 0; idxq2 < NQ; ++idxq2) {
//             // std::size_t ovidx = overallidxfinal(idx, idxq2);
//             // tmp += evectors(ovidx, i) * vec_lagderiv[idxq2][idxq];
//             tmp += vec_eigval[idxe][idxq2][idx] * vec_lag_deriv[idxq][idxq2];
//             // tmp += vec_lag_deriv[idxq][idxq2];
//           }
//           tmp *= 2.0 / elem_width;
//           vec_eigderiv[idxe][idxq][idx] = tmp;
//         }
//       }
//     }
//   }

//   // // calculate augmented
//   // this->augment_basis_calculate();
//   // vec_fullbasis = stdvvvec(_mesh.NumberOfElements(),
//   //                          stdvvec(NQ, stdvec(nummodes + numaug,
//   //  0.0)));
//   // vec_fullderiv = stdvvvec(_mesh.NumberOfElements(),
//   //                          stdvvec(NQ, stdvec(nummodes + numaug,
//   //  0.0)));
//   // for (int idxe = 0; idxe < _mesh.NumberOfElements(); ++idxe) {
//   //   for (int idxq = 0; idxq < NQ; ++idxq) {
//   //     for (int idx = 0; idx < evectors_gen.cols(); ++idx) {
//   //     }
//   //   }
//   // }
// };

// void
// sem::CalculateEigenfrequenciesSeparate(int N, double sigshift) {
//   // set number of modes
//   _num_modes_l = N;
//   evalues_gen_l.resize(_lmax);
//   evectors_gen_l.resize(_lmax);
//   vec_eigval_l.resize(_lmax);
//   vec_eigderiv_l.resize(_lmax);
//   vec_all_l = std::vector<stdvvvec>(
//       _lmax,
//       stdvvvec(_mesh.NE(), stdvvec(_mesh.NN(), stdvec(N + numaug, 0.0))));
//   vec_allderiv_l = std::vector<stdvvvec>(
//       _lmax,
//       stdvvvec(_mesh.NE(), stdvvec(_mesh.NN(), stdvec(N + numaug, 0.0))));
//   // initiate matrix multiplicaiton wrapper
//   // std::cout << "Eig Check 1\n";
//   // SparseGenMatProd<double> op_seig(mat_seig);
//   // std::cout << "Eig Check 2\n";
//   using OpType = SymShiftInvert<double, Eigen::Sparse, Eigen::Sparse>;
//   using BOpType = SparseSymMatProd<double>;
//   for (int idxl = 0; idxl < _lmax; ++idxl) {
//     OpType op_mult(mat_l_ke[idxl], mat_l_inertia[idxl]);
//     BOpType op_ke(mat_l_ke[idxl]), op_inertia(mat_l_inertia[idxl]);
//     // std::cout << "Eig Check 3\n";

//     // std::cout << "Eig Check 4\n";

//     // get nc
//     long int maxn, maxn2;
//     if ((5 * N) > mat_l_ke[idxl].rows()) {
//       maxn = mat_l_ke[idxl].rows();
//     } else {
//       maxn = 5 * N;
//     }
//     if ((3 * N) > mat_l_ke[idxl].rows()) {
//       maxn2 = mat_l_ke[idxl].rows();
//     } else {
//       maxn2 = 3 * N;
//     }

//     // initiate eigensolver
//     SymGEigsShiftSolver<OpType, BOpType, GEigsMode::ShiftInvert> eig_gen(
//         op_mult, op_inertia, N, maxn2, sigshift);
//     eig_gen.init();
//     int nconv_gen = eig_gen.compute(SortRule::LargestMagn);

//     if (eig_gen.info() == CompInfo::Successful) {
//       // std::cout << "Eigenvalues from generalised: \n"
//       //           << eig_gen.eigenvalues() * _freq_norm * _freq_norm <<
//       // "\n\n";

//       evalues_gen_l[idxl] = eig_gen.eigenvalues();
//       // std::cout << "\n\n" << evalues_gen_l[idxl] << "\n\n";
//       Eigen::MatrixXcd normval = eig_gen.eigenvectors().transpose() *
//                                  mat_l_inertia[idxl] *
//                                  eig_gen.eigenvectors();
//       // std::cout << normval << "\n\n";

//       // check if we are computing for the IC
//       if ((_has_fluid && (_il == 0)) || (!_has_fluid)) {
//         std::size_t nrow = eig_gen.eigenvectors().rows() + 1;
//         std::size_t ncol = eig_gen.eigenvectors().cols();
//         evectors_gen_l[idxl].resize(nrow, ncol);
//         evectors_gen_l[idxl].block(0, 0, 1, ncol) =
//             Eigen::MatrixXcd::Zero(1, ncol);
//         evectors_gen_l[idxl].block(1, 0, nrow - 1, ncol) =
//             eig_gen.eigenvectors();

//       } else {
//         std::size_t nrow = eig_gen.eigenvectors().rows();
//         std::size_t ncol = eig_gen.eigenvectors().cols();
//         evectors_gen_l[idxl].resize(nrow, ncol);
//         evectors_gen_l[idxl] = eig_gen.eigenvectors();
//       }

//       // normalise:
//       for (int i = 0; i < evectors_gen_l[idxl].rows(); ++i) {
//         for (int j = 0; j < evectors_gen_l[idxl].cols(); ++j) {
//           evectors_gen_l[idxl](i, j) *=
//               1.0 / (normint * sqrt(evalues_gen_l[idxl](j).real()));
//         }
//       }
//       _calc_gen = true;
//       // std::cout << "POST\n";
//     } else {
//       std::cout << "Unsuccessful generalised\n";
//     }

//     if (_calc_gen) {
//       // std::cout << "Saving eigenvectors for l = " << idxl + 1 << "\n";
//       // save the eigenvectors in the "standard format":
//       int NQ = _mesh.NN();
//       for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
//         stdvvec vec_out;
//         if ((idxe < _el) || (idxe > _eu - 1)) {
//           vec_eigval_l[idxl].push_back(
//               stdvvec(NQ, stdvec(evectors_gen_l[idxl].cols(), 0.0)));
//         } else {
//           for (int idxq = 0; idxq < NQ; ++idxq) {
//             std::vector<double> vec_x;
//             std::size_t ovidx = (idxe - _el) * (NQ - 1) + idxq;
//             for (int idx = 0; idx < evectors_gen_l[idxl].cols(); ++idx) {
//               std::size_t colidx = evectors_gen_l[idxl].cols() - 1 - idx;
//               std::complex<double> eigint = evectors_gen_l[idxl](ovidx,
//               colidx); vec_x.push_back(eigint.real());
//               vec_all_l[idxl][idxe][idxq][idx] = eigint.real();
//             }
//             vec_out.push_back(vec_x);
//           }
//           vec_eigval_l[idxl].push_back(vec_out);
//         }
//       }
//       // std::cout << "Check 2\n";
//       // calculate the derivatives of the eigenvectors:
//       vec_eigderiv_l[idxl] = stdvvvec(
//           _mesh.NE(), stdvvec(NQ, stdvec(evectors_gen_l[idxl].cols(), 0.0)));
//       for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
//         stdvvec vec_out;
//         double elem_width = _mesh.EW(idxe);
//         for (int idxq = 0; idxq < NQ; ++idxq) {
//           for (int idx = 0; idx < evectors_gen_l[idxl].cols(); ++idx) {
//             double tmp = 0.0;
//             for (int idxq2 = 0; idxq2 < NQ; ++idxq2) {
//               // std::size_t ovidx = overallidxfinal(idx, idxq2);
//               // tmp += evectors(ovidx, i) * vec_lagderiv[idxq2][idxq];
//               tmp += vec_eigval_l[idxl][idxe][idxq2][idx] *
//                      vec_lag_deriv[idxq][idxq2];
//               // tmp += vec_lag_deriv[idxq][idxq2];
//             }
//             tmp *= 2.0 / elem_width;
//             vec_eigderiv_l[idxl][idxe][idxq][idx] = tmp;
//             vec_allderiv_l[idxl][idxe][idxq][idx] = tmp;
//           }
//         }
//       }
//     }
//   }

//   // add to all
//   // vec_all_l = vec_eigval_l;
//   // vec_allderiv_l = vec_eigderiv_l;
// };

// auto
// sem::efrequencies_gen() const {
//   Eigen::VectorXd ret_eig(evalues_gen.rows());
//   for (int i = 0; i < evalues_gen.rows(); ++i) {
//     ret_eig(evalues_gen.rows() - 1 - i) =
//         std::sqrt(std::abs(evalues_gen(i))) * _freq_norm;
//   }
//   return ret_eig;
// };

// auto
// sem::efrequencies_gen_separate(int idxl) const {
//   Eigen::VectorXd ret_eig(evalues_gen_l[idxl - 1].rows());
//   for (int i = 0; i < evalues_gen_l[idxl - 1].rows(); ++i) {
//     ret_eig(evalues_gen_l[idxl - 1].rows() - 1 - i) =
//         std::sqrt(std::abs(evalues_gen_l[idxl - 1](i))) * _freq_norm;
//   }
//   return ret_eig;
// };

// auto
// sem::efunctions(int idxl) const {
//   return vec_eigval_l[idxl - 1];
// };

// Eigen::VectorXcd
// sem::CalculateForce(SourceInfo::EarthquakeCMT &cmt) {

//   // resize force
//   vec_force.resize(totlen);
//   vec_force = Eigen::VectorXcd::Zero(totlen);

//   // find element within which the source sits
//   double depth = cmt.Depth();
//   double rad_source = _mesh.PR() - 1000.0 * depth / _length_norm;
//   // std::cout << "Source radius: " << rad_source << "\n";
//   int idxsource = 0;

//   int NQ = _mesh.NN();

//   // get the y0-, y0+ values etc at the source location
//   double theta_s = (90.0 - cmt.Latitude()) * EIGEN_PI / (180.0);
//   double phi_s = cmt.Longitude() * EIGEN_PI / (180.0);

//   // wigner d matrix
//   auto wigdmat =
//       GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
//                        GSHTrans::Single, GSHTrans::ColumnMajor>(_lmax, _lmax,
//                        2,
//                                                                 theta_s);
//   auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
//     auto dl = wigdmat[N];
//     auto tmp = dl[l, m];
//     auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
//     return ylm;
//   };

//   // std::cout << "Source location: " << theta_s << " " << phi_s << "\n";
//   // int m = 0, _l = 2;
//   // parameters
//   double invsqrt2 = 1.0 / std::sqrt(2.0);
//   std::complex<double> isq2 = std::complex<double>(0.0, invsqrt2);

//   // loop through the elements
//   //   to find the element that contains the source
//   for (int idx = _el; idx < _eu; ++idx) {
//     // std::cout << _mesh.ELR(idx) << " " << _mesh.EUR(idx) << " " <<
//     // rad_source
//     //           << "\n";
//     if ((_mesh.ELR(idx) <= rad_source) && (_mesh.EUR(idx) >= rad_source)) {
//       // std::cout << idx << " " << rad_source << " " << _mesh.ELR(idx) << "
//       "
//       //           << _mesh.EUR(idx) << "\n";
//       std::vector<double> vec_nodes(NQ, 0.0);
//       for (int idxn = 0; idxn < NQ; ++idxn) {
//         vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
//       }
//       // std::cout << "Check 1\n";
//       auto pleg =
//           Interpolation::LagrangePolynomial(vec_nodes.begin(),
//           vec_nodes.end());
//       // std::cout << "Check 2\n";
//       for (int idxq = 0; idxq < NQ; ++idxq) {
//         auto w_val = pleg(idxq, rad_source) / rad_source;
//         auto w_prefactor = pleg.Derivative(idxq, rad_source) -
//                            w_val;   // prefactor for first term
//         for (int idxl = 1; idxl < _lmax + 1; ++idxl) {
//           double omegal2 = (idxl + 2) * (idxl - 1) / 2.0;
//           double lprefac =
//               std::exp(-2.0 * 3.141592653589793 * (idxl + 1) / (1 + 0.5));
//           // lprefac = 1.0;
//           for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
//             // std::cout << "Check 3\n";
//             // spherical harmonic
//             Complex ymc = std::conj(ylmn(idxl, idxm, -1, phi_s));
//             // std::cout << "Check 4\n";
//             Complex ypc = std::conj(ylmn(idxl, idxm, 1, phi_s));
//             // std::cout << "Check 5\n";
//             Complex ymmc = 0.0, yppc = 0.0;
//             if (idxl > 1) {
//               ymmc = std::conj(ylmn(idxl, idxm, -2, phi_s));
//               // std::cout << "Check 6\n";
//               yppc = std::conj(ylmn(idxl, idxm, 2, phi_s));
//               //  std::cout << "Check 7\n";
//             }

//             // get the index of the spherical harmonic
//             std::size_t ovidx = mat_indices.mat_index(idx, idxq, idxl, idxm);

//             // std::cout << "Size: " << vec_force.rows() << ", ovidx: " <<
//             // ovidx
//             //           << "\n";
//             Complex tmp = w_prefactor * (cmt.MC0m() * ymc - cmt.MC0p() *
//             ypc); tmp += w_val * omegal2 * (cmt.MCmm() * ymmc - cmt.MCpp() *
//             yppc); tmp *= isq2;

//             vec_force(ovidx) = tmp * lprefac;
//           };
//         };
//       };
//     };
//   };
//   return vec_force;
// };

// Eigen::MatrixXcd
// sem::CalculateForce(SourceInfo::EarthquakeCMT &cmt, int idxl) {

//   std::size_t flen = mat_indices.mat_size_l();
//   std::size_t fcols = 2 * idxl + 1;
//   Eigen::MatrixXcd vec_lforce = Eigen::MatrixXcd::Zero(flen, fcols);

//   // resize force
//   // vec_force.resize(totlen);
//   // vec_force = Eigen::VectorXcd::Zero(totlen);
//   // auto ylmn = [](int l, int m, int N, double theta, double phi) {
//   //   auto wigtemp = GSHTrans::Wigner(l, l, N, theta);   // Wigner D-matrix
//   //   auto dl = wigtemp(l);
//   //   auto tmp = dl(m);
//   //   auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
//   //   return ylm;
//   // };

//   // find element within which the source sits
//   double depth = cmt.Depth();
//   double rad_source = _mesh.PR() - 1000.0 * depth / _length_norm;
//   // std::cout << "\nSource radius: " << rad_source << "\n";
//   int idxsource = 0;

//   int NQ = _mesh.NN();

//   // get the y0-, y0+ values etc at the source location
//   double theta_s = (90.0 - cmt.Latitude()) * EIGEN_PI / (180.0);
//   double phi_s = cmt.Longitude() * EIGEN_PI / (180.0);

//   auto wigdmat =
//       GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
//                        GSHTrans::Single, GSHTrans::ColumnMajor>(_lmax, _lmax,
//                        2,
//                                                                 theta_s);
//   auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
//     auto dl = wigdmat[N];
//     auto tmp = dl[l, m];
//     auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
//     return ylm;
//   };
//   // std::cout << "\nSource location: " << theta_s << " " << phi_s << "\n";
//   // int m = 0, _l = 2;
//   // parameters
//   double invsqrt2 = 1.0 / std::sqrt(2.0);
//   std::complex<double> isq2 = std::complex<double>(0.0, invsqrt2);
//   double omegal2 = std::sqrt((idxl + 2) * (idxl - 1) / 2.0);
//   double lprefac = std::exp(-2.0 * 3.141592653589793 * (idxl + 1) / (1 +
//   0.5));
//   // lprefac = 1.0;
//   // loop through the elements
//   //   to find the element that contains the source
//   for (int idx = _el; idx < _eu; ++idx) {
//     // std::cout << _mesh.ELR(idx) << " " << _mesh.EUR(idx) << " " <<
//     // rad_source
//     //           << "\n";
//     if ((_mesh.ELR(idx) <= rad_source) && (_mesh.EUR(idx) >= rad_source)) {
//       // std::cout << idx << " " << rad_source << " " << _mesh.ELR(idx) << "
//       "
//       //           << _mesh.EUR(idx) << "\n";
//       std::vector<double> vec_nodes(NQ, 0.0);
//       for (int idxn = 0; idxn < NQ; ++idxn) {
//         vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
//       }

//       auto pleg =
//           Interpolation::LagrangePolynomial(vec_nodes.begin(),
//           vec_nodes.end());

//       for (int idxq = 0; idxq < NQ; ++idxq) {
//         auto w_val = pleg(idxq, rad_source) / rad_source;
//         auto w_prefactor = pleg.Derivative(idxq, rad_source) -
//                            w_val;   // prefactor for first term
//         // for (int idxl = 1; idxl < _lmax + 1; ++idxl) {

//         // get the index of the spherical harmonic
//         std::size_t ridx = mat_indices.mat_index_l(idx, idxq);
//         for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
//           // spherical harmonic
//           Complex ymc = std::conj(ylmn(idxl, idxm, -1, phi_s));
//           Complex ypc = std::conj(ylmn(idxl, idxm, 1, phi_s));
//           Complex ymmc = 0.0, yppc = 0.0;

//           if (idxl > 1) {
//             ymmc = std::conj(ylmn(idxl, idxm, -2, phi_s));
//             yppc = std::conj(ylmn(idxl, idxm, 2, phi_s));
//           }

//           // std::cout << "Size: " << vec_force.rows() << ", ovidx: " <<
//           ovidx
//           //           << "\n";
//           Complex tmp = w_prefactor * (cmt.MC0m() * ymc - cmt.MC0p() * ypc);
//           tmp += w_val * omegal2 * (cmt.MCmm() * ymmc - cmt.MCpp() * yppc);
//           tmp *= isq2;

//           vec_lforce(ridx, idxm + idxl) = tmp * lprefac;
//         };
//         // };
//       };
//     };
//   };
//   return vec_lforce;
// };

// Eigen::VectorXcd
// sem::ReceiverVectorTheta(double rad_r, double theta_r, double phi_r) {
//   // create the receiver vector
//   Eigen::VectorXcd vec_receiver(totlen);
//   vec_receiver = Eigen::VectorXcd::Zero(totlen);
//   auto wigdmat =
//       GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
//                        GSHTrans::Single, GSHTrans::ColumnMajor>(_lmax, _lmax,
//                        2,
//                                                                 theta_r);
//   auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
//     auto dl = wigdmat[N];
//     auto tmp = dl[l, m];
//     auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
//     return ylm;
//   };
//   Complex i1 = std::complex<double>(0.0, 1.0);
//   // loop through the elements
//   for (int idx = _el; idx < _eu; ++idx) {
//     if ((_mesh.ELR(idx) <= rad_r) && (_mesh.EUR(idx) >= rad_r)) {
//       std::vector<double> vec_nodes(_mesh.NN(), 0.0);
//       for (int idxn = 0; idxn < _mesh.NN(); ++idxn) {
//         vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
//       }
//       auto pleg =
//           Interpolation::LagrangePolynomial(vec_nodes.begin(),
//           vec_nodes.end());
//       for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
//         for (int idxl = 1; idxl < _lmax + 1; ++idxl) {
//           for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
//             Complex ylm = ylmn(idxl, idxm, -1, phi_r);
//             Complex ylp = ylmn(idxl, idxm, 1, phi_r);
//             std::size_t ovidx = mat_indices.mat_index(idx, idxq, idxl, idxm);
//             vec_receiver(ovidx) = -i1 * pleg(idxq, rad_r) * (ylm + ylp)
//             / 2.0;
//           }
//         }
//       }
//     }
//   }
//   return vec_receiver;
// };

// Eigen::VectorXcd
// sem::ReceiverVectorThetaSurface(double rad_r, double theta_r, double phi_r) {
//   // create the receiver vector
//   Eigen::VectorXcd vec_receiver(totlen);
//   vec_receiver = Eigen::VectorXcd::Zero(totlen);
//   auto wigdmat =
//       GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
//                        GSHTrans::Single, GSHTrans::ColumnMajor>(_lmax, _lmax,
//                        2,
//                                                                 theta_r);
//   auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
//     auto dl = wigdmat[N];
//     auto tmp = dl[l, m];
//     auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
//     return ylm;
//   };
//   Complex i1 = std::complex<double>(0.0, 1.0);

//   int idx = _eu - 1;
//   int idxq = _mesh.NN() - 1;
//   for (int idxl = 1; idxl < _lmax + 1; ++idxl) {
//     for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
//       Complex ylm = ylmn(idxl, idxm, -1, phi_r);
//       Complex ylp = ylmn(idxl, idxm, 1, phi_r);
//       std::size_t ovidx = mat_indices.mat_index(idx, idxq, idxl, idxm);
//       vec_receiver(ovidx) = -i1 * (ylm + ylp) / 2.0;
//     }
//   }
//   return vec_receiver;
// };

// Eigen::MatrixXcd
// sem::ReceiverVectorThetaSurfaceL(double rad_r, double theta_r, double phi_r,
//                                  int idxl) {

//   std::size_t flen = mat_indices.mat_size_l();
//   std::size_t fcols = 2 * idxl + 1;
//   // create the receiver vector
//   Eigen::MatrixXcd vec_receiver(flen, fcols);
//   vec_receiver = Eigen::MatrixXcd::Zero(flen, fcols);
//   auto wigdmat =
//       GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
//                        GSHTrans::Single, GSHTrans::ColumnMajor>(_lmax, _lmax,
//                        2,
//                                                                 theta_r);
//   auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
//     auto dl = wigdmat[N];
//     auto tmp = dl[l, m];
//     auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
//     return ylm;
//   };
//   Complex i1 = std::complex<double>(0.0, 1.0);

//   int idx = _eu - 1;
//   int idxq = _mesh.NN() - 1;
//   std::size_t ovidx = mat_indices.mat_index_l(idx, idxq);
//   for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
//     Complex ylm = ylmn(idxl, idxm, -1, phi_r);
//     Complex ylp = ylmn(idxl, idxm, 1, phi_r);
//     vec_receiver(ovidx, idxm + idxl) = -i1 * (ylm + ylp) / 2.0;
//   }
//   return vec_receiver;
// };

// Eigen::VectorXcd
// sem::ReceiverVectorPhi(double rad_r, double theta_r, double phi_r) {
//   // create the receiver vector
//   Eigen::VectorXcd vec_receiver(totlen);
//   vec_receiver = Eigen::VectorXcd::Zero(totlen);
//   auto wigdmat =
//       GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
//                        GSHTrans::Single, GSHTrans::ColumnMajor>(_lmax, _lmax,
//                        2,
//                                                                 theta_r);
//   auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
//     auto dl = wigdmat[N];
//     auto tmp = dl[l, m];
//     auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
//     return ylm;
//   };
//   Complex i1 = std::complex<double>(0.0, 1.0);
//   // loop through the elements
//   for (int idx = _el; idx < _eu; ++idx) {
//     if ((_mesh.ELR(idx) <= rad_r) && (_mesh.EUR(idx) >= rad_r)) {
//       std::vector<double> vec_nodes(_mesh.NN(), 0.0);
//       for (int idxn = 0; idxn < _mesh.NN(); ++idxn) {
//         vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
//       }
//       auto pleg =
//           Interpolation::LagrangePolynomial(vec_nodes.begin(),
//           vec_nodes.end());
//       for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
//         for (int idxl = 1; idxl < _lmax + 1; ++idxl) {
//           for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
//             Complex ylm = ylmn(idxl, idxm, -1, phi_r);
//             Complex ylp = ylmn(idxl, idxm, 1, phi_r);
//             std::size_t ovidx = mat_indices.mat_index(idx, idxq, idxl, idxm);
//             vec_receiver(ovidx) = pleg(idxq, rad_r) * (ylp - ylm) / 2.0;
//           }
//         }
//       }
//     }
//   }
//   return vec_receiver;
// };

// Eigen::VectorXcd
// sem::ReceiverVectorPhiSurface(double rad_r, double theta_r, double phi_r) {
//   // create the receiver vector
//   Eigen::VectorXcd vec_receiver(totlen);
//   vec_receiver = Eigen::VectorXcd::Zero(totlen);
//   auto wigdmat =
//       GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
//                        GSHTrans::Single, GSHTrans::ColumnMajor>(_lmax, _lmax,
//                        2,
//                                                                 theta_r);
//   auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
//     auto dl = wigdmat[N];
//     auto tmp = dl[l, m];
//     auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
//     return ylm;
//   };
//   Complex i1 = std::complex<double>(0.0, 1.0);

//   int idx = _eu - 1;
//   int idxq = _mesh.NN() - 1;
//   for (int idxl = 1; idxl < _lmax + 1; ++idxl) {
//     for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
//       Complex ylm = ylmn(idxl, idxm, -1, phi_r);
//       Complex ylp = ylmn(idxl, idxm, 1, phi_r);
//       std::size_t ovidx = mat_indices.mat_index(idx, idxq, idxl, idxm);
//       vec_receiver(ovidx) = (ylp - ylm) / 2.0;
//     }
//   }
//   return vec_receiver;
// };

// Eigen::MatrixXcd
// sem::ReceiverVectorPhiSurfaceL(double rad_r, double theta_r, double phi_r,
//                                int idxl) {
//   std::size_t flen = mat_indices.mat_size_l();
//   std::size_t fcols = 2 * idxl + 1;
//   // create the receiver vector
//   Eigen::MatrixXcd vec_receiver(flen, fcols);
//   vec_receiver = Eigen::MatrixXcd::Zero(flen, fcols);
//   auto wigdmat =
//       GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
//                        GSHTrans::Single, GSHTrans::ColumnMajor>(_lmax, _lmax,
//                        2,
//                                                                 theta_r);
//   auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
//     auto dl = wigdmat[N];
//     auto tmp = dl[l, m];
//     auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
//     return ylm;
//   };
//   Complex i1 = std::complex<double>(0.0, 1.0);

//   int idx = _eu - 1;
//   int idxq = _mesh.NN() - 1;
//   std::size_t ovidx = mat_indices.mat_index_l(idx, idxq);
//   // for (int idxl = 1; idxl < _lmax + 1; ++idxl) {
//   for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
//     Complex ylm = ylmn(idxl, idxm, -1, phi_r);
//     Complex ylp = ylmn(idxl, idxm, 1, phi_r);

//     vec_receiver(ovidx, idxm + idxl) = (ylp - ylm) / 2.0;
//   }
//   // }
//   return vec_receiver;
// };

// auto
// sem::PrintModesUpToFreq(double fmax) const {
//   for (int idxl = 0; idxl < _lmax; ++idxl) {
//     Eigen::VectorXcd eigvals = evalues_gen_l[idxl];
//     for (int i = 0; i < eigvals.size(); ++i) {
//       double freq = std::sqrt(std::abs(eigvals(eigvals.size() - 1 - i))) *
//                     _freq_norm / (2.0 * EIGEN_PI);
//       if (freq <= fmax) {
//         std::cout << "l = " << (idxl + 1) << " , overtone = " << i
//                   << ", freq = " << freq << "\n";
//       }
//     }
//   }
// };
// // ...existing code...

// auto
// sem::FindModesForCoupling(double fmax) {
//   for (int idxl = 0; idxl < _lmax; ++idxl) {
//     Eigen::VectorXcd eigvals = evalues_gen_l[idxl];
//     for (int i = 0; i < eigvals.size(); ++i) {
//       double freq = std::sqrt(std::abs(eigvals(eigvals.size() - 1 - i))) *
//                     _freq_norm / (2.0 * EIGEN_PI);
//       if (freq <= fmax) {
//         if (!((idxl == 0) && (i == 0))) {
//           nmodes += 2 * (idxl + 1) + 1;   // number of modes for this l
//           vec_store.push_back(TStore(idxl + 1, i, freq));
//         }
//       }
//     }
//   }
//   std::sort(vec_store.begin(), vec_store.end(),
//             [](const TStore &a, const TStore &b) {
//               return a.GetFreq() < b.GetFreq();
//             });
//   for (auto idx : vec_store) {
//     // std::cout << "l = " << idx.GetL() << ", n = " << idx.GetN()
//     //           << ", freq = " << idx.GetFreq() << "\n";
//   }
//   vec_indices.push_back(0);   // first index is 0
//   for (int idx = 0; idx < vec_store.size() - 1; ++idx) {
//     int idxl = vec_store[idx].GetL();
//     int tmp = vec_indices[idx] + (2 * idxl + 1);
//     vec_indices.push_back(tmp);   // number of modes for this l
//   }
// };

// // function that gets the matrices for modes up to a certain frequency
// template <class model1d>
// Eigen::MatrixXcd
// sem::NMC_INERTIA(const model1d &mod_new, bool aug_inertia) {

//   Eigen::MatrixXcd retmat;

//   // total size
//   auto totnum = nmodes;
//   if (aug_inertia) {
//     totnum += ((_lmax + 1) * (_lmax + 1) - 1) * numaug;
//   }

//   // resize the matrix
//   retmat.resize(totnum, totnum);
//   retmat.setZero();

//   // fill the matrix
//   // inertia matrix
//   std::size_t ovidx = 0;
//   auto q = _mesh.GLL();
//   // couple normal modes together
//   for (int is = 0; is < vec_store.size(); ++is) {
//     int idxl = vec_store[is].GetL();
//     int idxn = vec_store[is].GetN();
//     for (int js = 0; js < vec_store.size(); ++js) {
//       int idxl2 = vec_store[js].GetL();
//       int idxn2 = vec_store[js].GetN();

//       // if the l values are the same, we can fill the matrix
//       if (idxl == idxl2) {
//         double tmp = 0.0;
//         for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
//           double elem_width = _mesh.EW(idxe);
//           double d_val = 2.0 / elem_width;
//           int laynum = _mesh.LayerNumber(idxe);
//           double tmp1 = 0.0;
//           for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
//             double cr = _mesh.NodeRadius(idxe, idxq);
//             double tmp2 = vec_eigval_l[idxl - 1][idxe][idxq][idxn] *
//                           vec_eigval_l[idxl - 1][idxe][idxq][idxn2];
//             tmp2 *= mod_new.Density(laynum)(cr) * cr * cr;
//             tmp1 += q.W(idxq) * tmp2;
//           }
//           tmp += tmp1 * elem_width * 0.5;
//         }

//         // go through m values
//         for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
//           int idx1 = vec_indices[is] + idxm + idxl;
//           int idx2 = vec_indices[js] + idxm + idxl;
//           if (idx1 != idx2) {
//             retmat(idx1, idx2) = tmp;
//             retmat(idx2, idx1) = tmp;
//           } else {
//             retmat(idx1, idx2) = tmp;
//           }
//         }
//       }
//     }
//   }

//   auto idxouter = vec_indices.back() + (2 * vec_store.back().GetL() + 1);
//   if (aug_inertia) {   // cross augmentation with modes
//     for (int idxl = 0; idxl < _lmax; ++idxl) {
//       for (int idxn = 0; idxn < numaug; ++idxn) {
//         auto lval = idxl + 1;
//         for (int is = 0; is < vec_store.size(); ++is) {
//           if (lval == vec_store[is].GetL()) {
//             double tmp = 0.0;
//             for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
//               double elem_width = _mesh.EW(idxe);
//               double d_val = 2.0 / elem_width;
//               int laynum = _mesh.LayerNumber(idxe);
//               double tmp1 = 0.0;
//               for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
//                 double cr = _mesh.NodeRadius(idxe, idxq);
//                 double tmp2 =
//                     vec_augval_l[idxl][idxe][idxq][idxn] *
//                     vec_eigval_l[idxl][idxe][idxq][vec_store[is].GetN()];
//                 tmp2 *= mod_new.Density(laynum)(cr) * cr * cr;
//                 tmp1 += q.W(idxq) * tmp2;
//               }
//               tmp += tmp1 * elem_width * 0.5;
//             }

//             // go through m values
//             for (int idxm = -lval; idxm < lval + 1; ++idxm) {
//               int idx1 = vec_indices[is] + idxm + lval;
//               int idx2 = idxouter;
//               idx2 += (lval * lval - 1) * numaug +
//                       (2 * lval + 1) * idxn;   // offset for previous l
//                       values
//               idx2 += idxm + lval;
//               if (idx1 != idx2) {
//                 retmat(idx1, idx2) = tmp;
//                 retmat(idx2, idx1) = tmp;
//               } else {
//                 retmat(idx1, idx2) = tmp;
//               }
//             }
//           }
//         }
//       }
//     }

//     // cross augmentation functions against each other
//     for (int idxl = 0; idxl < _lmax; ++idxl) {
//       for (int idxn = 0; idxn < numaug; ++idxn) {
//         for (int idxn2 = 0; idxn2 < numaug; ++idxn2) {
//           auto lval = idxl + 1;

//           double tmp = 0.0;
//           for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
//             double elem_width = _mesh.EW(idxe);
//             double d_val = 2.0 / elem_width;
//             int laynum = _mesh.LayerNumber(idxe);
//             double tmp1 = 0.0;
//             for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
//               double cr = _mesh.NodeRadius(idxe, idxq);
//               double tmp2 = vec_augval_l[idxl][idxe][idxq][idxn] *
//                             vec_augval_l[idxl][idxe][idxq][idxn2];
//               tmp2 *= mod_new.Density(laynum)(cr) * cr * cr;
//               tmp1 += q.W(idxq) * tmp2;
//             }
//             tmp += tmp1 * elem_width * 0.5;
//           }

//           // go through m values
//           for (int idxm = -lval; idxm < lval + 1; ++idxm) {
//             int idx1 = idxouter;
//             idx1 += (lval * lval - 1) * numaug +
//                     (2 * lval + 1) * idxn;   // offset for previous l values
//             idx1 += idxm + lval;
//             int idx2 = idxouter;
//             idx2 += (lval * lval - 1) * numaug +
//                     (2 * lval + 1) * idxn2;   // offset for previous l values
//             idx2 += idxm + lval;
//             if (idx1 != idx2) {
//               retmat(idx1, idx2) = tmp;
//               retmat(idx2, idx1) = tmp;
//             } else {
//               retmat(idx1, idx2) = tmp;
//             }
//           }
//         }
//       }
//     }
//   }

//   return retmat;
// };

// // function that gets the matrices for modes up to a certain frequency
// template <class model1d>
// Eigen::MatrixXcd
// sem::NMC_KE(const model1d &mod_new, bool aug_ke) {

//   Eigen::MatrixXcd retmat;

//   // total size
//   auto totnum = nmodes;
//   if (aug_ke) {
//     totnum += ((_lmax + 1) * (_lmax + 1) - 1) * numaug;
//   }

//   // resize the matrix
//   retmat.resize(totnum, totnum);
//   retmat.setZero();

//   // fill the matrix
//   // ke matrix
//   std::size_t ovidx = 0;
//   auto q = _mesh.GLL();
//   for (int is = 0; is < vec_store.size(); ++is) {
//     int idxl = vec_store[is].GetL();
//     int idxn = vec_store[is].GetN();
//     for (int js = 0; js < vec_store.size(); ++js) {
//       int idxl2 = vec_store[js].GetL();
//       int idxn2 = vec_store[js].GetN();

//       // if the l values are the same, we can fill the matrix
//       if (idxl == idxl2) {
//         double tmp = 0.0;
//         auto k2 = (idxl * (idxl + 1));
//         for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
//           double elem_width = _mesh.EW(idxe);
//           double d_val = 2.0 / elem_width;
//           int laynum = _mesh.LayerNumber(idxe);
//           double tmp1 = 0.0;
//           for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
//             double cr = _mesh.NodeRadius(idxe, idxq);

//             // first integral
//             double tmp2 = cr * vec_eigderiv_l[idxl - 1][idxe][idxq][idxn] -
//                           vec_eigval_l[idxl - 1][idxe][idxq][idxn];
//             double tmp3 = cr * vec_eigderiv_l[idxl - 1][idxe][idxq][idxn2] -
//                           vec_eigval_l[idxl - 1][idxe][idxq][idxn2];
//             double tmp4 = tmp2 * tmp3 * mod_new.L(laynum)(cr);

//             // second integral
//             double tmp5 = (k2 - 2.0) *
//                           vec_eigval_l[idxl - 1][idxe][idxq][idxn] *
//                           vec_eigval_l[idxl - 1][idxe][idxq][idxn2];
//             tmp5 *= mod_new.N(laynum)(cr);

//             // sum
//             tmp1 += q.W(idxq) * (tmp4 + tmp5);
//           }
//           tmp += tmp1 * elem_width * 0.5;
//         }

//         // go through m values
//         for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
//           int idx1 = vec_indices[is] + idxm + idxl;
//           int idx2 = vec_indices[js] + idxm + idxl;
//           if (idx1 != idx2) {
//             retmat(idx1, idx2) = tmp;
//             retmat(idx2, idx1) = tmp;
//           } else {
//             retmat(idx1, idx2) = tmp;
//           }
//         }
//       }
//     }
//   }

//   auto idxouter = vec_indices.back() + (2 * vec_store.back().GetL() + 1);
//   if (aug_ke) {   // cross augmentation with modes
//     for (int idxl = 0; idxl < _lmax; ++idxl) {
//       for (int idxn = 0; idxn < numaug; ++idxn) {
//         auto lval = idxl + 1;
//         for (int is = 0; is < vec_store.size(); ++is) {
//           int idxn2 = vec_store[is].GetN();
//           if (lval == vec_store[is].GetL()) {

//             double tmp = 0.0;
//             auto k2 = (lval * (lval + 1));
//             for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
//               double elem_width = _mesh.EW(idxe);
//               double d_val = 2.0 / elem_width;
//               int laynum = _mesh.LayerNumber(idxe);
//               double tmp1 = 0.0;
//               for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
//                 double cr = _mesh.NodeRadius(idxe, idxq);

//                 // first integral
//                 double tmp2 = cr * vec_augderiv_l[idxl][idxe][idxq][idxn] -
//                               vec_augval_l[idxl][idxe][idxq][idxn];
//                 double tmp3 = cr * vec_eigderiv_l[idxl][idxe][idxq][idxn2] -
//                               vec_eigval_l[idxl][idxe][idxq][idxn2];
//                 double tmp4 = tmp2 * tmp3 * mod_new.L(laynum)(cr);

//                 // second integral
//                 double tmp5 = (k2 - 2.0) *
//                               vec_augval_l[idxl][idxe][idxq][idxn] *
//                               vec_eigval_l[idxl][idxe][idxq][idxn2];
//                 tmp5 *= mod_new.N(laynum)(cr);

//                 // sum
//                 tmp1 += q.W(idxq) * (tmp4 + tmp5);
//               }
//               tmp += tmp1 * elem_width * 0.5;
//             }

//             // go through m values
//             for (int idxm = -lval; idxm < lval + 1; ++idxm) {
//               int idx1 = vec_indices[is] + idxm + lval;
//               int idx2 = idxouter;
//               idx2 += (lval * lval - 1) * numaug +
//                       (2 * lval + 1) * idxn;   // offset for previous l
//                       values
//               idx2 += idxm + lval;
//               if (idx1 != idx2) {
//                 retmat(idx1, idx2) = tmp;
//                 retmat(idx2, idx1) = tmp;
//               } else {
//                 retmat(idx1, idx2) = tmp;
//               }
//             }
//           }
//         }
//       }
//     }

//     // cross augmentation functions against each other
//     for (int idxl = 0; idxl < _lmax; ++idxl) {
//       for (int idxn = 0; idxn < numaug; ++idxn) {
//         for (int idxn2 = 0; idxn2 < numaug; ++idxn2) {
//           auto lval = idxl + 1;

//           double tmp = 0.0;
//           auto k2 = (lval * (lval + 1));
//           for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
//             double elem_width = _mesh.EW(idxe);
//             double d_val = 2.0 / elem_width;
//             int laynum = _mesh.LayerNumber(idxe);
//             double tmp1 = 0.0;
//             for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
//               double cr = _mesh.NodeRadius(idxe, idxq);

//               // first integral
//               double tmp2 = cr * vec_augderiv_l[idxl][idxe][idxq][idxn] -
//                             vec_augval_l[idxl][idxe][idxq][idxn];
//               double tmp3 = cr * vec_augderiv_l[idxl][idxe][idxq][idxn2] -
//                             vec_augval_l[idxl][idxe][idxq][idxn2];
//               double tmp4 = tmp2 * tmp3 * mod_new.L(laynum)(cr);

//               // second integral
//               double tmp5 = (k2 - 2.0) * vec_augval_l[idxl][idxe][idxq][idxn]
//               *
//                             vec_augval_l[idxl][idxe][idxq][idxn2];
//               tmp5 *= mod_new.N(laynum)(cr);

//               // sum
//               tmp1 += q.W(idxq) * (tmp4 + tmp5);
//             }
//             tmp += tmp1 * elem_width * 0.5;
//           }

//           // go through m values
//           for (int idxm = -lval; idxm < lval + 1; ++idxm) {
//             int idx1 = idxouter;
//             idx1 += (lval * lval - 1) * numaug +
//                     (2 * lval + 1) * idxn;   // offset for previous l values
//             idx1 += idxm + lval;
//             int idx2 = idxouter;
//             idx2 += (lval * lval - 1) * numaug +
//                     (2 * lval + 1) * idxn2;   // offset for previous l values
//             idx2 += idxm + lval;
//             if (idx1 != idx2) {
//               retmat(idx1, idx2) = tmp;
//               retmat(idx2, idx1) = tmp;
//             } else {
//               retmat(idx1, idx2) = tmp;
//             }
//           }
//         }
//       }
//     }
//   }

//   return retmat;
// };

// auto
// sem::NMC_FORCE(SourceInfo::EarthquakeCMT &cmt, bool aug_force) {

//   auto totnum = nmodes;
//   if (aug_force) {
//     totnum += ((_lmax + 1) * (_lmax + 1) - 1) * numaug;
//   }

//   Eigen::VectorXcd vec_force = Eigen::VectorXcd::Zero(totnum);

//   // resize force
//   // vec_force.resize(totlen);
//   // vec_force = Eigen::VectorXcd::Zero(totlen);

//   // find element within which the source sits
//   double depth = cmt.Depth();
//   double rad_source = _mesh.PR() - 1000.0 * depth / _length_norm;
//   // std::cout << "Source radius: " << rad_source << "\n";
//   int idxsource = 0;

//   int NQ = _mesh.NN();

//   // get the y0-, y0+ values etc at the source location
//   double theta_s = (90.0 - cmt.Latitude()) * EIGEN_PI / (180.0);
//   double phi_s = cmt.Longitude() * EIGEN_PI / (180.0);

//   auto wigdmat =
//       GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
//                        GSHTrans::Single, GSHTrans::ColumnMajor>(_lmax, _lmax,
//                        2,
//                                                                 theta_s);
//   auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
//     auto dl = wigdmat[N];
//     auto tmp = dl[l, m];
//     auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
//     return ylm;
//   };
//   // std::cout << "Source location: " << theta_s << " " << phi_s << "\n";
//   // int m = 0, _l = 2;
//   // parameters
//   double invsqrt2 = 1.0 / std::sqrt(2.0);
//   std::complex<double> isq2 = std::complex<double>(0.0, invsqrt2);

//   // loop through the elements
//   //   to find the element that contains the source
//   for (int is = 0; is < vec_store.size(); ++is) {
//     int idxl = vec_store[is].GetL();
//     int idxn = vec_store[is].GetN();
//     for (int idx = _el; idx < _eu; ++idx) {
//       // std::cout << _mesh.ELR(idx) << " " << _mesh.EUR(idx) << " " <<
//       // rad_source
//       //           << "\n";
//       if ((_mesh.ELR(idx) <= rad_source) && (_mesh.EUR(idx) >= rad_source)) {
//         // std::cout << idx << " " << rad_source << " " << _mesh.ELR(idx) <<
//         "
//         // "
//         //           << _mesh.EUR(idx) << "\n";
//         std::vector<double> vec_nodes(NQ, 0.0);
//         for (int idxv = 0; idxv < NQ; ++idxv) {
//           vec_nodes[idxv] = _mesh.NodeRadius(idx, idxv);
//         }

//         auto pleg = Interpolation::LagrangePolynomial(vec_nodes.begin(),
//                                                       vec_nodes.end());

//         auto w_val = 0.0;
//         auto w_deriv = 0.0;
//         for (int idxq = 0; idxq < NQ; ++idxq) {
//           w_val +=
//               vec_eigval_l[idxl - 1][idx][idxq][idxn] * pleg(idxq,
//               rad_source);
//           w_deriv += vec_eigval_l[idxl - 1][idx][idxq][idxn] *
//                      pleg.Derivative(idxq, rad_source);
//         };
//         auto w_vald = w_val / rad_source;
//         auto w_prefactor = w_deriv - w_vald;
//         double omegal2 = (idxl + 2) * (idxl - 1) / 2.0;
//         double lprefac =
//             std::exp(-2.0 * 3.141592653589793 * (idxl + 1) / (1 + 0.5));
//         // lprefac = 1.0;
//         for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
//           // spherical harmonic
//           Complex ymc = std::conj(ylmn(idxl, idxm, -1, phi_s));
//           Complex ypc = std::conj(ylmn(idxl, idxm, 1, phi_s));
//           Complex ymmc = 0.0, yppc = 0.0;

//           if (idxl > 1) {
//             ymmc = std::conj(ylmn(idxl, idxm, -2, phi_s));
//             yppc = std::conj(ylmn(idxl, idxm, 2, phi_s));
//           }

//           Complex tmp = w_prefactor * (cmt.MC0m() * ymc - cmt.MC0p() * ypc);
//           tmp += w_vald * omegal2 * (cmt.MCmm() * ymmc - cmt.MCpp() * yppc);
//           tmp *= isq2 * lprefac;

//           vec_force(vec_indices[is] + idxm + idxl) = tmp;
//         };
//       };
//     };
//   };

//   auto idxouter = vec_indices.back() + (2 * vec_store.back().GetL() + 1);
//   if (aug_force) {   // cross augmentation with modes

//     for (int idxl = 0; idxl < _lmax; ++idxl) {
//       for (int idxn = 0; idxn < numaug; ++idxn) {
//         auto lval = idxl + 1;
//         for (int idx = _el; idx < _eu; ++idx) {
//           // std::cout << _mesh.ELR(idx) << " " << _mesh.EUR(idx) << " " <<
//           // rad_source
//           //           << "\n";
//           if ((_mesh.ELR(idx) <= rad_source) &&
//               (_mesh.EUR(idx) >= rad_source)) {
//             // std::cout << idx << " " << rad_source << " " << _mesh.ELR(idx)
//             // << "
//             // "
//             //           << _mesh.EUR(idx) << "\n";
//             std::vector<double> vec_nodes(NQ, 0.0);
//             for (int idxv = 0; idxv < NQ; ++idxv) {
//               vec_nodes[idxv] = _mesh.NodeRadius(idx, idxv);
//             }

//             auto pleg = Interpolation::LagrangePolynomial(vec_nodes.begin(),
//                                                           vec_nodes.end());

//             auto w_val = 0.0;
//             auto w_deriv = 0.0;
//             for (int idxq = 0; idxq < NQ; ++idxq) {
//               w_val +=
//                   vec_augval_l[idxl][idx][idxq][idxn] * pleg(idxq,
//                   rad_source);
//               w_deriv += vec_augval_l[idxl][idx][idxq][idxn] *
//                          pleg.Derivative(idxq, rad_source);
//             };
//             auto w_vald = w_val / rad_source;
//             auto w_prefactor = w_deriv - w_vald;
//             double omegal2 = (lval + 2) * (lval - 1) / 2.0;
//             double lprefac =
//                 std::exp(-2.0 * 3.141592653589793 * (idxl + 1) / (1 + 0.5));
//             // lprefac = 1.0;
//             for (int idxm = -lval; idxm < lval + 1; ++idxm) {
//               // std::cout << "l = " << lval << ", n = " << idxn
//               //           << ", m = " << idxm << "\n";
//               // spherical harmonic
//               Complex ymc = std::conj(ylmn(lval, idxm, -1, phi_s));
//               Complex ypc = std::conj(ylmn(lval, idxm, 1, phi_s));
//               Complex ymmc = 0.0, yppc = 0.0;

//               if (lval > 1) {
//                 ymmc = std::conj(ylmn(lval, idxm, -2, phi_s));
//                 yppc = std::conj(ylmn(lval, idxm, 2, phi_s));
//               }

//               Complex tmp = w_prefactor * (cmt.MC0m() * ymc - cmt.MC0p() *
//               ypc); tmp += w_vald * omegal2 * (cmt.MCmm() * ymmc - cmt.MCpp()
//               * yppc); tmp *= isq2;

//               int idx1 = idxouter;
//               idx1 += (lval * lval - 1) * numaug +
//                       (2 * lval + 1) * idxn;   // offset for previous l
//                       values
//               idx1 += idxm + lval;
//               vec_force(idx1) = tmp * lprefac;
//             };
//           };
//         }
//       };
//     };
//   }

//   return vec_force;
// };

// Eigen::VectorXcd
// sem::ReceiverVectorThetaSurfaceCoupling(double rad_r, double theta_r,
//                                         double phi_r, bool toaug) {
//   // std::size_t flen = mat_indices.mat_size_l();
//   // std::size_t fcols = 2 * idxl + 1;
//   // create the receiver vector
//   auto totnum = nmodes;
//   if (toaug) {
//     totnum += ((_lmax + 1) * (_lmax + 1) - 1) * numaug;
//   }
//   Eigen::VectorXcd vec_receiver(totnum);
//   vec_receiver = Eigen::VectorXcd::Zero(totnum);
//   auto wigdmat =
//       GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
//                        GSHTrans::Single, GSHTrans::ColumnMajor>(_lmax, _lmax,
//                        2,
//                                                                 theta_r);
//   auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
//     auto dl = wigdmat[N];
//     auto tmp = dl[l, m];
//     auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
//     return ylm;
//   };
//   Complex i1 = std::complex<double>(0.0, 1.0);

//   int idx = _eu - 1;
//   int idxq = _mesh.NN() - 1;
//   // std::size_t ovidx = mat_indices.mat_index_l(idx, idxq);
//   for (int is = 0; is < vec_store.size(); ++is) {
//     int idxl = vec_store[is].GetL();
//     int idxn = vec_store[is].GetN();
//     // for (int idxl = 1; idxl < _lmax + 1; ++idxl) {
//     for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
//       Complex ylm = ylmn(idxl, idxm, -1, phi_r);
//       Complex ylp = ylmn(idxl, idxm, 1, phi_r);

//       vec_receiver(vec_indices[is] + idxm + idxl) =
//           vec_eigval_l[idxl - 1][idx][idxq][idxn] * -i1 * (ylm + ylp) / 2.0;
//     }
//     // }
//   }
//   auto idxouter = vec_indices.back() + (2 * vec_store.back().GetL() + 1);
//   if (toaug) {   // cross augmentation with modes

//     for (int idxl = 0; idxl < _lmax; ++idxl) {
//       for (int idxn = 0; idxn < numaug; ++idxn) {
//         auto lval = idxl + 1;
//         for (int idxm = -lval; idxm < lval + 1; ++idxm) {
//           Complex ylm = ylmn(lval, idxm, -1, phi_r);
//           Complex ylp = ylmn(lval, idxm, 1, phi_r);

//           int idx1 = idxouter;
//           idx1 += (lval * lval - 1) * numaug +
//                   (2 * lval + 1) * idxn;   // offset for previous l values
//           idx1 += idxm + lval;
//           vec_receiver(idx1) =
//               vec_augval_l[idxl][idx][idxq][idxn] * -i1 * (ylm + ylp) / 2.0;
//         }
//       }
//     }
//   }
//   return vec_receiver;
// };

// Eigen::VectorXcd
// sem::ReceiverVectorPhiSurfaceCoupling(double rad_r, double theta_r,
//                                       double phi_r, bool toaug) {
//   // std::size_t flen = mat_indices.mat_size_l();
//   // std::size_t fcols = 2 * idxl + 1;
//   // create the receiver vector
//   auto totnum = nmodes;
//   if (toaug) {
//     totnum += ((_lmax + 1) * (_lmax + 1) - 1) * numaug;
//   }

//   Eigen::VectorXcd vec_receiver(totnum);
//   vec_receiver = Eigen::VectorXcd::Zero(totnum);
//   auto wigdmat =
//       GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
//                        GSHTrans::Single, GSHTrans::ColumnMajor>(_lmax, _lmax,
//                        2,
//                                                                 theta_r);
//   auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
//     auto dl = wigdmat[N];
//     auto tmp = dl[l, m];
//     auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
//     return ylm;
//   };
//   Complex i1 = std::complex<double>(0.0, 1.0);

//   int idx = _eu - 1;
//   int idxq = _mesh.NN() - 1;
//   // std::size_t ovidx = mat_indices.mat_index_l(idx, idxq);
//   for (int is = 0; is < vec_store.size(); ++is) {
//     int idxl = vec_store[is].GetL();
//     int idxn = vec_store[is].GetN();
//     // for (int idxl = 1; idxl < _lmax + 1; ++idxl) {
//     for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
//       Complex ylm = ylmn(idxl, idxm, -1, phi_r);
//       Complex ylp = ylmn(idxl, idxm, 1, phi_r);

//       vec_receiver(vec_indices[is] + idxm + idxl) =
//           vec_eigval_l[idxl - 1][idx][idxq][idxn] * (ylp - ylm) / 2.0;
//     }
//     // }
//   }

//   auto idxouter = vec_indices.back() + (2 * vec_store.back().GetL() + 1);
//   if (toaug) {   // cross augmentation with modes
//     for (int idxl = 0; idxl < _lmax; ++idxl) {
//       for (int idxn = 0; idxn < numaug; ++idxn) {
//         auto lval = idxl + 1;
//         for (int idxm = -lval; idxm < lval + 1; ++idxm) {
//           Complex ylm = ylmn(lval, idxm, -1, phi_r);
//           Complex ylp = ylmn(lval, idxm, 1, phi_r);

//           int idx1 = idxouter;
//           idx1 += (lval * lval - 1) * numaug +
//                   (2 * lval + 1) * idxn;   // offset for previous l values
//           idx1 += idxm + lval;
//           vec_receiver(idx1) =
//               vec_augval_l[idxl][idx][idxq][idxn] * (ylp - ylm) / 2.0;
//         }
//       }
//     }
//   }
//   return vec_receiver;
// };

// void
// sem::augment_basis_calculate() {
//   std::cout << "Calculating augmented basis\n";

//   std::size_t ncols = _mesh.NL() + 1;
//   _augbasis_l.resize(_lmax);
//   vec_augval_l = std::vector<stdvvvec>(_lmax);
//   vec_augderiv_l = std::vector<stdvvvec>(_lmax);
//   _calc_l = std::vector<bool>(_lmax, false);
//   int lowidx = 0;
//   if ((_has_fluid && (_il == 0)) || (!_has_fluid)) {
//     lowidx = 1;
//   } else {
//     lowidx = 0;
//   }
//   // std::cout << "Check 1\n";
//   // number of layers within layer
//   // int num_aug = _mesh.LayerNumber(_eu - 1) - _mesh.LayerNumber(_el) + 2;
//   // numaug = num_aug;
//   // std::cout << "\nNumber of augmentation functions: " << num_aug <<
//   "\n\n";
//   // std::cout << "Check 2\n";
//   for (int idx = 0; idx < _lmax; ++idx) {
//     // std::cout << "Check 3, l = " << (idx + 1) << "\n";
//     auto mat_keuse = mat_l_ke[idx];
//     Eigen::SparseLU<Eigen::SparseMatrix<double>> _solver;
//     _solver.analyzePattern(mat_keuse);
//     _solver.factorize(mat_keuse);
//     std::size_t nrows = lowidx + mat_keuse.rows();
//     _augbasis_l[idx] = Eigen::MatrixXd::Zero(nrows, numaug);
//     // need to account for fact that in solid we only find the result in one
//     // solid
//     // section
//     // whether or not we're looking at IC or purely solid planet

//     // std::cout << "Check 4\n";
//     // force matrix, initialise and set (0,0) element (for first boundary)
//     Eigen::MatrixXd mat_force = Eigen::MatrixXd::Zero(mat_keuse.cols(),
//     numaug);
//     // std::cout << "mat_force size: " << mat_force.rows() << " "
//     //           << mat_force.cols() << "\n\n";
//     // mat_force(0, 0) = -1;
//     int colnum = 0;
//     auto idxlow = _el;
//     if (idxlow == 0) {
//       idxlow += 1;
//     }
//     for (int idxe = idxlow; idxe < _eu; ++idxe) {
//       if ((_mesh.LayerNumber(idxe) - _mesh.LayerNumber(idxe - 1)) == 1) {
//         std::size_t ovidx;
//         ovidx = mat_indices.mat_index_l(idxe, 0);
//         // std::cout << "Setting force at element " << idxe << " row " <<
//         ovidx
//         //           << " col " << colnum << "\n";
//         mat_force(ovidx, colnum) = -1.0;
//         colnum += 1;
//       }
//     }
//     // mat_force(mat_ke.cols() - 1, colnum) = -1;
//     // std::cout << "Check 5\n";
//     // std::cout << "augbasis size: " << _augbasis_l[idx].rows() << " "
//     //           << _augbasis_l[idx].cols() << "\n";
//     // std::cout << "Block trying to make: " << lowidx << " " <<
//     // mat_keuse.cols()
//     //           << " " << num_aug << "\n";
//     // calculate basis
//     _augbasis_l[idx].block(lowidx, 0, mat_keuse.cols(), numaug) =
//         _solver.solve(mat_force);
//     _calc_l[idx] = true;
//     // Eigen::MatrixXcd normval =
//     //     _augbasis.block(lowidx, 0, mat_ke.cols(), num_aug).transpose() *
//     //     mat_inertia * _augbasis.block(lowidx, 0, mat_ke.cols(), num_aug);

//     // std::cout << std::setprecision(4) << normval << "\n";
//     // fill out augval and augderiv
//     //  save the augmented basis in the "standard format":
//     // std::cout << "Check 6\n";
//     int NQ = _mesh.NN();
//     // std::cout << "num_aug: " << num_aug << "\n";
//     // std::cout << "NQ: " << NQ << "\n";
//     // std::cout << "_mesh.NE(): " << _mesh.NE() << "\n";
//     vec_augval_l[idx] = stdvvvec(_mesh.NE(), stdvvec(NQ, stdvec(numaug,
//     0.0)));
//     // std::cout << "Check 6.1\n";
//     for (int idxe = _el; idxe < _eu; ++idxe) {
//       for (int idxq = 0; idxq < NQ; ++idxq) {
//         std::size_t ovidx = (idxe - _el) * (NQ - 1) + idxq;
//         for (int idxa = 0; idxa < numaug; ++idxa) {
//           // std::cout << "idxe, idxq, idxa: " << idxe << " " << idxq << " "
//           //           << idxa << "\n";
//           vec_augval_l[idx][idxe][idxq][idxa] = _augbasis_l[idx](ovidx,
//           idxa);
//         }
//       }
//     }
//     // calculate the derivatives of the augmented basis:
//     vec_augderiv_l[idx] =
//         stdvvvec(_mesh.NE(), stdvvec(NQ, stdvec(numaug, 0.0)));
//     // std::cout << "Check 7\n";
//     for (int idxe = _el; idxe < _eu; ++idxe) {
//       double elem_width = _mesh.EW(idxe);
//       for (int idxq = 0; idxq < NQ; ++idxq) {
//         for (int idxa = 0; idxa < numaug; ++idxa) {
//           double tmp = 0.0;
//           for (int idxq2 = 0; idxq2 < NQ; ++idxq2) {
//             tmp += vec_augval_l[idx][idxe][idxq2][idxa] *
//                    vec_lag_deriv[idxq][idxq2];
//           }
//           tmp *= 2.0 / elem_width;
//           vec_augderiv_l[idx][idxe][idxq][idxa] = tmp;
//         }
//       }
//     }
//   }

//   // we now need to add the augmentation to the eigenvector matrix
//   for (int idxl = 0; idxl < _lmax; ++idxl) {
//     for (int idxe = _el; idxe < _eu; ++idxe) {
//       int NQ = _mesh.NN();
//       for (int idxq = 0; idxq < NQ; ++idxq) {
//         for (int idxa = 0; idxa < numaug; ++idxa) {
//           vec_all_l[idxl][idxe][idxq][idxa + _num_modes_l] =
//               vec_augval_l[idxl][idxe][idxq][idxa];
//           vec_allderiv_l[idxl][idxe][idxq][idxa + _num_modes_l] =
//               vec_augderiv_l[idxl][idxe][idxq][idxa];
//         }
//       }
//     }
//   }
// };

// template <class model1d>
// Eigen::MatrixXcd
// sem::NMC_INERTIA(const model1d &mod_new, int idxl, bool toaug) {
//   Eigen::MatrixXcd retmat;

//   // total size
//   auto totnum = _num_modes_l;
//   if (toaug) {
//     totnum += numaug;
//   }

//   // resize the matrix
//   retmat.resize(totnum, totnum);
//   retmat.setZero();

//   // fill the matrix
//   // inertia matrix
//   std::size_t ovidx = 0;
//   auto q = _mesh.GLL();
//   // couple normal modes together

//   if (!toaug) {
//     for (int idxn = 0; idxn < _num_modes_l; ++idxn) {

//       for (int idxn2 = 0; idxn2 < idxn + 1; ++idxn2) {

//         double tmp = 0.0;
//         for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
//           double elem_width = _mesh.EW(idxe);
//           double d_val = 2.0 / elem_width;
//           int laynum = _mesh.LayerNumber(idxe);
//           double tmp1 = 0.0;
//           for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
//             double cr = _mesh.NodeRadius(idxe, idxq);
//             double tmp2 = vec_eigval_l[idxl - 1][idxe][idxq][idxn] *
//                           vec_eigval_l[idxl - 1][idxe][idxq][idxn2];
//             tmp2 *= mod_new.Density(laynum)(cr) * cr * cr;
//             tmp1 += q.W(idxq) * tmp2;
//           }
//           tmp += tmp1 * elem_width * 0.5;
//         }

//         // go through m values

//         if (idxn != idxn2) {
//           retmat(idxn, idxn2) = tmp;
//           retmat(idxn2, idxn) = tmp;
//         } else {
//           retmat(idxn, idxn2) = tmp;
//         }
//       }
//     }
//   } else if (toaug) {
//     for (int idxn = 0; idxn < _num_modes_l + numaug; ++idxn) {

//       for (int idxn2 = 0; idxn2 < idxn + 1; ++idxn2) {

//         double tmp = 0.0;
//         for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
//           double elem_width = _mesh.EW(idxe);
//           double d_val = 2.0 / elem_width;
//           int laynum = _mesh.LayerNumber(idxe);
//           double tmp1 = 0.0;
//           for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
//             double cr = _mesh.NodeRadius(idxe, idxq);
//             double tmp2 = vec_all_l[idxl - 1][idxe][idxq][idxn] *
//                           vec_all_l[idxl - 1][idxe][idxq][idxn2];
//             tmp2 *= mod_new.Density(laynum)(cr) * cr * cr;
//             tmp1 += q.W(idxq) * tmp2;
//           }
//           tmp += tmp1 * elem_width * 0.5;
//         }

//         // go through m values

//         if (idxn != idxn2) {
//           retmat(idxn, idxn2) = tmp;
//           retmat(idxn2, idxn) = tmp;
//         } else {
//           retmat(idxn, idxn2) = tmp;
//         }
//       }
//     }
//   }
//   return retmat;
// };

// template <class model1d>
// Eigen::MatrixXcd
// sem::NMC_KE(const model1d &mod_new, int idxl, bool toaug) {

//   // total size
//   auto totnum = _num_modes_l;
//   if (toaug) {
//     totnum += numaug;
//   }

//   // resize the matrix
//   Eigen::MatrixXcd retmat = Eigen::MatrixXcd::Zero(totnum, totnum);
//   // retmat.setZero();

//   // fill the matrix
//   // inertia matrix
//   std::size_t ovidx = 0;
//   auto q = _mesh.GLL();

//   // couple normal modes together
//   if (!toaug) {
//     std::cout << "\n Not augmenting\n\n";
//     for (int idxn = 0; idxn < _num_modes_l; ++idxn) {
//       for (int idxn2 = 0; idxn2 < idxn + 1; ++idxn2) {

//         double tmp = 0.0;
//         auto k2 = (idxl * (idxl + 1));
//         for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
//           double elem_width = _mesh.EW(idxe);
//           double d_val = 2.0 / elem_width;
//           int laynum = _mesh.LayerNumber(idxe);
//           double tmp1 = 0.0;
//           for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
//             double cr = _mesh.NodeRadius(idxe, idxq);

//             // first integral
//             double tmp2 = cr * vec_eigderiv_l[idxl - 1][idxe][idxq][idxn] -
//                           vec_eigval_l[idxl - 1][idxe][idxq][idxn];
//             double tmp3 = cr * vec_eigderiv_l[idxl - 1][idxe][idxq][idxn2] -
//                           vec_eigval_l[idxl - 1][idxe][idxq][idxn2];
//             double tmp4 = tmp2 * tmp3 * mod_new.L(laynum)(cr);

//             // second integral
//             double tmp5 = (k2 - 2.0) *
//                           vec_eigval_l[idxl - 1][idxe][idxq][idxn] *
//                           vec_eigval_l[idxl - 1][idxe][idxq][idxn2];
//             tmp5 *= mod_new.N(laynum)(cr);

//             // sum
//             tmp1 += q.W(idxq) * (tmp4 + tmp5);
//           }
//           tmp += tmp1 * elem_width * 0.5;
//         }

//         // go through m values

//         if (idxn != idxn2) {
//           retmat(idxn, idxn2) = tmp;
//           retmat(idxn2, idxn) = tmp;
//         } else {
//           retmat(idxn, idxn2) = tmp;
//         }

//         // }
//         // }
//       }
//     }
//   } else if (toaug) {
//     std::cout << "\n Augmenting\n\n";
//     for (int idxn = 0; idxn < totnum; ++idxn) {
//       for (int idxn2 = 0; idxn2 < idxn + 1; ++idxn2) {

//         double tmp = 0.0;
//         auto k2 = (idxl * (idxl + 1));
//         for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
//           double elem_width = _mesh.EW(idxe);
//           double d_val = 2.0 / elem_width;
//           int laynum = _mesh.LayerNumber(idxe);
//           double tmp1 = 0.0;
//           for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
//             double cr = _mesh.NodeRadius(idxe, idxq);

//             // first integral
//             double tmp2 = cr * vec_allderiv_l[idxl - 1][idxe][idxq][idxn] -
//                           vec_all_l[idxl - 1][idxe][idxq][idxn];
//             double tmp3 = cr * vec_allderiv_l[idxl - 1][idxe][idxq][idxn2] -
//                           vec_all_l[idxl - 1][idxe][idxq][idxn2];
//             double tmp4 = tmp2 * tmp3 * mod_new.L(laynum)(cr);

//             // second integral
//             double tmp5 = (k2 - 2.0) * vec_all_l[idxl - 1][idxe][idxq][idxn]
//             *
//                           vec_all_l[idxl - 1][idxe][idxq][idxn2];
//             tmp5 *= mod_new.N(laynum)(cr);

//             // sum
//             tmp1 += q.W(idxq) * (tmp4 + tmp5);
//           }
//           tmp += tmp1 * elem_width * 0.5;
//         }

//         // go through m values

//         if (idxn != idxn2) {
//           retmat(idxn, idxn2) = tmp;
//           retmat(idxn2, idxn) = tmp;
//         } else {
//           retmat(idxn, idxn2) = tmp;
//         }

//         // }
//         // }
//       }
//     }
//   }
//   return retmat;
// };

// Eigen::MatrixXcd
// sem::ReceiverVectorThetaSurface_NMCL(double rad_r, double theta_r, double
// phi_r,
//                                      int idxl, bool toaug) {

//   auto totnum = _num_modes_l;
//   if (toaug) {
//     totnum += numaug;
//   }
//   std::size_t flen = _num_modes_l;
//   std::size_t fcols = 2 * idxl + 1;
//   // create the receiver vector
//   Eigen::MatrixXcd vec_receiver(totnum, fcols);
//   vec_receiver = Eigen::MatrixXcd::Zero(totnum, fcols);
//   auto wigdmat =
//       GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
//                        GSHTrans::Single, GSHTrans::ColumnMajor>(_lmax, _lmax,
//                        2,
//                                                                 theta_r);
//   auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
//     auto dl = wigdmat[N];
//     auto tmp = dl[l, m];
//     auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
//     return ylm;
//   };
//   Complex i1 = std::complex<double>(0.0, 1.0);

//   int idx = _eu - 1;
//   int idxq = _mesh.NN() - 1;
//   if (!toaug) {
//     for (int idxn = 0; idxn < _num_modes_l; ++idxn) {
//       for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
//         Complex ylm = ylmn(idxl, idxm, -1, phi_r);
//         Complex ylp = ylmn(idxl, idxm, 1, phi_r);
//         vec_receiver(idxn, idxm + idxl) =
//             vec_eigval_l[idxl - 1][idx][idxq][idxn] * -i1 * (ylm + ylp)
//             / 2.0;
//       }
//     }
//   } else if (toaug) {
//     for (int idxn = 0; idxn < totnum; ++idxn) {
//       for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
//         Complex ylm = ylmn(idxl, idxm, -1, phi_r);
//         Complex ylp = ylmn(idxl, idxm, 1, phi_r);
//         vec_receiver(idxn, idxm + idxl) =
//             vec_all_l[idxl - 1][idx][idxq][idxn] * -i1 * (ylm + ylp) / 2.0;
//       }
//     }
//   }
//   return vec_receiver;
// };

// Eigen::MatrixXcd
// sem::ReceiverVectorPhiSurface_NMCL(double rad_r, double theta_r, double
// phi_r,
//                                    int idxl, bool toaug) {
//   auto totnum = _num_modes_l;
//   if (toaug) {
//     totnum += numaug;
//   }
//   // std::size_t flen = _num_modes_l;
//   std::size_t fcols = 2 * idxl + 1;
//   // create the receiver vector
//   Eigen::MatrixXcd vec_receiver(totnum, fcols);
//   vec_receiver = Eigen::MatrixXcd::Zero(totnum, fcols);
//   auto wigdmat =
//       GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
//                        GSHTrans::Single, GSHTrans::ColumnMajor>(_lmax, _lmax,
//                        2,
//                                                                 theta_r);
//   auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
//     auto dl = wigdmat[N];
//     auto tmp = dl[l, m];
//     auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
//     return ylm;
//   };
//   Complex i1 = std::complex<double>(0.0, 1.0);

//   int idx = _eu - 1;
//   int idxq = _mesh.NN() - 1;
//   if (!toaug) {
//     for (int idxn = 0; idxn < _num_modes_l; ++idxn) {
//       for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
//         Complex ylm = ylmn(idxl, idxm, -1, phi_r);
//         Complex ylp = ylmn(idxl, idxm, 1, phi_r);
//         vec_receiver(idxn, idxm + idxl) =
//             vec_eigval_l[idxl - 1][idx][idxq][idxn] * (ylp - ylm) / 2.0;
//       }
//     }
//   } else if (toaug) {
//     for (int idxn = 0; idxn < totnum; ++idxn) {
//       for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
//         Complex ylm = ylmn(idxl, idxm, -1, phi_r);
//         Complex ylp = ylmn(idxl, idxm, 1, phi_r);
//         vec_receiver(idxn, idxm + idxl) =
//             vec_all_l[idxl - 1][idx][idxq][idxn] * (ylp - ylm) / 2.0;
//       }
//     }
//   }
//   return vec_receiver;
// };

// Eigen::MatrixXcd
// sem::CalculateForceNMC(SourceInfo::EarthquakeCMT &cmt, int idxl, bool toaug)
// {
//   auto totnum = _num_modes_l;
//   if (toaug) {
//     totnum += numaug;
//   }
//   // std::size_t flen = _num_modes_l;
//   std::size_t fcols = 2 * idxl + 1;
//   Eigen::MatrixXcd vec_lforce = Eigen::MatrixXcd::Zero(totnum, fcols);

//   // resize force
//   // vec_force.resize(totlen);
//   // vec_force = Eigen::VectorXcd::Zero(totlen);

//   // find element within which the source sits
//   double depth = cmt.Depth();
//   double rad_source = _mesh.PR() - 1000.0 * depth / _length_norm;
//   // std::cout << "\nSource radius: " << rad_source << "\n";
//   int idxsource = 0;

//   int NQ = _mesh.NN();

//   // get the y0-, y0+ values etc at the source location
//   double theta_s = (90.0 - cmt.Latitude()) * EIGEN_PI / (180.0);
//   double phi_s = cmt.Longitude() * EIGEN_PI / (180.0);
//   auto wigdmat =
//       GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
//                        GSHTrans::Single, GSHTrans::ColumnMajor>(_lmax, _lmax,
//                        2,
//                                                                 theta_s);
//   auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
//     auto dl = wigdmat[N];
//     auto tmp = dl[l, m];
//     auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
//     return ylm;
//   };
//   // std::cout << "\nSource location: " << theta_s << " " << phi_s << "\n";
//   // int m = 0, _l = 2;
//   // parameters
//   double invsqrt2 = 1.0 / std::sqrt(2.0);
//   std::complex<double> isq2 = std::complex<double>(0.0, invsqrt2);
//   double omegal2 = std::sqrt((idxl + 2) * (idxl - 1) / 2.0);
//   double lprefac = std::exp(-2.0 * 3.141592653589793 * (idxl + 1) / (1 +
//   0.5));
//   // loop through the elements
//   //   to find the element that contains the source
//   for (int idx = _el; idx < _eu; ++idx) {
//     // std::cout << _mesh.ELR(idx) << " " << _mesh.EUR(idx) << " " <<
//     // rad_source
//     //           << "\n";
//     if ((_mesh.ELR(idx) <= rad_source) && (_mesh.EUR(idx) >= rad_source)) {
//       // std::cout << idx << " " << rad_source << " " << _mesh.ELR(idx) << "
//       "
//       //           << _mesh.EUR(idx) << "\n";
//       std::vector<double> vec_nodes(NQ, 0.0);
//       for (int idxi = 0; idxi < NQ; ++idxi) {
//         vec_nodes[idxi] = _mesh.NodeRadius(idx, idxi);
//       }
//       auto pleg =
//           Interpolation::LagrangePolynomial(vec_nodes.begin(),
//           vec_nodes.end());
//       for (int idxn = 0; idxn < totnum; ++idxn) {

//         auto w_val = 0.0;
//         auto w_deriv = 0.0;
//         for (int idxq = 0; idxq < NQ; ++idxq) {
//           w_val +=
//               vec_all_l[idxl - 1][idx][idxq][idxn] * pleg(idxq, rad_source);
//           w_deriv += vec_all_l[idxl - 1][idx][idxq][idxn] *
//                      pleg.Derivative(idxq, rad_source);
//           if (idxn == 0) {
//             // std::cout << "vec_eigval_l: " << "idxl, idx, idxq, idxn: " <<
//             // idxl
//             //           << " " << idx << " " << idxq << " " << idxn << " "
//             //           << vec_eigval_l[idxl - 1][idx][idxq][idxn] << "\n";
//           }
//         };
//         auto w_vald = w_val / rad_source;
//         auto w_prefactor = w_deriv - w_vald;
//         // std::cout << "w_vald, w_deriv, w_prefactor: " << w_vald << " "
//         //           << w_deriv << " " << w_prefactor << "\n";
//         // double omegal2 = (idxl + 2) * (idxl - 1) / 2.0;

//         // get the index of the spherical harmonic
//         for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
//           // spherical harmonic
//           Complex ymc = std::conj(ylmn(idxl, idxm, -1, phi_s));
//           Complex ypc = std::conj(ylmn(idxl, idxm, 1, phi_s));
//           Complex ymmc = 0.0, yppc = 0.0;

//           if (idxl > 1) {
//             ymmc = std::conj(ylmn(idxl, idxm, -2, phi_s));
//             yppc = std::conj(ylmn(idxl, idxm, 2, phi_s));
//           }

//           Complex tmp = w_prefactor * (cmt.MC0m() * ymc - cmt.MC0p() * ypc);
//           tmp += w_vald * omegal2 * (cmt.MCmm() * ymmc - cmt.MCpp() * yppc);
//           tmp *= isq2 * lprefac;

//           vec_lforce(idxn, idxm + idxl) = tmp;
//           // };
//           // };
//         };
//       };
//     };
//   }
//   return vec_lforce;
// };

#endif