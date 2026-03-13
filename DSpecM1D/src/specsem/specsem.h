#ifndef SPECSEM_H
#define SPECSEM_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <Eigen/Core>
#include <GSHTrans/Core>
#include <Interpolation/Lagrange>
#include "../SourceInfo.h"
#include <EarthMesh/All>
#include "../InputParser.h"
#include "../MeshModel.h"
#include "../NormClass.h"

namespace Full1D {

// spectral element solver
class specsem {
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
  int _lmax, _il, _el = 0, _eu = 0, _en = 0, numlen = 0, _k2 = 0, _solint = 0,
                  _NQ;

  // dealing with fluid regions
  bool _has_fluid = false;
  std::vector<int> _vec_fluid, fsb, vec_offset{0};
  std::vector<bool> _vec_dof;

  // matrix parameters
  std::size_t mlen, totlen;

  // derivatives of basis functions
  stdvvec vec_lag_deriv, vec_delta;

  // normalisation factors
  double densitynorm = 5515.0;
  double pi_db = EIGEN_PI;
  double bigg_nd = 6.6723 * std::pow(10.0, -11.0);
  double bigg_db;
  double frequencynorm = std::sqrt(pi_db * bigg_nd * densitynorm);
  double normint;
  double _freq_norm, _length_norm, _moment_norm;

  // base sparse matrices
  SMAT mat_inertia_0, mat_ke_0, mat_ke_0_atten, mat_in_t_base;

  // private helpers
  double _SourceRadius(const SourceInfo::EarthquakeCMT &cmt) const {
    return _mesh.PR() - 1000.0 * cmt.Depth() / _length_norm;
  }
  double _ReceiverRadius(const InputParameters &param) const {
    return _mesh.PR() - 1000.0 * param.receiver_depth() / _length_norm;
  }
  std::vector<SMAT> vec_ke_t_base, vec_ke_t_atten, vec_ke_s_base,
      vec_ke_s_atten, vec_in_s_base;

public:
  specsem() {};
  template <class model1d> specsem(const model1d &, double, int, int);

  // local-to-global maps
  auto LtG_S(int, int, int) const;
  auto LtG_T(int, int) const;
  auto LtG_R(int, int, int) const;
  auto EL() const { return _el; };
  auto EU() const { return _eu; };

  // source / receiver element queries
  auto Receiver_Elements(InputParameters &) const;
  auto Source_Element(SourceInfo::EarthquakeCMT &) const;

  // accessors
  const EarthMesh::RadialMesh &mesh() const { return _mesh; };
  const MeshModel &mesh_model() const { return _mesh_model; };

  // spheroidal force vectors
  Eigen::MatrixXcd CalculateForce(SourceInfo::EarthquakeCMT &, int);
  Eigen::MatrixXcd CalculateForce_All(SourceInfo::EarthquakeCMT &, int);
  Eigen::MatrixXcd CalculateForce_Coefficients(SourceInfo::EarthquakeCMT &,
                                               int);
  Eigen::MatrixXcd CalculateForce_RED_Coefficients(SourceInfo::EarthquakeCMT &,
                                                   int, double);

  // toroidal force vectors
  Eigen::MatrixXcd CalculateForce_T(SourceInfo::EarthquakeCMT &, int);
  Eigen::MatrixXcd CalculateForce_All_T(SourceInfo::EarthquakeCMT &, int);
  Eigen::MatrixXcd CalculateForce_Coefficients_T(SourceInfo::EarthquakeCMT &,
                                                 int);
  Eigen::MatrixXcd
  CalculateForce_RED_Coefficients_T(SourceInfo::EarthquakeCMT &, int, double);

  // radial force vectors
  Eigen::MatrixXcd CalculateForce_R(SourceInfo::EarthquakeCMT &);
  Eigen::MatrixXcd CalculateForce_Red_R(SourceInfo::EarthquakeCMT &);

  // system matrices
  SMAT H_TA(int) const;
  SMAT H_SA(int) const;
  SMAT H_S(int) const;
  SMAT P_S(int) const;
  SMAT H_TK(int) const;
  SMAT P_TK(int) const;
  SMAT H_R() const { return mat_ke_0; };
  SMAT H_RA() const { return mat_ke_0_atten; };
  SMAT P_R() const { return mat_inertia_0; };

  // spheroidal receiver vectors
  Eigen::MatrixXcd RV_FULL(InputParameters &, int);
  Eigen::MatrixXcd RV_FULL_T(InputParameters &, int);
  Eigen::MatrixXcd RV_BASE_Z(InputParameters &, int, int);
  Eigen::MatrixXcd RV_VAL_Z(InputParameters &, int, int);
  Eigen::MatrixXcd RV_BASE_THETA(InputParameters &, int, int);
  Eigen::MatrixXcd RV_VAL_THETA(InputParameters &, int, int);
  Eigen::MatrixXcd RV_BASE_THETA_T(InputParameters &, int, int);
  Eigen::MatrixXcd RV_VAL_THETA_T(InputParameters &, int, int);
  Eigen::MatrixXcd RV_VAL_PHI(InputParameters &, int, int);
  Eigen::MatrixXcd RV_BASE_PHI_T(InputParameters &, int, int);
  Eigen::MatrixXcd RV_VAL_PHI_T(InputParameters &, int, int);
  Eigen::MatrixXcd RV_BASE_FULL(InputParameters &, int);
  Eigen::MatrixXcd RV_BASE_FULL_T(InputParameters &, int);

  // toroidal receiver vectors
  Eigen::MatrixXcd RV_THETA_T(InputParameters &, int, int);
  Eigen::MatrixXcd RV_PHI_T(InputParameters &, int, int);

  // radial receiver vectors
  Eigen::MatrixXcd RV_Z_R(InputParameters &, int);
  Eigen::MatrixXcd RV_RED_Z_R(InputParameters &);
};

}   // namespace Full1D

#include "specsem_ltg.h"
#include "specsem_constructor.h"
#include "specsem_matrices.h"
#include "specsem_force_spheroidal.h"
#include "specsem_force_toroidal.h"
#include "specsem_force_radial.h"
#include "specsem_receivers.h"

#endif   // SPECSEM_H
