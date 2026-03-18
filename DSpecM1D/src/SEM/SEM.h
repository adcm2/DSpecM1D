#ifndef SEM_H
#define SEM_H

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

class InputParametersNew;

namespace Full1D {

// spectral element solver
class SEM {
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
  SEM() {};
  explicit SEM(const InputParametersNew &);
  template <class model1d> SEM(const model1d &, double, int, int);

  // local-to-global maps
  auto ltgS(int, int, int) const;
  auto ltgT(int, int) const;
  auto ltgR(int, int, int) const;
  auto el() const { return _el; };
  auto eu() const { return _eu; };

  // source / receiver element queries
  auto receiverElements(InputParameters &) const;
  auto sourceElement(SourceInfo::EarthquakeCMT &) const;

  // accessors
  const EarthMesh::RadialMesh &mesh() const { return _mesh; };
  const MeshModel &meshModel() const { return _mesh_model; };

  // spheroidal force vectors
  Eigen::MatrixXcd calculateForce(SourceInfo::EarthquakeCMT &, int);
  Eigen::MatrixXcd calculateForceAll(SourceInfo::EarthquakeCMT &, int);
  Eigen::MatrixXcd calculateForceCoefficients(SourceInfo::EarthquakeCMT &, int);
  Eigen::MatrixXcd calculateForceRedCoefficients(SourceInfo::EarthquakeCMT &,
                                                 int, double);

  // toroidal force vectors
  Eigen::MatrixXcd calculateForceT(SourceInfo::EarthquakeCMT &, int);
  Eigen::MatrixXcd calculateForceAllT(SourceInfo::EarthquakeCMT &, int);
  Eigen::MatrixXcd calculateForceCoefficientsT(SourceInfo::EarthquakeCMT &,
                                               int);
  Eigen::MatrixXcd calculateForceRedCoefficientsT(SourceInfo::EarthquakeCMT &,
                                                  int, double);

  // radial force vectors
  Eigen::MatrixXcd calculateForceR(SourceInfo::EarthquakeCMT &);
  Eigen::MatrixXcd calculateForceRedR(SourceInfo::EarthquakeCMT &);

  // system matrices
  SMAT hTa(int) const;
  SMAT hSa(int) const;
  SMAT hS(int) const;
  SMAT pS(int) const;
  SMAT hTk(int) const;
  SMAT pTk(int) const;
  SMAT hR() const { return mat_ke_0; };
  SMAT hRa() const { return mat_ke_0_atten; };
  SMAT pR() const { return mat_inertia_0; };

  // spheroidal receiver vectors
  Eigen::MatrixXcd rvFull(InputParameters &, int);
  Eigen::MatrixXcd rvFullT(InputParameters &, int);
  Eigen::MatrixXcd rvBaseZ(InputParameters &, int, int);
  Eigen::MatrixXcd rvValZ(InputParameters &, int, int);
  Eigen::MatrixXcd rvBaseTheta(InputParameters &, int, int);
  Eigen::MatrixXcd rvValTheta(InputParameters &, int, int);
  Eigen::MatrixXcd rvBaseThetaT(InputParameters &, int, int);
  Eigen::MatrixXcd rvValThetaT(InputParameters &, int, int);
  Eigen::MatrixXcd rvValPhi(InputParameters &, int, int);
  Eigen::MatrixXcd rvBasePhiT(InputParameters &, int, int);
  Eigen::MatrixXcd rvValPhiT(InputParameters &, int, int);
  Eigen::MatrixXcd rvBaseFull(InputParameters &, int);
  Eigen::MatrixXcd rvBaseFullT(InputParameters &, int);

  // toroidal receiver vectors
  Eigen::MatrixXcd rvThetaT(InputParameters &, int, int);
  Eigen::MatrixXcd rvPhiT(InputParameters &, int, int);

  // radial receiver vectors
  Eigen::MatrixXcd rvZR(InputParameters &, int);
  Eigen::MatrixXcd rvRedZR(InputParameters &);
};

}   // namespace Full1D

#include "SEMLG.h"
#include "SEMConstructor.h"
#include "SEMMatrices.h"
#include "SEMForceSpheroidal.h"
#include "SEMForceToroidal.h"
#include "SEMForceRadial.h"
#include "SEMReceivers.h"

#endif   // SEM_H
