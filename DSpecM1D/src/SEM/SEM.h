#ifndef DSPECM1D_SEM_H
#define DSPECM1D_SEM_H

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

/**
 * @brief Spectral element model builder and operator factory for 1D problems.
 *
 * `SEM` owns the radial mesh, material interpolation, local-to-global maps,
 * and sparse matrices/vectors required by the frequency-domain solver. It is
 * part of the public API, but most release users should access it through the
 * higher-level `SparseFSpec` and `InputParametersNew` workflow unless they
 * need low-level matrix or vector access.
 */
class SEM {
private:
  // abbreviations
  using stdvec = std::vector<double>;
  using stdvvec = std::vector<stdvec>;
  using stdvvvec = std::vector<stdvvec>;
  using Complex = std::complex<double>;
  using SMAT = Eigen::SparseMatrix<double>;

  // mesh
  EarthMesh::RadialMesh m_mesh;
  MeshModel m_meshModel;

  // integers for calc
  int m_lmax, _il, m_el = 0, m_eu = 0, _en = 0, numlen = 0, _k2 = 0, _solint = 0,
                  m_nq;

  // dealing with fluid regions
  bool m_hasFluid = false;
  std::vector<int> m_vecFluid, m_fsb, m_vecOffset{0};
  std::vector<bool> m_vecDof;

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
  double m_frequencyNorm, m_lengthNorm, m_momentNorm;

  // base sparse matrices
  SMAT m_matInertia0, m_matKe0, m_matKe0Atten, m_matInTBase;

  // private helpers
  double _SourceRadius(const SourceInfo::EarthquakeCMT &cmt) const {
    return m_mesh.PR() - 1000.0 * cmt.Depth() / m_lengthNorm;
  }
  double _ReceiverRadius(const InputParameters &param) const {
    return m_mesh.PR() - 1000.0 * param.receiver_depth() / m_lengthNorm;
  }
  std::vector<SMAT> m_vecKeTBase, m_vecKeTAtten, m_vecKeSBase,
      m_vecKeSAtten, m_vecInSBase;

public:
  SEM() {};
  /// Constructs the SEM directly from the preferred workflow aggregate.
  explicit SEM(const InputParametersNew &);
  /// Legacy low-level constructor retained for compatibility.
  template <class model1d> SEM(const model1d &, double, int, int);

  /// Local-to-global map for spheroidal unknowns.
  auto ltgS(int, int, int) const;
  /// Local-to-global map for toroidal unknowns.
  auto ltgT(int, int) const;
  /// Local-to-global map for radial unknowns.
  auto ltgR(int, int, int) const;
  auto el() const { return m_el; };
  auto eu() const { return m_eu; };

  /// Returns the element indices containing all receivers.
  auto receiverElements(InputParameters &) const;
  /// Returns the source element index.
  auto sourceElement(SourceInfo::EarthquakeCMT &) const;

  /// Returns the SEM radial mesh.
  const EarthMesh::RadialMesh &mesh() const { return m_mesh; };
  /// Returns the mesh-interpolated material model.
  const MeshModel &meshModel() const { return m_meshModel; };

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
  SMAT hR() const { return m_matKe0; };
  SMAT hRa() const { return m_matKe0Atten; };
  SMAT pR() const { return m_matInertia0; };

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

#endif   // DSPECM1D_SEM_H
