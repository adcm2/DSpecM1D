#ifndef DSPECM1D_SEM_SOURCE_TREATMENT_H
#define DSPECM1D_SEM_SOURCE_TREATMENT_H

#include "SEM.h"

namespace Full1D {

namespace source_treatment_detail {
inline Eigen::MatrixXcd radialElementOperator(const SEM &sem, int element,
                                              std::complex<double> z,
                                              std::complex<double> attenuationFactor = {}) {
  const auto &mesh = sem.mesh(); const auto &model = sem.meshModel();
  const int nq = mesh.NN();
  const double pi = static_cast<double>(EIGEN_PI), width = mesh.EW(element);
  const auto gll = const_cast<EarthMesh::RadialMesh &>(mesh).GLL();
  const auto basis = Interpolation::LagrangeBasis(gll.Points());
  Eigen::MatrixXcd local = Eigen::MatrixXcd::Zero(2 * nq, 2 * nq);
  for (int q = 0; q < nq; ++q) {
    const double r = mesh.NodeRadius(element, q), rho = model.Density(element, q);
    const double gravity = model.Gravity(element, q), dr = 0.5 * width * gll.W(q);
    const auto C = std::complex<double>(model.C(element, q)) +
                   attenuationFactor * model.C_atten(element, q);
    const auto F = std::complex<double>(model.F(element, q)) +
                   attenuationFactor * model.F_atten(element, q);
    const auto A = std::complex<double>(model.A(element, q)) +
                   attenuationFactor * model.A_atten(element, q);
    const auto N = std::complex<double>(model.N(element, q)) +
                   attenuationFactor * model.N_atten(element, q);
    Eigen::Matrix2cd P = Eigen::Matrix2cd::Zero(), Q = Eigen::Matrix2cd::Zero();
    Eigen::Matrix2cd R = Eigen::Matrix2cd::Zero(), M = Eigen::Matrix2cd::Zero();
    P(0,0) = r*r*C; P(1,1) = r*r/(4.0*pi*sem.gravityCoefficient());
    Q(0,0) = 2.0*r*F; Q(1,0) = rho*r*r;
    R(0,0) = 4.0*(A-N) + 4.0*rho*(pi*sem.gravityCoefficient()*rho*r-gravity)*r;
    M(0,0) = rho*r*r;
    for (int i=0;i<nq;++i) {
      const double ni=basis.Evaluate(i,gll.X(q));
      const double dni=2.0/width*basis.template Evaluate<1>(i,gll.X(q));
      for (int j=0;j<nq;++j) {
        const double nj=basis.Evaluate(j,gll.X(q));
        const double dnj=2.0/width*basis.template Evaluate<1>(j,gll.X(q));
        for (int a=0;a<2;++a) for (int b=0;b<2;++b)
          local(a*nq+i,b*nq+j)+=dr*(dni*P(a,b)*dnj+dni*Q(a,b)*nj+
              ni*Q(b,a)*dnj+ni*(R(a,b)-z*z*M(a,b))*nj);
      }
    }
  }
  return local;
}

inline Eigen::MatrixXcd toroidalElementOperator(const SEM &sem, int element,
                                                int degree, std::complex<double> z,
                                                std::complex<double> attenuationFactor = {}) {
  const auto &mesh=sem.mesh(); const auto &model=sem.meshModel();
  const int nq=mesh.NN();
  const double lambda=static_cast<double>(degree)*(degree+1.0), width=mesh.EW(element);
  const auto gll=const_cast<EarthMesh::RadialMesh &>(mesh).GLL();
  const auto basis=Interpolation::LagrangeBasis(gll.Points());
  Eigen::MatrixXcd local=Eigen::MatrixXcd::Zero(nq,nq);
  for(int q=0;q<nq;++q) {
    const double r=mesh.NodeRadius(element,q), dr=0.5*width*gll.W(q);
    const auto L=std::complex<double>(model.L(element,q))+
                 attenuationFactor*model.L_atten(element,q);
    const auto N=std::complex<double>(model.N(element,q))+
                 attenuationFactor*model.N_atten(element,q);
    const double rho=model.Density(element,q);
    const auto P=lambda*r*r*L, Q=-lambda*r*L;
    const auto R=lambda*(L-2.0*N)+lambda*lambda*N;
    const double M=lambda*rho*r*r;
    for(int i=0;i<nq;++i) {
      const double ni=basis.Evaluate(i,gll.X(q));
      const double dni=2.0/width*basis.template Evaluate<1>(i,gll.X(q));
      for(int j=0;j<nq;++j) {
        const double nj=basis.Evaluate(j,gll.X(q));
        const double dnj=2.0/width*basis.template Evaluate<1>(j,gll.X(q));
        local(i,j)+=dr*(dni*P*dnj+dni*Q*nj+ni*Q*dnj+ni*(R-z*z*M)*nj);
      }
    }
  }
  return local;
}

inline Eigen::MatrixXcd spheroidalElementOperator(
    const SEM &sem, int element, int degree, std::complex<double> z,
    std::complex<double> attenuationFactor, bool cowling) {
  const auto &mesh = sem.mesh();
  const auto &model = sem.meshModel();
  const int nq = mesh.NN(), fields = cowling ? 2 : 3;
  const double lambda = static_cast<double>(degree) * (degree + 1.0);
  const double pi = static_cast<double>(EIGEN_PI);
  const double width = mesh.EW(element), jac = 2.0 / width;
  const auto gll = const_cast<EarthMesh::RadialMesh &>(mesh).GLL();
  const auto basis = Interpolation::LagrangeBasis(gll.Points());
  const auto effective = [&](int node, double (MeshModel::*base)(int, int) const,
                             double (MeshModel::*atten)(int, int) const) {
    return std::complex<double>((model.*base)(element, node), 0.0) +
           attenuationFactor * (model.*atten)(element, node);
  };
  Eigen::MatrixXcd local = Eigen::MatrixXcd::Zero(fields * nq, fields * nq);
  for (int q = 0; q < nq; ++q) {
    const double r = mesh.NodeRadius(element, q);
    const double rho = model.Density(element, q), gravity = model.Gravity(element, q);
    const double dr = 0.5 * width * gll.W(q);
    const auto A = effective(q, &MeshModel::A, &MeshModel::A_atten);
    const auto C = effective(q, &MeshModel::C, &MeshModel::C_atten);
    const auto F = effective(q, &MeshModel::F, &MeshModel::F_atten);
    const auto L = effective(q, &MeshModel::L, &MeshModel::L_atten);
    const auto N = effective(q, &MeshModel::N, &MeshModel::N_atten);
    Eigen::Matrix3cd P = Eigen::Matrix3cd::Zero();
    Eigen::Matrix3cd Q = Eigen::Matrix3cd::Zero();
    Eigen::Matrix3cd R = Eigen::Matrix3cd::Zero();
    Eigen::Matrix3cd M = Eigen::Matrix3cd::Zero();
    P(0, 0) = r * r * C;
    P(1, 1) = lambda * r * r * L;
    Q(0, 0) = 2.0 * r * F; Q(0, 1) = -lambda * r * F;
    Q(1, 0) = lambda * r * L; Q(1, 1) = -lambda * r * L;
    if (!cowling) {
      P(2, 2) = r * r / (4.0 * pi * sem.gravityCoefficient());
      Q(2, 0) = rho * r * r;
      R(2, 2) = lambda / (4.0 * pi * sem.gravityCoefficient());
      R(1, 2) = R(2, 1) = lambda * rho * r;
    }
    R(0, 0) = 4.0 * (A - N) +
        4.0 * rho * (pi * sem.gravityCoefficient() * rho * r - gravity) * r +
        lambda * L;
    if (cowling) R(0, 0) -= 4.0 * pi * sem.gravityCoefficient() * rho * rho * r * r;
    R(0, 1) = R(1, 0) = lambda *
        (rho * gravity * r - L - 2.0 * (A - N));
    R(1, 1) = lambda * (L - 2.0 * N) + lambda * lambda * A;
    M(0, 0) = rho * r * r; M(1, 1) = lambda * rho * r * r;
    for (int i = 0; i < nq; ++i) {
      const double ni = basis.Evaluate(i, gll.X(q));
      const double dni = jac * basis.template Evaluate<1>(i, gll.X(q));
      for (int j = 0; j < nq; ++j) {
        const double nj = basis.Evaluate(j, gll.X(q));
        const double dnj = jac * basis.template Evaluate<1>(j, gll.X(q));
        for (int a = 0; a < fields; ++a)
          for (int b = 0; b < fields; ++b)
            local(a * nq + i, b * nq + j) += dr * (
                dni * P(a, b) * dnj + dni * Q(a, b) * nj +
                ni * Q(b, a) * dnj + ni * (R(a, b) - z * z * M(a, b)) * nj);
      }
    }
  }
  return local;
}
} // namespace source_treatment_detail

inline SEM::SourceTreatmentRadial SEM::sourceTreatedRadial(
    SourceInfo::EarthquakeCMT &cmt, std::complex<double> z,
    std::complex<double> attenuationFactor) {
  if (!sourceMatches(cmt)) throw std::invalid_argument(
      "Source treatment requires an SEM aligned to the current source radius.");
  const int e=sourceElement(cmt), nq=m_mesh.NN();
  const int ndof=ltgR(1,m_mesh.NE()-1,nq-1)+1;
  const double r=m_sourceRadius;
  auto wigner=GSHTrans::Wigner<double,GSHTrans::All,GSHTrans::All,GSHTrans::Single>(0,0,0,0.0);
  const auto y0=wigner[0][0,0];
  const std::complex<double> b=cmt.MC00()*y0/m_momentNorm;
  const std::complex<double> a=-2.0*cmt.MCmp()*y0/(r*m_momentNorm);
  const std::complex<double> C=m_meshModel.C(e,0)+
      attenuationFactor*m_meshModel.C_atten(e,0);
  const std::complex<double> F=m_meshModel.F(e,0)+
      attenuationFactor*m_meshModel.F_atten(e,0);
  const std::complex<double> j=b/(r*r*C);
  SourceTreatmentRadial result;
  result.force=Eigen::MatrixXcd::Zero(ndof,1);
  result.lift=Eigen::MatrixXcd::Zero(ndof,1); result.sourceElement=e;
  const int source=ltgR(0,e,0);
  result.force(source,0)=a-2.0*r*F*j;
  result.lift(source,0)=j;
  const auto local=source_treatment_detail::radialElementOperator(
      *this,e,z,attenuationFactor);
  Eigen::VectorXcd g=Eigen::VectorXcd::Zero(2*nq); g[0]=j;
  const Eigen::VectorXcd action=local*g;
  for(int field=0;field<2;++field) for(int i=0;i<nq;++i)
    result.force(ltgR(field,e,i),0)-=action[field*nq+i];
  return result;
}

inline SEM::SourceTreatmentToroidal SEM::sourceTreatedToroidal(
    SourceInfo::EarthquakeCMT &cmt, int degree, std::complex<double> z,
    std::complex<double> attenuationFactor) {
  if (!sourceMatches(cmt)) throw std::invalid_argument(
      "Source treatment requires an SEM aligned to the current source radius.");
  if(degree<1 || degree>m_lmax) throw std::invalid_argument("Toroidal source degree is invalid.");
  const int e=sourceElement(cmt), nq=m_mesh.NN();
  if (e < m_el || e >= m_eu)
    throw std::invalid_argument(
        "Toroidal source treatment is supported only in the connected solid mantle.");
  const int ndof=ltgT(m_eu-1,nq-1)+1;
  const double r=m_sourceRadius, lambda=static_cast<double>(degree)*(degree+1.0);
  const double k=std::sqrt(lambda);
  const std::complex<double> L=m_meshModel.L(e,0)+
      attenuationFactor*m_meshModel.L_atten(e,0);
  const std::complex<double> j=k/(lambda*r*r*L);
  SourceTreatmentToroidal result;
  result.force=Eigen::MatrixXcd::Zero(ndof,2);
  result.lift=Eigen::MatrixXcd::Zero(ndof,2); result.sourceElement=e;
  const int source=ltgT(e,0);
  result.force(source,0)=k/r;
  result.force(source,1)=k/r;
  result.lift(source,0)=j;
  const auto local=source_treatment_detail::toroidalElementOperator(
      *this,e,degree,z,attenuationFactor);
  Eigen::VectorXcd g=Eigen::VectorXcd::Zero(nq); g[0]=j;
  const Eigen::VectorXcd action=local*g;
  for(int i=0;i<nq;++i) result.force(ltgT(e,i),0)-=action[i];
  return result;
}

/** Build the exact discrete spheroidal source load and its one-sided lift. */
inline SEM::SourceTreatmentSpheroidal SEM::sourceTreatedSpheroidal(
    SourceInfo::EarthquakeCMT &cmt, int degree, std::complex<double> z,
    std::complex<double> attenuationFactor, bool cowling) {
  if (!sourceMatches(cmt))
    throw std::invalid_argument(
        "Source treatment requires an SEM aligned to the current source radius.");
  if (degree < 1 || degree > m_lmax)
    throw std::invalid_argument("Source treatment requires 1 <= degree <= lmax.");

  const int element = sourceElement(cmt);
  const int nq = m_mesh.NN();
  const int fields = cowling ? 2 : 3;
  const int ndof = cowling ? ltgSC(1, m_mesh.NE() - 1, nq - 1) + 1
                           : ltgS(2, m_mesh.NE() - 1, nq - 1) + 1;
  const double lambda = static_cast<double>(degree) *
                        static_cast<double>(degree + 1);
  const double angularSourceScale = std::sqrt(static_cast<double>(lambda) / 2.0);
  const double radius = m_sourceRadius;
  const auto dof = [this, cowling, element](int field, int node) {
    return cowling ? ltgSC(field, element, node)
                   : ltgS(field, element, node);
  };

  SourceTreatmentSpheroidal result;
  result.force = Eigen::MatrixXcd::Zero(ndof, 4);
  result.lift = Eigen::MatrixXcd::Zero(ndof, 4);
  result.sourceElement = element;

  const auto effective = [this, element, attenuationFactor](
                             int node, double (MeshModel::*base)(int, int) const,
                             double (MeshModel::*atten)(int, int) const) {
    return std::complex<double>((m_meshModel.*base)(element, node), 0.0) +
           attenuationFactor * (m_meshModel.*atten)(element, node);
  };
  const int sourceNode = 0;
  const auto Csource = effective(sourceNode, &MeshModel::C,
                                 &MeshModel::C_atten);
  const auto Lsource = effective(sourceNode, &MeshModel::L,
                                 &MeshModel::L_atten);
  // The source coefficients are material-scaled by the sampled operator. For
  // attenuating solves use the same complex P as H(z), including its phase.
  const std::complex<double> jumpU = 1.0 / (radius * radius * Csource);
  const std::complex<double> jumpV =
      angularSourceScale /
      (lambda * radius * radius * Lsource);
  result.force(dof(0, sourceNode), 1) = 1.0 / radius;
  result.force(dof(1, sourceNode), 3) = angularSourceScale / radius;

  const auto Fsource = effective(sourceNode, &MeshModel::F,
                                 &MeshModel::F_atten);
  const auto LsourceEff = Lsource;
  const std::complex<double> qUFromU = 2.0 * radius * Fsource * jumpU;
  const std::complex<double> qVFromU = -lambda * radius * Fsource * jumpU;
  const std::complex<double> qUFromV =
      lambda * radius * LsourceEff * jumpV;
  const std::complex<double> qVFromV =
      -lambda * radius * LsourceEff * jumpV;
  result.force(dof(0, sourceNode), 0) -= qUFromU;
  result.force(dof(1, sourceNode), 0) -= qVFromU;
  result.force(dof(0, sourceNode), 2) -= qUFromV;
  result.force(dof(1, sourceNode), 2) -= qVFromV;
  result.lift(dof(0, sourceNode), 0) = jumpU;
  result.lift(dof(1, sourceNode), 2) = jumpV;

  // Apply H_loc to the element-local cardinal lift.
  const Eigen::MatrixXcd localOperator = source_treatment_detail::spheroidalElementOperator(
      *this, element, degree, z, attenuationFactor, cowling);
  Eigen::MatrixXcd localLift = Eigen::MatrixXcd::Zero(fields * nq, 4);
  for (int column : {0, 2}) {
    const int field = column == 0 ? 0 : 1;
    localLift(field * nq + sourceNode, column) =
        result.lift(dof(field, sourceNode), column);
  }
  for (const int sourceColumn : {0, 2}) {
    const Eigen::VectorXcd action =
        localOperator * localLift.col(sourceColumn);
    for (int testField = 0; testField < fields; ++testField)
      for (int testNode = 0; testNode < nq; ++testNode)
        result.force(dof(testField, testNode), sourceColumn) -=
            action[testField * nq + testNode];
  }

  if (!result.force.allFinite() || !result.lift.allFinite())
    throw std::runtime_error("Source lift produced non-finite values.");
  return result;
}

inline Eigen::MatrixXcd SEM::sourceLiftReceiver(
    InputParameters &params, int degree, const Eigen::MatrixXcd &lift,
    int sourceElement, bool cowling) {
  const auto receivers = receiverElements(params);
  const int low = cowling ? ltgSC(0, receivers.front(), 0)
                          : ltgS(0, receivers.front(), 0);
  const int high = cowling ? ltgSC(1, receivers.back(), m_nq - 1)
                           : ltgS(1, receivers.back(), m_nq - 1);
  Eigen::MatrixXcd result = rvBaseFull(params, degree, cowling) *
                            lift.block(low, 0, high - low + 1, lift.cols());
  const double sourceRadius = m_sourceRadius;
  const double sourcePatchTop = m_mesh.EUR(sourceElement);
  for (int receiver = 0; receiver < params.num_receivers(); ++receiver) {
    const double radius = _ReceiverRadius(params);
    const bool upperTrace = radius >= sourceRadius && radius <= sourcePatchTop;
    if (!upperTrace) result.block(3 * receiver, 0, 3, result.cols()).setZero();
  }
  return result;
}

inline Eigen::MatrixXcd SEM::sourceLiftReceiverRadial(
    InputParameters &params, const Eigen::MatrixXcd &lift, int sourceElement) {
  Eigen::MatrixXcd result=Eigen::MatrixXcd::Zero(3*params.num_receivers(),lift.cols());
  const double patchTop=m_mesh.EUR(sourceElement), source=m_sourceRadius;
  for(int r=0;r<params.num_receivers();++r) {
    if(_ReceiverRadius(params)<source || _ReceiverRadius(params)>patchTop) continue;
    result(3*r, Eigen::all)=rvZR(params,r).transpose()*lift;
  }
  return result;
}

inline Eigen::MatrixXcd SEM::sourceLiftReceiverToroidal(
    InputParameters &params,int degree,const Eigen::MatrixXcd &lift,int sourceElement) {
  if (sourceElement < m_el || sourceElement >= m_eu)
    throw std::invalid_argument(
        "Toroidal source lift is supported only in the connected solid mantle.");
  const auto elems=receiverElements(params,true);
  const int low=ltgT(elems.front(),0), high=ltgT(elems.back(),m_nq-1);
  Eigen::MatrixXcd result=rvBaseFullT(params,degree)*lift.block(low,0,high-low+1,lift.cols());
  const double patchTop=m_mesh.EUR(sourceElement), source=m_sourceRadius;
  for(int r=0;r<params.num_receivers();++r)
    if(_ReceiverRadius(params)<source || _ReceiverRadius(params)>patchTop)
      result.block(3*r,0,3,result.cols()).setZero();
  return result;
}

} // namespace Full1D

#endif // DSPECM1D_SEM_SOURCE_TREATMENT_H
