#ifndef DSPECM1D_PHASE_C_RECONSTRUCTION_HPP
#define DSPECM1D_PHASE_C_RECONSTRUCTION_HPP

#include "phase_a_reconstruction.hpp"

#include <algorithm>
#include <bit>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <tuple>
#include <vector>

namespace phasec {

using Complex = std::complex<double>;
using ModeKey = std::tuple<int,int,int,int>; // element, local node, role, receiver

struct NativeModeTable {
  std::unique_ptr<phasea::StrainTable> table;
  double radiusNormalized=0.0;
  bool fluid=false;
};

struct AngularKey {
  int role=0, receiver=-1;
  std::size_t angle=0;
  friend bool operator<(const AngularKey &a,const AngularKey &b) noexcept {
    return std::tie(a.role,a.receiver,a.angle)<std::tie(b.role,b.receiver,b.angle);
  }
};

struct AngularEval {
  bool isotropic=false;
  double inputTheta=0.0,inputPhi=0.0,phi=0.0,a=0.0,b=0.0;
  std::unique_ptr<phasea::LegendreTable> legendre;
};

struct AngularPlanCache {
  using Axes=std::array<std::array<double,3>,3>;
  int lmax=-1;
  std::size_t angleCount=0;
  std::map<std::pair<int,int>,Axes> frames;
  std::map<AngularKey,AngularEval> plans;
  std::map<std::size_t,std::pair<double,double>> inputAngles;
};

inline SPARSESPEC::SpinStrainField toSpin(const phasea::Strain &v) {
  return {v[0],v[1],v[2],v[3],v[4],v[5]};
}

/** Reconstruct all requested native one-sided nodes from the accepted Phase A
 * radial modes. Every call owns its mode tables, so frequency workers do not
 * retain modal histories or share mutable reconstruction state. */
inline std::vector<SPARSESPEC::KernelField3DSample> sampleNativeSides(
    SPARSESPEC::Windowed3DFrequencyWorkspace &workspace,
    const std::vector<SPARSESPEC::Field3DPoint> &points,
    const std::vector<int> &localNodes,std::size_t angleCount,
    AngularPlanCache &angularCache) {
  if(localNodes.empty()||angleCount*localNodes.size()!=points.size())
    throw std::invalid_argument("Phase C native-side point inventory is inconsistent.");
  const int lmax=workspace.parameters().lmax();
  double minimumRadius=std::numeric_limits<double>::infinity();
  for(std::size_t n=0;n<localNodes.size();++n)
    minimumRadius=std::min(minimumRadius,points[n].radiusMetres);
  std::map<ModeKey,NativeModeTable> tables;
  std::map<std::pair<int,int>,AngularPlanCache::Axes> frameAxes;
  mathprofile::Scope radialModes("phase_c_stream_projected_radial_modes_to_strain_tables");
  workspace.visitProjectedRadialModes(
      [&](SPARSESPEC::SensitivitySolveFamily,int degree,bool adjoint,int receiver,
          int m,int element,int localNode,double radius,bool fluid,bool,
          const std::array<std::array<double,3>,3>&axes,const Eigen::MatrixXcd &values) {
        if(values.rows()!=6||values.cols()!=1)
          throw std::runtime_error("Phase C projected radial mode shape changed.");
        if(degree<0||degree>lmax)throw std::runtime_error("Phase C mode degree exceeds table.");
        const int role=adjoint?1:0;
        const auto frameKey=std::make_pair(role,receiver);
        const auto [frameIt,frameInserted]=frameAxes.emplace(frameKey,axes);
        if(!frameInserted&&frameIt->second!=axes)
          throw std::runtime_error("Phase C visitor changed axes within one semantic frame.");
        const ModeKey key{element,localNode,role,receiver};
        auto it=tables.find(key);
        if(it==tables.end()){
          NativeModeTable item;
          item.table=std::make_unique<phasea::StrainTable>(1,adjoint?1:2,lmax);
          item.radiusNormalized=radius;item.fluid=fluid;
          it=tables.emplace(key,std::move(item)).first;
        }else if(it->second.radiusNormalized!=radius||it->second.fluid!=fluid){
          throw std::runtime_error("Phase C radial modes disagree on native-side metadata.");
        }
        std::array<Complex,6> radial{};
        for(int c=0;c<6;++c)radial[static_cast<std::size_t>(c)]=values(c,0);
        const auto strain=phasea::strainOf(radial,degree,radius,fluid);
        // Spheroidal and toroidal modes occupy the same degree/order slot and
        // their strain contributions add linearly.
        for(std::size_t c=0;c<phasea::Components;++c)
          (*it->second.table)(0,c,m,degree)+=strain[c];
      },minimumRadius);
  radialModes.stop();

  std::map<std::pair<int,int>,std::vector<ModeKey>> byNode;
  std::set<std::pair<int,int>> frames;
  for(const auto &[key,value]:tables){
    (void)value;
    const int role=std::get<2>(key),receiver=std::get<3>(key);
    byNode[{std::get<0>(key),std::get<1>(key)}].push_back(key);
    frames.emplace(role,receiver);
  }
  if(angularCache.lmax<0){
    angularCache.lmax=lmax;angularCache.angleCount=angleCount;angularCache.frames=frameAxes;
    for(std::size_t a=0;a<angleCount;++a){
      const auto &point=points[a*localNodes.size()];
      angularCache.inputAngles.emplace(a,std::make_pair(point.thetaRadians,point.phiRadians));
    }
  }else if(angularCache.lmax!=lmax||angularCache.angleCount!=angleCount||
           angularCache.frames!=frameAxes){
    throw std::runtime_error("Phase C worker angular-plan cache semantic key changed.");
  }
  for(const auto &[role,receiver]:frames)for(std::size_t a=0;a<angleCount;++a){
    const auto &point=points[a*localNodes.size()];
    const AngularKey key{role,receiver,a};
    const auto prior=angularCache.plans.find(key);
    if(prior!=angularCache.plans.end()){
      if(prior->second.inputTheta!=point.thetaRadians||prior->second.inputPhi!=point.phiRadians||
         angularCache.inputAngles.at(a)!=std::make_pair(point.thetaRadians,point.phiRadians))
        throw std::runtime_error("Phase C worker angular-plan angle key changed.");
      mathprofile::count("phase_c_worker_angular_plan_cache_hits");
      continue;
    }
    AngularEval eval;eval.inputTheta=point.thetaRadians;eval.inputPhi=point.phiRadians;
    eval.isotropic=(role==1&&receiver<0);
    double theta=point.thetaRadians;
    if(eval.isotropic){eval.phi=point.phiRadians;}
    else {
      const auto coordinates=workspace.phaseAngularCoordinates(role==1,receiver,
          point.thetaRadians,point.phiRadians);
      theta=coordinates[0];eval.phi=coordinates[1];eval.a=coordinates[2];eval.b=coordinates[3];
    }
    eval.legendre=std::make_unique<phasea::LegendreTable>(lmax,role?1:2,
        std::vector<double>{theta});
    angularCache.plans.emplace(key,std::move(eval));
    mathprofile::count("phase_c_worker_angular_plan_cache_misses");
  }

  std::vector<SPARSESPEC::KernelField3DSample> result(points.size());
  mathprofile::Scope angularReconstruction("phase_c_degree_inner_legendre_reconstruction");
  for(std::size_t a=0;a<angleCount;++a)for(std::size_t n=0;n<localNodes.size();++n){
    const std::size_t p=a*localNodes.size()+n;
    const auto &point=points[p];const int localNode=localNodes[n];
    auto &sample=result[p];sample.radialElement=point.radialElement;
    sample.fluid=workspace.sem().elementIsFluid(point.radialElement);
    phasea::Strain forward{},adjoint{};
    const auto found=byNode.find({point.radialElement,localNode});
    if(found==byNode.end())throw std::runtime_error("Phase C lacks a projected mode at a requested native side.");
    const double expectedRadius=point.radiusMetres/workspace.sem().lengthScaleMetres();
    for(const auto &key:found->second){
      const auto &mode=tables.at(key);
      const double tolerance=8.0*std::numeric_limits<double>::epsilon()*
          std::max({1.0,std::abs(expectedRadius),std::abs(mode.radiusNormalized)});
      if(std::abs(mode.radiusNormalized-expectedRadius)>tolerance||mode.fluid!=sample.fluid)
        throw std::runtime_error("Phase C mode table does not match the requested one-sided native radius.");
    }
    for(const auto &key:found->second){
      const int role=std::get<2>(key),receiver=std::get<3>(key);
      const auto &eval=angularCache.plans.at(AngularKey{role,receiver,a});
      auto value=phasea::evaluate(*tables.at(key).table,0,*eval.legendre,0,eval.phi);
      if(!eval.isotropic)value=phasea::rotateToGlobal(value,eval.a,eval.b);
      auto &target=role?adjoint:forward;
      for(std::size_t c=0;c<phasea::Components;++c)target[c]+=value[c];
    }
    if(sample.fluid){
      forward=phasea::inheritedFluidProjection(forward);
      adjoint=phasea::inheritedFluidProjection(adjoint);
    }
    sample.forwardRawStrain=toSpin(forward);
    sample.adjointDualStrain=toSpin(adjoint);
  }
  angularReconstruction.stop();
  return result;
}

} // namespace phasec

#endif
