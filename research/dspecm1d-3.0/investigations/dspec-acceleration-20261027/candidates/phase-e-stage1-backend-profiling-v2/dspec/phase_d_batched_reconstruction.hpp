#ifndef DSPECM1D_PHASE_D_BATCHED_RECONSTRUCTION_HPP
#define DSPECM1D_PHASE_D_BATCHED_RECONSTRUCTION_HPP

#include "phase_c_reconstruction.hpp"
#include <Eigen/Dense>

namespace phased {
using Complex=phasea::Complex;
using ModeKey=phasec::ModeKey;
using MatrixXcd=Eigen::Matrix<Complex,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;

struct WeightKey {
  int role=0,receiver=-1,spin=0,m=0;
  friend bool operator<(const WeightKey &a,const WeightKey &b) noexcept {
    return std::tie(a.role,a.receiver,a.spin,a.m)<std::tie(b.role,b.receiver,b.spin,b.m);
  }
};
struct FrameGeometry {
  bool isotropic=false;
  std::vector<double> phi,rotA,rotB;
  std::unique_ptr<phasea::LegendreTable> legendre;
};
struct BatchedAngularPlanCache {
  using Axes=phasec::AngularPlanCache::Axes;
  int lmax=-1;
  std::size_t angleCount=0,retainedBytes=0,maxTemporaryBytes=0;
  std::map<std::pair<int,int>,Axes> frames;
  std::map<std::size_t,std::pair<double,double>> inputAngles;
  std::map<std::pair<int,int>,FrameGeometry> geometry;
  std::map<WeightKey,std::vector<Complex>> phaseRotation;
};

inline Complex globalRotationFactor(int spin,double a,double b) {
  const Complex c(a,-b);
  switch(spin){
    case 0:return 1.0;
    case 1:return c;
    case -1:return std::conj(c);
    case 2:return c*c;
    case -2:return std::conj(c*c);
    default:throw std::invalid_argument("invalid spin in Phase D rotation cache");
  }
}

inline void prepareWeights(SPARSESPEC::Windowed3DFrequencyWorkspace &workspace,
    const std::vector<SPARSESPEC::Field3DPoint> &points,
    std::size_t angleCount,int lmax,
    const std::map<std::pair<int,int>,BatchedAngularPlanCache::Axes> &frameAxes,
    BatchedAngularPlanCache &cache) {
  mathprofile::Scope angularPlanAndWeightPreparation(
      "phase_e_angular_plan_and_weight_preparation");
  const std::size_t nNodes=points.size()/angleCount;
  std::vector<double> canonicalTheta(angleCount),canonicalPhi(angleCount);
  for(std::size_t a=0;a<angleCount;++a){
    canonicalTheta[a]=points[a*nNodes].thetaRadians;
    canonicalPhi[a]=points[a*nNodes].phiRadians;
    for(std::size_t n=1;n<nNodes;++n)
      if(points[a*nNodes+n].thetaRadians!=canonicalTheta[a]||points[a*nNodes+n].phiRadians!=canonicalPhi[a])
        throw std::runtime_error("Phase D angle tile contains inconsistent node coordinates.");
  }
  if(cache.lmax>=0){
    if(cache.lmax!=lmax||cache.angleCount!=angleCount||cache.frames!=frameAxes)
      throw std::runtime_error("Phase D worker angular cache semantic key changed.");
    for(std::size_t a=0;a<angleCount;++a)
      if(cache.inputAngles.at(a)!=std::make_pair(canonicalTheta[a],canonicalPhi[a]))
        throw std::runtime_error("Phase D worker angular cache angle key changed.");
    mathprofile::count("phase_d_batched_angular_plan_cache_hits");return;
  }
  cache.lmax=lmax;cache.angleCount=angleCount;cache.frames=frameAxes;
  for(std::size_t a=0;a<angleCount;++a)cache.inputAngles.emplace(a,std::make_pair(canonicalTheta[a],canonicalPhi[a]));
  for(const auto &[frame,axes]:frameAxes){(void)axes;
    const int role=frame.first,receiver=frame.second,orders=role?1:2;
    FrameGeometry item;item.isotropic=(role==1&&receiver<0);
    item.phi.resize(angleCount);item.rotA.resize(angleCount);item.rotB.resize(angleCount);
    std::vector<double> theta(angleCount);
    for(std::size_t a=0;a<angleCount;++a){
      if(item.isotropic){theta[a]=canonicalTheta[a];item.phi[a]=canonicalPhi[a];item.rotA[a]=1;item.rotB[a]=0;}
      else {const auto coordinate=workspace.phaseAngularCoordinates(role==1,receiver,canonicalTheta[a],canonicalPhi[a]);
        theta[a]=coordinate[0];item.phi[a]=coordinate[1];item.rotA[a]=coordinate[2];item.rotB[a]=coordinate[3];}
    }
    item.legendre=std::make_unique<phasea::LegendreTable>(lmax,orders,theta);
    cache.retainedBytes+=angleCount*5u*(2u*static_cast<std::size_t>(orders)+1u)*
        static_cast<std::size_t>(lmax+1)*sizeof(double)+angleCount*3u*sizeof(double);
    for(int spin=-2;spin<=2;++spin)for(int m=-orders;m<=orders;++m){
      auto &factor=cache.phaseRotation[WeightKey{role,receiver,spin,m}];factor.resize(angleCount);
      for(std::size_t a=0;a<angleCount;++a){
        const Complex phase=std::polar(1.0,m*item.phi[a]);
        const Complex rotation=item.isotropic?Complex(1):globalRotationFactor(spin,item.rotA[a],item.rotB[a]);
        factor[a]=phase*rotation;
      }
      cache.retainedBytes+=angleCount*sizeof(Complex);
    }
    cache.geometry.emplace(frame,std::move(item));
  }
  constexpr std::size_t WorkerBound=32ull<<20;
  if(cache.retainedBytes>(WorkerBound/2))
    throw std::runtime_error("Phase D cached Legendre geometry exceeds the 16 MiB per-worker geometry reservation.");
  mathprofile::count("phase_d_batched_angular_plan_cache_misses");
  mathprofile::count("phase_d_batched_angular_plan_bytes",cache.retainedBytes);
}

/** Matrix-batched angular reconstruction over contiguous degree rows. Products
 * use ordinary transpose (never conjugate transpose): each output is
 * sum_l radial[l]*angularWeight[l,a], with phase and frame rotation in weight. */
template<class Sink>
inline void visitNativeSidesBatched(
    SPARSESPEC::Windowed3DFrequencyWorkspace &workspace,
    const std::vector<SPARSESPEC::Field3DPoint> &points,
    const std::vector<int> &localNodes,std::size_t angleCount,
    BatchedAngularPlanCache &angularCache,Sink &&sink) {
  if(localNodes.empty()||angleCount==0||points.size()%angleCount!=0||
     angleCount*localNodes.size()!=points.size())
    throw std::invalid_argument("Phase D native-side point inventory is inconsistent.");
  const int lmax=workspace.parameters().lmax();
  const std::size_t nNodes=localNodes.size();
  double minimumRadius=std::numeric_limits<double>::infinity();
  for(std::size_t n=0;n<nNodes;++n)minimumRadius=std::min(minimumRadius,points[n].radiusMetres);
  std::map<ModeKey,phasec::NativeModeTable> tables;
  std::map<std::pair<int,int>,BatchedAngularPlanCache::Axes> frameAxes;
  mathprofile::Scope radialModes("phase_d_stream_projected_radial_modes_to_strain_tables");
  mathprofile::Scope radialInterpolationAndStrain(
      "phase_e_radial_mode_interpolation_and_strain_table_build");
  workspace.visitProjectedRadialModes(
    [&](SPARSESPEC::SensitivitySolveFamily,int degree,bool adjoint,int receiver,
        int m,int element,int localNode,double radius,bool fluid,bool,
        const std::array<std::array<double,3>,3>&axes,const Eigen::MatrixXcd &values) {
      if(values.rows()!=6||values.cols()!=1||degree<0||degree>lmax)
        throw std::runtime_error("Phase D projected radial mode shape/degree changed.");
      const int role=adjoint?1:0;const auto frame=std::make_pair(role,receiver);
      const auto [fit,inserted]=frameAxes.emplace(frame,axes);
      if(!inserted&&fit->second!=axes)throw std::runtime_error("Phase D visitor changed axes within a frame.");
      const ModeKey key{element,localNode,role,receiver};auto it=tables.find(key);
      if(it==tables.end()){
        phasec::NativeModeTable item;item.table=std::make_unique<phasea::StrainTable>(1,adjoint?1:2,lmax);
        item.radiusNormalized=radius;item.fluid=fluid;it=tables.emplace(key,std::move(item)).first;
      }else if(it->second.radiusNormalized!=radius||it->second.fluid!=fluid)
        throw std::runtime_error("Phase D modes disagree on native-side metadata.");
      std::array<Complex,6> radial{};for(int c=0;c<6;++c)radial[static_cast<std::size_t>(c)]=values(c,0);
      const auto strain=phasea::strainOf(radial,degree,radius,fluid);
      for(std::size_t c=0;c<phasea::Components;++c)(*it->second.table)(0,c,m,degree)+=strain[c];
    },minimumRadius);
  radialInterpolationAndStrain.stop();
  radialModes.stop();
  prepareWeights(workspace,points,angleCount,lmax,frameAxes,angularCache);

  std::map<std::tuple<int,int,int,int>,ModeKey> lookup;
  mathprofile::Scope angularLookupMaterialization(
      "phase_e_angular_mode_lookup_materialization");
  for(const auto &[key,value]:tables){(void)value;
    const int element=std::get<0>(key),node=std::get<1>(key),role=std::get<2>(key),receiver=std::get<3>(key);
    lookup.emplace(std::make_tuple(element,node,role,receiver),key);
  }
  angularLookupMaterialization.stop();
  constexpr Eigen::Index AngleTile=16,RadiusTile=32;
  const std::size_t matrixScratch=static_cast<std::size_t>(RadiusTile)*(static_cast<std::size_t>(lmax+1)+AngleTile)*sizeof(Complex);
  const std::size_t accumulationScratch=static_cast<std::size_t>(AngleTile)*nNodes*2u*phasea::Components*sizeof(Complex);
  const std::size_t maxWeightTile=5u*5u*static_cast<std::size_t>(lmax+1)*static_cast<std::size_t>(AngleTile)*sizeof(Complex);
  const std::size_t wignerBuildScratch=25u*static_cast<std::size_t>(lmax+1)*sizeof(double);
  constexpr std::size_t EigenGemmPackingReserve=2ull<<20;
  const std::size_t totalTemporary=angularCache.retainedBytes+matrixScratch+accumulationScratch+maxWeightTile+wignerBuildScratch+EigenGemmPackingReserve;
  if(totalTemporary>(32ull<<20))throw std::runtime_error("Phase D angular cache and tiled temporary budget exceeds 32 MiB per worker.");
  angularCache.maxTemporaryBytes=std::max(angularCache.maxTemporaryBytes,totalTemporary);
  const std::size_t tileScratchBytes=matrixScratch+accumulationScratch;
  if(tileScratchBytes>(32ull<<20))throw std::runtime_error("Phase D angular tile scratch exceeds 32 MiB per worker.");
  mathprofile::Scope angularReconstruction("phase_d_tiled_eigen_angular_reconstruction");
  for(std::size_t angleStart=0;angleStart<angleCount;angleStart+=static_cast<std::size_t>(AngleTile)){
    const Eigen::Index tile=static_cast<Eigen::Index>(std::min<std::size_t>(AngleTile,angleCount-angleStart));
    std::vector<phasea::Strain> forward(static_cast<std::size_t>(tile)*nNodes),adjoint(static_cast<std::size_t>(tile)*nNodes);
    for(const auto &[frame,axes]:frameAxes){(void)axes;
      const int role=frame.first,receiver=frame.second,orders=role?1:2;
      const auto &geometry=angularCache.geometry.at(frame);
      std::map<WeightKey,MatrixXcd> tileWeights;
      mathprofile::Scope tileWeightPacking(
          "phase_e_angular_tile_weight_packing");
      mathprofile::count("phase_e_angular_weight_tiles");
      for(int spin=-2;spin<=2;++spin)for(int m=-orders;m<=orders;++m){
        const WeightKey key{role,receiver,spin,m};MatrixXcd weight(lmax+1,tile);
        const auto &factor=angularCache.phaseRotation.at(key);
        for(Eigen::Index a=0;a<tile;++a){
          const std::size_t angle=angleStart+static_cast<std::size_t>(a);
          const double *legendre=geometry.legendre->row(angle,spin,m);
          for(int l=0;l<=lmax;++l)weight(l,a)=legendre[l]*factor[angle];
        }
        tileWeights.emplace(key,std::move(weight));
      }
      tileWeightPacking.stop();
      for(std::size_t radiusStart=0;radiusStart<nNodes;radiusStart+=static_cast<std::size_t>(RadiusTile)){
       const Eigen::Index radiusCount=static_cast<Eigen::Index>(std::min<std::size_t>(RadiusTile,nNodes-radiusStart));
       MatrixXcd coefficients(radiusCount,lmax+1),values(radiusCount,tile);
       for(Eigen::Index r=0;r<radiusCount;++r){
        const std::size_t n=radiusStart+static_cast<std::size_t>(r);
        const int element=points[n].radialElement,node=localNodes[n];
        const auto found=lookup.find(std::make_tuple(element,node,role,receiver));
        if(found==lookup.end())coefficients.row(r).setZero();
        else if(tables.at(found->second).radiusNormalized*workspace.sem().lengthScaleMetres()!=points[n].radiusMetres)
          throw std::runtime_error("Phase D radial table and requested native radius differ.");
       }
       for(int c=0;c<static_cast<int>(phasea::Components);++c)for(int m=-orders;m<=orders;++m){
        mathprofile::Scope coefficientRows("phase_e_angular_coefficient_row_lookup_and_packing");
        for(Eigen::Index r=0;r<radiusCount;++r){
          const std::size_t n=radiusStart+static_cast<std::size_t>(r);
          const int element=points[n].radialElement,node=localNodes[n];
          const auto found=lookup.find(std::make_tuple(element,node,role,receiver));
          if(found==lookup.end())continue;
          const auto &table=*tables.at(found->second).table;
          Eigen::Map<const Eigen::RowVectorXcd> radial(table.row(0,static_cast<std::size_t>(c),m),lmax+1);
          coefficients.row(r)=radial;
        }
       coefficientRows.stop();
       mathprofile::count("phase_e_angular_coefficient_row_batches");
        const MatrixXcd &weight=tileWeights.at(WeightKey{role,receiver,phasea::Spin[static_cast<std::size_t>(c)],m});
        mathprofile::Scope angularGemm("phase_e_angular_eigen_gemm");
        values.noalias()=coefficients*weight;
        angularGemm.stop();
        mathprofile::count("phase_e_angular_eigen_gemm_calls");
        auto &target=role?adjoint:forward;
        mathprofile::Scope angularScatter("phase_e_angular_tile_scatter");
        for(Eigen::Index r=0;r<radiusCount;++r)for(Eigen::Index a=0;a<tile;++a)
          target[static_cast<std::size_t>(a)*nNodes+radiusStart+static_cast<std::size_t>(r)][static_cast<std::size_t>(c)]+=values(r,a);
        angularScatter.stop();
       }
      }
    }
    mathprofile::Scope materialEvaluationAndContraction(
        "phase_e_angular_material_evaluation_and_contraction");
    for(Eigen::Index a=0;a<tile;++a)for(std::size_t n=0;n<nNodes;++n){
      const std::size_t p=(angleStart+static_cast<std::size_t>(a))*nNodes+n;
      SPARSESPEC::KernelField3DSample sample;sample.radialElement=points[p].radialElement;
      sample.fluid=workspace.sem().elementIsFluid(points[p].radialElement);
      auto f=forward[static_cast<std::size_t>(a)*nNodes+n],d=adjoint[static_cast<std::size_t>(a)*nNodes+n];
      if(sample.fluid){f=phasea::inheritedFluidProjection(f);d=phasea::inheritedFluidProjection(d);}
      sample.forwardRawStrain=phasec::toSpin(f);sample.adjointDualStrain=phasec::toSpin(d);
      sink(p,sample);
    }
    materialEvaluationAndContraction.stop();
    mathprofile::count("phase_e_angular_material_contraction_batches");
  }
  angularReconstruction.stop();
}

inline std::vector<SPARSESPEC::KernelField3DSample> sampleNativeSidesBatched(
    SPARSESPEC::Windowed3DFrequencyWorkspace &workspace,
    const std::vector<SPARSESPEC::Field3DPoint> &points,
    const std::vector<int> &localNodes,std::size_t angleCount,
    BatchedAngularPlanCache &angularCache) {
  std::vector<SPARSESPEC::KernelField3DSample> result(points.size());
  visitNativeSidesBatched(workspace,points,localNodes,angleCount,angularCache,
      [&](std::size_t index,const auto &sample){result[index]=sample;});
  return result;
}
}
#endif
