#ifndef DSPECM1D_PHASE_D_BATCHED_RECONSTRUCTION_HPP
#define DSPECM1D_PHASE_D_BATCHED_RECONSTRUCTION_HPP

#include "phase_c_reconstruction.hpp"
#include <Eigen/Dense>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <string>

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
    BatchedAngularPlanCache &angularCache,Sink &&sink,int frequencyBin=-1,
    bool allowRadialOracle=true) {
  if(localNodes.empty()||angleCount==0||points.size()%angleCount!=0||
     angleCount*localNodes.size()!=points.size())
    throw std::invalid_argument("Phase D native-side point inventory is inconsistent.");
  const int lmax=workspace.parameters().lmax();
  const std::size_t nNodes=localNodes.size();
  double minimumRadius=std::numeric_limits<double>::infinity();
  for(std::size_t n=0;n<nNodes;++n)minimumRadius=std::min(minimumRadius,points[n].radiusMetres);
  const char *radialRouteEnv=std::getenv("DSFAST_RADIAL_ROUTE");
  const bool useGenericRadial=radialRouteEnv&&std::string(radialRouteEnv)=="generic";
  const char *oracleEnv=std::getenv("DSFAST_RADIAL_ORACLE");
  const bool checkRadialOracle=allowRadialOracle&&oracleEnv&&
      std::string(oracleEnv)=="1"&&!useGenericRadial;
  struct RadialOracleRecord {
    std::array<Complex,6> values{};
    double radius=0.0;
    bool fluid=false,liftApplied=false;
    std::array<std::array<double,3>,3> axes{};
  };
  using RadialOracleKey=std::tuple<int,int,int,int,int,int,int>;
  std::map<RadialOracleKey,RadialOracleRecord> radialOracle;
  std::array<long double,8> oracleError2{},oracleReference2{};
  std::size_t radialOracleCompared=0;
  std::array<std::size_t,4> oracleCategoryCounts{};
  std::size_t radialOracleConservativeBytes=0;
  constexpr std::size_t OracleMemoryBound=256ull<<20;
  if(checkRadialOracle){
    workspace.visitProjectedRadialModes(
        [&](SPARSESPEC::SensitivitySolveFamily family,int degree,bool adjoint,int receiver,
            int m,int element,int localNode,double radius,bool fluid,bool liftApplied,
            const std::array<std::array<double,3>,3>&axes,const Eigen::MatrixXcd &values){
          if(values.rows()!=6||values.cols()!=1)
            throw std::runtime_error("Stage 2 radial oracle shape changed.");
          RadialOracleRecord record;record.radius=radius;record.fluid=fluid;
          record.liftApplied=liftApplied;record.axes=axes;
          for(int c=0;c<6;++c)record.values[static_cast<std::size_t>(c)]=values(c,0);
          const RadialOracleKey key{static_cast<int>(family),degree,adjoint?1:0,
              receiver,m,element,localNode};
          if(!radialOracle.emplace(key,std::move(record)).second)
            throw std::runtime_error("Stage 2 radial oracle produced a duplicate mode key.");
          radialOracleConservativeBytes=radialOracle.size()*
              (sizeof(RadialOracleRecord)+sizeof(RadialOracleKey)+128u);
          if(radialOracleConservativeBytes>OracleMemoryBound)
            throw std::runtime_error("Stage 2 radial oracle capture exceeds its 256 MiB per-workspace bound.");
        },minimumRadius);
  }
  const int meshElements=workspace.sem().mesh().NE();
  const int meshNodes=workspace.sem().mesh().NN();
  const std::size_t frameSlots=static_cast<std::size_t>(workspace.parameters().num_receivers())+2u;
  std::vector<std::unique_ptr<phasec::NativeModeTable>> tables;
  if(!useGenericRadial)tables.resize(
      static_cast<std::size_t>(meshElements)*static_cast<std::size_t>(meshNodes)*frameSlots);
  std::map<ModeKey,phasec::NativeModeTable> baselineTables;
  const auto tableIndex=[&](int element,int localNode,int role,int receiver){
    if(element<0||element>=meshElements||localNode<0||localNode>=meshNodes)
      throw std::runtime_error("Phase D radial table key is outside the native mesh.");
    std::size_t frameSlot=0;
    if(role==0&&receiver==-1)frameSlot=0;
    else if(role==1&&receiver==-1)frameSlot=1;
    else if(role==1&&receiver>=0)frameSlot=static_cast<std::size_t>(receiver)+2u;
    else throw std::runtime_error("Phase D radial table role/receiver key is invalid.");
    if(frameSlot>=frameSlots)
      throw std::runtime_error("Phase D radial table receiver key is outside the workspace.");
    return (static_cast<std::size_t>(element)*static_cast<std::size_t>(meshNodes)+
        static_cast<std::size_t>(localNode))*frameSlots+frameSlot;
  };
  std::map<std::pair<int,int>,BatchedAngularPlanCache::Axes> frameAxes;
  mathprofile::Scope radialModes("phase_d_stream_projected_radial_modes_to_strain_tables");
  mathprofile::Scope radialInterpolationAndStrain(
      "phase_e_radial_mode_interpolation_and_strain_table_build");
  const auto visitRadial=[&](const auto &visitor,double minimumRadius){
    if(useGenericRadial)workspace.visitProjectedRadialModes(visitor,minimumRadius);
    else workspace.visitProjectedRadialModesBatched(visitor,minimumRadius);
  };
  visitRadial(
    [&](SPARSESPEC::SensitivitySolveFamily family,int degree,bool adjoint,int receiver,
        int m,int element,int localNode,double radius,bool fluid,bool liftApplied,
        const std::array<std::array<double,3>,3>&axes,const Eigen::MatrixXcd &values) {
      if(values.rows()!=6||values.cols()!=1||degree<0||degree>lmax)
        throw std::runtime_error("Phase D projected radial mode shape/degree changed.");
      std::array<Complex,6> oracleValues{};
      if(checkRadialOracle){
        const RadialOracleKey key{static_cast<int>(family),degree,adjoint?1:0,
            receiver,m,element,localNode};
        const auto found=radialOracle.find(key);
        if(found==radialOracle.end())
          throw std::runtime_error("Stage 2 fast radial mode is missing from the generic oracle.");
        const auto &reference=found->second;
        oracleValues=reference.values;
        if(reference.radius!=radius||reference.fluid!=fluid||
           reference.liftApplied!=liftApplied||reference.axes!=axes)
          throw std::runtime_error("Stage 2 radial metadata differs from generic oracle.");
        long double error2=0.0L,reference2=0.0L;
        for(int c=0;c<6;++c){
          const Complex delta=values(c,0)-reference.values[static_cast<std::size_t>(c)];
          error2+=std::norm(delta);reference2+=std::norm(reference.values[static_cast<std::size_t>(c)]);
          oracleError2[0]+=std::norm(delta);oracleReference2[0]+=std::norm(reference.values[static_cast<std::size_t>(c)]);
          const std::size_t group=(c%2==0)?1u:2u;
          oracleError2[group]+=std::norm(delta);oracleReference2[group]+=std::norm(reference.values[static_cast<std::size_t>(c)]);
          const std::size_t roleGroup=adjoint?4u:3u;
          oracleError2[roleGroup]+=std::norm(delta);oracleReference2[roleGroup]+=std::norm(reference.values[static_cast<std::size_t>(c)]);
          const std::size_t fluidGroup=fluid?5u:6u;
          oracleError2[fluidGroup]+=std::norm(delta);oracleReference2[fluidGroup]+=std::norm(reference.values[static_cast<std::size_t>(c)]);
        }
        const double relative=std::sqrt(error2/std::max(reference2,1e-300L));
        if(!std::isfinite(relative)||relative>1e-11)
          throw std::runtime_error("Stage 2 batched radial interpolation differs from generic oracle above 1e-11.");
        radialOracle.erase(found);
        ++radialOracleCompared;
        ++oracleCategoryCounts[adjoint?1u:0u];
        ++oracleCategoryCounts[fluid?2u:3u];
      }
      const int role=adjoint?1:0;const auto frame=std::make_pair(role,receiver);
      const auto [fit,inserted]=frameAxes.emplace(frame,axes);
      if(!inserted&&fit->second!=axes)throw std::runtime_error("Phase D visitor changed axes within a frame.");
      const ModeKey key{element,localNode,role,receiver};
      phasec::NativeModeTable *table=nullptr;
      if(useGenericRadial){
        auto it=baselineTables.find(key);
        if(it==baselineTables.end()){
          phasec::NativeModeTable fresh;
          fresh.table=std::make_unique<phasea::StrainTable>(1,adjoint?1:2,lmax);
          fresh.radiusNormalized=radius;fresh.fluid=fluid;
          it=baselineTables.emplace(key,std::move(fresh)).first;
        }
        table=&it->second;
        if(!table->table){
          table->table=std::make_unique<phasea::StrainTable>(1,adjoint?1:2,lmax);
          table->radiusNormalized=radius;table->fluid=fluid;
        }
      }else{
        auto &slot=tables[tableIndex(element,localNode,role,receiver)];
        if(!slot){
          slot=std::make_unique<phasec::NativeModeTable>();
          slot->table=std::make_unique<phasea::StrainTable>(1,adjoint?1:2,lmax);
          slot->radiusNormalized=radius;slot->fluid=fluid;
        }
        table=slot.get();
      }
      if(table->radiusNormalized!=radius||table->fluid!=fluid)
        throw std::runtime_error("Phase D modes disagree on native-side metadata.");
      std::array<Complex,6> radial{};for(int c=0;c<6;++c)radial[static_cast<std::size_t>(c)]=values(c,0);
      const auto strain=phasea::strainOf(radial,degree,radius,fluid);
      if(checkRadialOracle){
        const auto oracleStrain=phasea::strainOf(oracleValues,degree,radius,fluid);
        long double error2=0.0L,reference2=0.0L;
        for(std::size_t c=0;c<phasea::Components;++c){
          error2+=std::norm(strain[c]-oracleStrain[c]);
          reference2+=std::norm(oracleStrain[c]);
          oracleError2[7]+=std::norm(strain[c]-oracleStrain[c]);
          oracleReference2[7]+=std::norm(oracleStrain[c]);
        }
        const double relative=std::sqrt(error2/std::max(reference2,1e-300L));
        if(!std::isfinite(relative)||relative>1e-11)
          throw std::runtime_error("Stage 2 batched radial strain differs from generic oracle above 1e-11.");
      }
      for(std::size_t c=0;c<phasea::Components;++c)(*table->table)(0,c,m,degree)+=strain[c];
    },minimumRadius);
  if(checkRadialOracle){
    if(!radialOracle.empty())
      throw std::runtime_error("Stage 2 generic radial oracle contains modes absent from batched path.");
    std::array<double,8> oracleRelative{};
    for(std::size_t i=0;i<oracleRelative.size();++i){
      oracleRelative[i]=std::sqrt(static_cast<double>(oracleError2[i]/
          std::max(oracleReference2[i],1e-300L)));
      if(!std::isfinite(oracleRelative[i])||oracleRelative[i]>1e-11)
        throw std::runtime_error("Stage 2 radial oracle category exceeds the 1e-11 gate.");
    }
    const char *proofPath=std::getenv("DSFAST_RADIAL_ORACLE_PROOF");
    if(!proofPath||frequencyBin<0)
      throw std::runtime_error("Stage 2 radial oracle proof requires a path and frequency bin.");
    const std::string file=std::string(proofPath)+".bin"+
        std::to_string(frequencyBin)+".csv";
    std::ofstream proof(file,std::ios::out|std::ios::trunc);
    if(!proof)throw std::runtime_error("Cannot write Stage 2 radial oracle proof row.");
    proof<<frequencyBin<<','<<radialOracleCompared<<','<<radialOracleConservativeBytes<<','<<oracleCategoryCounts[0]<<','
         <<oracleCategoryCounts[1]<<','<<oracleCategoryCounts[2]<<','
         <<oracleCategoryCounts[3]<<std::setprecision(17);
    for(double relative:oracleRelative)proof<<','<<relative;
    proof<<",generic-and-batched-keyed-modes\n";
  }
  radialInterpolationAndStrain.stop();
  radialModes.stop();
  prepareWeights(workspace,points,angleCount,lmax,frameAxes,angularCache);

  std::map<std::tuple<int,int,int,int>,ModeKey> baselineLookup;
  if(useGenericRadial){
    mathprofile::Scope angularLookupMaterialization(
        "phase_e_angular_mode_lookup_materialization");
    for(const auto &[key,value]:baselineTables){(void)value;
      baselineLookup.emplace(std::make_tuple(std::get<0>(key),std::get<1>(key),
          std::get<2>(key),std::get<3>(key)),key);
    }
    angularLookupMaterialization.stop();
  }
  const auto getTable=[&](int element,int node,int role,int receiver)->const phasec::NativeModeTable*{
    if(useGenericRadial){
      const auto found=baselineLookup.find(std::make_tuple(element,node,role,receiver));
      return found==baselineLookup.end()?nullptr:&baselineTables.at(found->second);
    }
    return tables[tableIndex(element,node,role,receiver)].get();
  };

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
          const auto *table=getTable(element,node,role,receiver);
          if(!table)coefficients.row(r).setZero();
          else if(table->radiusNormalized*workspace.sem().lengthScaleMetres()!=points[n].radiusMetres)
          throw std::runtime_error("Phase D radial table and requested native radius differ.");
       }
       for(int c=0;c<static_cast<int>(phasea::Components);++c)for(int m=-orders;m<=orders;++m){
        mathprofile::Scope coefficientRows("phase_e_angular_coefficient_row_lookup_and_packing");
        for(Eigen::Index r=0;r<radiusCount;++r){
          const std::size_t n=radiusStart+static_cast<std::size_t>(r);
          const int element=points[n].radialElement,node=localNodes[n];
          const auto *table=getTable(element,node,role,receiver);
          if(!table)continue;
          Eigen::Map<const Eigen::RowVectorXcd> radial(table->table->row(0,static_cast<std::size_t>(c),m),lmax+1);
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
    BatchedAngularPlanCache &angularCache,int frequencyBin=-1) {
  std::vector<SPARSESPEC::KernelField3DSample> result(points.size());
  visitNativeSidesBatched(workspace,points,localNodes,angleCount,angularCache,
      [&](std::size_t index,const auto &sample){result[index]=sample;},frequencyBin,false);
  return result;
}
}
#endif
