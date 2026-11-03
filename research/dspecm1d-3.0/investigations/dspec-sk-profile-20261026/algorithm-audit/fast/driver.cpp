#include "fast_engine.hpp"
#include "engine_headers/SensitivityKernels.h"
#include <DSpecM1D/All>
#include <DSpecM1D/ModelInput>
#include "engine_headers/Windowed3DKernels.h"
#include <DSpecM1D/src/InputParser.h>
#include <DSpecM1D/src/NormClass.h>
#include <DSpecM1D/src/SourceInfo.h>

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>

namespace fs = std::filesystem;
using Complex = std::complex<double>;
using Clock = std::chrono::steady_clock;
constexpr int nt = 128, bins = nt / 2 + 1;

std::string paramsText(const fs::path &model, const fs::path &output, int lmax) {
  return "\"" + (output / "unused").string() + "\"\n\"" + model.string() +
      "\"\n4\n0\n2\n0\n0\n1e-5\n0\n" + std::to_string(lmax) +
      "\n0.1\n0.8\n30\n20\n0\n0\n0\n0\n50\n10\n20\n1e19\n2e18\n-3e18\n-5e18\n4e18\n-5e18\n3\n1\n45 90\n";
}

std::vector<SPARSESPEC::Field3DPoint> readPoints(const fs::path &path) {
  std::ifstream file(path);
  if (!file) throw std::runtime_error("cannot read sample CSV");
  std::string line; std::getline(file, line);
  std::vector<SPARSESPEC::Field3DPoint> out;
  while (std::getline(file, line)) {
    std::istringstream row(line); std::array<std::string, 11> v;
    for (auto &x : v) if (!std::getline(row, x, ',')) throw std::runtime_error("invalid sample row");
    if (v[0] == "J_Z" && v[1] == "A")
      out.push_back({std::stod(v[2]), std::stod(v[3]), std::stod(v[4]), std::stoi(v[5])});
  }
  if (out.empty()) throw std::runtime_error("no J_Z,A sample rows found");
  return out;
}

Eigen::VectorXd window() {
  Eigen::VectorXd w = Eigen::VectorXd::Zero(nt);
  for (int n = 0; n <= 90; ++n) w[n] = 0.5 - 0.5 * std::cos(2 * std::numbers::pi * n / 90.0);
  return w;
}

Eigen::MatrixXcd weights(int objective) {
  Eigen::MatrixXcd q = Eigen::MatrixXcd::Zero(3, bins);
  q(objective % 3, 1) = Complex(0.7, 0.2);
  q(objective % 3, 2) = Complex(-0.25, 0.45);
  return q;
}

double relNorm(const Eigen::MatrixXcd &a, const Eigen::MatrixXcd &b) {
  const double den=b.norm(),num=(a-b).norm();
  return den==0.0?(num==0.0?0.0:std::numeric_limits<double>::infinity()):num/den;
}
void addKernel(SPARSESPEC::ModulusKernel3D &a,const SPARSESPEC::ModulusKernel3D &b){
  a.A+=b.A;a.C+=b.C;a.F+=b.F;a.L+=b.L;a.N+=b.N;a.kappa+=b.kappa;a.fluid=b.fluid;
}
struct KernelError { double maxL2=0,maxPointwise=0; };
KernelError kernelError(const std::vector<SPARSESPEC::ModulusKernel3D> &a,
                        const std::vector<SPARSESPEC::ModulusKernel3D> &b) {
  KernelError result;
  for(int modulus=0;modulus<6;++modulus){
    std::vector<double> diff(a.size()),ref(b.size());double peak=0,normRef=0,normDiff=0;
    for(std::size_t i=0;i<a.size();++i){
      const std::array<double,6> av{a[i].A,a[i].C,a[i].F,a[i].L,a[i].N,a[i].kappa};
      const std::array<double,6> bv{b[i].A,b[i].C,b[i].F,b[i].L,b[i].N,b[i].kappa};
      diff[i]=av[modulus]-bv[modulus];ref[i]=bv[modulus];peak=std::max(peak,std::abs(ref[i]));
      normDiff+=diff[i]*diff[i];normRef+=ref[i]*ref[i];
    }
    result.maxL2=std::max(result.maxL2,normRef==0.0?(normDiff==0.0?0.0:std::numeric_limits<double>::infinity()):std::sqrt(normDiff/normRef));
    const double floor=0.01*peak;
    for(std::size_t i=0;i<a.size();++i){
      const double den=std::max(std::abs(ref[i]),floor);
      if(den==0.0){if(diff[i]!=0.0)result.maxPointwise=std::numeric_limits<double>::infinity();}
      else result.maxPointwise=std::max(result.maxPointwise,std::abs(diff[i])/den);
    }
  }
  return result;
}

template<class T> void writeScalar(std::ofstream &out, const T &x) {
  out.write(reinterpret_cast<const char*>(&x), sizeof(T));
}
void writeMatrix(std::ofstream &out, const Eigen::MatrixXcd &a) {
  const std::uint32_t rows=static_cast<std::uint32_t>(a.rows()), cols=static_cast<std::uint32_t>(a.cols());
  writeScalar(out,rows); writeScalar(out,cols);
  for(Eigen::Index r=0;r<a.rows();++r)for(Eigen::Index c=0;c<a.cols();++c){
    const double re=a(r,c).real(),im=a(r,c).imag();writeScalar(out,re);writeScalar(out,im);
  }
}
void writeEvaluation(std::ofstream &out,const SPARSESPEC::Windowed3DEvaluationResult &r){
  writeMatrix(out,r.rawDisplacementSpectrum);writeMatrix(out,r.processedSISpectrum);
  const std::uint64_t nk=r.pointKernels.size(),nc=r.sourceContact.size();writeScalar(out,nk);
  for(const auto &k:r.pointKernels){const std::array<double,6> v{k.A,k.C,k.F,k.L,k.N,k.kappa};
    for(double x:v)writeScalar(out,x);}
  writeScalar(out,nc);
  for(const auto &c:r.sourceContact){for(double x:{c.A,c.C,c.F,c.L,c.N})writeScalar(out,x);}
  writeScalar(out,r.objective);
}
std::array<std::array<double,3>,3> frameAxes(double latitudeDeg,double longitudeDeg){
  const double lat=latitudeDeg*std::numbers::pi/180.0,lon=longitudeDeg*std::numbers::pi/180.0;
  return {{{std::sin(lat)*std::cos(lon),std::sin(lat)*std::sin(lon),-std::cos(lat)},
           {-std::sin(lon),std::cos(lon),0.0},
           {std::cos(lat)*std::cos(lon),std::cos(lat)*std::sin(lon),std::sin(lat)}}};
}

SPARSESPEC::Windowed3DEvaluationResult candidate(
    Full1D::SEM &sem, SourceInfo::EarthquakeCMT &source, InputParameters &parameters,
    const SpectraSolver::FreqFull &frequency, const Eigen::VectorXd &timeWindow,
    const Eigen::MatrixXcd &weightsSI, int endBin,
    const std::vector<SPARSESPEC::Field3DPoint> &points,
    const std::vector<std::array<double,2>> &angles, double &seconds,
    std::ofstream *packet, std::uint32_t objective,int backend,
    SPARSESPEC::Windowed3DFieldEngine::SamplingPlanCache &plan,
    dssolver::BatchResponseCache *batchResponses,
    bool tensorMode,
    bool tensorRealMode,
    std::array<double,5> &stages) {
  const auto t0=Clock::now();
  SPARSESPEC::WindowedDisplacementMap measurement(frequency,sem.lengthScaleMetres(),timeWindow);
  auto stage0=Clock::now();
  const Eigen::MatrixXcd qhat=measurement.adjoint(weightsSI);
  stages[0]+=std::chrono::duration<double>(Clock::now()-stage0).count();
  SPARSESPEC::Windowed3DEvaluationResult out;
  out.rawDisplacementSpectrum=Eigen::MatrixXcd::Zero(weightsSI.rows(),weightsSI.cols());
  out.pointKernels.resize(points.size()); out.sourceContact.resize(angles.size());
  std::vector<SPARSESPEC::Field3DPoint> receiverPoints;
  const auto receiverElements=sem.receiverElements(parameters,false);
  const int e=receiverElements.at(0);
  const double rad=std::numbers::pi/180.0;
  for(int rec=0;rec<parameters.num_receivers();++rec){
    const auto station=parameters.receivers().at(rec);
    const double r=(sem.mesh().PR()-1000.0*parameters.receiver_depth()/sem.lengthScaleMetres())*sem.lengthScaleMetres();
    receiverPoints.push_back({r,(90.0-station.first)*rad,station.second*rad,e});
  }
  std::vector<SPARSESPEC::Field3DPoint> all=points;
  const std::size_t receiverBegin=all.size(); all.insert(all.end(),receiverPoints.begin(),receiverPoints.end());
  const auto w=frequency.w();
  const char *modalEnv=std::getenv("DSFAST_EXPORT_MODAL");
  const bool exportModal=!modalEnv||std::string(modalEnv)!="0";
  for(int bin=1;bin<endBin;++bin){
    const Complex z=w.at(bin)+Complex(0,-frequency.ep());
    stage0=Clock::now();
    auto workspace=SPARSESPEC::Windowed3DFieldEngine::prepareFrequencyWorkspace(
        sem,source,parameters,z,qhat.col(bin),batchResponses);
    stages[1]+=std::chrono::duration<double>(Clock::now()-stage0).count();
    std::vector<SPARSESPEC::PoleAngularCoefficientRecord3D> poleRecords;
    std::set<std::tuple<std::uint64_t,bool,int,std::uint64_t>> packetTraceKeys;
    std::function<void(const SPARSESPEC::PoleAngularCoefficientRecord3D&)> visitor;
    if(packet)visitor=[&](const auto &record){
      const auto key=std::tuple{record.synthesisSet,record.adjoint,record.element,
          std::bit_cast<std::uint64_t>(record.radiusModel)};
      if(packetTraceKeys.insert(key).second)poleRecords.push_back(record);
    };
    stage0=Clock::now();
    const auto fields=tensorMode?
        workspace.sampleForKernelsTensor(all,receiverBegin,&plan,backend,visitor,tensorRealMode):
        workspace.sampleForKernelsPacked(all,receiverBegin,&plan,backend,visitor);
    stages[2]+=std::chrono::duration<double>(Clock::now()-stage0).count();
    const std::uint32_t binId=static_cast<std::uint32_t>(bin);
    const std::uint64_t pointCount=static_cast<std::uint64_t>(points.size());
    const std::uint32_t blockCount=exportModal?static_cast<std::uint32_t>(workspace.modalBlocks().size()):0;
    const std::uint32_t poleCount=static_cast<std::uint32_t>(poleRecords.size());
    if(packet){
    writeScalar(*packet,objective);writeScalar(*packet,binId);writeScalar(*packet,pointCount);writeScalar(*packet,blockCount);writeScalar(*packet,poleCount);
    for(std::size_t p=0;p<points.size();++p){
      writeScalar(*packet,points[p].radiusMetres);writeScalar(*packet,points[p].thetaRadians);
      writeScalar(*packet,points[p].phiRadians);writeScalar(*packet,points[p].radialElement);
      const auto &f=fields[p].forwardRawStrain;const auto &a=fields[p].adjointDualStrain;
      const std::array<Complex,12> values{f.rr,f.horizontalTrace,f.radialPlus,f.radialMinus,
          f.horizontalPlus,f.horizontalMinus,a.rr,a.horizontalTrace,a.radialPlus,a.radialMinus,
          a.horizontalPlus,a.horizontalMinus};
      for(const auto &v:values){const double re=v.real(),im=v.imag();writeScalar(*packet,re);writeScalar(*packet,im);}
    }
    for(const auto &record:poleRecords){
      writeScalar(*packet,record.synthesisSet);writeScalar(*packet,record.adjoint);
      writeScalar(*packet,record.element);writeScalar(*packet,record.radiusModel);
      for(const auto &axis:record.frame.axes)for(double x:axis)writeScalar(*packet,x);
      const std::uint64_t n=record.labels.size();writeScalar(*packet,n);
      for(std::size_t q=0;q<record.labels.size();++q){
        writeScalar(*packet,record.labels[q].first);writeScalar(*packet,record.labels[q].second);
        for(const auto &v:record.rawStrain[q]){const double re=v.real(),im=v.imag();writeScalar(*packet,re);writeScalar(*packet,im);}
      }
    }
    if(exportModal)for(const auto &b:workspace.modalBlocks()){
      const std::int32_t family=static_cast<std::int32_t>(b.family),degree=b.degree,sourceElement=b.sourceElement;
      writeScalar(*packet,family);writeScalar(*packet,degree);writeScalar(*packet,sourceElement);
      writeMatrix(*packet,b.forward);writeMatrix(*packet,b.adjoint);writeMatrix(*packet,b.sourceLift);
      writeMatrix(*packet,b.sourceCoefficients);writeMatrix(*packet,b.modeReceiver);
    }
    }
    stage0=Clock::now();
    for(std::size_t p=0;p<points.size();++p)
      addKernel(out.pointKernels[p],SPARSESPEC::Windowed3DKernelEngine::pointKernel(sem,fields[p]));
    stages[3]+=std::chrono::duration<double>(Clock::now()-stage0).count();
    stage0=Clock::now();
    const auto contacts=SPARSESPEC::Windowed3DKernelEngine::sourceShellContact(workspace,source,angles);
    stages[4]+=std::chrono::duration<double>(Clock::now()-stage0).count();
    for(std::size_t p=0;p<contacts.size();++p){
      out.sourceContact[p].A+=contacts[p].A;out.sourceContact[p].C+=contacts[p].C;
      out.sourceContact[p].F+=contacts[p].F;out.sourceContact[p].L+=contacts[p].L;
      out.sourceContact[p].N+=contacts[p].N;
    }
    const Complex fz=-Complex(0,1)/z;
    for(int rec=0;rec<parameters.num_receivers();++rec)for(int c=0;c<3;++c)
      out.rawDisplacementSpectrum(3*rec+c,bin)=fz*fields[receiverBegin+rec].forwardRawDisplacement[c];
  }
  out.processedSISpectrum=measurement.forward(out.rawDisplacementSpectrum);
  out.objective=weightsSI.conjugate().cwiseProduct(out.processedSISpectrum).sum().real();
  seconds=std::chrono::duration<double>(Clock::now()-t0).count();
  return out;
}

int main(int argc,char **argv) try {
  if(argc!=10) throw std::runtime_error("usage: fast_driver MODEL SAMPLES EMPTY_OUTPUT LMAX NF OBJECTIVES BACKEND RADIAL_MODE SOLVER_BACKEND");
  const fs::path model=fs::absolute(argv[1]), samples=fs::absolute(argv[2]), output=fs::absolute(argv[3]);
  const int lmax=std::stoi(argv[4]), nfreq=std::stoi(argv[5]), nobj=std::stoi(argv[6]);
  const std::string backend=argv[7], radialMode=argv[8], solverBackend=argv[9];
  if((lmax!=3&&lmax!=6&&lmax!=12)||nfreq<1||nfreq>48||nobj<1||nobj>3||
      (backend!="scalar"&&backend!="eigen"&&backend!="openblas"&&backend!="tensor"&&backend!="tensor-real")||
      (radialMode!="direct"&&radialMode!="native")||
      (solverBackend!="equilibrated"&&solverBackend!="band"&&solverBackend!="natural"&&solverBackend!="condense"))
    throw std::runtime_error("invalid bounded case selection");
  const int solverBackendId=solverBackend=="equilibrated"?0:solverBackend=="band"?1:solverBackend=="natural"?2:3;
  SPARSESPEC::FixedMeshAdjointSolver::setDiagnosticBackend(solverBackendId);
  if(!fs::is_directory(output)||fs::directory_iterator(output)!=fs::directory_iterator{})
    throw std::runtime_error("output must be an existing empty directory");
  const auto paramPath=output/"parameters.txt";
  {std::ofstream f(paramPath);f<<paramsText(model,output,lmax);}
  InputParametersNew input(paramPath.string(),8,1,0.05,1.0,0.1,0.0,1);
  auto &parameters=input.inputParameters();auto &frequency=input.freqFull();auto &source=input.cmt();auto &base=input.earthModel();
  Full1D::SEM sem(base,0.05,8,lmax,source);
  if(frequency.nt()!=nt||nfreq>frequency.nt()/2) throw std::runtime_error("resolved frequency grid mismatch");
  auto points=readPoints(samples);
  if(radialMode=="native"){
    std::map<std::pair<std::uint64_t,std::uint64_t>,std::array<double,2>> uniqueAngles;
    for(const auto &p:points) uniqueAngles.try_emplace(
        std::pair{std::bit_cast<std::uint64_t>(p.thetaRadians),std::bit_cast<std::uint64_t>(p.phiRadians)},
        std::array<double,2>{p.thetaRadians,p.phiRadians});
    if(uniqueAngles.empty())throw std::runtime_error("native grid requires the sample angular set");
    std::vector<std::array<double,2>> angleSet;for(const auto &[k,v]:uniqueAngles){(void)k;angleSet.push_back(v);}
    points.clear();
    for(int e=0;e<sem.mesh().NE();++e)for(int q=0;q<sem.mesh().NN();++q){
      const double radius=sem.mesh().NodeRadius(e,q)*sem.lengthScaleMetres();
      if(!(radius>0.0))continue; // Raw point kernels are undefined at the center.
      for(const auto &a:angleSet)points.push_back({radius,a[0],a[1],e});
    }
  }
  const auto timeWindow=window();
  const bool tensorMode=backend=="tensor"||backend=="tensor-real";
  const bool tensorRealMode=backend=="tensor-real";
  const int backendId=backend=="scalar"?1:backend=="eigen"?2:3;
  const double thetaSource=(90.0-source.Latitude())*std::numbers::pi/180.0;
  const double phiSource=source.Longitude()*std::numbers::pi/180.0;
  const std::vector<std::array<double,2>> angles{{thetaSource,phiSource},
      {std::numbers::pi-thetaSource,phiSource+std::numbers::pi},{0.0,0.3},
      {std::numbers::pi,1.2},{std::numbers::pi/3,65.0*std::numbers::pi/180.0},
      {0.73,1.17},{std::nextafter(0.73,1.0),1.17}};
  double maxFieldErr=0,maxKernelL2Err=0,maxKernelPointErr=0,maxObjectiveErr=0,maxContactErr=0;
  std::vector<double> candidateSeconds,referenceSeconds;
  std::vector<std::array<double,5>> candidateStages;
  std::ofstream candidateArrays(output/"arrays_candidate.bin",std::ios::binary);
  std::ofstream referenceArrays(output/"arrays_reference.bin",std::ios::binary);
  std::ofstream coordinates(output/"coordinates.bin",std::ios::binary);
  if(!candidateArrays||!referenceArrays||!coordinates) throw std::runtime_error("cannot create binary arrays");
  const std::uint64_t coordCount=points.size();writeScalar(coordinates,coordCount);
  for(const auto &p:points){writeScalar(coordinates,p.radiusMetres);writeScalar(coordinates,p.thetaRadians);
    writeScalar(coordinates,p.phiRadians);writeScalar(coordinates,p.radialElement);}
  const bool exportPacket=std::getenv("DSFAST_EXPORT_PACKET")&&std::string(std::getenv("DSFAST_EXPORT_PACKET"))=="1";
  const char *modalEnv=std::getenv("DSFAST_EXPORT_MODAL");
  const bool exportModal=!modalEnv||std::string(modalEnv)!="0";
  const bool batchResponsesEnabled=std::getenv("DSFAST_BATCH_RESPONSES")&&
      std::string(std::getenv("DSFAST_BATCH_RESPONSES"))=="1";
  std::ofstream packetFile;
  std::ofstream *packet=nullptr;
  if(exportPacket){packetFile.open(output/"solver_packet.bin",std::ios::binary);if(!packetFile) throw std::runtime_error("cannot create solver_packet.bin");packet=&packetFile;
  packetFile.write("DSFASTPK",8);const std::uint32_t schema=1;writeScalar(packetFile,schema);
  const std::uint32_t objectiveCount=static_cast<std::uint32_t>(nobj),frequencyCount=static_cast<std::uint32_t>(nfreq);
  const std::uint64_t packetPoints=static_cast<std::uint64_t>(points.size());
  writeScalar(packetFile,objectiveCount);writeScalar(packetFile,frequencyCount);writeScalar(packetFile,packetPoints);}
  SPARSESPEC::Windowed3DFieldEngine::SamplingPlanCache sharedPlan;
  dssolver::BatchResponseCache batchResponses;
  for(int obj=0;obj<nobj;++obj){
    const auto q=weights(obj);
    double candS=0;std::array<double,5> stages{};
    auto fast=candidate(sem,source,parameters,frequency,timeWindow,q,nfreq+1,points,angles,
        candS,packet,static_cast<std::uint32_t>(obj),backendId,sharedPlan,
        batchResponsesEnabled?&batchResponses:nullptr,tensorMode,tensorRealMode,stages);
    candidateStages.push_back(stages);
    const auto t0=Clock::now();
    const auto ref=SPARSESPEC::Windowed3DKernelEngine::evaluate(sem,source,parameters,frequency,
        timeWindow,q,1,nfreq+1,points,angles);
    const double refS=std::chrono::duration<double>(Clock::now()-t0).count();
    candidateSeconds.push_back(candS);referenceSeconds.push_back(refS);
    maxFieldErr=std::max(maxFieldErr,relNorm(fast.rawDisplacementSpectrum,ref.rawDisplacementSpectrum));
    const auto ke=kernelError(fast.pointKernels,ref.pointKernels);
    maxKernelL2Err=std::max(maxKernelL2Err,ke.maxL2);maxKernelPointErr=std::max(maxKernelPointErr,ke.maxPointwise);
    const double od=std::abs(ref.objective),on=std::abs(fast.objective-ref.objective);
    const double oe=od==0.0?(on==0.0?0.0:std::numeric_limits<double>::infinity()):on/od;
    maxObjectiveErr=std::max(maxObjectiveErr,oe);
    for(std::size_t c=0;c<std::min(fast.sourceContact.size(),ref.sourceContact.size());++c){
      const auto &a=fast.sourceContact[c],&b=ref.sourceContact[c];
      for(const auto pair:{std::pair{a.A,b.A},std::pair{a.C,b.C},std::pair{a.F,b.F},
                           std::pair{a.L,b.L},std::pair{a.N,b.N}}){
        const double den=std::abs(pair.second),num=std::abs(pair.first-pair.second);
        const double err=den==0.0?(num==0.0?0.0:std::numeric_limits<double>::infinity()):num/den;
        maxContactErr=std::max(maxContactErr,err);
      }
    }
    writeEvaluation(candidateArrays,fast);writeEvaluation(referenceArrays,ref);
  }
  candidateArrays.close();referenceArrays.close();coordinates.close();if(packet)packetFile.close(); fs::remove(paramPath);
  std::ofstream meta(output/"result.json");
  meta<<std::setprecision(17)<<"{\"degree\":"<<lmax<<",\"frequency_bins\":"<<nfreq
      <<",\"objectives\":"<<nobj<<",\"backend\":\""<<backend<<"\",\"radial_mode\":\""<<radialMode<<"\",\"solver_backend\":\""<<solverBackend
      <<"\",\"sample_count\":"<<points.size()<<",\"max_field_relative_error\":"<<maxFieldErr
      <<",\"max_kernel_l2_relative_error\":"<<maxKernelL2Err<<",\"max_kernel_pointwise_scaled_error\":"<<maxKernelPointErr<<",\"max_contact_relative_error\":"<<maxContactErr<<",\"max_objective_relative_error\":"<<maxObjectiveErr
      <<",\"arrays\":{\"candidate\":\"arrays_candidate.bin\",\"reference\":\"arrays_reference.bin\",\"each_objective_order\":\"raw_complex_matrix,processed_complex_matrix,kernel_count,kernels[A,C,F,L,N,kappa],contact_count,contacts[A,C,F,L,N],objective\"},\"coordinates\":\"coordinates.bin: u64 count then radius_m,theta_rad,phi_rad,i32 element per point\",\"frames\":{\"axis_vectors_global_xyz\":[\"theta\",\"east\",\"radial\"],\"source\":";
  const auto sf=frameAxes(source.Latitude(),source.Longitude());meta<<"[";
  for(int i=0;i<3;++i){meta<<(i?",":"")<<"[";for(int j=0;j<3;++j)meta<<(j?",":"")<<sf[i][j];meta<<"]";}meta<<"],\"receivers\":[";
  for(int rec=0;rec<parameters.num_receivers();++rec){const auto station=parameters.receivers()[rec];const auto rf=frameAxes(station.first,station.second);meta<<(rec?",":"")<<"[";
    for(int i=0;i<3;++i){meta<<(i?",":"")<<"[";for(int j=0;j<3;++j)meta<<(j?",":"")<<rf[i][j];meta<<"]";}meta<<"]";}meta<<"]},\"packet\":{\"file\":\"solver_packet.bin\",\"magic\":\"DSFASTPK\",\"schema\":1,\"endianness\":\"little\",\"record_order\":\"objective,bin,point_samples,pole_radial_coefficients,modal_blocks\",\"record_header\":\"u32 objective,u32 bin,u64 point_count,u32 modal_block_count,u32 pole_coefficient_count\",\"point_value_order\":\"radius_m,theta_rad,phi_rad,element,forward_strain[rr,trace,+1,-1,+2,-2],adjoint_strain[same]\",\"pole_coeff_order\":\"u64 synthesisSet,bool adjoint,i32 element,f64 radiusModel,frame_axes[3x3],u64 count,(i32 m,i32 l,6 complex raw strain)[]\",\"matrix_order\":\"row_major_complex_interleaved_f64\",\"modal_matrix_order\":\"forward,adjoint,sourceLift,sourceCoefficients,modeReceiver\",\"azimuthal_phase\":\"exp(+i*m*phi)\",\"angular_convention\":\"GSHTrans Wigner d^l_{n,m}, unnormalized solver convention; direct samples are global-frame values\"},\"candidate_s\":[";
  for(std::size_t i=0;i<candidateSeconds.size();++i)meta<<(i?",":"")<<candidateSeconds[i];
  meta<<"],\"reference_s\":[";
  for(std::size_t i=0;i<referenceSeconds.size();++i)meta<<(i?",":"")<<referenceSeconds[i];
  meta<<"],\"candidate_stage_s\":[";
  for(std::size_t i=0;i<candidateStages.size();++i){if(i)meta<<",";meta<<"{\"measurement_adjoint\":"<<candidateStages[i][0]<<",\"modal_solve_and_basis\":"<<candidateStages[i][1]<<",\"packed_sampling\":"<<candidateStages[i][2]<<",\"final_kernel_contraction\":"<<candidateStages[i][3]<<",\"source_contact\":"<<candidateStages[i][4]<<"}";}
  meta<<"],\"packet_export_enabled\":"<<(exportPacket?"true":"false")<<",\"modal_export_enabled\":"<<(exportModal?"true":"false")<<",\"batch_responses_enabled\":"<<(batchResponsesEnabled?"true":"false")<<",\"batch_cache_entries\":"<<batchResponses.entryCount()<<",\"batch_cache_bytes\":"<<batchResponses.retainedBytes()<<",\"batch_cache_hits\":"<<batchResponses.cacheHits()<<",\"batch_cache_misses\":"<<batchResponses.cacheMisses()<<",\"batch_max_basis_reconstruction_relative\":"<<batchResponses.maximumBasisReconstructionRelative()<<",\"contacts_enabled\":"<<(!(std::getenv("DSFAST_CONTACTS")&&std::string(std::getenv("DSFAST_CONTACTS"))=="0")?"true":"false")<<",\"rstrain_bridge_enabled\":"<<(std::getenv("DSFAST_RSTRAIN_BRIDGE")&&std::string(std::getenv("DSFAST_RSTRAIN_BRIDGE"))=="1"?"true":"false")<<",\"acceptance\":{\"field_limit\":1e-12,\"kernel_l2_limit\":1e-12,\"pointwise_scaled_limit\":1e-10,\"contact_limit\":1e-12,\"objective_limit\":1e-12}}\n";
  std::cout<<"{\"status\":\"done\",\"degree\":"<<lmax<<",\"nfreq\":"<<nfreq
      <<",\"objectives\":"<<nobj<<",\"max_field_relative_error\":"<<maxFieldErr
      <<",\"max_kernel_l2_relative_error\":"<<maxKernelL2Err<<",\"max_kernel_pointwise_scaled_error\":"<<maxKernelPointErr<<",\"max_contact_relative_error\":"<<maxContactErr<<",\"max_objective_relative_error\":"<<maxObjectiveErr<<"}\n";
  return 0;
}catch(const std::exception &e){std::cerr<<"fast diagnostic: "<<e.what()<<'\n';return 1;}
