#include "engine_headers/Windowed3DKernels.h"
#include "profile.hpp"
#include <DSpecM1D/All>
#include <DSpecM1D/ModelInput>
#include <DSpecM1D/src/InputParser.h>
#include <DSpecM1D/src/NormClass.h>
#include <DSpecM1D/src/SourceInfo.h>
#include <GSHTrans/GaussLegendreGrid.hpp>
#include <GSHTrans/Indexing.hpp>
#include <GSHTrans/Wigner.hpp>
#include <GaussQuad/All>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numbers>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/sysinfo.h>
#include <vector>

namespace fs = std::filesystem;
using Complex = std::complex<double>;
using Clock = std::chrono::steady_clock;
constexpr double kPi = std::numbers::pi;

namespace {
struct Options {
  fs::path model, output;
  int lmax=12, nq=6, qex=1, frequencyBin=2;
  std::size_t radialSideIndex=0;
  double radialRefinement=0.05, maxMemoryGiB=8.0, maxOutputGiB=8.0;
  bool reference=false;
};

Options parse(int argc,char **argv) {
  if(argc<3)throw std::runtime_error("usage: physical_global_gsh MODEL EMPTY_OUTPUT --lmax L --nq N --radial-side-index I [--radial-refinement R] [--qex Q] [--frequency-bin K] [--reference on|off] [--max-memory-gib G] [--max-output-gib G]");
  Options o;o.model=fs::absolute(argv[1]);o.output=fs::absolute(argv[2]);
  for(int i=3;i<argc;++i){
    if(i+1>=argc)throw std::runtime_error("option is missing its value");
    const std::string key=argv[i++],value=argv[i];
    if(key=="--lmax")o.lmax=std::stoi(value);
    else if(key=="--nq")o.nq=std::stoi(value);
    else if(key=="--qex")o.qex=std::stoi(value);
    else if(key=="--frequency-bin")o.frequencyBin=std::stoi(value);
    else if(key=="--radial-side-index")o.radialSideIndex=static_cast<std::size_t>(std::stoull(value));
    else if(key=="--radial-refinement")o.radialRefinement=std::stod(value);
    else if(key=="--max-memory-gib")o.maxMemoryGiB=std::stod(value);
    else if(key=="--max-output-gib")o.maxOutputGiB=std::stod(value);
    else if(key=="--reference"){if(value!="on"&&value!="off")throw std::runtime_error("--reference must be on or off");o.reference=value=="on";}
    else throw std::runtime_error("unknown option: "+key);
  }
  if(!fs::exists(o.model))throw std::runtime_error("MODEL path does not exist");
  if(o.lmax<1||o.lmax>1536||o.nq<2||o.qex<1||o.frequencyBin<0||o.radialRefinement<=0.0||
     !(std::isfinite(o.maxMemoryGiB)&&o.maxMemoryGiB>0.0)||!(std::isfinite(o.maxOutputGiB)&&o.maxOutputGiB>0.0))
    throw std::runtime_error("degree, SEM points, frequency, and resource limits are invalid");
  fs::create_directories(o.output);
  if(fs::directory_iterator(o.output)!=fs::directory_iterator{})throw std::runtime_error("EMPTY_OUTPUT must be empty");
  return o;
}

struct StageTimes {
  std::map<std::string,double> seconds;
  template<class Fn> auto time(const std::string &name,Fn &&fn) {
    const auto start=Clock::now();
    if constexpr(std::is_void_v<std::invoke_result_t<Fn>>) {
      fn();seconds[name]+=std::chrono::duration<double>(Clock::now()-start).count();
    } else {
      auto result=fn();seconds[name]+=std::chrono::duration<double>(Clock::now()-start).count();return result;
    }
  }
};

struct OptionsParams {
  static constexpr double sourceLat=10.0,sourceLon=20.0,receiverLat=45.0,receiverLon=75.0;
  static constexpr double sourceDepthKm=50.0,receiverDepthKm=3.0;
  static constexpr double dt=5.0,duration=3600.0;
};

std::string paramsText(const Options &o) {
  std::ostringstream s;
  s<<'"'<<(o.output/"unused").string()<<"\"\n\""<<o.model.string()<<"\"\n"
   <<"4\n0\n2\n0\n0\n1e-5\n0\n"<<o.lmax<<"\n0\n10\n"
   <<OptionsParams::duration/60.0<<'\n'<<OptionsParams::dt<<"\n0\n0\n0\n0\n"
   <<OptionsParams::sourceDepthKm<<'\n'<<OptionsParams::sourceLat<<'\n'<<OptionsParams::sourceLon<<'\n'
   <<"1e19\n2e18\n-3e18\n-5e18\n4e18\n-5e18\n"
   <<OptionsParams::receiverDepthKm<<"\n1\n"<<OptionsParams::receiverLat<<' '<<OptionsParams::receiverLon<<"\n";
  return s.str();
}

struct NativeSide { int element=-1,node=-1;double radius=0.0,width=0.0;bool fluid=false; };
std::vector<NativeSide> nativeSides(const Full1D::SEM &sem) {
  std::vector<NativeSide> out;
  for(int e=0;e<sem.mesh().NE();++e)for(int q=0;q<sem.mesh().NN();++q){
    const double r=sem.mesh().NodeRadius(e,q)*sem.lengthScaleMetres();
    if(r>0.0)out.push_back({e,q,r,(sem.mesh().EUR(e)-sem.mesh().ELR(e))*sem.lengthScaleMetres(),sem.elementIsFluid(e)});
  }
  return out;
}

class NativeComplete2LGrid final : public GSHTrans::SphericalGrid<double,GSHTrans::All,GSHTrans::All> {
  using Base=GSHTrans::SphericalGrid<double,GSHTrans::All,GSHTrans::All>;
  static std::pair<std::vector<double>,std::vector<double>> nodes(int lmax) {
    auto q=GaussQuad::Quadrature1D<double>(GaussQuad::LegendrePolynomial<double>{}.GaussQuadrature(2*lmax+1));
    q.Transform([](double x){return std::acos(-x);},[](double){return 1.0;});
    return {q.Points(),q.Weights()};
  }
 public:
  explicit NativeComplete2LGrid(int lmax)
      :NativeComplete2LGrid(lmax,nodes(lmax)){}
 private:
  NativeComplete2LGrid(int lmax,std::pair<std::vector<double>,std::vector<double>> q)
      :Base(lmax,2,std::move(q.first),std::move(q.second),FFTWpp::Measure,
            GSHTrans::Chunking::Automatic(),GSHTrans::WignerValues::Generated(),
            GSHTrans::TransformKernel::Loop(),GSHTrans::FastFFTSize(4*lmax+1)){}
};

struct Channel { int record=0,component=0,spin=0,slot=0;std::vector<Complex> coefficients; };
constexpr std::array<int,6> kSpin{0,0,1,-1,2,-2};
constexpr std::array<const char*,12> kChannelNames{
  "forward_rr","forward_horizontal_trace","forward_radial_plus","forward_radial_minus","forward_horizontal_plus","forward_horizontal_minus",
  "adjoint_rr","adjoint_horizontal_trace","adjoint_radial_plus","adjoint_radial_minus","adjoint_horizontal_plus","adjoint_horizontal_minus"};

struct RecordPair {
  SPARSESPEC::PoleAngularCoefficientRecord3D pole[2];
  bool hasPole[2]{false,false};
  Complex radialL0[2][2]{};
  bool hasRadialL0[2]{false,false};
};

std::pair<double,double> frameAngles(const SPARSESPEC::PoleAngularCoefficientRecord3D &r) {
  const auto &radial=r.frame.axes[2];
  return {std::acos(std::clamp(radial[2],-1.0,1.0)),std::atan2(radial[1],radial[0])};
}

std::vector<Channel> rotatePoleCoefficients(const RecordPair &records,int lmax,
                                             std::array<GSHTrans::Wigner<double,GSHTrans::All,GSHTrans::All,GSHTrans::Single>,2> &rotations,
                                             StageTimes &times) {
  constexpr std::array<int,2> scopes{0,1};
  std::vector<Channel> channels;channels.reserve(12);
  const auto start=Clock::now();
  for(int ri:scopes){
    const auto &rec=records.pole[ri];
    std::map<std::pair<int,int>,std::array<Complex,6>> local;
    for(std::size_t i=0;i<rec.labels.size();++i)local.emplace(rec.labels[i],rec.rawStrain[i]);
    const auto [theta0,longitude]=frameAngles(rec);(void)theta0;
    for(int comp=0;comp<6;++comp){
      const int spin=kSpin[static_cast<std::size_t>(comp)];
      GSHTrans::GSHIndices<GSHTrans::All> indices(lmax,lmax,spin);
      Channel ch;ch.record=ri;ch.component=comp;ch.spin=spin;ch.slot=ri*6+comp;
      ch.coefficients.assign(static_cast<std::size_t>(indices.Size()),Complex{});
      for(int l=std::abs(spin);l<=lmax;++l){
        const double norm=std::sqrt((2.0*l+1.0)/(4.0*kPi));
        for(int m=-l;m<=l;++m){
          Complex sum{};
          for(int mp=-std::min(2,l);mp<=std::min(2,l);++mp){
            const auto it=local.find({mp,l});if(it==local.end())continue;
            // The Wigner upper index is the local harmonic order m-prime,
            // not the strain spin. The lower index is the global order m.
            sum+=rotations[ri][mp][l][m]*it->second[static_cast<std::size_t>(comp)]/norm;
          }
          ch.coefficients[static_cast<std::size_t>(indices.Index(l,m))]=
              std::exp(Complex(0.0,-m*longitude))*sum;
        }
      }
      channels.push_back(std::move(ch));
    }
  }
  times.seconds["coefficient_rotation_and_global_layout"]+=std::chrono::duration<double>(Clock::now()-start).count();
  return channels;
}

std::uint64_t availableRamBytes() {
  std::ifstream f("/proc/meminfo");std::string key,unit;std::uint64_t value=0;
  while(f>>key>>value>>unit)if(key=="MemAvailable:")return value*1024;
  struct sysinfo info{};if(sysinfo(&info)!=0)return 0;
  return static_cast<std::uint64_t>(info.freeram)*info.mem_unit;
}

void writeRaw(std::fstream &f,std::uint64_t offset,const Complex *data,std::size_t n) {
  f.seekp(static_cast<std::streamoff>(offset));
  f.write(reinterpret_cast<const char*>(data),static_cast<std::streamsize>(n*sizeof(Complex)));
  if(!f)throw std::runtime_error("binary field write failed");
}

double relativeL2(const std::vector<Complex> &a,const std::vector<Complex> &b) {
  long double num=0.0L,den=0.0L;
  for(std::size_t i=0;i<a.size();++i){num+=std::norm(a[i]-b[i]);den+=std::norm(b[i]);}
  return std::sqrt(static_cast<double>(num/std::max(den,std::numeric_limits<long double>::min())));
}
}

int main(int argc,char **argv) try {
  mathprofile::Root top("global_gsh_probe",-1);
  const auto started=Clock::now();StageTimes times;
  const Options o=parse(argc,argv);
  const std::size_t nTheta=static_cast<std::size_t>(2*o.lmax+1);
  const std::size_t nPhi=static_cast<std::size_t>(GSHTrans::FastFFTSize(4*o.lmax+1));
  if(nTheta>std::numeric_limits<std::size_t>::max()/nPhi)throw std::runtime_error("angular field size overflows size_t");
  const std::size_t fieldSize=nTheta*nPhi;
  const std::uint64_t fieldBytes=static_cast<std::uint64_t>(fieldSize)*sizeof(Complex);
  const std::uint64_t gshOutputBytes=12*fieldBytes;
  const std::uint64_t referenceBytes=o.reference?12*fieldBytes:0;
  const std::uint64_t coeffSize=static_cast<std::uint64_t>((o.lmax+1)*(o.lmax+1));
  const std::uint64_t rotationBytes=2*5*coeffSize*sizeof(double);
  const std::uint64_t transformBatchFields=o.lmax<=32?4:1;
  const std::uint64_t knownArrayEstimate=rotationBytes+transformBatchFields*fieldBytes+
      12*coeffSize*sizeof(Complex)+(o.reference?fieldSize*(sizeof(SPARSESPEC::KernelField3DSample)+sizeof(SPARSESPEC::Field3DPoint)):0);
  const double gib=1024.0*1024.0*1024.0;
  const auto availRam=availableRamBytes();
  if(knownArrayEstimate>o.maxMemoryGiB*gib|| (availRam&&knownArrayEstimate>availRam))
    throw std::runtime_error("memory preflight rejected known prototype arrays; lower L/reference or raise --max-memory-gib after review");
  if(gshOutputBytes+referenceBytes>o.maxOutputGiB*gib)
    throw std::runtime_error("disk preflight rejected planned global field output; raise --max-output-gib after review");

  const auto freeBefore=fs::space(o.output).available;
  if(gshOutputBytes+referenceBytes>freeBefore)
    throw std::runtime_error("disk preflight rejected global field output against available bytes");

  fs::path params=o.output/"parameters.txt";
  std::unique_ptr<InputParametersNew> input;
  std::unique_ptr<Full1D::SEM> sem;
  std::vector<NativeSide> sides;
  times.time("DSpec_model_frequency_and_SEM_setup",[&]{
    {std::ofstream f(params);f<<paramsText(o);}
    input=std::make_unique<InputParametersNew>(params.string(),o.nq,1,o.radialRefinement,1.0,0.1,0.0,o.qex);
    auto &base=input->earthModel();auto &source=input->cmt();
    sem=std::make_unique<Full1D::SEM>(base,o.radialRefinement,o.nq,o.lmax,source);
    sides=nativeSides(*sem);
  });
  if(o.frequencyBin>=input->freqFull().nt()/2+1)throw std::runtime_error("frequency bin is outside the resolved one-sided spectrum");
  if(o.radialSideIndex>=sides.size())throw std::runtime_error("radial-side index exceeds available positive-radius native SEM sides");
  const NativeSide side=sides[o.radialSideIndex];

  NativeComplete2LGrid grid(o.lmax);
  bool parityChecked=false;double maxThetaError=0.0,maxWeightError=0.0;bool parityCompatible=false;
  if(o.lmax<=64){
    auto referenceGrid=GSHTrans::GaussLegendreGrid<double,GSHTrans::All,GSHTrans::All>::ForBand(
        o.lmax,2,2.0,FFTWpp::Estimate,GSHTrans::Chunking::Automatic(),
        GSHTrans::WignerValues::Generated(),GSHTrans::TransformKernel::Loop());
    parityCompatible=grid.NumberOfCoLatitudes()==referenceGrid.NumberOfCoLatitudes()&&
        grid.NumberOfLongitudes()==referenceGrid.NumberOfLongitudes()&&grid.FieldSize()==referenceGrid.FieldSize();
    auto a=grid.CoLatitudes().begin(),b=referenceGrid.CoLatitudes().begin();
    auto aw=grid.CoLatitudeWeights().begin(),bw=referenceGrid.CoLatitudeWeights().begin();
    for(;a!=grid.CoLatitudes().end()&&b!=referenceGrid.CoLatitudes().end();++a,++b,++aw,++bw){
      maxThetaError=std::max(maxThetaError,std::abs(*a-*b));maxWeightError=std::max(maxWeightError,std::abs(*aw-*bw));
    }
    parityCompatible=parityCompatible&&maxThetaError==0.0&&maxWeightError==0.0;
    parityChecked=true;
    if(!parityCompatible)throw std::runtime_error("custom band-L grid points/weights do not exactly match ForBand(L,2,2.0)");
  }

  const auto &freq=input->freqFull();
  const Complex z=freq.w().at(o.frequencyBin)+Complex(0.0,-freq.ep());
  Eigen::VectorXcd weights=Eigen::VectorXcd::Zero(3*input->inputParameters().num_receivers());
  weights[0]=Complex(0.7,0.2);
  std::optional<SPARSESPEC::Windowed3DFrequencyWorkspace> workspace;
  times.time("frequency_workspace_prepare_modal_solves",[&]{
    workspace.emplace(SPARSESPEC::Windowed3DFieldEngine::prepareFrequencyWorkspace(
        *sem,input->cmt(),input->inputParameters(),z,weights,nullptr,true,false));
  });

  const double thetaProbe=static_cast<double>(*grid.CoLatitudes().begin());
  const SPARSESPEC::Field3DPoint representative{side.radius,thetaProbe,0.0,side.element};
  RecordPair records;
  SPARSESPEC::Windowed3DFieldEngine::SamplingPlanCache samplePlan;
  times.time("pole_coefficient_record_capture",[&]{
    const std::vector<SPARSESPEC::Field3DPoint> one{representative};
    (void)workspace->sampleForKernelsPacked(one,0,&samplePlan,2,[&](const SPARSESPEC::PoleAngularCoefficientRecord3D &r){
      if(r.synthesisSet==0&&!r.adjoint){records.pole[0]=r;records.hasPole[0]=true;}
      else if(r.synthesisSet==1&&r.adjoint){records.pole[1]=r;records.hasPole[1]=true;}
      else if(r.synthesisSet==std::numeric_limits<std::uint64_t>::max()&&!r.adjoint){
        records.radialL0[0][0]=r.rawStrain.front()[0];records.radialL0[0][1]=r.rawStrain.front()[1];records.hasRadialL0[0]=true;
      }else if(r.synthesisSet==std::numeric_limits<std::uint64_t>::max()-1&&r.adjoint){
        records.radialL0[1][0]=r.rawStrain.front()[0];records.radialL0[1][1]=r.rawStrain.front()[1];records.hasRadialL0[1]=true;
      }
    });
  });
  if(!records.hasPole[0]||!records.hasPole[1])throw std::runtime_error("expected source and receiver pole coefficient records were not captured");

  std::array<GSHTrans::Wigner<double,GSHTrans::All,GSHTrans::All,GSHTrans::Single>,2> rotations;
  times.time("rotation_wigner_table_generation",[&]{
    for(int ri=0;ri<2;++ri){const auto [theta,phi]=frameAngles(records.pole[ri]);(void)phi;
      rotations[ri]=GSHTrans::Wigner<double,GSHTrans::All,GSHTrans::All,GSHTrans::Single>(o.lmax,o.lmax,2,theta);}
  });
  auto channels=rotatePoleCoefficients(records,o.lmax,rotations,times);

  std::vector<SPARSESPEC::KernelField3DSample> reference;
  std::vector<SPARSESPEC::Field3DPoint> points;
  if(o.reference){
    times.time("reference_native_point_generation",[&]{
      points.reserve(fieldSize);
      const double dPhi=2.0*kPi/nPhi;
      for(std::size_t it=0;it<nTheta;++it){const double theta=grid.CoLatitudes()[static_cast<std::ptrdiff_t>(it)];
        for(std::size_t ip=0;ip<nPhi;++ip)points.push_back({side.radius,theta,dPhi*ip,side.element});}
    });
    const auto receiverElement=sem->receiverElements(input->inputParameters(),false).at(0);
    const auto station=input->inputParameters().receivers().at(0);
    const double receiverRadius=(sem->mesh().PR()-1000.0*input->inputParameters().receiver_depth()/sem->lengthScaleMetres())*sem->lengthScaleMetres();
    points.push_back({receiverRadius,(90.0-station.first)*kPi/180.0,station.second*kPi/180.0,receiverElement});
    auto all=times.time("accepted_pole_tensor_reference_sampling",[&]{return workspace->sampleForKernelsTensor(points,fieldSize,&samplePlan,2,{},false);});
    all.resize(fieldSize);reference=std::move(all);
  }

  fs::path fieldPath=o.output/"gsh_fields.bin";
  {std::ofstream f(fieldPath,std::ios::binary|std::ios::trunc);f.seekp(static_cast<std::streamoff>(gshOutputBytes-1));f.put('\0');if(!f)throw std::runtime_error("could not size GSH field output");}
  std::fstream fieldOut(fieldPath,std::ios::binary|std::ios::in|std::ios::out);
  std::fstream referenceOut;
  if(o.reference){
    fs::path p=o.output/"reference.bin";
    {std::ofstream f(p,std::ios::binary|std::ios::trunc);f.seekp(static_cast<std::streamoff>(referenceBytes-1));f.put('\0');if(!f)throw std::runtime_error("could not size reference output");}
    referenceOut.open(p,std::ios::binary|std::ios::in|std::ios::out);
  }
  std::array<double,12> channelRel{};
  std::array<double,12> channelMaxError{},channelMaxReference{};
  std::array<std::array<double,6>,2> kernelMaxError{};
  std::array<std::array<double,6>,2> kernelMaxReference{};
  const bool batchLowL=o.lmax<=32;
  times.time("native_global_gsh_inverse_transforms",[&]{
    std::map<int,std::vector<const Channel*>> groups;
    for(const auto &ch:channels)groups[ch.spin].push_back(&ch);
    for(const auto &[spin,group]:groups){
      const std::size_t batchLimit=batchLowL?group.size():1;
      const std::size_t coefficientSize=group.front()->coefficients.size();
      for(std::size_t begin=0;begin<group.size();begin+=batchLimit){
        const std::size_t count=std::min(batchLimit,group.size()-begin);
        std::vector<Complex> coeff(count*coefficientSize),fields(count*fieldSize);
        for(std::size_t i=0;i<count;++i)
          std::copy(group[begin+i]->coefficients.begin(),group[begin+i]->coefficients.end(),coeff.begin()+static_cast<std::ptrdiff_t>(i*coefficientSize));
        grid.InverseTransformation(o.lmax,spin,coeff,GSHTrans::Batch::Contiguous(static_cast<std::ptrdiff_t>(count),static_cast<std::ptrdiff_t>(coefficientSize)),
            fields,GSHTrans::Batch::Contiguous(static_cast<std::ptrdiff_t>(count),static_cast<std::ptrdiff_t>(fieldSize)),GSHTrans::Execution::Sequential());
        for(std::size_t i=0;i<count;++i){
          const Channel &ch=*group[begin+i];
          Complex *values=fields.data()+static_cast<std::ptrdiff_t>(i*fieldSize);
          if(ch.component<2&&records.hasRadialL0[ch.record])
            for(std::size_t p=0;p<fieldSize;++p)values[p]+=records.radialL0[ch.record][ch.component];
          writeRaw(fieldOut,static_cast<std::uint64_t>(ch.slot)*fieldBytes,values,fieldSize);
          if(o.reference){
            std::vector<Complex> ref(fieldSize);
            for(std::size_t p=0;p<fieldSize;++p){
              const auto &s=reference[p];const auto &strain=ch.record==0?s.forwardRawStrain:s.adjointDualStrain;
              const std::array<Complex,6> x{strain.rr,strain.horizontalTrace,strain.radialPlus,strain.radialMinus,strain.horizontalPlus,strain.horizontalMinus};
              ref[p]=x[static_cast<std::size_t>(ch.component)];
            }
            writeRaw(referenceOut,static_cast<std::uint64_t>(ch.slot)*fieldBytes,ref.data(),fieldSize);
            channelRel[static_cast<std::size_t>(ch.slot)]=relativeL2(std::vector<Complex>(values,values+fieldSize),ref);
            for(std::size_t p=0;p<fieldSize;++p){
              channelMaxError[static_cast<std::size_t>(ch.slot)]=std::max(channelMaxError[static_cast<std::size_t>(ch.slot)],std::abs(values[p]-ref[p]));
              channelMaxReference[static_cast<std::size_t>(ch.slot)]=std::max(channelMaxReference[static_cast<std::size_t>(ch.slot)],std::abs(ref[p]));
            }
          }
        }
      }
    }
  });
  fieldOut.close();
  if(o.reference)referenceOut.close();

  if(o.reference){
    std::ifstream candidateInput(fieldPath,std::ios::binary);
    if(!candidateInput)throw std::runtime_error("cannot reopen GSH fields for pointwise kernel comparison");
    times.time("pointwise_modulus_kernel_contraction_and_comparison",[&]{
      for(std::size_t p=0;p<fieldSize;++p){
        auto candidate=reference[p];
        std::array<Complex,12> v{};
        // Read the 12 already-written candidate channels. This keeps only one point in memory.
        for(int slot=0;slot<12;++slot){candidateInput.seekg(static_cast<std::streamoff>(static_cast<std::uint64_t>(slot)*fieldBytes+p*sizeof(Complex)));candidateInput.read(reinterpret_cast<char*>(&v[static_cast<std::size_t>(slot)]),sizeof(Complex));if(!candidateInput)throw std::runtime_error("failed reading GSH point for kernel comparison");}
        candidate.forwardRawStrain={v[0],v[1],v[2],v[3],v[4],v[5]};
        candidate.adjointDualStrain={v[6],v[7],v[8],v[9],v[10],v[11]};
        const auto r=SPARSESPEC::Windowed3DKernelEngine::pointKernel(*sem,reference[p]);
        const auto c=SPARSESPEC::Windowed3DKernelEngine::pointKernel(*sem,candidate);
        const std::array<double,6> rv{r.A,r.C,r.F,r.L,r.N,r.kappa},cv{c.A,c.C,c.F,c.L,c.N,c.kappa};
        for(std::size_t j=0;j<6;++j){kernelMaxError[0][j]=std::max(kernelMaxError[0][j],std::abs(cv[j]-rv[j]));kernelMaxReference[0][j]=std::max(kernelMaxReference[0][j],std::abs(rv[j]));}
      }
    });
  }

  {std::ofstream f(o.output/"result.json");f<<std::setprecision(17);
   f<<"{\"mode\":\"exact_native_complete2L_points_bandL_GSH_probe\",\"degree\":"<<o.lmax<<",\"nq_gll_points_per_element\":"<<o.nq<<",\"polynomial_degree\":"<<o.nq-1;
   f<<",\"qex\":"<<o.qex<<",\"frequency_bin\":"<<o.frequencyBin<<",\"z_real\":"<<z.real()<<",\"z_imag\":"<<z.imag();
   f<<",\"adjoint_weight_vertical\":[0.7,0.2],\"source_lat_lon_depth_m\":[10,20,50000],\"receiver_lat_lon_depth_m\":[45,75,3000]";
   f<<",\"radial_side_index\":"<<o.radialSideIndex<<",\"radial_side_count\":"<<sides.size()<<",\"radial_element\":"<<side.element<<",\"radial_node\":"<<side.node<<",\"radius_m\":"<<side.radius<<",\"fluid\":"<<(side.fluid?"true":"false");
   f<<",\"n_theta\":"<<nTheta<<",\"n_phi\":"<<nPhi<<",\"field_points\":"<<fieldSize<<",\"native_grid_contract\":\"same theta nodes/weights and theta-major uniform longitudes as ForBand(L,2,2.0); transform declares band L using Generated Wigner and Loop kernel\"";
   f<<",\"native_grid_parity_checked\":"<<(parityChecked?"true":"false")<<",\"native_grid_parity_passed\":"<<(parityCompatible?"true":"false")<<",\"max_theta_abs_error\":"<<maxThetaError<<",\"max_weight_abs_error\":"<<maxWeightError;
   f<<",\"batch_low_degree\":"<<(batchLowL?"true":"false")<<",\"gsh_fields_file\":\"gsh_fields.bin; 12 channel-major complex128 arrays, source six then receiver-adjoint six, theta-major points\",\"reference_enabled\":"<<(o.reference?"true":"false")<<",\"reference_file\":\"reference.bin; same layout, accepted pole/tensor all-point result\"";
   f<<",\"gsh_output_bytes\":"<<gshOutputBytes<<",\"reference_output_bytes\":"<<referenceBytes<<",\"planned_output_bytes\":"<<gshOutputBytes+referenceBytes<<",\"output_cap_bytes\":"<<static_cast<std::uint64_t>(o.maxOutputGiB*gib)<<",\"available_disk_bytes_at_preflight\":"<<freeBefore;
   f<<",\"known_driver_array_estimate_bytes\":"<<knownArrayEstimate<<",\"max_memory_gib\":"<<o.maxMemoryGiB<<",\"available_ram_bytes_at_preflight\":"<<availRam<<",\"memory_estimate_note\":\"known arrays only; excludes SEM/modal solver, FFTW, and GSH internal scratch; not peak RSS\"";
   f<<",\"rotation_formula\":\"global(l,m)=exp(-i*m*longitude)*sum_local_mprime(d_l[mprime,m]*local(l,mprime)); d is obtained by dividing GSHTrans Wigner values by sqrt((2l+1)/(4pi))\",\"radial_l0_included\":true,\"polar_conventions_runtime_validated\":false";
   f<<",\"channel_relative_l2_errors\":[";
   for(int i=0;i<12;++i){if(i)f<<',';if(o.reference)f<<channelRel[static_cast<std::size_t>(i)];else f<<"null";}
   f<<"],\"pointwise_kernel_max_abs_error\":[";for(int i=0;i<6;++i){if(i)f<<',';if(o.reference)f<<kernelMaxError[0][static_cast<std::size_t>(i)];else f<<"null";}
   f<<"],\"pointwise_kernel_max_abs_reference\":[";for(int i=0;i<6;++i){if(i)f<<',';if(o.reference)f<<kernelMaxReference[0][static_cast<std::size_t>(i)];else f<<"null";}
   f<<"],\"channel_max_abs_error\":[";for(int i=0;i<12;++i){if(i)f<<',';if(o.reference)f<<channelMaxError[static_cast<std::size_t>(i)];else f<<"null";}
   f<<"],\"channel_max_abs_reference\":[";for(int i=0;i<12;++i){if(i)f<<',';if(o.reference)f<<channelMaxReference[static_cast<std::size_t>(i)];else f<<"null";}
   f<<"],\"stage_timings_s\":{";bool first=true;for(const auto &[name,value]:times.seconds){f<<(first?"":",")<<'\"'<<name<<"\":"<<value;first=false;}
   f<<"},\"elapsed_s\":"<<std::chrono::duration<double>(Clock::now()-started).count()<<"}\n";}
  fs::remove(params);
  return 0;
}catch(const std::exception &e){std::cerr<<"FAIL: "<<e.what()<<'\n';return 1;}
