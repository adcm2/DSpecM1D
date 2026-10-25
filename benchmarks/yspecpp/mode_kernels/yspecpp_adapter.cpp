#include <YSpec/Kernels.hpp>
#include <YSpec/Synthetics.hpp>
#include "Spectra.hpp"
#include "TimeSeries.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs=std::filesystem;
using namespace YSpec;
constexpr double duration=72*3600.0, dt=60.0, sourceDepth=50000.0, receiverDepth=3000.0;
Model alter(const Model& model,double h,const std::string& property) {
  std::vector<Layer> layers; double mass=0;
  for(const auto& old:model.Layers()) {
    std::vector<double> radii(old.Knots().begin(),old.Knots().end());
    std::vector<Properties> values; values.reserve(radii.size());
    for(double r:radii) {
      const auto m=old(r); Properties p;
      p[Property::Rho]=m.rho*(property=="rho"?1+h:1);
      p[Property::A]=m.A*(property=="A"?1+h:1); p[Property::C]=m.C*(property=="C"?1+h:1);
      p[Property::F]=m.F*(property=="F"?1+h:1);
      p[Property::L]=m.L*((property=="L"||property=="LN")?1+h:1);
      p[Property::N]=m.N*((property=="N"||property=="LN")?1+h:1);
      p[Property::Kappa]=(p[Property::C]+4*p[Property::A]-4*p[Property::N]+4*p[Property::F])/9;
      p[Property::Mu]=(p[Property::C]+p[Property::A]+6*p[Property::L]+5*p[Property::N]-2*p[Property::F])/15;
      p[Property::QKappa]=p[Property::QMu]=p[Property::QA]=p[Property::QF]=0;
      values.push_back(p);
    }
    layers.emplace_back(std::move(radii),values,old.GetPhase(),mass);
    mass=old.MassBelowTop()*(property=="rho"?1+h:1);
  }
  return Model(model.GetScales(),model.ReferenceFrequency(),std::move(layers));
}
void guard(const fs::path& p) {
  if(fs::exists(p)&&(!fs::is_directory(p)||fs::directory_iterator(p)!=fs::directory_iterator{}))
    throw std::runtime_error("output directory must be empty: "+p.string());
  fs::create_directories(p);
}
std::map<std::string,std::vector<std::array<double,3>>> readQ(const fs::path& p) {
  std::ifstream in(p); if(!in) throw std::runtime_error("cannot read q weights");
  std::string line; std::getline(in,line); std::map<std::string,std::vector<std::array<double,3>>> q;
  while(std::getline(in,line)) { std::replace(line.begin(),line.end(),',',' '); std::istringstream row(line); std::string name; int k; double re,im; if(row>>name>>k>>re>>im) q[name].push_back({double(k),re,im}); }
  return q;
}
int main(int argc,char**argv) try {
  if(argc!=10&&argc!=12) throw std::runtime_error("usage: yspecpp_adapter ACTION MODEL OUT FMIN_MHZ FMAX_MHZ DELTA TOLERANCE QCSV MESH_M [MODE PROPERTY]");
  const std::string action=argv[1]; const fs::path modelPath=fs::absolute(argv[2]),out=fs::absolute(argv[3]);
  const double fmin=std::stod(argv[4])*1e-3, fmax=std::stod(argv[5])*1e-3;
  const double delta=std::stod(argv[6]), tolerance=std::stod(argv[7]), maxElement=std::stod(argv[9]);
  const std::string mode=argc==12?argv[10]:"toroidal", property=argc==12?argv[11]:"LN";
  const bool spheroidal=mode=="spheroidal"; if(!spheroidal&&mode!="toroidal") throw std::runtime_error("mode must be toroidal or spheroidal");
  guard(out);
  auto model=alter(Model::ReadDeck(modelPath),delta,property);
  Calculation calc; calc.gravity=Gravity::Full; calc.attenuation=false; calc.quantity=Quantity::Displacement;
  calc.degrees={2,2}; calc.band={fmin,fmax}; calc.sampling={duration,dt};
  calc.solver.tolerance=tolerance; calc.solver.startDecay=11; calc.solver.stepFactor=0.25;
  const Source source{{0,0},sourceDepth,spheroidal?MomentTensor{1e19,0,0,0,0,0}:MomentTensor{0,1e19,0,0,0,0},{}};
  const Receiver receiver{{0,90},receiverDepth,{}};
  const auto grid=FrequencyGrid(fmin*model.GetScales().Time(),fmax*model.GetScales().Time(),
                                dt/model.GetScales().Time(),duration/model.GetScales().Time());
  const std::size_t expectedFirst=spheroidal?146:116, expectedLast=spheroidal?177:158, expectedCount=spheroidal?32:43;
  if(grid.Samples()!=8192 || grid.First()!=expectedFirst || grid.Last()!=expectedLast || grid.Count()!=expectedCount ||
     std::abs(grid.TimeStep()*model.GetScales().Time()-dt)>1e-12 ||
     std::abs(grid.Spacing()/model.GetScales().Time()-1.0/(8192.0*dt))>1e-15 ||
     std::abs(grid.Damping()/model.GetScales().Time()-5.0/duration)>1e-15)
    throw std::runtime_error("unexpected YSpec grid for selected modal band");
  const auto raw=SeismogramSpectra(model,calc,std::span{&source,1},std::span{&receiver,1});
  std::vector<Complex> north(grid.Count());
  const Component component=spheroidal?Vertical:North;
  for(std::size_t j=0;j<grid.Count();++j) north[j]=raw(0,0,component,j);
  const auto time=InverseTransform(grid,std::nullopt,north,1);
  constexpr std::size_t windowSamples=static_cast<std::size_t>(duration/dt)+1;
  auto processed=[&](std::size_t k) {
    Complex total{};
    for(std::size_t n=0;n<windowSamples;++n) {
      const double w=0.5-0.5*std::cos(2*std::numbers::pi*n/(windowSamples-1));
      const double phase=2*std::numbers::pi*double(k*n)/double(grid.Samples());
      total += time[n]*w*Complex(std::cos(phase),-std::sin(phase));
    }
    return total*dt*model.GetScales().Length();
  };
  std::ofstream meta(out/"meta.json");
  meta<<std::setprecision(17)<<"{\"solver\":\"YSpecpp\",\"nt\":"<<grid.Samples()
      <<",\"dt_s\":"<<grid.TimeStep()*model.GetScales().Time()<<",\"df_hz\":"<<grid.Spacing()/model.GetScales().Time()
      <<",\"damping_s\":"<<grid.Damping()/model.GetScales().Time()<<",\"first\":"<<grid.First()
      <<",\"last\":"<<grid.Last()<<",\"source_radius_m\":"<<model.SurfaceRadius()*model.GetScales().Length()-sourceDepth
      <<",\"tolerance\":"<<tolerance<<",\"delta\":"<<delta<<",\"mesh_m\":"<<maxElement
      <<",\"mode\":\""<<mode<<"\",\"property\":\""<<property<<"\"}\n";
  if(action=="spectrum") {
    std::ofstream spec(out/"spectrum.csv"); spec<<"bin,frequency_hz,raw_real,raw_imag,processed_real,processed_imag\n"<<std::setprecision(17);
    for(std::size_t j=0;j<grid.Count();++j) {
      const auto k=grid.First()+j; const auto value=raw(0,0,component,j)*model.GetScales().Length()*model.GetScales().Time();
      const auto proc=processed(k);
      spec<<k<<','<<grid.Frequency(j)/model.GetScales().Time()<<','<<value.real()<<','<<value.imag()<<','<<proc.real()<<','<<proc.imag()<<'\n';
    }
    return 0;
  }
  if(action!="kernel") throw std::runtime_error("unknown action");
  const auto q=readQ(argv[8]); if(q.empty()) throw std::runtime_error("q file has no weights");
  std::vector<Measurement> measurements; std::vector<std::string> names;
  for(const auto& [name,rows]:q) {
    std::vector<double> adjoint(windowSamples,0.0);
    for(std::size_t n=0;n<windowSamples;++n) {
      const double w=0.5-0.5*std::cos(2*std::numbers::pi*n/(windowSamples-1));
      for(const auto& row:rows) { const auto k=static_cast<std::size_t>(row[0]);
        const double phase=2*std::numbers::pi*double(k*n)/double(grid.Samples());
        adjoint[n]+=w*(row[1]*std::cos(phase)-row[2]*std::sin(phase)); }
    }
    Direction direction; if(spheroidal) direction.up=1; else direction.north=1;
    measurements.push_back({receiver,direction,std::move(adjoint)}); names.push_back(name);
  }
  const auto kernels=ComputeRadialKernels(model,calc,source,measurements,RadialMesh{5,maxElement});
  std::ofstream summary(out/"kernel_summary.csv");
  summary<<"basis,distributed_A,distributed_C,distributed_F,distributed_L,distributed_N,distributed_rho,source_A,source_C,source_F,source_L,source_N,total\n"<<std::setprecision(17);
  for(std::size_t j=0;j<kernels.size();++j) {
    const auto& kernel=kernels[j]; double distributed[6]={};
    std::ofstream curve(out/("kernel_"+names[j]+".csv"));
    curve<<"radius_m,radial_A,radial_C,radial_F,radial_L,radial_N,radial_rho,weight_m3\n"<<std::setprecision(17);
    const auto radii=kernel.Radii(), weights=kernel.Weights();
    for(std::size_t i=0;i<kernel.Size();++i) {
      const double r=radii[i], rn=r/model.GetScales().Length();
      const auto material=model(model.LayerIndex(rn,Side::Above),rn);
      const double L=material.L*model.GetScales().Modulus(),N=material.N*model.GetScales().Modulus();
      const double A=material.A*model.GetScales().Modulus(),C=material.C*model.GetScales().Modulus(),F=material.F*model.GetScales().Modulus(),rho=material.rho*model.GetScales().Density();
      const double values[]={kernel[RadialParameter::A][i]*A*r*r,kernel[RadialParameter::C][i]*C*r*r,
        kernel[RadialParameter::F][i]*F*r*r,kernel[RadialParameter::L][i]*L*r*r,
        kernel[RadialParameter::N][i]*N*r*r,kernel[RadialParameter::Rho][i]*rho*r*r};
      const RadialParameter params[]={RadialParameter::A,RadialParameter::C,RadialParameter::F,
        RadialParameter::L,RadialParameter::N,RadialParameter::Rho};
      const double materials[]={A,C,F,L,N,rho};
      for(int p=0;p<6;++p) distributed[p]+=kernel[params[p]][i]*materials[p]*weights[i];
      curve<<r;for(double v:values)curve<<','<<v;curve<<','<<weights[i]<<'\n';
    }
    const auto at=model(model.LayerIndex(kernel.SourceRadius()/model.GetScales().Length(),Side::Above),
                        kernel.SourceRadius()/model.GetScales().Length());
    const double sourceL=kernel.AtSource(RadialParameter::L)*at.L*model.GetScales().Modulus();
    const double sourceN=kernel.AtSource(RadialParameter::N)*at.N*model.GetScales().Modulus();
    const double srcA=kernel.AtSource(RadialParameter::A)*at.A*model.GetScales().Modulus();
    const double srcC=kernel.AtSource(RadialParameter::C)*at.C*model.GetScales().Modulus();
    const double srcF=kernel.AtSource(RadialParameter::F)*at.F*model.GetScales().Modulus();
    summary<<names[j];for(double v:distributed)summary<<','<<v;
    summary<<','<<srcA<<','<<srcC<<','<<srcF<<','<<sourceL<<','<<sourceN<<',';
    const double total=spheroidal ? distributed[0]+distributed[1]+distributed[2]+
      distributed[3]+distributed[4]+distributed[5]+srcA+srcC+srcF+sourceL+sourceN
      : distributed[3]+distributed[4]+sourceL+sourceN;
    summary<<total<<'\n';
  }
  return 0;
} catch(const std::exception&e) { std::cerr<<"YSpec modal adapter: "<<e.what()<<'\n'; return 1; }
