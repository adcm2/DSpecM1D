#include <YSpec/YSpec.hpp>
#include <YSpec/Kernels.hpp>
#include <YSpec/Model.hpp>
#include <YSpec/Seismograms.hpp>
#include "Spectra.hpp"
#include "TimeSeries.hpp"
#include <YSpec/Source.hpp>
#include <YSpec/Receiver.hpp>
#include "TimeSeries.hpp"
#include <YSpec/Units.hpp>
#include "Coefficients.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <sstream>
#include <string>
#include <vector>

namespace fs=std::filesystem;
using namespace YSpec;
constexpr int nt=128, bins=65, nSamples=91;
constexpr double dt=20.0, duration=1800.0;
double geographic(double gc){const double x=gc*std::numbers::pi/180.0; return std::atan(std::tan(x)/0.993277)*180.0/std::numbers::pi;}
double geocentric(double geo){const double x=geo*std::numbers::pi/180.0; return std::atan(0.993277*std::tan(x))*180.0/std::numbers::pi;}
std::vector<double> readRadial(const fs::path&p,std::vector<int>&layers,std::vector<double>&drs,const Model&m){
  std::ifstream f(p); if(!f)throw std::runtime_error("cannot read DSpec radial_nodes.csv"); std::string line; std::getline(f,line); std::vector<double> r;
  while(std::getline(f,line)){std::array<double,4>x{}; char c; std::istringstream s(line); if(!(s>>x[0]>>c>>x[1]>>c>>x[2]>>c>>x[3]))throw std::runtime_error("invalid radial quadrature row"); r.push_back(x[2]); drs.push_back(x[3]); layers.push_back(m.LayerIndex((x[2]/m.GetScales().Length()),Side::Above));}
  if(r.empty())throw std::runtime_error("empty shared DSpec radial quadrature"); return r;
}
struct Angle {double theta,phi,dOmega,chi;};
std::vector<Angle> readAngles(const fs::path&p){std::ifstream f(p);if(!f)throw std::runtime_error("cannot read DSpec angles.csv");std::string line;std::getline(f,line);std::vector<Angle>a;while(std::getline(f,line)){std::replace(line.begin(),line.end(),',',' ');Angle x{};std::istringstream s(line);if(!(s>>x.theta>>x.phi>>x.dOmega>>x.chi))throw std::runtime_error("invalid shared angle row");a.push_back(x);}if(a.size()!=192)throw std::runtime_error("shared angle grid must contain 192 points");return a;}
Measurement measurement(const Receiver&r,int component){std::vector<double> j(nSamples);for(int n=0;n<nSamples;++n){const double w=0.5-0.5*std::cos(2.0*std::numbers::pi*n/(nSamples-1));auto term=[&](int k,double re,double im){const double a=2.0*std::numbers::pi*k*n/nt;return re*std::cos(a)-im*std::sin(a);};j[n]=w*(term(1,.7,.2)+term(2,-.25,.45));}return {r,component==0?Direction{.up=1}:Direction{.north=1},std::move(j)};}
int main(int argc,char**argv) try{
  if(argc!=5)throw std::runtime_error("usage: phase_d_adapter MODEL OUTPUT D_SPEC_RADIAL_NODES D_SPEC_ANGLES");
  const fs::path modelPath=fs::absolute(argv[1]),out=fs::absolute(argv[2]);if(!fs::exists(out)||!fs::is_directory(out)||fs::directory_iterator(out)!=fs::directory_iterator{})throw std::runtime_error("adapter output must be a newly created empty directory");
  const auto model=Model::ReadDeck(modelPath); Calculation calc;calc.gravity=Gravity::Full;calc.attenuation=false;calc.quantity=Quantity::Displacement;calc.degrees={0,3};calc.band={0.1e-3,0.8e-3};calc.sampling={duration,dt};calc.solver.tolerance=1e-9;calc.solver.startDecay=11;calc.solver.stepFactor=.25;
  const auto source=Source{{geographic(10),20},50000,{1e19,2e18,-3e18,-5e18,4e18,-5e18},{}};const auto receiver=Receiver{{geographic(45),90},3000,{}};
  const auto grid=MakeGrid(calc,model.GetScales());if(grid.Samples()!=nt||grid.First()!=1||grid.Last()!=3||std::abs(grid.Damping()/model.GetScales().Time()-1.0/360.0)>1e-14)throw std::runtime_error("resolved YSpec frequency grid differs from frozen grid");
  std::vector<int> layers;std::vector<double> drs;auto radii=readRadial(argv[3],layers,drs,model);auto angles=readAngles(argv[4]);KernelSampling sampling;for(std::size_t i=0;i<radii.size();++i)sampling.nodes.push_back({layers[i],radii[i]});for(const auto&a:angles)sampling.places.push_back({a.theta,a.phi});std::array<KernelSampling,1> samplings{std::move(sampling)};
  std::array<Measurement,2> measurements{measurement(receiver,0),measurement(receiver,1)};
  auto kernels=ComputeSampledKernels(model,calc,source,measurements,ParameterSet::Love,Scaling::Absolute,samplings);
  const auto raw=SeismogramSpectra(model,calc,std::span{&source,1},std::span{&receiver,1});const auto seismograms=ComputeSeismograms(model,calc,std::span{&source,1},std::span{&receiver,1});
  std::ofstream spectra(out/"spectra.csv"),processed(out/"processed.csv"),profiles(out/"profiles.csv"),meta(out/"meta.json");spectra<<"objective,component,bin,frequency_hz,raw_real_m_s,raw_imag_m_s\n"<<std::setprecision(17);processed<<"objective,component,bin,real_m_s,imag_m_s\n"<<std::setprecision(17);profiles<<"objective,modulus,radius_m,theta_rad,phi_rad,chi,w,dV_m3,deltaX_Pa,kernel_per_Pa_m3\n"<<std::setprecision(17);
  constexpr std::array<const char*,3> comps{"Z","N","E"};constexpr std::array<const char*,5> mods{"A","C","F","L","N"};constexpr std::array<double,5>x0{352960000000.0,352960000000.0,129602500000.0,111678750000.0,111678750000.0};
  std::array<std::array<Complex,bins>,3> rawFull{};for(int c=0;c<3;++c)for(int k=0;k<bins;++k)rawFull[c][k]=Complex{};
  for(std::size_t j=0;j<grid.Count();++j){const int k=static_cast<int>(grid.First()+j);for(int c=0;c<3;++c){rawFull[c][k]=raw(0,0,c,j)*model.GetScales().Length()*model.GetScales().Time();}}
  for(int c=0;c<3;++c){const auto tr=seismograms.Trace(0,0,static_cast<Component>(c));for(int k=0;k<bins;++k){Complex value{};for(int n=0;n<nSamples;++n){const double w=.5-.5*std::cos(2.0*std::numbers::pi*n/(nSamples-1));value+=tr[n]*w*std::exp(Complex(0,-2.0*std::numbers::pi*k*n/nt))*dt;}for(int obs=0;obs<2;++obs)processed<<(obs==0?"J_Z":"J_N")<<','<<comps[c]<<','<<k<<','<<value.real()<<','<<value.imag()<<'\n';}}
  for(int obs=0;obs<2;++obs)for(int c=0;c<3;++c)for(int k=0;k<bins;++k)spectra<<(obs==0?"J_Z":"J_N")<<','<<comps[c]<<','<<k<<','<<static_cast<double>(k)*grid.Spacing()/model.GetScales().Time()<<','<<rawFull[c][k].real()<<','<<rawFull[c][k].imag()<<'\n';
  for(int obs=0;obs<2;++obs){const auto&sampled=kernels[obs][0];for(std::size_t ri=0;ri<radii.size();++ri)for(std::size_t ai=0;ai<angles.size();++ai){const double r=radii[ri],theta=angles[ai].theta,phi=angles[ai].phi,wr=std::pow(std::sin(std::numbers::pi*(r-4200000.0)/2000000.0),2);const double dV=r*r*drs[ri]*angles[ai].dOmega;for(int m=0;m<5;++m){auto idx=std::find(sampled.Names().begin(),sampled.Names().end(),mods[m])-sampled.Names().begin();if(idx>=static_cast<std::ptrdiff_t>(sampled.Names().size()))throw std::runtime_error("YSpec sampled kernel missing Love modulus");const double kval=sampled.Values(idx,ri)[ai];profiles<<(obs==0?"J_Z":"J_N")<<','<<mods[m]<<','<<r<<','<<theta<<','<<phi<<','<<angles[ai].chi<<','<<wr<<','<<dV<<','<<x0[m]*wr*angles[ai].chi<<','<<kval<<'\n';}}
  }
  const double sourceRadius=model.Layers().back().Knots().back()*model.GetScales().Length()-source.depth;
  const auto sourceLayer=model.LayerIndex(sourceRadius/model.GetScales().Length(),Side::Above);
  const auto supportInnerLayer=model.LayerIndex(4200000.0/model.GetScales().Length(),Side::Above);
  const auto supportOuterLayer=model.LayerIndex(6200000.0/model.GetScales().Length(),Side::Above);
  meta<<std::setprecision(17)<<"{\"solver\":\"YSpecpp\",\"nt\":"<<grid.Samples()<<",\"dt_s\":"<<grid.TimeStep()*model.GetScales().Time()<<",\"df_hz\":"<<grid.Spacing()/model.GetScales().Time()<<",\"epsilon_s_inv\":"<<grid.Damping()/model.GetScales().Time()<<",\"physical_G\":"<<GravitationalConstant<<",\"radius_m\":"<<model.Layers().back().Knots().back()*model.GetScales().Length()<<",\"source_radius_m\":"<<sourceRadius<<",\"source_layer_upper\":"<<sourceLayer<<",\"regular_support_inner_layer\":"<<supportInnerLayer<<",\"regular_support_outer_layer\":"<<supportOuterLayer<<",\"source_geocentric\":["<<geocentric(source.location.latitude)<<",20,50],\"receiver_geocentric\":["<<geocentric(receiver.location.latitude)<<",90,3],\"moment_tensor_Nm\":[1e19,2e18,-3e18,-5e18,4e18,-5e18],\"gravity\":\"full\",\"attenuation\":false,\"degrees\":[0,3],\"solver_bins\":[1,2,3],\"tolerance\":1e-9,\"start_decay\":11,\"step_factor\":0.25,\"kernel_sample_count\":"<<radii.size()<<"}\n";
  return 0;
}catch(const std::exception&e){std::cerr<<"Phase D YSpec adapter: "<<e.what()<<'\n';return 1;}
