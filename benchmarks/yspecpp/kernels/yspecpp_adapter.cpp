#include <YSpec/Kernels.hpp>
#include <YSpec/Synthetics.hpp>
#include "Spectra.hpp"
#include "TimeSeries.hpp"

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>

namespace fs=std::filesystem;
using namespace YSpec;
enum class Mod { A,C,F,L,N };
constexpr std::array<const char*,5> names{"A","C","F","L","N"};
constexpr std::array<Mod,5> mods{Mod::A,Mod::C,Mod::F,Mod::L,Mod::N};
constexpr std::array<double,6> extendedSteps{0.08,0.04,0.02,0.01,0.005,0.0025};
constexpr int firstBin=1,endBin=4,nSamples=91;
constexpr double dt=20.0, duration=1800.0, r0=2000000.0, r1=4000000.0;
constexpr double q1r=0.7,q1i=0.2,q2r=-0.25,q2i=0.45;
double geo(double geocentric) { const double x=geocentric*std::numbers::pi/180.0; return std::atan(std::tan(x)/0.993277)*180.0/std::numbers::pi; }
double geocentric(double geographic) { const double x=geographic*std::numbers::pi/180.0; return std::atan(0.993277*std::tan(x))*180.0/std::numbers::pi; }
Model alter(const Model& model,Mod selected,std::size_t selectedLayer,double amount) {
  std::vector<Layer> layers; double mass=0;
  for(std::size_t li=0;li<model.Layers().size();++li) {
    const auto& old=model.Layers()[li]; std::vector<double> radii(old.Knots().begin(),old.Knots().end());
    std::vector<Properties> values; values.reserve(radii.size());
    for(double r:radii) {
      const auto m=old(r); Properties p;
      p[Property::Rho]=m.rho; p[Property::A]=m.A*(li==selectedLayer&&selected==Mod::A?1+amount:1);
      p[Property::C]=m.C*(li==selectedLayer&&selected==Mod::C?1+amount:1);
      p[Property::F]=m.F*(li==selectedLayer&&selected==Mod::F?1+amount:1);
      p[Property::L]=m.L*(li==selectedLayer&&selected==Mod::L?1+amount:1);
      p[Property::N]=m.N*(li==selectedLayer&&selected==Mod::N?1+amount:1);
      p[Property::Kappa]=(p[Property::C]+4*p[Property::A]-4*p[Property::N]+4*p[Property::F])/9;
      p[Property::Mu]=(p[Property::C]+p[Property::A]+6*p[Property::L]+5*p[Property::N]-2*p[Property::F])/15;
      p[Property::QKappa]=p[Property::QMu]=p[Property::QA]=p[Property::QF]=0;
      values.push_back(p);
    }
    layers.emplace_back(std::move(radii),values,old.GetPhase(),mass);
    mass=old.MassBelowTop();
  }
  return Model(model.GetScales(),model.ReferenceFrequency(),std::move(layers));
}
Calculation calc(double tol) {
  Calculation c; c.gravity=Gravity::Full; c.attenuation=false; c.quantity=Quantity::Displacement;
  c.degrees={0,2}; c.band={0.1e-3,0.8e-3}; c.sampling={duration,dt};
  c.solver.tolerance=tol; c.solver.startDecay=11; return c;
}
Source source() { return {{geo(10),20},50000,{1e19,2e18,-3e18,-5e18,4e18,-5e18},{}}; }
Receiver receiver() { return {{geo(45),90},3000,{}}; }
Direction direction(int o) { return o==0?Direction{.up=1}:Direction{.north=1}; }
Measurement measurement(const Receiver& r,int obs) {
  std::vector<double> j(nSamples);
  for(int n=0;n<nSamples;++n) {
    const double w=0.5-0.5*std::cos(2*std::numbers::pi*n/(nSamples-1));
    const auto term=[&](int k,double re,double im) { const double a=2*std::numbers::pi*k*n/128.0; return re*std::cos(a)-im*std::sin(a); };
    j[n]=w*(term(1,q1r,q1i)+term(2,q2r,q2i));
  }
  return {r,direction(obs),std::move(j)};
}
double objective(const Model& m,const Calculation& c,const Source& s,const Receiver& r,const Measurement& j) {
  const auto u=ComputeSeismograms(m,c,std::span{&s,1},std::span{&r,1});
  double value=0; const auto d=j.direction;
  for(int n=0;n<nSamples;++n) value+=j.adjointSource[n]*(d.up*u(0,0,Vertical,n)+d.north*u(0,0,North,n)+d.east*u(0,0,East,n))*dt;
  return value;
}
int main(int argc,char**argv) try {
  if(argc<6 || argc>7) throw std::runtime_error("usage: kernel_adapter MODEL OUTPUT TOLERANCE LEVEL SUPPORT_LAYER [source-local-other|source-local-za]");
  const fs::path modelPath=fs::absolute(argv[1]),out=fs::absolute(argv[2]);
  const auto supportLayer=static_cast<std::size_t>(std::stoul(argv[5]));
  const std::string selector=argc==7?argv[6]:"all";
  if(selector!="all" && selector!="source-local-other" && selector!="source-local-za")
    throw std::runtime_error("unknown direction selector: "+selector);
  if(fs::exists(out)&&(!fs::is_directory(out)||fs::directory_iterator(out)!=fs::directory_iterator{})) throw std::runtime_error("output directory must be empty");
  fs::create_directories(out);
  const auto model=Model::ReadDeck(modelPath); const auto c=calc(std::stod(argv[3])); const auto s=source(); const auto r=receiver();
  const auto grid=FrequencyGrid(c.band.min*model.GetScales().Time(),c.band.max*model.GetScales().Time(),dt/model.GetScales().Time(),duration/model.GetScales().Time());
  if(grid.Samples()!=128||grid.First()>firstBin||grid.Last()<endBin-1) throw std::runtime_error("unexpected YSpec frequency grid");
  const auto spectra=SeismogramSpectra(model,c,std::span{&s,1},std::span{&r,1});
  const auto shellIndex=model.LayerIndex(3000e3/model.GetScales().Length());
  const auto shell=model(shellIndex,3000e3/model.GetScales().Length());
  const auto& support=model.GetLayer(supportLayer);
  const double supportLower=support.Lower()*model.GetScales().Length();
  const double supportUpper=support.Upper()*model.GetScales().Length();
  const double supportProbe=(support.Lower()+support.Upper())/2.0;
  const auto supportMaterial=model(supportLayer,supportProbe);
  const std::array<double,5> shellPa{supportMaterial.A,supportMaterial.C,supportMaterial.F,supportMaterial.L,supportMaterial.N};
  const std::array<double,3> sourceGeocentric{
      geocentric(s.location.latitude),s.location.longitude,s.depth};
  const std::array<double,3> receiverGeocentric{
      geocentric(r.location.latitude),r.location.longitude,r.depth};
  const std::array<double,6> tensor{s.moment.rr,s.moment.rt,s.moment.rp,
                                    s.moment.tt,s.moment.tp,s.moment.pp};
  const double densitySI=shell.rho*model.GetScales().Density();
  const double vpSI=std::sqrt(shell.C/shell.rho)*model.GetScales().Velocity();
  const double vsSI=std::sqrt(shell.L/shell.rho)*model.GetScales().Velocity();
  const double sourceRadiusM=model.Layers().back().Knots().back()*model.GetScales().Length()-s.depth;
  const double sourceRadius=sourceRadiusM/model.GetScales().Length();
  const double receiverRadiusM=model.Layers().back().Knots().back()*model.GetScales().Length()-r.depth;
  const double receiverRadius=receiverRadiusM/model.GetScales().Length();
  const auto sourceLayer=model.LayerIndex(sourceRadius,Side::Above);
  const auto receiverLayer=model.LayerIndex(receiverRadius,Side::Above);
  const bool sourceLocal=(sourceLayer==supportLayer && receiverLayer==supportLayer);
  const std::vector<double> steps=sourceLocal?std::vector<double>(extendedSteps.begin(),extendedSteps.end()):std::vector<double>{0.01,0.005,0.0025};
  const double kernelMaxElement=sourceLocal?1000.0:250e3;
  std::ofstream meta(out/"meta.json");
  meta<<std::setprecision(17)<<"{\"solver\":\"YSpecpp\",\"nt\":"<<grid.Samples()
      <<",\"dt_s\":"<<grid.TimeStep()*model.GetScales().Time()
      <<",\"df_hz\":"<<grid.Spacing()/model.GetScales().Time()
      <<",\"damping_s\":"<<grid.Damping()/model.GetScales().Time()
      <<",\"pi\":"<<std::numbers::pi
      <<",\"omega_bin1\":"<<grid.Omega(1-grid.First()).real()/model.GetScales().Time()
      <<",\"tolerance\":"<<c.solver.tolerance<<",\"physical_G\":"<<GravitationalConstant
      <<",\"density_kg_m3\":"<<densitySI<<",\"vp_m_s\":"<<vpSI<<",\"vs_m_s\":"<<vsSI
      <<",\"radius_m\":"<<model.Layers().back().Knots().back()*model.GetScales().Length()
      <<",\"requested_band_hz\":["<<c.band.min<<','<<c.band.max<<']'
      <<",\"shell_A_C_F_L_N_Pa\":[";
  for(std::size_t i=0;i<shellPa.size();++i)meta<<shellPa[i]*model.GetScales().Modulus()<<(i+1<shellPa.size()?",":"");
  meta<<"],\"gravity\":\"full\",\"attenuation\":false,\"degrees\":[0,2]"
      <<",\"q_bins\":[1,2],\"q_values\":["<<q1r<<','<<q1i<<','<<q2r<<','<<q2i
      <<"],\"window_samples\":"<<nSamples<<",\"window\":\"Hann\""
      <<",\"source_geocentric\":["<<sourceGeocentric[0]<<','<<sourceGeocentric[1]<<','<<sourceGeocentric[2]<<']'
      <<",\"receiver_geocentric\":["<<receiverGeocentric[0]<<','<<receiverGeocentric[1]<<','<<receiverGeocentric[2]<<']'
      <<",\"source_radius_m\":"<<model.Layers().back().Knots().back()*model.GetScales().Length()-s.depth
      <<",\"source_layer\":"<<sourceLayer
      <<",\"receiver_radius_m\":"<<receiverRadiusM<<",\"receiver_layer\":"<<receiverLayer
      <<",\"support_layer\":"<<supportLayer
      <<",\"direction_selector\":\""<<selector<<"\""
      <<",\"support_bounds_m\":["<<supportLower<<','<<supportUpper<<']'
      <<",\"source_local_support_case\":"<<(sourceLocal?"true":"false")
      <<",\"kernel_mesh_order\":5,\"kernel_mesh_max_element_m\":"<<kernelMaxElement
      <<",\"kernel_source_receiver_same_layer\":"<<(sourceLayer==receiverLayer?"true":"false")
      <<",\"source_geographic\":["<<s.location.latitude<<','<<s.location.longitude<<','<<s.depth<<']'
      <<",\"receiver_geographic\":["<<r.location.latitude<<','<<r.location.longitude<<','<<r.depth<<']'
      <<",\"moment_tensor_Nm\":["<<tensor[0]<<','<<tensor[1]<<','<<tensor[2]<<','<<tensor[3]<<','<<tensor[4]<<','<<tensor[5]<<"]}\n";
  std::ofstream spec(out/"spectra.csv"),curves(out/"samples.csv"),dirs(out/"derivatives.csv");
  spec<<"component,bin,frequency_hz,omega_real,omega_imag,real_m_s,imag_m_s\n"<<std::setprecision(17);
  for(int k=firstBin;k<endBin;++k) for(int comp=0;comp<3;++comp) {
    const auto x=spectra(0,0,comp,k-grid.First())*model.GetScales().Length()*model.GetScales().Time();
    spec<<(comp==0?"Z":comp==1?"N":"E")<<','<<k<<','<<grid.Frequency(k-grid.First())/model.GetScales().Time()<<','<<grid.Omega(k-grid.First()).real()/model.GetScales().Time()<<','<<grid.Omega(k-grid.First()).imag()/model.GetScales().Time()<<','<<x.real()<<','<<x.imag()<<'\n';
  }
  curves<<"observable,modulus,radius_m,kernel_dr_per_Pa_m,weight_m3,delta_modulus_Pa\n"<<std::setprecision(17);
  dirs<<"observable,modulus,objective,kernel,volume_kernel,positive_volume_kernel,negative_volume_kernel,fd_h_0_08,fd_h_0_04,fd_h_0_02,fd_h_0_01,fd_h_0_005,fd_h_0_0025,r4_0_08_0_04,r4_0_04_0_02,r4_0_02_0_01,r4_0_01_0_005,r4_0_005_0_0025,r6_0_08_0_04_0_02,r6_0_04_0_02_0_01,r6_0_02_0_01_0_005,r6_0_01_0_005_0_0025,richardson,relative_kernel_fd,source_kernel,source_delta,source_contribution,total_kernel\n"<<std::setprecision(17);
  const RadialMesh mesh{.order=5,.maxElement=kernelMaxElement};
  const auto sourceMaterial=model(sourceLayer,sourceRadius);
  for(int obs=0;obs<2;++obs) {
    const auto j=measurement(r,obs); const auto val=objective(model,c,s,r,j);
    const auto k=ComputeRadialKernels(model,c,s,j,mesh);
    for(std::size_t mi=0;mi<mods.size();++mi) {
      if((selector=="source-local-other" && obs==0 && mi==0) ||
         (selector=="source-local-za" && (obs!=0 || mi!=0))) continue;
      RadialParameter p=static_cast<RadialParameter>(static_cast<int>(RadialParameter::A)+static_cast<int>(mi));
      double predicted=0,positiveVolume=0,negativeVolume=0;
      for(std::size_t i=0;i<k.Size();++i) {
        const double rr=k.Radii()[i], xrad=rr/model.GetScales().Length();
        const auto layer=model.LayerIndex(xrad,Side::Above); const auto material=model(layer,xrad);
        const double v=mi==0?material.A:mi==1?material.C:mi==2?material.F:mi==3?material.L:material.N;
        const bool active=layer==supportLayer;
        if(active) {
          const double contribution=k[p][i]*k.Weights()[i]*v*model.GetScales().Modulus();
          predicted+=contribution;
          if(contribution>0) positiveVolume+=contribution;
          else negativeVolume+=contribution;
        }
        curves<<(obs==0?"Z":"N")<<','<<names[mi]<<','<<rr<<','<<k[p][i]*rr*rr<<','<<k.Weights()[i]<<','<<(active?v*model.GetScales().Modulus():0.0)<<'\n';
      }
      const double sourceValue=mi==0?sourceMaterial.A:mi==1?sourceMaterial.C:
                               mi==2?sourceMaterial.F:mi==3?sourceMaterial.L:sourceMaterial.N;
      const double sourceDelta=sourceLayer==supportLayer?sourceValue*model.GetScales().Modulus():0.0;
      const double sourceKernel=k.AtSource(p);
      const double sourceContribution=sourceKernel*sourceDelta;
      std::vector<double> fd(steps.size());
      for(std::size_t h=0;h<steps.size();++h) {
        const auto plus=alter(model,mods[mi],supportLayer,steps[h]),minus=alter(model,mods[mi],supportLayer,-steps[h]);
        fd[h]=(objective(plus,c,s,r,j)-objective(minus,c,s,r,j))/(2*steps[h]);
      }
      std::vector<double> r4, r6;
      for(std::size_t h=0;h+1<fd.size();++h) r4.push_back((4*fd[h+1]-fd[h])/3);
      for(std::size_t h=0;h+1<r4.size();++h) r6.push_back((16*r4[h+1]-r4[h])/15);
      const double rich=sourceLocal?r4.back():(4*fd[2]-fd[1])/3;
      const double total=predicted+sourceContribution;
      const double rel=std::abs(total-rich)/std::max({std::abs(total),std::abs(rich),1e-300});
      dirs<<(obs==0?"Z":"N")<<','<<names[mi]<<','<<val<<','<<total<<','<<predicted<<','<<positiveVolume<<','<<negativeVolume;
      if(sourceLocal) { for(const double v:fd) dirs<<','<<v; }
      else { dirs<<",,,,"<<fd[0]<<','<<fd[1]<<','<<fd[2]; }
      if(sourceLocal) { for(const double v:r4) dirs<<','<<v; for(const double v:r6) dirs<<','<<v; }
      else dirs<<",,,,,,,,,";
      dirs<<','<<rich<<','<<rel<<','<<sourceKernel<<','<<sourceDelta<<','<<sourceContribution<<','<<total<<'\n';
    }
    (void)val;
  }
  return 0;
} catch(const std::exception&e) { std::cerr<<"YSpec kernel adapter: "<<e.what()<<'\n'; return 1; }
