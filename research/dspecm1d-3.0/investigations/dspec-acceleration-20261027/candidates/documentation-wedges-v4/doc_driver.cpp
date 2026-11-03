// Matched documentation-case adapter: common model/measurement, independent solvers.
#include <YSpec/YSpec.hpp>
#include "Spectra.hpp"
#include "KernelParts.hpp"
#include "Sphere.hpp"
#include "matched_mesh_model.hpp"
#include "engine_headers/Windowed3DKernels.h"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <sstream>
#include <omp.h>
namespace fs=std::filesystem;
using Complex=std::complex<double>;
using Clock=std::chrono::steady_clock;
constexpr double pi=std::numbers::pi;
double seconds(Clock::time_point t){return std::chrono::duration<double>(Clock::now()-t).count();}
double profile(double r){return r>5.3e6&&r<5.8e6?std::pow(std::sin(pi*(r-5.3e6)/.5e6),2):0.;}
struct CommonModel {
 const YSpec::Model &model;double epsilon=0;
 auto values(int l,double r) const {
  auto m=model(l,r);double d=m.mu*std::expm1(2*epsilon*profile(r*LengthNorm()));
  m.F-=2*d;m.L+=d;m.N+=d;m.kappa-=4*d/3;m.mu+=d;return m;
 }
 double LengthNorm() const{return model.GetScales().Length();}
 double TimeNorm() const{return model.GetScales().Time();}
 double DensityNorm() const{return model.GetScales().Density();}
 double MassNorm() const{return DensityNorm()*std::pow(LengthNorm(),3);}
 double GravitationalConstant() const{return 1/(DensityNorm()*TimeNorm()*TimeNorm());}
 int NumberOfLayers() const{return model.NumberOfLayers();}
 double OuterRadius() const{return model.SurfaceRadius();}
 double LowerRadius(int l) const{return model.GetLayer(l).Lower();}
 double UpperRadius(int l) const{return model.GetLayer(l).Upper();}
 bool IsFluid(int l) const{return model.GetLayer(l).IsFluid();}
 bool IsSolid(int l) const{return !IsFluid(l);}
 auto LayerRadii(int l) const {auto r=model.GetLayer(l).Knots();return std::vector<double>(r.begin(),r.end());}
 #define PROPERTY(name,expression) auto name(int l) const{return [this,l](double r){auto m=values(l,r);return (expression);};}
 PROPERTY(Density,m.rho) PROPERTY(A,m.A) PROPERTY(C,m.C) PROPERTY(F,m.F) PROPERTY(L,m.L) PROPERTY(N,m.N)
 PROPERTY(Kappa,m.kappa) PROPERTY(Mu,m.mu)
 PROPERTY(VPV,std::sqrt(m.C/m.rho)) PROPERTY(VPH,std::sqrt(m.A/m.rho))
 PROPERTY(VSV,std::sqrt(std::max(0.,m.L/m.rho))) PROPERTY(VSH,std::sqrt(std::max(0.,m.N/m.rho)))
 PROPERTY(VP,std::sqrt((m.kappa+4*m.mu/3)/m.rho)) PROPERTY(VS,std::sqrt(std::max(0.,m.mu/m.rho)))
 PROPERTY(QKappa,m.qKappa>0?1/m.qKappa:0.) PROPERTY(QMu,m.qMu>0?1/m.qMu:0.)
 PROPERTY(Eta,m.F/(m.A-2*m.L)) PROPERTY(AAtten,m.A*m.qA) PROPERTY(CAtten,m.C*m.qA)
 PROPERTY(FAtten,m.F*m.qF) PROPERTY(LAtten,m.L*m.qMu) PROPERTY(NAtten,m.N*m.qMu)
 #undef PROPERTY
};
struct Node {int element,node,layer;double radius,dr;bool solid;};
std::string params(const YSpec::KernelInput &in,int degree){
 const auto &s=in.source;const auto &r=in.receiver;const auto &m=s.moment;
 auto lat=[](double x){return 90-YSpec::Sphere::Colatitude(x)*180/pi;};
 std::ostringstream o;o<<std::setprecision(17)<<"\"unused\"\n\""<<in.model.string()<<"\"\n4\n"<<(in.calculation.attenuation?1:0)<<"\n2\n0\n0\n1e-5\n0\n"<<degree<<"\n2\n10\n40\n2\n2\n4\n8\n10\n"<<s.depth/1000<<'\n'<<lat(s.location.latitude)<<'\n'<<s.location.longitude<<'\n'<<m.rr<<'\n'<<m.rt<<'\n'<<m.rp<<'\n'<<m.tt<<'\n'<<m.tp<<'\n'<<m.pp<<'\n'<<r.depth/1000<<"\n1\n"<<lat(r.location.latitude)<<' '<<r.location.longitude<<'\n';return o.str();
}
struct Measured {YSpec::Measurement measurement;std::vector<double> velocity,rate,window;double norm=0,delayAudit=0;};
Measured measure(const YSpec::KernelInput &in,const YSpec::Model &model,const YSpec::FrequencyGrid &grid,const std::vector<Complex> &raw){
 Measured out;auto rate=raw;for(std::size_t k=0;k<rate.size();++k)rate[k]*=Complex(0,1)*grid.Omega(k);
 out.velocity=YSpec::InverseTransform(grid,YSpec::MakeTaper(in.calculation,model.GetScales()),raw,1);
 out.rate=YSpec::InverseTransform(grid,YSpec::MakeTaper(in.calculation,model.GetScales()),rate,1);
 const auto &s=model.GetScales();double dt=in.calculation.sampling.timeStep;
 std::size_t samples=YSpec::SampleCount(grid,in.calculation.sampling);
 out.velocity.resize(samples);out.rate.resize(samples);out.measurement={in.receiver,in.measurements[0].direction,{}};
 for(std::size_t n=0;n<samples;++n){out.velocity[n]*=s.Velocity();out.rate[n]*=s.Velocity()/s.Time();double w=in.measurements[0].window(n*dt);out.window.push_back(w);out.norm+=w*out.rate[n]*out.rate[n]*dt;out.measurement.adjointSource.push_back(-w*out.rate[n]);}
 if(!(out.norm>0))throw std::runtime_error("zero traveltime normalisation");
 for(auto &j:out.measurement.adjointSource)j/=out.norm;
 auto adjoint=YSpec::MakeAdjointSpectrum(model,in.calculation,out.measurement);
 for(std::size_t k=0;k<raw.size();++k){Complex z=grid.Omega(k),q=std::conj(adjoint.weight[k]*adjoint.force[k]);out.delayAudit+=std::real(std::conj(q)*Complex(0,-1)/z*(-Complex(0,1)*z*raw[k]/s.Time()));}
 if(std::abs(out.delayAudit-1)>1e-10)throw std::runtime_error("discrete time-shift adjoint identity failed");
 return out;
}
// Complex weak-form contractions, with fixed effective Q, before real projection.
double kernelVs(const Full1D::SEM &sem,const SPARSESPEC::KernelField3DSample &sample,const YSpec::Material &m,Complex dispersion){
 if(sample.fluid)return 0;
 const auto &f=sample.forwardRawStrain;const auto &a=sample.adjointDualStrain;
 Complex F=-(std::conj(a.rr)*f.horizontalTrace+std::conj(a.horizontalTrace)*f.rr);
 Complex L=-2.*(std::conj(a.radialPlus)*f.radialPlus+std::conj(a.radialMinus)*f.radialMinus);
 Complex N=-(.5*(std::conj(a.horizontalPlus)*f.horizontalPlus+std::conj(a.horizontalMinus)*f.horizontalMinus)-std::conj(a.horizontalTrace)*f.horizontalTrace);
 double length=sem.lengthScaleMetres();
 return 2*m.mu*std::real(-2.*F*(1.+m.qF*dispersion)+(L+N)*(1.+m.qMu*dispersion))/std::pow(length,3);
}
int main(int argc,char **argv)try{
 if(argc!=10)throw std::runtime_error("usage: doc-driver dspec|yspec|dspec-forward INPUT EMPTY_OUTPUT L NQ REFINEMENT WORKERS sections|audit EPSILON");
 auto start=Clock::now();std::string backend=argv[1],gridKind=argv[8];fs::path out=fs::absolute(argv[3]);
 int degree=std::stoi(argv[4]),nq=std::stoi(argv[5]),workers=std::stoi(argv[7]);double refinement=std::stod(argv[6]),epsilon=std::stod(argv[9]);
 if(workers<1||workers>8||degree<1||degree>200||nq<6||nq>16||!(refinement>0)||std::abs(epsilon)>.001)throw std::runtime_error("bounded documentation scope exceeded");
 fs::create_directories(out);if(fs::directory_iterator(out)!=fs::directory_iterator{})throw std::runtime_error("output must be empty");
 Eigen::setNbThreads(1);omp_set_dynamic(0);omp_set_num_threads(1);
 mathprofile::Root root("documentation-wedge",-1);
 auto in=YSpec::ReadKernelInput(argv[2]);in.calculation.degrees.max=degree;
 if(in.calculation.quantity!=YSpec::Quantity::Velocity||in.measurements.size()!=1||in.measurements[0].component!="T")throw std::runtime_error("expected documentation transverse-velocity case");
 auto model=YSpec::Model::ReadDeck(in.model);CommonModel common{model,epsilon};
 const auto &scales=model.GetScales();auto grid=YSpec::MakeGrid(in.calculation,scales);
 auto geoLat=[](double lat){return 90-YSpec::Sphere::Colatitude(lat)*180/pi;};const auto &moment=in.source.moment;
 SourceInfo::EarthquakeCMT source(geoLat(in.source.location.latitude),in.source.location.longitude,in.source.depth/1000,moment.rr,moment.rt,moment.rp,moment.tt,moment.tp,moment.pp);
 Full1D::SEM sem(common,refinement,nq,degree,source);
 std::ofstream pf(out/"parameters.txt");pf<<params(in,degree);pf.close();InputParameters parameters((out/"parameters.txt").string());
 std::vector<Node> nodes;auto gll=GaussQuad::GaussLobattoLegendreQuadrature1D<double>(nq);double gravityError=0;
 for(int e=0;e<sem.mesh().NE();++e)for(int q=0;q<nq;++q){double r=sem.mesh().NodeRadius(e,q);int layer=sem.mesh().LayerNumber(e);gravityError=std::max(gravityError,std::abs(sem.meshModel().Gravity(e,q)-model(layer,r).g));if(r*scales.Length()>=3e6)nodes.push_back({e,q,layer,r*scales.Length(),.5*sem.mesh().EW(e)*gll.W(q)*scales.Length(),!sem.elementIsFluid(e)});}
 if(gravityError>1e-10)throw std::runtime_error("common-model gravity mismatch");
 std::vector<YSpec::KernelSampling> samplings(2);std::vector<double> angularWeights;
 const double spacing=360./(4*degree+1);
 if(gridKind=="sections"){
  YSpec::SectionOptions options{0,4,spacing,500e3};samplings[0]=YSpec::PathSection(model,in.source.location,in.receiver.location,options);samplings[1]=YSpec::CrossSection(model,in.source.location,in.receiver.location,options);
 }else if(gridKind=="audit"){
  auto quad=GaussQuad::Quadrature1D<double>(GaussQuad::LegendrePolynomial<double>{}.GaussQuadrature(2*degree+1));int np=4*degree+1;
  for(int t=0;t<quad.N();++t)for(int p=0;p<np;++p){samplings[0].places.push_back({std::acos(quad.X(t)),2*pi*p/np,0,0});angularWeights.push_back(quad.W(t)*2*pi/np);}samplings.resize(1);
 }else throw std::runtime_error("unknown sampling kind");
 for(auto &s:samplings){s.nodes.clear();for(auto n:nodes)s.nodes.push_back({n.layer,n.radius});}
 const auto dir=in.measurements[0].direction;Eigen::Vector3cd direction;direction<<dir.up,-dir.north,dir.east;
 int receiverElement=sem.receiverElements(parameters,true).front();double receiverRadius=scales.Length()*model.SurfaceRadius()-in.receiver.depth;
 std::vector<SPARSESPEC::Field3DPoint> receiver{{receiverRadius,YSpec::Sphere::Colatitude(in.receiver.location.latitude),in.receiver.location.longitude*pi/180,receiverElement}};
 std::vector<Complex> raw(grid.Count());double forwardSeconds=0,kernelSeconds=0,measurementSeconds=0;
 dssolver::BatchResponseCache cache({static_cast<std::size_t>(2*degree)*grid.Count()+1,16ull<<30,1e-12});
 auto dispersion=[&](Complex z){return in.calculation.attenuation?2/pi*std::log(Complex(0,1)*z/model.ReferenceFrequency()):Complex{};};
 auto timer=Clock::now();
 if(backend=="yspec"){
  if(epsilon!=0)throw std::runtime_error("YSpec perturbation not requested");omp_set_num_threads(workers);
  auto spectra=YSpec::SeismogramSpectra(model,in.calculation,std::span{&in.source,1},std::span{&in.receiver,1});
  for(std::size_t k=0;k<grid.Count();++k)raw[k]=dir.up*spectra(0,0,YSpec::Vertical,k)+dir.north*spectra(0,0,YSpec::North,k)+dir.east*spectra(0,0,YSpec::East,k);
 }else if(backend=="dspec"||backend=="dspec-forward"){
  for(std::size_t k=0;k<grid.Count();++k){auto z=grid.Omega(k);auto workspace=SPARSESPEC::Windowed3DFieldEngine::prepareFrequencyWorkspace(sem,source,parameters,z,Eigen::Vector3cd::Zero(),&cache,true,false,nullptr,workers,dispersion(z));auto field=workspace.sampleForKernelsPacked(receiver,0,nullptr,2).front();for(int c=0;c<3;++c)raw[k]+=direction[c]*field.forwardRawDisplacement[c];}
 }else throw std::runtime_error("unknown backend");
 forwardSeconds=seconds(timer);timer=Clock::now();auto measured=measure(in,model,grid,raw);measurementSeconds=seconds(timer);
 {std::ofstream f(out/"trace.csv");f<<"time_s,velocity_m_s,rate_m_s2,window,adjoint_source\n"<<std::setprecision(17);for(std::size_t n=0;n<measured.velocity.size();++n)f<<n*in.calculation.sampling.timeStep<<','<<measured.velocity[n]<<','<<measured.rate[n]<<','<<measured.window[n]<<','<<measured.measurement.adjointSource[n]<<'\n';}
 {std::ofstream f(out/"spectrum.csv");f<<"bin,frequency_hz,real_m,imag_m\n"<<std::setprecision(17);for(std::size_t k=0;k<raw.size();++k)f<<grid.First()+k<<','<<grid.Frequency(k)/scales.Time()<<','<<raw[k].real()*scales.Length()<<','<<raw[k].imag()*scales.Length()<<'\n';}
 {std::ofstream f(out/"nodes.csv");f<<"element,node,layer,radius_m,dr_m,solid\n"<<std::setprecision(17);for(auto n:nodes)f<<n.element<<','<<n.node<<','<<n.layer<<','<<n.radius<<','<<n.dr<<','<<n.solid<<'\n';}
 std::vector<std::vector<double>> values;for(auto &s:samplings)values.emplace_back(s.nodes.size()*s.places.size());
 timer=Clock::now();
 if(backend=="yspec"){
  auto kernels=YSpec::ComputeSampledKernels(model,in.calculation,in.source,std::span{&measured.measurement,1},YSpec::ParameterSet::Velocities,YSpec::Scaling::Relative,samplings,[](std::size_t n,std::size_t total){std::cout<<"YSpec frequencies "<<n<<'/'<<total<<std::endl;});
  for(std::size_t s=0;s<samplings.size();++s){auto &k=kernels[0][s];auto names=k.Names();auto it=std::find(names.begin(),names.end(),"vs");if(it==names.end())throw std::runtime_error("missing vs");auto parameter=it-names.begin();for(std::size_t n=0;n<nodes.size();++n){auto v=k.Values(parameter,n);std::copy(v.begin(),v.end(),values[s].begin()+n*v.size());}}
 }else if(backend=="dspec"){
  auto adjoint=YSpec::MakeAdjointSpectrum(model,in.calculation,measured.measurement);
  std::vector<SPARSESPEC::Windowed3DFieldEngine::SamplingPlanCache> plans(samplings.size());
  for(std::size_t k=0;k<grid.Count();++k){auto z=grid.Omega(k),d=dispersion(z);auto weights=direction*std::conj(adjoint.weight[k]*adjoint.force[k]);auto workspace=SPARSESPEC::Windowed3DFieldEngine::prepareFrequencyWorkspace(sem,source,parameters,z,weights,&cache,true,false,nullptr,workers,d);
   for(std::size_t s=0;s<samplings.size();++s){const auto &sampling=samplings[s];std::vector<SPARSESPEC::Field3DPoint> points;for(auto p:sampling.places)for(auto n:nodes)points.push_back({n.radius,p.colatitude,p.longitude,n.element});auto fields=workspace.sampleForKernelsTensor(points,points.size(),&plans[s],2,{},false);for(std::size_t a=0;a<sampling.places.size();++a)for(std::size_t n=0;n<nodes.size();++n)values[s][n*sampling.places.size()+a]+=kernelVs(sem,fields[a*nodes.size()+n],model(nodes[n].layer,nodes[n].radius/scales.Length()),d);}
   std::cout<<"DSpec frequencies "<<k+1<<'/'<<grid.Count()<<std::endl;
  }
 }
 kernelSeconds=seconds(timer);
 for(std::size_t s=0;s<samplings.size();++s){std::string name=s==0?"along":"across";std::ofstream f(out/(name+".bin"),std::ios::binary);f.write(reinterpret_cast<char*>(values[s].data()),values[s].size()*sizeof(double));std::ofstream a(out/(name+"-angles.csv"));a<<"colatitude_rad,longitude_rad,along_deg,across_deg,angular_weight_sr\n"<<std::setprecision(17);for(std::size_t p=0;p<samplings[s].places.size();++p){auto x=samplings[s].places[p];a<<x.colatitude<<','<<x.longitude<<','<<x.along<<','<<x.across<<','<<(gridKind=="audit"?angularWeights[p]:0)<<'\n';}}
 double integrated=0;if(gridKind=="audit")for(std::size_t n=0;n<nodes.size();++n)for(std::size_t a=0;a<angularWeights.size();++a)integrated+=values[0][n*angularWeights.size()+a]*profile(nodes[n].radius)*nodes[n].radius*nodes[n].radius*nodes[n].dr*angularWeights[a];
 {std::ofstream f(out/"result.json");f<<std::setprecision(17)<<"{\"backend\":\""<<backend<<"\",\"degree\":"<<degree<<",\"nq\":"<<nq<<",\"refinement\":"<<refinement<<",\"workers\":"<<workers<<",\"attenuation\":"<<(in.calculation.attenuation?"true":"false")<<",\"epsilon_perturbation\":"<<epsilon<<",\"frequency_count\":"<<grid.Count()<<",\"first_bin\":"<<grid.First()<<",\"last_bin\":"<<grid.Last()<<",\"fft_samples\":"<<grid.Samples()<<",\"fft_period_s\":"<<grid.Samples()*in.calculation.sampling.timeStep<<",\"damping_si\":"<<grid.Damping()/scales.Time()<<",\"native_side_count\":"<<nodes.size()<<",\"along_angles\":"<<samplings[0].places.size()<<",\"across_angles\":"<<(samplings.size()>1?samplings[1].places.size():0)<<",\"gravity_max_abs_nd_error\":"<<gravityError<<",\"time_shift_adjoint_identity\":"<<measured.delayAudit<<",\"measurement_rate_norm\":"<<measured.norm<<",\"integrated_profile_kernel\":"<<integrated<<",\"forward_s\":"<<forwardSeconds<<",\"measurement_s\":"<<measurementSeconds<<",\"kernel_s\":"<<kernelSeconds<<",\"cache_bytes\":"<<cache.retainedBytes()<<",\"elapsed_s\":"<<seconds(start)<<",\"kernel_parameter\":\"ln_vs_at_fixed_vp_density_and_effective_Q\",\"kernel_units\":\"s m^-3\",\"layout\":\"native-side-major then angle, float64\",\"source_shell_contacts\":false}\n";}
 return 0;
}catch(const std::exception &e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}
