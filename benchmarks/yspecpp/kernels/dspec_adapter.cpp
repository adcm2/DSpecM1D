#include <DSpecM1D/All>
#include <DSpecM1D/SensitivityKernels>
#include <DSpecM1D/ModelInput>
#include <DSpecM1D/src/InputParser.h>
#include <DSpecM1D/src/NormClass.h>
#include <DSpecM1D/src/SourceInfo.h>

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>

namespace fs = std::filesystem;
enum class Modulus { A, C, F, L, N };
template<class Base, class Norm> class Perturbed final : public Base {
  using Field = std::function<double(double)>;
  template<class F> auto change(F f, int layer, Modulus which, Modulus selected,
                                int selectedLayer, double amount) const {
    const bool active = layer == selectedLayer && which == selected;
    return [f=std::move(f), active, amount](double r) {
      return f(r) * (1.0 + (active ? amount : 0.0));
    };
  }
public:
  Perturbed(const std::string& path, const Norm& norm, Modulus selected,
            int selectedLayer, double amount)
    : Base(path, norm), selected_(selected), selectedLayer_(selectedLayer), amount_(amount) {}
  auto A(int l) const { return change(Base::A(l),l,Modulus::A,selected_,selectedLayer_,amount_); }
  auto C(int l) const { return change(Base::C(l),l,Modulus::C,selected_,selectedLayer_,amount_); }
  auto F(int l) const { return change(Base::F(l),l,Modulus::F,selected_,selectedLayer_,amount_); }
  auto L(int l) const { return change(Base::L(l),l,Modulus::L,selected_,selectedLayer_,amount_); }
  auto N(int l) const { return change(Base::N(l),l,Modulus::N,selected_,selectedLayer_,amount_); }
private: Modulus selected_; int selectedLayer_; double amount_;
};

constexpr std::array<const char*,5> names{"A","C","F","L","N"};
constexpr std::array<Modulus,5> mods{Modulus::A,Modulus::C,Modulus::F,Modulus::L,Modulus::N};
constexpr std::array<double,6> extendedSteps{0.08,0.04,0.02,0.01,0.005,0.0025};
constexpr int firstBin=1, endBin=4;
constexpr double radius=6371000.0, interfaceRadius=3200000.0;
std::string paramsText(const fs::path& model, const fs::path& out) {
  return "\""+(out/"unused").string()+"\"\n\""+model.string()+"\"\n4\n0\n2\n0\n0\n1e-5\n0\n2\n0.1\n0.8\n30\n20\n0\n0\n0\n0\n50\n10\n20\n1e19\n2e18\n-3e18\n-5e18\n4e18\n-5e18\n3\n1\n45 90\n";
}
Eigen::VectorXd window(const SpectraSolver::FreqFull& f) {
  Eigen::VectorXd w=Eigen::VectorXd::Zero(f.nt());
  constexpr int n=91;
  for(int i=0;i<n;++i) w[i]=0.5-0.5*std::cos(2.0*3.14159265358979323846*i/(n-1));
  return w;
}
Eigen::MatrixXcd qweights(int component,int bins) {
  Eigen::MatrixXcd q=Eigen::MatrixXcd::Zero(3,bins);
  q(component,1)={0.7,0.2}; q(component,2)={-0.25,0.45};
  return q;
}
int main(int argc,char**argv) try {
  if(argc<6 || argc>7) throw std::runtime_error("usage: kernel_adapter MODEL OUTPUT MAXSTEP NQ SUPPORT_LAYER [source-local-other|source-local-za]");
  const fs::path modelPath=fs::absolute(argv[1]), out=fs::absolute(argv[2]);
  const int supportLayer=std::stoi(argv[5]);
  const std::string selector=argc==7?argv[6]:"all";
  if(selector!="all" && selector!="source-local-other" && selector!="source-local-za")
    throw std::runtime_error("unknown direction selector: "+selector);
  if(fs::exists(out) && (!fs::is_directory(out)||fs::directory_iterator(out)!=fs::directory_iterator{}))
    throw std::runtime_error("output directory must be empty: "+out.string());
  fs::create_directories(out);
  prem_norm<double> norm;
  using Base=EarthModels::ModelInput<double>;
  Base base(modelPath.string(),norm);
  auto pfile=out/"params.txt"; {std::ofstream p(pfile); p<<paramsText(modelPath,out);}
  InputParametersNew p(pfile.string(),std::stoi(argv[4]),1,std::stod(argv[3]),1.0,0.1,0.0,1);
  const auto&freq=p.freqFull(); const auto win=window(freq);
  if(freq.nt()!=128 || freq.i1()>firstBin || freq.i2()<endBin) throw std::runtime_error("unexpected DSpec frequency grid");
  auto& parameters=p.inputParameters();
  SourceInfo::EarthquakeCMT source(parameters);
  Full1D::SEM sem(base,std::stod(argv[3]),std::stoi(argv[4]),2,source);
  const double supportProbe=(base.LowerRadius(supportLayer)+base.UpperRadius(supportLayer))/2.0;
  const double supportLowerM=base.LowerRadius(supportLayer)*sem.lengthScaleMetres();
  const double supportUpperM=base.UpperRadius(supportLayer)*sem.lengthScaleMetres();
  std::array<double,5> shellPa{base.A(supportLayer)(supportProbe),base.C(supportLayer)(supportProbe),
      base.F(supportLayer)(supportProbe),base.L(supportLayer)(supportProbe),base.N(supportLayer)(supportProbe)};
  for(auto&v:shellPa)v*=base.StressNorm();
  const int sourceElem=sem.sourceElement(source);
  const auto receiverElems=sem.receiverElements(parameters);
  const bool sourceLocal=(sem.mesh().LayerNumber(sourceElem)==supportLayer &&
                          sem.mesh().LayerNumber(receiverElems.front())==supportLayer);
  const std::vector<double> steps=sourceLocal?std::vector<double>(extendedSteps.begin(),extendedSteps.end()):std::vector<double>{0.01,0.005,0.0025};
  const double gPhysical=sem.gravityCoefficient()*base.GravitationalConstant();
  const auto receiver=parameters.receivers().front();
  const std::array<double,3> sourceLocation{source.Latitude(),source.Longitude(),source.Depth()*1000.0};
  const std::array<double,3> receiverLocation{receiver.first,receiver.second,parameters.receiver_depth()*1000.0};
  const std::array<double,6> tensor{source.Mrr(),source.Mrt(),source.Mrp(),source.Mtt(),source.Mtp(),source.Mpp()};
  const auto qmeta=qweights(0,freq.nt()/2+1);
  int windowSamples=0;
  for(int i=0;i<win.size();++i) if(win[i]!=0.0) windowSamples=i+2;
  const double densitySI=base.Density(1)(3000)*base.DensityNorm();
  const double vpSI=base.VPV(1)(3000)*base.VelocityNorm();
  const double vsSI=base.VSV(1)(3000)*base.VelocityNorm();
  std::ofstream meta(out/"meta.json");
  meta<<std::setprecision(17)<<"{\"solver\":\"DSpecM1D\",\"nt\":"<<freq.nt()
      <<",\"dt_s\":"<<freq.dt()*freq.timeNorm()<<",\"df_hz\":"<<freq.df()/freq.timeNorm()
      <<",\"damping_s\":"<<freq.ep()/freq.timeNorm()<<",\"pi\":3.1415926535"
      <<",\"omega_bin1\":"<<2*3.1415926535*freq.df()/freq.timeNorm()
      <<",\"maxstep_fraction\":"<<argv[3]<<",\"nq\":"<<argv[4]
      <<",\"physical_G\":"<<gPhysical<<",\"density_kg_m3\":"<<densitySI
      <<",\"vp_m_s\":"<<vpSI<<",\"vs_m_s\":"<<vsSI
      <<",\"radius_m\":"<<sem.mesh().PR()*sem.lengthScaleMetres()
      <<",\"source_element\":"<<sourceElem<<",\"source_element_bounds_m\":["
      <<sem.mesh().ELR(sourceElem)*sem.lengthScaleMetres()<<','
      <<sem.mesh().EUR(sourceElem)*sem.lengthScaleMetres()<<']'
      <<",\"source_aligned\":"<<(sem.sourceAligned()?"true":"false")
      <<",\"source_radius_m\":"<<sem.sourceRadius()*sem.lengthScaleMetres()
      <<",\"source_layer\":"<<sem.mesh().LayerNumber(sourceElem)
      <<",\"support_layer\":"<<supportLayer
      <<",\"direction_selector\":\""<<selector<<"\""
      <<",\"support_bounds_m\":["<<supportLowerM<<','<<supportUpperM<<']'
      <<",\"receiver_element\":"<<receiverElems.front()<<",\"receiver_element_bounds_m\":["
      <<sem.mesh().ELR(receiverElems.front())*sem.lengthScaleMetres()<<','
      <<sem.mesh().EUR(receiverElems.front())*sem.lengthScaleMetres()<<']'
      <<",\"shell_A_C_F_L_N_Pa\":["<<shellPa[0]<<','<<shellPa[1]<<','<<shellPa[2]<<','<<shellPa[3]<<','<<shellPa[4]<<']'
      <<",\"gravity\":\""<<(parameters.gravitation()==2?"full":"other")<<"\",\"attenuation\":"
      <<(parameters.attenuation() ? "true":"false")<<",\"degrees\":["<<parameters.lmin()<<','<<parameters.lmax()<<']'
      <<",\"q_bins\":[1,2],\"q_values\":["<<qmeta(0,1).real()<<','<<qmeta(0,1).imag()<<','
      <<qmeta(0,2).real()<<','<<qmeta(0,2).imag()<<"],\"window_samples\":"<<windowSamples
      <<",\"window\":\"Hann\",\"source_geocentric\":["<<sourceLocation[0]<<','<<sourceLocation[1]<<','<<sourceLocation[2]<<']'
      <<",\"receiver_geocentric\":["<<receiverLocation[0]<<','<<receiverLocation[1]<<','<<receiverLocation[2]<<']'
      <<",\"moment_tensor_Nm\":["<<tensor[0]<<','<<tensor[1]<<','<<tensor[2]<<','<<tensor[3]<<','<<tensor[4]<<','<<tensor[5]<<"]}\n";
  std::ofstream spectra(out/"spectra.csv"), derivatives(out/"derivatives.csv");
  spectra<<"component,bin,frequency_hz,omega_real,omega_imag,real_m_s,imag_m_s\n"<<std::setprecision(17);
  std::ofstream curves(out/"samples.csv");
  curves<<"observable,modulus,radius_m,kernel_dr_per_Pa_m,weight_m,delta_modulus_Pa\n"<<std::setprecision(17);
  derivatives<<"observable,modulus,objective,kernel,volume_kernel,positive_volume_kernel,negative_volume_kernel,fd_h_0_08,fd_h_0_04,fd_h_0_02,fd_h_0_01,fd_h_0_005,fd_h_0_0025,r4_0_08_0_04,r4_0_04_0_02,r4_0_02_0_01,r4_0_01_0_005,r4_0_005_0_0025,r6_0_08_0_04_0_02,r6_0_04_0_02_0_01,r6_0_02_0_01_0_005,r6_0_01_0_005_0_0025,richardson,relative_kernel_fd,source_kernel,source_delta,source_contribution,total_kernel\n"<<std::setprecision(17);
  for(int obs=0;obs<2;++obs) {
    auto q=qweights(obs==0?0:1,freq.nt()/2+1);
    const auto baseResult=SPARSESPEC::WindowedSensitivityEngine::evaluate(sem,source,parameters,freq,win,q,firstBin,endBin);
    if(obs==0) for(int b=firstBin;b<endBin;++b) for(int c=0;c<3;++c) {
      const auto v=baseResult.rawDisplacementSpectrum(c,b)*sem.lengthScaleMetres()*freq.timeNorm();
      spectra<<(c==0?"Z":c==1?"N":"E")<<','<<b<<','<<freq.f(b)/freq.timeNorm()<<','<<freq.w(b)/freq.timeNorm()<<','<<-freq.ep()/freq.timeNorm()<<','<<v.real()<<','<<v.imag()<<'\n';
    }
    const auto&k=baseResult.kernels;
    const std::array<const Eigen::MatrixXd*,5> kernel{&k.A,&k.C,&k.F,&k.L,&k.N};
    for(std::size_t mi=0;mi<mods.size();++mi) {
      if((selector=="source-local-other" && obs==0 && mi==0) ||
         (selector=="source-local-za" && (obs!=0 || mi!=0))) continue;
      double volume=0;
      double positiveVolume=0, negativeVolume=0;
      for(int e=0;e<sem.mesh().NE();++e)
        for(int i=0;i<sem.mesh().NN();++i) {
          const double value=(*kernel[mi])(e,i);
          const bool active=sem.mesh().LayerNumber(e)==supportLayer;
          const double delta=active?shellPa[mi]:0.0;
          if(active) {
            const double contribution=value*k.drWeightsMetres(e,i)*delta;
            volume+=contribution;
            if(contribution>0) positiveVolume+=contribution;
            else negativeVolume+=contribution;
          }
          curves<<(obs==0?"Z":"N")<<','<<names[mi]<<','<<sem.mesh().NodeRadius(e,i)*sem.lengthScaleMetres()<<','<<value<<','<<k.drWeightsMetres(e,i)<<','<<delta<<'\n';
        }
      const auto& sourceTerms=baseResult.sourcePointKernels;
      const double sourceKernel=mi==0?sourceTerms.A:mi==1?sourceTerms.C:mi==2?sourceTerms.F:mi==3?sourceTerms.L:sourceTerms.N;
      const int sourceLayer=sem.mesh().LayerNumber(sourceElem);
      const double sourceDelta=sourceLayer==supportLayer?shellPa[mi]:0.0;
      const double sourceContribution=sourceKernel*sourceDelta;
      const double predicted=volume+sourceContribution;
      std::vector<double> fd(steps.size());
      for(std::size_t h=0;h<steps.size();++h) {
        auto objective=[&](double amount) {
          Perturbed<Base,prem_norm<double>> changed(modelPath.string(),norm,mods[mi],supportLayer,amount);
          Full1D::SEM changedSem(changed,std::stod(argv[3]),std::stoi(argv[4]),2,source);
          auto val=SPARSESPEC::WindowedSensitivityEngine::evaluate(changedSem,source,parameters,freq,win,q,firstBin,endBin).objective;
          return val;
        };
        fd[h]=(objective(steps[h])-objective(-steps[h]))/(2*steps[h]);
      }
      std::vector<double> r4, r6;
      for(std::size_t h=0;h+1<fd.size();++h) r4.push_back((4*fd[h+1]-fd[h])/3);
      for(std::size_t h=0;h+1<r4.size();++h) r6.push_back((16*r4[h+1]-r4[h])/15);
      const double rich=sourceLocal?r4.back():(4*fd[2]-fd[1])/3;
      const double rel=std::abs(predicted-rich)/std::max({std::abs(predicted),std::abs(rich),1e-300});
      derivatives<<(obs==0?"Z":"N")<<','<<names[mi]<<','<<baseResult.objective<<','<<predicted<<','<<volume<<','<<positiveVolume<<','<<negativeVolume;
      if(sourceLocal) { for(const double v:fd) derivatives<<','<<v; }
      else { derivatives<<",,,,"<<fd[0]<<','<<fd[1]<<','<<fd[2]; }
      if(sourceLocal) { for(const double v:r4) derivatives<<','<<v; for(const double v:r6) derivatives<<','<<v; }
      else derivatives<<",,,,,,,,,";
      derivatives<<','<<rich<<','<<rel<<','<<sourceKernel<<','<<sourceDelta<<','<<sourceContribution<<','<<predicted<<'\n';
    }
  }
  return 0;
} catch(const std::exception&e) { std::cerr<<"DSpec kernel adapter: "<<e.what()<<'\n'; return 1; }
