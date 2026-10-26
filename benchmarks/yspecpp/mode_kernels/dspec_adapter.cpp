#include <DSpecM1D/All>
#include <DSpecM1D/SensitivityKernels>
#include <DSpecM1D/ModelInput>
#include <DSpecM1D/src/InputParser.h>
#include <DSpecM1D/src/NormClass.h>
#include <DSpecM1D/src/SourceInfo.h>
#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>

namespace fs = std::filesystem;
template<class Base, class Norm> class Perturbed final : public Base {
public:
  Perturbed(const std::string& path,const Norm& norm,double amount,std::string property)
      : Base(path,norm),amount_(amount),property_(std::move(property)) {}
  auto Density(int layer) const { return scaled(Base::Density(layer),"rho"); }
  auto A(int layer) const { return scaled(Base::A(layer),"A"); }
  auto C(int layer) const { return scaled(Base::C(layer),"C"); }
  auto F(int layer) const { return scaled(Base::F(layer),"F"); }
  auto L(int layer) const { return scaled(Base::L(layer),"L"); }
  auto N(int layer) const { return scaled(Base::N(layer),"N"); }
private:
  template<class Field> auto scaled(Field field,const std::string& name) const {
    const bool selected=property_==name || (property_=="LN"&&(name=="L"||name=="N"));
    const double factor=1.0+(selected?amount_:0.0);
    return [field=std::move(field),factor](double r){return field(r)*factor;};
  }
  double amount_; std::string property_;
};
std::string paramsText(const fs::path& model,const fs::path& out,int type,double fmin,double fmax,const std::string& geometry) {
  std::ostringstream p;
  p<<std::setprecision(17)<<'"'<<(out/"unused").string()<<"\n\""<<model.string()<<"\"\n"
   <<type<<"\n0\n2\n0\n0\n1e-8\n2\n2\n"<<fmin<<'\n'<<fmax<<"\n4320\n60\n0\n0\n0\n0\n50\n0\n0\n";
  if(geometry=="toroidal") p<<"0\n1e19\n0\n0\n0\n0\n";
  else p<<"1e19\n0\n0\n0\n0\n0\n";
  p<<"3\n1\n0 90\n";
  return p.str();
}
Eigen::VectorXd window(const SpectraSolver::FreqFull& f) {
  constexpr int samples=4321; Eigen::VectorXd w=Eigen::VectorXd::Zero(f.nt());
  for(int i=0;i<samples&&i<f.nt();++i) w[i]=0.5-0.5*std::cos(2.0*3.14159265358979323846*i/(samples-1));
  return w;
}
void guardOutput(const fs::path& out) {
  if(fs::exists(out)&&(!fs::is_directory(out)||fs::directory_iterator(out)!=fs::directory_iterator{}))
    throw std::runtime_error("output directory must be empty: "+out.string());
  fs::create_directories(out);
}
int main(int argc,char**argv) try {
  if(argc<11||argc>14) throw std::runtime_error("usage: kernel_adapter ACTION MODEL OUT MAXSTEP NQ FMIN_MHZ FMAX_MHZ DELTA TYPE [PROPERTY QCSV [COMPONENT [GEOMETRY]]]");
  const std::string action=argv[1]; const fs::path modelPath=fs::absolute(argv[2]),out=fs::absolute(argv[3]);
  const double maxstep=std::stod(argv[4]),delta=std::stod(argv[8]);
  const int nq=std::stoi(argv[5]),type=std::stoi(argv[9]);
  const double fmin=std::stod(argv[6]),fmax=std::stod(argv[7]);
  const std::string property=argc>=12?argv[10]:"LN";
  const fs::path qpath=argc>=12?argv[11]:argv[10];
  const std::string component=argc>=13?argv[12]:"north";
  const std::string geometry=argc==14?argv[13]:"toroidal";
  if(geometry!="toroidal"&&geometry!="spheroidal") throw std::runtime_error("geometry must be toroidal or spheroidal");
  const int outputRow=component=="up"?0:component=="north"?1:-1;
  if(outputRow<0) throw std::runtime_error("component must be up or north");
  guardOutput(out);
  prem_norm<double> norm; using Base=EarthModels::ModelInput<double>;
  Base base(modelPath.string(),norm);
  {std::ofstream p(out/"params.txt");p<<paramsText(modelPath,out,type,fmin,fmax,geometry);}
  InputParametersNew input((out/"params.txt").string(),nq,1,maxstep,1.0,0.1,0.0,1);
  auto& parameters=input.inputParameters(); SourceInfo::EarthquakeCMT source(parameters);
  Perturbed<Base,prem_norm<double>> changed(modelPath.string(),norm,delta,property);
  Full1D::SEM sem(changed,maxstep,nq,2,source);
  const auto& freq=input.freqFull(); const auto win=window(freq);
  const int expectedFirst=fmin>=0.28?146:116, expectedEnd=fmin>=0.28?178:159;
  if(freq.nt()!=8192||std::abs(freq.dt()*freq.timeNorm()-60.0)>1e-12||
     std::abs(freq.ep()/freq.timeNorm()-5.0/259200.0)>1e-15||
     std::abs(freq.df()/freq.timeNorm()-1.0/(8192.0*60.0))>1e-15||freq.i1()!=expectedFirst||freq.i2()!=expectedEnd)
    throw std::runtime_error("unexpected modal FFT grid for the selected frequency band");
  {std::ofstream m(out/"meta.json");m<<std::setprecision(17)<<"{\"solver\":\"DSpecM1D\",\"nt\":"<<freq.nt()
    <<",\"dt_s\":"<<freq.dt()*freq.timeNorm()<<",\"df_hz\":"<<freq.df()/freq.timeNorm()
    <<",\"damping_s\":"<<freq.ep()/freq.timeNorm()<<",\"first\":"<<freq.i1()<<",\"end\":"<<freq.i2()
    <<",\"source_radius_m\":"<<sem.sourceRadius()*sem.lengthScaleMetres()<<",\"type\":"<<type
    <<",\"component\":\""<<component<<"\",\"geometry\":\""<<geometry<<"\",\"property\":\""<<property<<"\",\"delta\":"<<delta<<"}\n";}
  if(action=="spectrum") {
    Eigen::MatrixXcd q=Eigen::MatrixXcd::Zero(3,freq.nt()/2+1);
    const auto result=SPARSESPEC::WindowedSensitivityEngine::evaluate(sem,source,parameters,freq,win,q,freq.i1(),freq.i2());
    std::ofstream spec(out/"spectrum.csv");spec<<"bin,frequency_hz,raw_real,raw_imag,processed_real,processed_imag\n"<<std::setprecision(17);
    for(int k=freq.i1();k<freq.i2();++k){auto raw=result.rawDisplacementSpectrum(outputRow,k)*sem.lengthScaleMetres()*freq.timeNorm();auto proc=result.processedSISpectrum(outputRow,k);spec<<k<<','<<freq.f(k)/freq.timeNorm()<<','<<raw.real()<<','<<raw.imag()<<','<<proc.real()<<','<<proc.imag()<<'\n';}
    return 0;
  }
  if(action!="kernel") throw std::runtime_error("unknown action");
  std::map<std::string,Eigen::MatrixXcd> groups;std::ifstream in(qpath);if(!in)throw std::runtime_error("cannot read q weights");
  std::string line;std::getline(in,line);while(std::getline(in,line)){std::replace(line.begin(),line.end(),',',' ');std::istringstream row(line);std::string name;int k;double re,im;if(row>>name>>k>>re>>im){if(!groups.contains(name))groups.emplace(name,Eigen::MatrixXcd::Zero(3,freq.nt()/2+1));groups.at(name)(outputRow,k)={re,im};}}
  if(groups.empty())throw std::runtime_error("q file has no weights");
  std::ofstream summary(out/"kernel_summary.csv");summary<<"basis,distributed_A,distributed_C,distributed_F,distributed_L,distributed_N,distributed_rho,source_A,source_C,source_F,source_L,source_N,total\n"<<std::setprecision(17);
  for(const auto&[name,q]:groups){
    const auto result=SPARSESPEC::WindowedSensitivityEngine::evaluate(sem,source,parameters,freq,win,q,freq.i1(),freq.i2());
    double sums[6]={};
    std::ofstream curve(out/("kernel_"+name+".csv"));
    curve<<"radius_m,radial_A,radial_C,radial_F,radial_L,radial_N,radial_rho,weight_m\n"<<std::setprecision(17);
    for(int e=0;e<sem.mesh().NE();++e) {
      for(int i=0;i<sem.mesh().NN();++i) {
        const int layer=sem.mesh().LayerNumber(e);
        const double r=sem.mesh().NodeRadius(e,i), w=result.kernels.drWeightsMetres(e,i);
        const double A=base.A(layer)(r)*sem.pressureScalePascals();
        const double C=base.C(layer)(r)*sem.pressureScalePascals();
        const double F=base.F(layer)(r)*sem.pressureScalePascals();
        const double L=base.L(layer)(r)*sem.pressureScalePascals();
        const double N=base.N(layer)(r)*sem.pressureScalePascals();
        const double rho=base.Density(layer)(r)*base.DensityNorm();
        const double values[]={result.kernels.A(e,i)*A,result.kernels.C(e,i)*C,
          result.kernels.F(e,i)*F,result.kernels.L(e,i)*L,result.kernels.N(e,i)*N,
          result.materialKernels.densityKgPerM3(e,i)*rho};
        for(int j=0;j<6;++j) sums[j]+=values[j]*w;
        curve<<r*sem.lengthScaleMetres();
        for(double v:values) curve<<','<<v;
        curve<<','<<w<<'\n';
      }
    }
    const double sourceRadius=sem.sourceRadius();
    const int sourceElement=sem.sourceElement(source);
    const int sourceLayer=sem.mesh().LayerNumber(sourceElement);
    const double sourceA=base.A(sourceLayer)(sourceRadius)*sem.pressureScalePascals();
    const double sourceC=base.C(sourceLayer)(sourceRadius)*sem.pressureScalePascals();
    const double sourceF=base.F(sourceLayer)(sourceRadius)*sem.pressureScalePascals();
    const double sourceL=base.L(sourceLayer)(sourceRadius)*sem.pressureScalePascals();
    const double sourceN=base.N(sourceLayer)(sourceRadius)*sem.pressureScalePascals();
    const auto& source=result.sourcePointKernels;
    const double point[5]={source.A*sourceA,source.C*sourceC,source.F*sourceF,
                           source.L*sourceL,source.N*sourceN};
    summary<<name;
    for(double value:sums) summary<<','<<value;
    for(double value:point) summary<<','<<value;
    const double total=geometry=="toroidal" ? sums[3]+sums[4]+point[3]+point[4]
      : sums[0]+sums[1]+sums[2]+sums[3]+sums[4]+sums[5]+
        point[0]+point[1]+point[2]+point[3]+point[4];
    summary<<','<<total<<'\n';
  }
  return 0;
} catch(const std::exception&e){std::cerr<<"DSpec spheroidal modal adapter: "<<e.what()<<'\n';return 1;}
