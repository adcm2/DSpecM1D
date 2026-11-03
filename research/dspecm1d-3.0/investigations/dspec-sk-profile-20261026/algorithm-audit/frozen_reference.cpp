#include <map>
#include <bit>
#include <cstdint>
#include <limits>
#include <DSpecM1D/SensitivityKernels>
#include <DSpecM1D/All>
#include <DSpecM1D/ModelInput>
#include <DSpecM1D/Windowed3DKernels>
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
#include <sstream>
#include <stdexcept>
#include <string>

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
  return {{{-std::sin(lat)*std::cos(lon),-std::sin(lat)*std::sin(lon),std::cos(lat)},
           {-std::sin(lon),std::cos(lon),0.0},
           {std::cos(lat)*std::cos(lon),std::cos(lat)*std::sin(lon),std::sin(lat)}}};
}

int main(int argc,char **argv) try {
 if(argc!=9)throw std::runtime_error("usage: frozen-reference MODEL SAMPLES OUT L NF OBJECTIVES BACKEND RADIAL");
 const fs::path model=fs::absolute(argv[1]),samples=fs::absolute(argv[2]),output=fs::absolute(argv[3]);
 const int lmax=std::stoi(argv[4]),nfreq=std::stoi(argv[5]),nobj=std::stoi(argv[6]);
 const std::string radial=argv[8];
 if(!fs::is_directory(output)||fs::directory_iterator(output)!=fs::directory_iterator{})throw std::runtime_error("nonempty output");
 const auto setupStart=Clock::now();const auto paramPath=output/"parameters.txt";
 {std::ofstream f(paramPath);f<<paramsText(model,output,lmax);}
 InputParametersNew input(paramPath.string(),8,1,0.05,1.0,0.1,0.0,1);
 auto &parameters=input.inputParameters();auto &frequency=input.freqFull();auto &source=input.cmt();auto &base=input.earthModel();
 Full1D::SEM sem(base,0.05,8,lmax,source);auto points=readPoints(samples);
 if(radial=="native"){
  std::map<std::pair<std::uint64_t,std::uint64_t>,std::array<double,2>> unique;
  for(const auto&p:points)unique.try_emplace(std::pair{std::bit_cast<std::uint64_t>(p.thetaRadians),std::bit_cast<std::uint64_t>(p.phiRadians)},std::array<double,2>{p.thetaRadians,p.phiRadians});
  points.clear();for(int e=0;e<sem.mesh().NE();++e)for(int q=0;q<sem.mesh().NN();++q){
   const double radius=sem.mesh().NodeRadius(e,q)*sem.lengthScaleMetres();if(radius<=0)continue;
   for(const auto &[key,a]:unique){(void)key;points.push_back({radius,a[0],a[1],e});}
  }
 }
 const auto timeWindow=window();const std::vector<std::array<double,2>> angles{{std::numbers::pi/4,std::numbers::pi/2}};
 const double setup=std::chrono::duration<double>(Clock::now()-setupStart).count();
 std::ofstream arrays(output/"arrays_reference.bin",std::ios::binary),coords(output/"coordinates.bin",std::ios::binary);
 const std::uint64_t count=points.size();writeScalar(coords,count);
 for(const auto&p:points){writeScalar(coords,p.radiusMetres);writeScalar(coords,p.thetaRadians);writeScalar(coords,p.phiRadians);writeScalar(coords,p.radialElement);}
 coords.close();std::vector<double> times;
 for(int obj=0;obj<nobj;++obj){const auto start=Clock::now();
  const auto result=SPARSESPEC::Windowed3DKernelEngine::evaluate(sem,source,parameters,frequency,timeWindow,weights(obj),1,nfreq+1,points,angles);
  times.push_back(std::chrono::duration<double>(Clock::now()-start).count());writeEvaluation(arrays,result);
 }
 arrays.close();fs::remove(paramPath);std::ofstream meta(output/"result.json");meta<<std::setprecision(17);
 meta<<"{\"reference\":\"accepted production Phase1\",\"degree\":"<<lmax<<",\"frequency_bins\":"<<nfreq<<",\"objectives\":"<<nobj<<",\"sample_count\":"<<count<<",\"setup_s\":"<<setup<<",\"api_s\":[";
 for(std::size_t i=0;i<times.size();++i)meta<<(i?",":"")<<times[i];meta<<"]}\n";
 std::cout<<"Reference "<<lmax<<" "<<nfreq<<" "<<count<<" complete\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}
