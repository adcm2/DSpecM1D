#include <DSpecM1D/All>
#include <DSpecM1D/ModelInput>
#include <DSpecM1D/Windowed3DKernels>
#include <DSpecM1D/src/InputParser.h>
#include <DSpecM1D/src/NormClass.h>
#include <DSpecM1D/src/SourceInfo.h>

#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/resource.h>
#include <vector>

namespace fs = std::filesystem;
using Complex = std::complex<double>;
constexpr int nt = 128;
constexpr int bins = nt / 2 + 1;
constexpr std::array<const char *, 5> moduli{"A", "C", "F", "L", "N"};

std::string paramsText(const fs::path &model, const fs::path &output, int lmax) {
  return "\"" + (output / "unused").string() + "\"\n\"" + model.string() +
      "\"\n4\n0\n2\n0\n0\n1e-5\n0\n" + std::to_string(lmax) +
      "\n0.1\n0.8\n30\n20\n0\n0\n0\n0\n50\n10\n20\n1e19\n2e18\n-3e18\n-5e18\n4e18\n-5e18\n3\n1\n45 90\n";
}

Eigen::VectorXd hannWindow() {
  Eigen::VectorXd window = Eigen::VectorXd::Zero(nt);
  for (int n = 0; n <= 90; ++n)
    window[n] = 0.5 - 0.5 * std::cos(2.0 * std::numbers::pi * n / 90.0);
  return window;
}

Eigen::MatrixXcd objectiveWeights(int objective) {
  Eigen::MatrixXcd q = Eigen::MatrixXcd::Zero(3, bins);
  q(objective, 1) = Complex(0.7, 0.2);
  q(objective, 2) = Complex(-0.25, 0.45);
  return q;
}

std::vector<SPARSESPEC::Field3DPoint> readPoints(const fs::path &path) {
  std::ifstream file(path);
  if (!file)
    throw std::runtime_error("cannot read shared tutorial samples CSV");
  std::string line;
  std::getline(file, line);
  std::vector<SPARSESPEC::Field3DPoint> points;
  while (std::getline(file, line)) {
    std::istringstream row(line);
    std::array<std::string, 11> values;
    for (auto &value : values)
      if (!std::getline(row, value, ','))
        throw std::runtime_error("invalid shared tutorial sample row");
    if (values[0] == "J_Z" && values[1] == "A")
      points.push_back({std::stod(values[2]), std::stod(values[3]),
                        std::stod(values[4]), std::stoi(values[5])});
  }
  if (points.empty())
    throw std::runtime_error("shared samples contain no A-modulus coordinates");
  return points;
}

long peakRssBytes() {
  rusage usage{};
  if (getrusage(RUSAGE_SELF, &usage) != 0)
    throw std::runtime_error("getrusage failed");
#if defined(__APPLE__)
  return usage.ru_maxrss;
#else
  return usage.ru_maxrss * 1024L;
#endif
}


int main(int argc,char**argv) try {
  if(argc!=4)throw std::runtime_error("MODEL SAMPLES OUTPUT");
  const fs::path output=fs::absolute(argv[3]);
  const auto paramPath=output/"parameters.txt";
  {std::ofstream f(paramPath);f<<paramsText(fs::absolute(argv[1]),output,12);}
  InputParametersNew input(paramPath.string(),8,1,0.05,1.0,0.1,0.0,1);
  auto &params=input.inputParameters();auto &freq=input.freqFull();auto &source=input.cmt();
  Full1D::SEM sem(input.earthModel(),0.05,8,12,source);
  const auto points=readPoints(fs::absolute(argv[2]));
  SPARSESPEC::WindowedDisplacementMap map(freq,sem.lengthScaleMetres(),hannWindow());
  const Eigen::MatrixXcd q=map.adjoint(objectiveWeights(0));
  const Complex z=freq.w().at(1)+Complex(0,-freq.ep());
  auto start=std::chrono::steady_clock::now();
  auto work=SPARSESPEC::Windowed3DFieldEngine::prepareFrequencyWorkspace(sem,source,params,z,q.col(1));
  std::cout<<std::setprecision(17)<<"{\"prepare_s\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<",\"points\":"<<points.size()<<",\"sample_bytes\":"<<sizeof(SPARSESPEC::Field3DSample)<<",\"runs\":[";
  std::vector<SPARSESPEC::ModulusKernel3D> reference;
  bool first=true;
  for(std::size_t chunk:{points.size(),std::size_t{4096},std::size_t{256},points.size()}){
    start=std::chrono::steady_clock::now();
    std::vector<SPARSESPEC::ModulusKernel3D> kernels; kernels.reserve(points.size());
    for(std::size_t a=0;a<points.size();a+=chunk){
      std::vector<SPARSESPEC::Field3DPoint> p(points.begin()+a,points.begin()+std::min(points.size(),a+chunk));
      auto k=SPARSESPEC::Windowed3DKernelEngine::evaluateFrequency(work,p);
      kernels.insert(kernels.end(),k.begin(),k.end());
    }
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    if(reference.empty())reference=kernels;
    double maxdiff=0,peak=0;
    for(std::size_t i=0;i<kernels.size();++i){
      const auto &a=kernels[i],&b=reference[i];
      for(auto pair:{std::pair{a.A,b.A},std::pair{a.C,b.C},std::pair{a.F,b.F},std::pair{a.L,b.L},std::pair{a.N,b.N}}){maxdiff=std::max(maxdiff,std::abs(pair.first-pair.second));peak=std::max(peak,std::abs(pair.second));}
    }
    std::cout<<(first?"":",")<<"{\"chunk\":"<<chunk<<",\"seconds\":"<<seconds<<",\"max_diff_over_peak\":"<<maxdiff/peak<<"}"<<std::flush;first=false;
  }
  std::cout<<"]}\n";
  return 0;
} catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
