from pathlib import Path
import json,subprocess
root=Path(__file__).resolve().parent
s=(root/'timing_current.cpp').read_text().split('int main(int argc, char **argv) try {')[0]
s+=r'''
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
'''
(root/'chunk_probe.cpp').write_text(s)
c=json.loads((root/'build_commands.json').read_text())['current']
c[c.index(str(root/'timing_current.cpp'))]=str(root/'chunk_probe.cpp');c[c.index('-o')+1]=str(root/'bin/chunk_probe')
(root/'chunk_build_command.json').write_text(json.dumps(c,indent=2)+'\n')
with (root/'chunk_build.log').open('w') as out:p=subprocess.run(c,stdout=out,stderr=out)
if p.returncode:raise RuntimeError((root/'chunk_build.log').read_text()[-4000:])
print('chunk driver built')
