from pathlib import Path
import json,subprocess
R=Path(__file__).resolve().parent;repo=Path('/home/adcm2/Documents/Research/Projects/DSpecM1D')
s=(R/'fast/driver.cpp').read_text();s=s[:s.index('SPARSESPEC::Windowed3DEvaluationResult candidate(')]
s=s.replace('#include "fast_engine.hpp"','#include <map>\n#include <bit>\n#include <cstdint>\n#include <limits>')
s=s.replace('#include "engine_headers/SensitivityKernels.h"','#include <DSpecM1D/SensitivityKernels>').replace('#include "engine_headers/Windowed3DKernels.h"','#include <DSpecM1D/Windowed3DKernels>')
s+='''int main(int argc,char **argv) try {
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
 meta<<"{\\"reference\\":\\"accepted production Phase1\\",\\"degree\\":"<<lmax<<",\\"frequency_bins\\":"<<nfreq<<",\\"objectives\\":"<<nobj<<",\\"sample_count\\":"<<count<<",\\"setup_s\\":"<<setup<<",\\"api_s\\":[";
 for(std::size_t i=0;i<times.size();++i)meta<<(i?",":"")<<times[i];meta<<"]}\\n";
 std::cout<<"Reference "<<lmax<<" "<<nfreq<<" "<<count<<" complete\\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<"\\n";return 1;}
'''
(R/'frozen_reference.cpp').write_text(s)
cmd=json.loads((R/'profile_build_command.json').read_text());cmd[cmd.index(str(R/'profile_current.cpp'))]=str(R/'frozen_reference.cpp');cmd[cmd.index('-o')+1]=str(R/'bin/frozen-reference')
(R/'reference_build_command.json').write_text(json.dumps(cmd,indent=2)+'\n')
with (R/'reference-build.log').open('w') as log:subprocess.run(cmd,stdout=log,stderr=log,check=True)
print('Frozen reference built',flush=True)
