from pathlib import Path
import subprocess,json
R=Path(__file__).resolve().parent;repo=Path('/home/adcm2/Documents/Research/Projects/DSpecM1D');src=repo/'benchmarks/yspecpp/kernels/3d/dspec_3d_adapter.cpp'
d=(R/'fast/v5/driver.cpp').read_text();includes=d[:d.index('namespace fs =')];body=d[d.index('namespace fs ='):d.index('int main(int argc,char **argv)')]
includes=includes.replace('"fast_engine.hpp"','"fast/v5/fast_engine.hpp"').replace('"engine_headers/','"fast/v5/engine_headers/')
(R/'born_candidate.hpp').write_text(includes+'\nnamespace auditcandidate {\n'+body+'\n}\n')
s=src.read_text().replace('#include <DSpecM1D/All>','#include "born_candidate.hpp"\n#include <DSpecM1D/All>')
s=s.replace('  for(int obs=0;obs<2;++obs){','  SPARSESPEC::Windowed3DFieldEngine::SamplingPlanCache auditPlan;\n  for(int obs=0;obs<2;++obs){')
a='''    const auto result=SPARSESPEC::Windowed3DKernelEngine::evaluate(
        sem,source,parameters,freq,window,q,firstBin,endBin,evaluationPoints,angles);'''
b='''    double auditSeconds=0;std::array<double,5> auditStages{};
    const auto result=auditcandidate::candidate(sem,source,parameters,freq,window,q,endBin,
        evaluationPoints,angles,auditSeconds,nullptr,obs,3,auditPlan,nullptr,true,true,auditStages);'''
assert a in s;s=s.replace(a,b);(R/'born_candidate.cpp').write_text(s)
cmd=json.loads((R/'reference_build_command.json').read_text());cmd[cmd.index(str(R/'frozen_reference.cpp'))]=str(R/'born_candidate.cpp');cmd[cmd.index('-o')+1]=str(R/'bin/born-candidate');cmd+=['-DDSFAST_HAVE_CBLAS','-DDSFAST_DISABLE_PROFILE','-lopenblas']
(R/'born_build_command.json').write_text(json.dumps(cmd,indent=2)+'\n')
with (R/'born-build.log').open('w') as f:subprocess.run(['python3',str(R/'run_compile.py'),'--']+cmd,stdout=f,stderr=f,check=True)
print('Born candidate built',flush=True)
