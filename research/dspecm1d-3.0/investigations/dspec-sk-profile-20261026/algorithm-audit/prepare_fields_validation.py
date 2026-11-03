from pathlib import Path
import json,subprocess
R=Path(__file__).resolve().parent;repo=Path('/home/adcm2/Documents/Research/Projects/DSpecM1D')
s=(repo/'build/kernel-reconstruction-pole-frames-20261026/fields.cpp').read_text()
s='#include "fast/fast_engine.hpp"\n#include "fast/engine_headers/SensitivityKernels.h"\n'+s.replace('#include <DSpecM1D/Windowed3DKernels>','#include "fast/engine_headers/Windowed3DKernels.h"')
s=s.replace('if(argc!=2)throw std::invalid_argument("contact output required");','if(argc!=3)throw std::invalid_argument("contact output and backend required");\n  SPARSESPEC::FixedMeshAdjointSolver::setDiagnosticBackend(std::stoi(argv[2]));')
a='          const auto fields=workspace.sampleWeightedChunk(points);'
b='''          const auto fields=workspace.sampleWeightedChunk(points);
          std::vector<SPARSESPEC::Field3DPoint> noncenter;
          std::vector<std::size_t> original;
          for(std::size_t i=0;i<points.size();++i)if(points[i].radiusMetres>0){noncenter.push_back(points[i]);original.push_back(i);}
          SPARSESPEC::Windowed3DFieldEngine::SamplingPlanCache packedPlan;
          const auto packed=workspace.sampleForKernelsPacked(noncenter,0,&packedPlan,3);
          std::array<double,3> numerator{},denominator{};
          for(std::size_t j=0;j<packed.size();++j){
            const auto &a=packed[j];const auto &b=fields[original[j]];
            for(int c=0;c<3;++c){numerator[0]+=std::norm(a.forwardRawDisplacement[c]-b.forwardRawDisplacement[c]);denominator[0]+=std::norm(b.forwardRawDisplacement[c]);}
            const auto av=a.forwardRawStrain,bv=b.forwardRawStrain,aa=a.adjointDualStrain,ba=b.adjointDualStrain;
            if(a.fluid){numerator[1]+=std::norm(av.rr+av.horizontalTrace-bv.rr-bv.horizontalTrace);denominator[1]+=std::norm(bv.rr+bv.horizontalTrace);
              numerator[2]+=std::norm(aa.rr+aa.horizontalTrace-ba.rr-ba.horizontalTrace);denominator[2]+=std::norm(ba.rr+ba.horizontalTrace);}
            else{const std::array<Complex,6> af{av.rr,av.horizontalTrace,av.radialPlus,av.radialMinus,av.horizontalPlus,av.horizontalMinus},bf{bv.rr,bv.horizontalTrace,bv.radialPlus,bv.radialMinus,bv.horizontalPlus,bv.horizontalMinus}, ax{aa.rr,aa.horizontalTrace,aa.radialPlus,aa.radialMinus,aa.horizontalPlus,aa.horizontalMinus},bx{ba.rr,ba.horizontalTrace,ba.radialPlus,ba.radialMinus,ba.horizontalPlus,ba.horizontalMinus};
              for(int c=0;c<6;++c){numerator[1]+=std::norm(af[c]-bf[c]);denominator[1]+=std::norm(bf[c]);numerator[2]+=std::norm(ax[c]-bx[c]);denominator[2]+=std::norm(bx[c]);}}
          }
          std::cerr<<std::setprecision(17)<<"PACKERR {\\"case\\":"<<caseIndex<<",\\"frequency\\":"<<frequency<<",\\"relative_norms\\":[";
          for(int c=0;c<3;++c)std::cerr<<(c?",":"")<<(denominator[c]==0?(numerator[c]==0?0:std::numeric_limits<double>::infinity()):std::sqrt(numerator[c]/denominator[c]));
          std::cerr<<"]}\\n";
'''
assert s.count(a)==1;s=s.replace(a,b)
(R/'fields_validation.cpp').write_text(s)
cmd=json.loads((R.parent/'phase1/attempt-02/build_commands.json').read_text())['fields']
cmd[cmd.index(str(repo/'build/kernel-reconstruction-pole-frames-20261026/fields.cpp'))]=str(R/'fields_validation.cpp')
cmd[cmd.index('-o')+1]=str(R/'bin/fields-validation')
cmd+=['-DDSFAST_HAVE_CBLAS','-DDSFAST_DISABLE_PROFILE','-lopenblas']
(R/'fields_build_command.json').write_text(json.dumps(cmd,indent=2)+'\n')
with (R/'fields-build.log').open('w') as log:subprocess.run(['python3',str(R/'run_compile.py'),'--']+cmd,stdout=log,stderr=log,check=True)
print('Packed/full PREM field validator built',flush=True)
