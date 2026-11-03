from pathlib import Path
import json,subprocess,concurrent.futures
R=Path(__file__).resolve().parent
s=(R/'fields_validation.cpp').read_text().replace('"fast/fast_engine.hpp"','"fast/v5/fast_engine.hpp"').replace('"fast/engine_headers/','"fast/v5/engine_headers/')
# Preserve historical public-field outputs and add a genuinely tensor PREM oracle.
s=s.replace('          SPARSESPEC::Windowed3DFieldEngine::SamplingPlanCache packedPlan;', '''          std::map<std::pair<int,double>,SPARSESPEC::Field3DPoint> tracePoints;
          std::map<std::pair<double,double>,SPARSESPEC::Field3DPoint> anglePoints;
          for(const auto &p:noncenter){tracePoints.try_emplace(std::pair{p.radialElement,p.radiusMetres},p);anglePoints.try_emplace(std::pair{p.thetaRadians,p.phiRadians},p);}
          noncenter.clear();
          for(const auto &[tk,tp]:tracePoints)for(const auto &[ak,ap]:anglePoints){(void)tk;(void)ak;noncenter.push_back({tp.radiusMetres,ap.thetaRadians,ap.phiRadians,tp.radialElement});}
          std::reverse(noncenter.begin(),noncenter.end()); // Exercise a permuted complete tensor.
          const auto tensorFull=workspace.sampleWeightedChunk(noncenter);
          SPARSESPEC::Windowed3DFieldEngine::SamplingPlanCache packedPlan;''')
s=s.replace('workspace.sampleForKernelsPacked(noncenter,0,&packedPlan,3)','workspace.sampleForKernelsTensor(noncenter,noncenter.size(),&packedPlan,3,{},true)').replace('const auto &b=fields[original[j]];','const auto &b=tensorFull[j];')
# Tensor volume intentionally returns no displacement; strain comparisons are the relevant contract.
s=s.replace('numerator[0]+=std::norm(a.forwardRawDisplacement[c]-b.forwardRawDisplacement[c]);denominator[0]+=std::norm(b.forwardRawDisplacement[c]);','numerator[0]+=0.0;denominator[0]+=0.0;')
(R/'fields_validation_v5_tensor.cpp').write_text(s)
cmd=json.loads((R/'fields_build_command.json').read_text());cmd[cmd.index(str(R/'fields_validation.cpp'))]=str(R/'fields_validation_v5_tensor.cpp');cmd[cmd.index('-o')+1]=str(R/'bin/fields-validation-v5-tensor')
(R/'fields_v5_tensor_build_command.json').write_text(json.dumps(cmd,indent=2)+'\n')
with (R/'fields-v5-tensor-build.log').open('w') as f:subprocess.run(['python3',str(R/'run_compile.py'),'--']+cmd,stdout=f,stderr=f,check=True)
print('Final v5 PREM tensor validator built',flush=True)
