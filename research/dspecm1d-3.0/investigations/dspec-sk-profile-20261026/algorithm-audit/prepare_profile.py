from pathlib import Path
import json,hashlib,re,subprocess,difflib,sys
ROOT=Path(__file__).resolve().parent
OLD=ROOT.parent
REPO=Path('/home/adcm2/Documents/Research/Projects/DSpecM1D')
SHADOW=ROOT/'profile_headers';SHADOW.mkdir(exist_ok=True)
(ROOT/'bin').mkdir(exist_ok=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
if not (ROOT/'before.json').exists():
 files=subprocess.check_output(['git','ls-files'],cwd=REPO,text=True).splitlines()
 manifest={f:sha(REPO/f) for f in files if (REPO/f).is_file()}
 extra={str(p.relative_to(REPO)):sha(p) for p in (REPO/'validation/phase0').rglob('*') if p.is_file()}
 (ROOT/'before.json').write_text(json.dumps({'repo':str(REPO),'head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip(),'status':subprocess.check_output(['git','status','--short'],cwd=REPO,text=True),'tracked_hashes':manifest,'phase0_hashes':extra},indent=2)+'\n')
 for f in ['DSpecM1D/src/Windowed3DFields.h','DSpecM1D/src/Windowed3DKernels.h','DSpecM1D/src/SensitivityKernels.h','tests/test_sensitivity_kernels.cpp','implementation_status.md','DSpecM1D/src/SEM/SEM.h','DSpecM1D/src/SEM/SEMConstructor.h']:
  dst=ROOT/'before'/f;dst.parent.mkdir(parents=True,exist_ok=True);dst.write_bytes((REPO/f).read_bytes())
 (ROOT/'before.patch').write_bytes(subprocess.check_output(['git','diff'],cwd=REPO))
(SHADOW/'profile.hpp').write_text('''#pragma once
#include <chrono>
#include <map>
#include <string>
#include <iostream>
namespace auditprofile {
struct Entry {double seconds=0; unsigned long long calls=0;};
inline std::map<std::string,Entry> data;
inline std::map<std::string,unsigned long long> counters;
struct Scope {const char* name;std::chrono::steady_clock::time_point start; Scope(const char* n):name(n),start(std::chrono::steady_clock::now()){} ~Scope(){auto& e=data[name];e.seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();++e.calls;}};
struct Dump {~Dump(){std::cerr<<"AUDITPROFILE {\\\"timers\\\":{";bool first=true;for(auto&[n,e]:data){std::cerr<<(first?"":",")<<"\\\""<<n<<"\\\":{\\\"seconds\\\":"<<e.seconds<<",\\\"calls\\\":"<<e.calls<<"}";first=false;}std::cerr<<"},\\\"counters\\\":{";first=true;for(auto&[n,v]:counters){std::cerr<<(first?"":",")<<"\\\""<<n<<"\\\":"<<v;first=false;}std::cerr<<"}}\\n";}};inline Dump dump;
}
''')
def inc(s):
 return re.sub(r'#include "([^"]+)"',lambda m:'#include "'+str(REPO/'DSpecM1D/src'/m[1])+'"',s)
def one(s,a,b):
 assert s.count(a)==1,(a,s.count(a));return s.replace(a,b)
s=inc((REPO/'DSpecM1D/src/SensitivityKernels.h').read_text())
s='#include "profile.hpp"\n'+s
s=s.replace('    solver.compute(h);','    {auditprofile::Scope t("factorization");solver.compute(h);}')
s=s.replace('    SparseComplex h = systemMatrix(sem, family, degree, z);','    SparseComplex h;{auditprofile::Scope t("system_matrix");h=systemMatrix(sem,family,degree,z);}')
s=s.replace('    Eigen::MatrixXcd fields = solver.solve(compactRhs);','    Eigen::MatrixXcd fields;{auditprofile::Scope t("forward_substitution");fields=solver.solve(compactRhs);}')
s=s.replace('    Eigen::MatrixXcd fields = solver.solveAdjoint(rightHandSides);','    Eigen::MatrixXcd fields;{auditprofile::Scope t("adjoint_substitution");fields=solver.solveAdjoint(rightHandSides);}')
s=s.replace('    fields = solver.solveAdjoint(rightHandSides);','    {auditprofile::Scope t("adjoint_global_substitution");fields=solver.solveAdjoint(rightHandSides);}')
s=s.replace('    basisFields = solver.solveAdjoint(basisRightHandSides);','    {auditprofile::Scope t("adjoint_basis_substitution");basisFields=solver.solveAdjoint(basisRightHandSides);}')
(SHADOW/'SensitivityKernels.h').write_text(s)
s=inc((REPO/'DSpecM1D/src/Windowed3DFields.h').read_text()).replace(str(REPO/'DSpecM1D/src/SensitivityKernels.h'),'SensitivityKernels.h')
s=one(s,'    if(blocks.empty())return;','    auditprofile::Scope whole("pole_synthesis");\n    auto geometryStart=std::chrono::steady_clock::now();\n    if(blocks.empty())return;')
s=one(s,'    for(std::size_t chunk=0;chunk<chunkCount;++chunk){','    auditprofile::data["pole_geometry_grouping"].seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-geometryStart).count();\n    ++auditprofile::data["pole_geometry_grouping"].calls;\n    for(std::size_t chunk=0;chunk<chunkCount;++chunk){')
s=s.replace('if(pointIt!=frameIt->second.end()){angle=pointIt->second;found=true;}','if(pointIt!=frameIt->second.end()){angle=pointIt->second;found=true;++auditprofile::counters["geometry_hit"];}')
s=s.replace('      if(!found){','      if(!found){++auditprofile::counters["geometry_miss"];')
s=s.replace('      if(!angularTable){','      if(!angularTable){auditprofile::Scope t("wigner_build");++auditprofile::counters["wigner_miss"];')
s=s.replace('if(tableIt!=samplingPlan->angularTables.end())angularTable=tableIt->second;','if(tableIt!=samplingPlan->angularTables.end()){angularTable=tableIt->second;++auditprofile::counters["wigner_hit"];}')
s=s.replace('      if(!angularCoefficients){','      if(!angularCoefficients){auditprofile::Scope t("angular_coefficient_build");++auditprofile::counters["angular_coefficient_miss"];')
s=s.replace('          angularCoefficients=coefficientsIt->second;','          {angularCoefficients=coefficientsIt->second;++auditprofile::counters["angular_coefficient_hit"];}')
s=s.replace('        const auto buildCoefficients=[&](PoleCoefficientRows &coefficients){','        const auto buildCoefficients=[&](PoleCoefficientRows &coefficients){auditprofile::Scope t("radial_coefficients_build");++auditprofile::counters["radial_coefficients_miss"];')
s=s.replace('if(cached!=frequencyCoefficientCache->rows.end())coefficientsPtr=&cached->second;','if(cached!=frequencyCoefficientCache->rows.end()){coefficientsPtr=&cached->second;++auditprofile::counters["radial_coefficients_hit"];}')
s=s.replace('        for(std::size_t j=0;j<group.size();++j){','        {auditprofile::Scope t("degree_sums_phases");\n        for(const auto& [m,row]:coefficients)for(std::size_t l=0;l<row.present.size();++l)if(row.present[l])auditprofile::counters["strain_component_macs"]+=6*group.size();\n        for(std::size_t j=0;j<group.size();++j){')
s=one(s,'      for(std::size_t i=0;i<chunkPoints.size();++i){\n        const auto p=chunkPoints[i];\n        const bool adjoint','      }\n      auditprofile::Scope rotationTimer("rotate_add_scratch");\n      for(std::size_t i=0;i<chunkPoints.size();++i){\n        const auto p=chunkPoints[i];\n        const bool adjoint')
s=one(s,'    if(receiverBegin>points.size())','    auditprofile::Scope t("sample_compact_total");\n    if(receiverBegin>points.size())')
s=one(s,'    std::vector<KernelField3DSample> compact(full.size());','    auditprofile::Scope copyTimer("compact_copy");\n    std::vector<KernelField3DSample> compact(full.size());')
s=one(s,'    std::vector<OwnedModalFrequencyBlock3D> blocks;','    auditprofile::Scope t("prepare_workspace");\n    std::vector<OwnedModalFrequencyBlock3D> blocks;')
s=one(s,'    std::vector<OwnedModalFrequencyBlock3D> radialBlocks;','    auditprofile::Scope t("sample_full_total");\n    std::vector<OwnedModalFrequencyBlock3D> radialBlocks;')
s=one(s,'    auto result=Windowed3DFieldEngine::sampleModalBlocks(','    auto radialStart=std::chrono::steady_clock::now();\n    auto result=Windowed3DFieldEngine::sampleModalBlocks(')
s=one(s,'    if(points.empty())return result;','    auditprofile::data["radial_sampling_and_point_validation"].seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-radialStart).count();\n    ++auditprofile::data["radial_sampling_and_point_validation"].calls;\n    if(points.empty())return result;')
(SHADOW/'Windowed3DFields.h').write_text(s)
s=inc((REPO/'DSpecM1D/src/Windowed3DKernels.h').read_text()).replace(str(REPO/'DSpecM1D/src/Windowed3DFields.h'),'Windowed3DFields.h')
s=one(s,'    const Eigen::MatrixXcd qhat=measurement.adjoint(weightsSI);','    Eigen::MatrixXcd qhat;{auditprofile::Scope t("measurement_adjoint");qhat=measurement.adjoint(weightsSI);}')
s=one(s,'        for(std::size_t p=0;p<fields.size();++p)\n          addPoint(result.pointKernels[begin+p],pointKernel(sem,fields[p]));','        {auditprofile::Scope t("final_kernel_contraction");for(std::size_t p=0;p<fields.size();++p)\n          addPoint(result.pointKernels[begin+p],pointKernel(sem,fields[p]));}')
(SHADOW/'Windowed3DKernels.h').write_text(s)
s=(OLD/'timing_current.cpp').read_text()
s=s.replace('#include <DSpecM1D/Windowed3DKernels>','#include "profile_headers/Windowed3DKernels.h"')
s='#include "profile_headers/SensitivityKernels.h"\n'+s
(ROOT/'profile_current.cpp').write_text(s)
cmd=json.loads((OLD/'phase1/attempt-02/build_commands.json').read_text())['timing']
cmd[cmd.index(str(OLD/'timing_current.cpp'))]=str(ROOT/'profile_current.cpp')
cmd[cmd.index('-o')+1]=str(ROOT/'bin/profile-current')
(ROOT/'profile_build_command.json').write_text(json.dumps(cmd,indent=2)+'\n')
if '--no-build' not in sys.argv:
 with (ROOT/'profile-build.log').open('w') as out:
  subprocess.run(cmd,stdout=out,stderr=out,check=True)
 print('Profile build complete',flush=True)
else: print('Baseline profile sources restored; original binary retained',flush=True)
