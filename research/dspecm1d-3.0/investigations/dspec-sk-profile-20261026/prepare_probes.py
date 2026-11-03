from pathlib import Path
import json, hashlib, subprocess
root=Path(__file__).resolve().parent
repo=Path('/home/adcm2/Documents/Research/Projects/DSpecM1D')
yrepo=repo.parent/'YSpecpp'
shadow=root/'instrumented_headers'; shadow.mkdir(exist_ok=True)
(root/'bin').mkdir(exist_ok=True)
# Timers aggregate wall time and call counts; nested labels are not additive.
(shadow/'profile.hpp').write_text('''#pragma once
#include <chrono>
#include <array>
#include <iostream>
namespace skprofile {
struct Entry { double seconds=0; unsigned long long calls=0; };
inline std::array<Entry,20> entries{};
inline const char* names[]={"api", "solve_modes", "prepare_workspace", "sample_fields", "pole_synthesis", "geometry_grouping", "angular_tables", "interpolation", "point_accumulation", "system_matrix", "factorization", "forward_solve", "adjoint_solve", "point_contraction", "measurement_map"};
struct Scope { int id; std::chrono::steady_clock::time_point start; Scope(int i):id(i),start(std::chrono::steady_clock::now()){} ~Scope(){entries[id].seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();++entries[id].calls;} };
struct Dump { ~Dump(){std::cerr<<"SKPROFILE {"; for(int i=0;i<15;++i)std::cerr<<(i?",":"")<<"\\\""<<names[i]<<"\\\":{\\\"seconds\\\":"<<entries[i].seconds<<",\\\"calls\\\":"<<entries[i].calls<<"}"; std::cerr<<"}\\n";} };
inline Dump dump;
}
''')
def original(name): return (repo/'DSpecM1D/src'/name).read_text()
def includes(t):
    import re
    return re.sub(r'#include "([^"]+)"',lambda m: '#include "'+str(repo/'DSpecM1D/src'/m[1])+'"',t)
def replace(t,old,new):
    if t.count(old)!=1: raise ValueError((old,t.count(old)))
    return t.replace(old,new)
s=includes(original('SensitivityKernels.h'))
s='#include "profile.hpp"\n'+s
s=s.replace('    SparseComplex h = systemMatrix(sem, family, degree, z);','    SparseComplex h = systemMatrix(sem, family, degree, z);')
s=s.replace('    solver.compute(h);','    { skprofile::Scope timer(10); solver.compute(h); }')
s=s.replace('    Eigen::MatrixXcd fields = solver.solve(compactRhs);','    Eigen::MatrixXcd fields; { skprofile::Scope timer(11); fields = solver.solve(compactRhs); }')
s=s.replace('    Eigen::MatrixXcd fields = solver.solveAdjoint(rightHandSides);','    Eigen::MatrixXcd fields; { skprofile::Scope timer(12); fields = solver.solveAdjoint(rightHandSides); }')
s=s.replace('    fields = solver.solveAdjoint(rightHandSides);','    { skprofile::Scope timer(12); fields = solver.solveAdjoint(rightHandSides); }')
s=s.replace('    basisFields = solver.solveAdjoint(basisRightHandSides);','    { skprofile::Scope timer(12); basisFields = solver.solveAdjoint(basisRightHandSides); }')
s=replace(s,'                                    Complex z) {\n','                                    Complex z) {\n    skprofile::Scope timer(9);\n')
(shadow/'SensitivityKernels.h').write_text(s)
s=includes(original('Windowed3DFields.h')).replace(str(repo/'DSpecM1D/src/SensitivityKernels.h'),'SensitivityKernels.h')
s=replace(s,'      int element, Value value) {','      int element, Value value) {') if False else s
s=replace(s,'      SensitivitySolveFamily family, double radius, int element, Value value) {','      SensitivitySolveFamily family, double radius, int element, Value value) {\n    skprofile::Scope timer(7);')
s=replace(s,'      PoleFrame *sourceFrame=nullptr) {','      PoleFrame *sourceFrame=nullptr) {\n    skprofile::Scope timer(1);')
s=replace(s,'    if(blocks.empty())return;\n    struct PointAngle','    skprofile::Scope timer(4);\n    if(blocks.empty())return;\n    auto geometryStart=std::chrono::steady_clock::now();\n    struct PointAngle')
s=replace(s,'    for(std::size_t chunk=0;chunk<chunkCount;++chunk){','    skprofile::entries[5].seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-geometryStart).count(); ++skprofile::entries[5].calls;\n    for(std::size_t chunk=0;chunk<chunkCount;++chunk){')
s=replace(s,'      AngularWigner w(maxDegree,std::min(2,maxDegree),std::min(2,maxDegree),chunks[chunk].theta);','      auto tableStart=std::chrono::steady_clock::now();\n      AngularWigner w(maxDegree,std::min(2,maxDegree),std::min(2,maxDegree),chunks[chunk].theta);\n      skprofile::entries[6].seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-tableStart).count(); ++skprofile::entries[6].calls;')
s=replace(s,'            for(const auto p:group){\n              const std::size_t thetaIndex=angleForPoint[p]-begin;','            skprofile::Scope pointTimer(8);\n            for(const auto p:group){\n              const std::size_t thetaIndex=angleForPoint[p]-begin;')
s=replace(s,'      const Eigen::VectorXcd &adjointWeights) {\n    sensitivity_detail::requireConstructedSEM(sem);','      const Eigen::VectorXcd &adjointWeights) {\n    skprofile::Scope timer(2);\n    sensitivity_detail::requireConstructedSEM(sem);')
s=replace(s,'      const std::vector<Field3DPoint> &points,bool allowCenter) const {','      const std::vector<Field3DPoint> &points,bool allowCenter) const {\n    skprofile::Scope timer(3);')
(shadow/'Windowed3DFields.h').write_text(s)
s=includes(original('Windowed3DKernels.h')).replace(str(repo/'DSpecM1D/src/Windowed3DFields.h'),'Windowed3DFields.h')
s=replace(s,'      const std::vector<std::array<double,2>> &angles = {}) {','      const std::vector<std::array<double,2>> &angles = {}) {\n    skprofile::Scope timer(0);')
s=replace(s,'      for(std::size_t p=0;p<points.size();++p){\n        const auto kernel=pointKernel','      { skprofile::Scope contractionTimer(13);\n      for(std::size_t p=0;p<points.size();++p){\n        const auto kernel=pointKernel')
s=replace(s,'      for(std::size_t p=0;p<fields.sourceContact.size();++p){','      }\n      for(std::size_t p=0;p<fields.sourceContact.size();++p){')
(shadow/'Windowed3DKernels.h').write_text(s)
# Current live source, remove the warm-up for a clean first-call comparison.
s=(repo/'tutorials/sensitivity_kernels_3d/kernel_timing.cpp').read_text()
a=s.index('  (void)evaluate(0);'); b=s.index('  for (int repeat',a)
s=s[:a]+'  const double warmupSeconds = 0.0;\n'+s[b:]
s=s.replace('"warmup_batches":1,"warmup_api_calls":','"warmup_batches":0,"warmup_api_calls":').replace('<< measurementCount\n            << R"(,"setup_s":)','<< 0\n            << R"(,"setup_s":)')
s=s.replace('<< R"(],"first_workload_api_s":)" << warmupSeconds','<< R"(],"first_workload_api_s":)" << apiSeconds.front()').replace('<< R"(,"first_model_to_arrays_s":)" << setupSeconds + warmupSeconds','<< R"(,"first_model_to_arrays_s":)" << setupSeconds + apiSeconds.front()')
(root/'timing_current.cpp').write_text(s)
(root/'timing_instrumented.cpp').write_text('#include "instrumented_headers/SensitivityKernels.h"\n'+s.replace('#include <DSpecM1D/Windowed3DKernels>','#include "instrumented_headers/Windowed3DKernels.h"'))
cmd=json.loads((repo/'build/kernel-reconstruction-pole-frames-20261026/timing_command.json').read_text())
commands={}
for kind in ['current','instrumented']:
    c=cmd.copy(); c[c.index('/tmp/dspec-angular-reuse-20261026/timing_probe.cpp')]=str(root/f'timing_{kind}.cpp'); c[c.index('-o')+1]=str(root/'bin'/kind)
    commands[kind]=c
(root/'build_commands.json').write_text(json.dumps(commands,indent=2)+'\n')
provenance={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [repo/'DSpecM1D/src/Windowed3DFields.h',repo/'DSpecM1D/src/SensitivityKernels.h',repo/'DSpecM1D/src/Windowed3DKernels.h',yrepo/'src/SampledKernels.cpp']}
(root/'source_hashes.json').write_text(json.dumps(provenance,indent=2)+'\n')
for kind,c in commands.items():
    with (root/f'build_{kind}.log').open('w') as out:
        p=subprocess.run(c,stdout=out,stderr=out)
    print(kind,p.returncode,flush=True)
    if p.returncode: raise RuntimeError((root/f'build_{kind}.log').read_text()[-5000:])
