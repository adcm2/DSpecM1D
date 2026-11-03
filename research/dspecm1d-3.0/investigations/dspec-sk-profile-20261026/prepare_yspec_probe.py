from pathlib import Path
import json,subprocess
root=Path(__file__).resolve().parent
repo=Path('/home/adcm2/Documents/Research/Projects/DSpecM1D'); yrepo=repo.parent/'YSpecpp'
cache=repo/'build/phase-d-build-cache/build-yspecpp'
def replace(s,a,b):
    if s.count(a)!=1:raise ValueError((a,s.count(a)))
    return s.replace(a,b)
s=(yrepo/'src/SampledKernels.cpp').read_text()
s='#include "instrumented_headers/profile.hpp"\n'+s
s=replace(s,'    const Progress& progress) {\n  const auto& scales', '    const Progress& progress) {\n  skprofile::Scope apiTimer(0);\n  const auto& scales')
s=replace(s,'          for (auto l = degrees.min; l <= lMax; ++l) {','          { skprofile::Scope radialTimer(1);\n          for (auto l = degrees.min; l <= lMax; ++l) {')
s=replace(s,'          // The factors of each measurement\'s adjoint field:', '          }\n          // The factors of each measurement\'s adjoint field:')
s=replace(s,'            const auto sourceTable =','            auto tableStart=std::chrono::steady_clock::now();\n            const auto sourceTable =')
s=replace(s,'            // [measurement][node][parameter][place]', '            skprofile::entries[6].seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-tableStart).count(); ++skprofile::entries[6].calls;\n            // [measurement][node][parameter][place]')
s=s.replace('            for (auto n : std::views::iota(std::size_t{0}, nodes.size())) {','            { skprofile::Scope contractionTimer(8);\n            for (auto n : std::views::iota(std::size_t{0}, nodes.size())) {',1)
s=replace(s,'            const auto lock = std::scoped_lock(','            }\n            const auto lock = std::scoped_lock(')
s=replace(s,'  auto names = Names(set);','  skprofile::Scope convertTimer(13);\n  auto names = Names(set);')
(root/'yspec_sampled_instrumented.cpp').write_text(s)
# Recompile the live sampled-kernel translation unit and link it ahead of the archive.
flags=['/usr/bin/g++','-O3','-DNDEBUG','-std=gnu++23','-fopenmp','-DGSHTRANS_HAVE_INTERPOLATION','-DPLANETMODEL_VERSION="1.2.5"','-DPLANETMODEL_WITH_NETCDF=1',f'-I{yrepo}/src',f'-I{yrepo}']
for dep,sub in [('numericconcepts','include'),('interpolation','include'),('gaussquad',''),('gshtrans',''),('fftwpp',''),('planetmodel','cpp/include')]:
    flags.append('-I'+str(cache/'_deps'/f'{dep}-src'/sub))
obj=root/'yspec_sampled_instrumented.o'
compile=flags+['-c',str(root/'yspec_sampled_instrumented.cpp'),'-o',str(obj)]
link=(cache/'CMakeFiles/yspec_kernel_timing.dir/link.txt').read_text().split()
link[3:3]=[str(obj)]
link[link.index('-o')+1]=str(root/'bin/yspec_instrumented')
# All relative build paths resolve in the existing build cache; no cache file is changed.
(root/'yspec_build_commands.json').write_text(json.dumps({'compile':compile,'link':link,'link_cwd':str(cache)},indent=2)+'\n')
with (root/'build_yspec_instrumented.log').open('w') as out:
    a=subprocess.run(compile,stdout=out,stderr=out)
    if a.returncode: raise RuntimeError((root/'build_yspec_instrumented.log').read_text()[-5000:])
    b=subprocess.run(link,cwd=cache,stdout=out,stderr=out)
    if b.returncode:raise RuntimeError((root/'build_yspec_instrumented.log').read_text()[-5000:])
print('yspec instrumented built',flush=True)
