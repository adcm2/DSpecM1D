from pathlib import Path
import json,subprocess
ROOT=Path(__file__).resolve().parent
REPO=Path('/home/adcm2/Documents/Research/Projects/DSpecM1D')
YREPO=REPO.parent/'YSpecpp';cache=REPO/'build/phase-d-build-cache/build-yspecpp'
s=(REPO/'benchmarks/yspecpp/kernels/3d/yspec_kernel_timing.cpp').read_text()
s=s.replace('#include <algorithm>','#include <algorithm>\n#include <cstdlib>')
s=s.replace('  const int repeats = std::stoi(argv[6]);','  const int repeats = std::stoi(argv[6]);\n  const int auditBins=std::getenv("SK_BINS")?std::stoi(std::getenv("SK_BINS")):3;\n  if(auditBins<2||auditBins>48)throw std::runtime_error("audit bins must be2..48; original YSpec grid requires at least2 bins");')
s=s.replace('(lmax != 3 && lmax != 12)','(lmax != 3 && lmax != 6 && lmax != 12)')
s=s.replace('  calculation.band = {0.1e-3, 0.8e-3};','  calculation.band = {0.1e-3, auditBins==3?0.8e-3:(auditBins-0.5)/(128.0*20.0)};')
s=s.replace('grid.Last() != 3','grid.Last() != static_cast<std::size_t>(auditBins)')
s=s.replace('  (void)evaluate();\n  const double warmupSeconds = std::chrono::duration<double>(\n      std::chrono::steady_clock::now() - warmupStart).count();','  const double warmupSeconds = 0.0;')
s=s.replace('R"(,"warmup_batches":1,"warmup_api_calls":1,"setup_s":)"','R"(,"warmup_batches":0,"warmup_api_calls":0,"setup_s":)"')
s=s.replace('<< warmupSeconds\n            << R"(,"first_model_to_arrays_s":)" << setupSeconds + warmupSeconds','<< apiSeconds.front()\n            << R"(,"first_model_to_arrays_s":)" << setupSeconds + apiSeconds.front()')
s=s.replace('<< R"(,"solver_bins_zero_based":[1,2,3],"band_hz":[)"','<< R"(,"solver_first_bin":1,"solver_last_bin":)" << grid.Last() << R"(,"band_hz":[)"')
(ROOT/'yspec_audit.cpp').write_text(s)
flags=['/usr/bin/g++','-O3','-DNDEBUG','-std=gnu++23','-fopenmp','-DGSHTRANS_HAVE_INTERPOLATION','-DPLANETMODEL_VERSION="1.2.5"','-DPLANETMODEL_WITH_NETCDF=1',f'-I{YREPO}/src',f'-I{YREPO}']
for dep,sub in [('numericconcepts','include'),('interpolation','include'),('gaussquad',''),('gshtrans',''),('fftwpp',''),('planetmodel','cpp/include')]:flags.append('-I'+str(cache/'_deps'/f'{dep}-src'/sub))
obj=ROOT/'yspec_audit.o';compile=flags+['-c',str(ROOT/'yspec_audit.cpp'),'-o',str(obj)]
link=(cache/'CMakeFiles/yspec_kernel_timing.dir/link.txt').read_text().split()
link[link.index('CMakeFiles/yspec_kernel_timing.dir/yspec_kernel_timing.cpp.o')]=str(obj)
link[link.index('-o')+1]=str(ROOT/'bin/yspec-audit')
(ROOT/'yspec_build_commands.json').write_text(json.dumps({'compile':compile,'link':link,'link_cwd':str(cache)},indent=2)+'\n')
with (ROOT/'yspec-build.log').open('w') as log:
 subprocess.run(compile,stdout=log,stderr=log,check=True)
 subprocess.run(link,cwd=cache,stdout=log,stderr=log,check=True)
print('YSpec audit driver built',flush=True)
