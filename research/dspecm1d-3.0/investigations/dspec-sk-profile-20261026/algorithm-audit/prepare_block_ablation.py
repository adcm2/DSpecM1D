from pathlib import Path
import json,re,subprocess
R=Path(__file__).resolve().parent;repo=Path('/home/adcm2/Documents/Research/Projects/DSpecM1D')
S=R/'block_headers';S.mkdir(exist_ok=True)
for name in ['SensitivityKernels.h','Windowed3DFields.h','Windowed3DKernels.h']:
 src=R/'before/DSpecM1D/src'/name
 text=src.read_text()
 text=re.sub(r'#include "([^"]+)"',lambda m:'#include "'+str(repo/'DSpecM1D/src'/m[1])+'"',text)
 for other in ['SensitivityKernels.h','Windowed3DFields.h']:
  text=text.replace(str(repo/'DSpecM1D/src'/other),other)
 if name=='Windowed3DKernels.h':
  text='#include <cstdlib>\n'+text
  text=text.replace('constexpr std::size_t fieldBlockSize=4096;','const std::size_t fieldBlockSize=std::getenv("SK_BLOCK")?std::stoul(std::getenv("SK_BLOCK")):4096;\n    if(fieldBlockSize==0)throw std::invalid_argument("zero block size");')
 if name=='Windowed3DFields.h':
  text=text.replace('static constexpr std::size_t byteBudget=32u*1024u*1024u;','static inline const std::size_t byteBudget=(std::getenv("SK_PLAN_MB")?std::stoul(std::getenv("SK_PLAN_MB")):32u)*1024u*1024u;',1)
  text='#include <cstdlib>\n'+text
 (S/name).write_text(text)
text=(R.parent/'timing_current.cpp').read_text().replace('#include <DSpecM1D/Windowed3DKernels>','#include "block_headers/Windowed3DKernels.h"')
text='#include "block_headers/SensitivityKernels.h"\n'+text
(R/'block_ablation.cpp').write_text(text)
cmd=json.loads((R/'profile_build_command.json').read_text())
cmd[cmd.index(str(R/'profile_current.cpp'))]=str(R/'block_ablation.cpp');cmd[cmd.index('-o')+1]=str(R/'bin/block-ablation')
(R/'block_build_command.json').write_text(json.dumps(cmd,indent=2)+'\n')
with (R/'block-build.log').open('w') as log:subprocess.run(cmd,stdout=log,stderr=log,check=True)
print('Native block ablation built',flush=True)
