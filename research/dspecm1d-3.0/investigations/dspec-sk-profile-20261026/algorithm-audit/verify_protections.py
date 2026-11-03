from pathlib import Path
import hashlib,json,subprocess
R=Path(__file__).resolve().parent;b=json.loads((R/'before.json').read_text());repo=Path(b['repo'])
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
tracked={str(p):sha(repo/p)==v for p,v in b['tracked_hashes'].items()};phase0={str(p):sha(repo/p)==v for p,v in b['phase0_hashes'].items()}
head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip();status=subprocess.check_output(['git','status','--short'],cwd=repo,text=True)
y=Path('/home/adcm2/Documents/Research/Projects/YSpecpp');ystatus=subprocess.check_output(['git','status','--short'],cwd=y,text=True);yhead=subprocess.check_output(['git','rev-parse','HEAD'],cwd=y,text=True).strip()
r={'tracked_count':len(tracked),'tracked_unchanged':all(tracked.values()),'changed_tracked':[p for p,x in tracked.items() if not x],'phase0_count':len(phase0),'phase0_unchanged':all(phase0.values()),'head_unchanged':head==b['head'],'status_unchanged':status==b['status'],'dspec_head':head,'dspec_status':status,'yspec_head':yhead,'yspec_clean':not ystatus,'yspec_expected_head':yhead=='280eab0200495391539c0f97e62f978083c307bc'}
(R/'final_protections.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2));assert all(r[k] for k in ['tracked_unchanged','phase0_unchanged','head_unchanged','status_unchanged','yspec_clean','yspec_expected_head'])
