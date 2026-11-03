"""Upload fresh archive parts, preserving failed transfer logs and artifacts."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import json, subprocess, time

here = Path(__file__).resolve().parent
manifest = json.loads((here/'manifest.json').read_text())
target = 'tunya:/space/adcm2/research/runs/dspec-yspec-kernels/2026-10-27-tunya/transfer/'
options = ['-B','-F','/tmp/ccdsm-tunya-20261027/ssh_config',
           '-o','ControlMaster=no','-o','ControlPath=none',
           '-o','UpdateHostKeys=no','-o','AddKeysToAgent=no',
           '-o','ProxyCommand=ssh -F /tmp/ccdsm-tunya-20261027/ssh_config -o ControlMaster=no -o ControlPath=none -o UpdateHostKeys=no -o AddKeysToAgent=no -W %h:%p ssh-west']
def upload(name):
    started = time.time()
    with (here/(name+'.upload.log')).open('xb') as log:
        result = subprocess.run(['scp',*options,str(here/name),target+name],stdout=log,stderr=log)
    record = {'name':name,'returncode':result.returncode,'started_unix':started,'finished_unix':time.time()}
    (here/(name+'.upload.json')).write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record),flush=True)
    return result.returncode

for name in ['manifest.json','directories.json']:
    if upload(name):raise SystemExit('metadata upload failed')
failed=[]
with ThreadPoolExecutor(max_workers=4) as pool:
    for result in as_completed([pool.submit(upload,part['name']) for part in manifest['parts']]):
        if result.result():failed.append(True)
if failed:raise SystemExit('part upload failed; preserve all artifacts and retry separately')
(here/'UPLOADED').write_text('All parts and metadata uploaded successfully.\n')
