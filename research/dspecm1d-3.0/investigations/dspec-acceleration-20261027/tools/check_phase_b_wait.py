import json,time,sys,pathlib,hashlib
root=pathlib.Path(sys.argv[1]);output=pathlib.Path(sys.argv[2]);deadline=time.monotonic()+36000
while time.monotonic()<deadline:
 path=root/'state.json'
 if path.exists():
  raw=path.read_bytes();state=json.loads(raw)
  if state['status']!='running':
   assert state['status']=='completed',state['status']
   assert all(j['status']=='completed' and j['returncode']==0 and j['inputs_unchanged'] and j['cleanup']['verified_gone'] for j in state['jobs'])
   with output.open('x') as f:json.dump({'dependency_status':'completed','dependency_task':str(root),'state_sha256':hashlib.sha256(raw).hexdigest()},f,indent=2)
   sys.exit(0)
 time.sleep(15)
raise TimeoutError('Phase B numerical queue did not finish within10hours')
