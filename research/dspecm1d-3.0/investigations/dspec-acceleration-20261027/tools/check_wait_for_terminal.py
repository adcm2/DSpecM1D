import json,time,sys,pathlib
state=pathlib.Path(sys.argv[1])/'state.json';output=pathlib.Path(sys.argv[2]);deadline=time.monotonic()+28800
while time.monotonic()<deadline:
 if state.exists():
  data=json.loads(state.read_text())
  if data['status'] in ('completed','failed','resource_limit','timed_out'):
   assert all(j.get('cleanup',{}).get('verified_gone') is True for j in data['jobs'])
   with output.open('x') as stream:json.dump({'dependency_status':data['status'],'dependency_task':str(state.parent),'meaning':'Terminal dependency only; no scientific acceptance inferred.'},stream,indent=2)
   sys.exit(0)
 time.sleep(15)
raise TimeoutError('global queue did not reach verified terminal state within 8 hours')
