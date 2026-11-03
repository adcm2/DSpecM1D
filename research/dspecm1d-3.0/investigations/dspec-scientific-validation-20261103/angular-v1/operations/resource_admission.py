"""Reviewed resource contract and remote functions; no action on import."""
from pathlib import Path
import hashlib
import json

ROOT=Path(__file__).resolve().parents[1]
RECEIPT_SHA="dc59cc0f82a84409e5338c2f7c73f37a683a82e4aaf6b060e2c0ce35f38503e1"
def resource_payload():
    p=ROOT/"resource-coordination/deployment_receipt.json"
    assert hashlib.sha256(p.read_bytes()).hexdigest()==RECEIPT_SHA
    assert json.loads(p.read_text())["status"]=="ADMITTED"
    return {"resource_admission_receipt_sha256":RECEIPT_SHA,
            "resource_helper_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}

REMOTE_RESOURCE_CODE=r'''
import importlib.util as ra_import, hashlib as ra_hash, json as ra_json, os as ra_os
from pathlib import Path as RaPath
import subprocess as ra_sub, time as ra_time, sys as ra_sys
ra_sys.dont_write_bytecode=True
RA_SP=RaPath('/space/adcm2/research/runs/CCDSM/extended_elastic_20261031_v1/sp12_20261101_v1')
RA_SOURCE=RA_SP/'controller/monitor_tunya_global_ram_v4.py'
RA_STATE=RA_SP/'resource_coordination_20261101_v1/global_ram_guard_v4/state.json'
RA_RECEIPT=RA_SP/'sp12_ellipticity_factorial_20261103_v1/resource_coordination_dspecm1d_angular_v1/deployment_receipt.json'
RA_SOURCE_SHA='03421ab7cb49f3e28cfa08c3907ff264cae34d3a4b6f9ea170120b1976989d51'
RA_SLICE='dspecscienceangular20261103.slice'
RA_PREFIX='dspec-science-angular-v1-'
RA_GUARD='ccdsm-node-memory-guard-v4.service'
RA_CAPS={'build':32*(1<<30),'run':16*(1<<30),'check':8*(1<<30)}
def ra_sha(path): return ra_hash.sha256(RaPath(path).read_bytes()).hexdigest()
def ra_props(unit):
    keys=['MainPID','ControlGroup','ActiveState','SubState','FreezerState','MemoryCurrent',
          'MemoryMax','MemorySwapMax','TasksCurrent','TasksMax','CPUQuotaPerSecUSec','CPUUsageNSec']
    args=['systemctl','--user','show',unit]
    for key in keys: args+=['-p',key]
    return dict(line.split('=',1) for line in ra_sub.check_output(args,text=True).splitlines() if '=' in line)
def ra_fresh_guard():
    assert ra_sha(RA_RECEIPT)==CFG['resource_admission_receipt_sha256'],'resource receipt changed'
    receipt=ra_json.loads(RA_RECEIPT.read_text())
    assert receipt['status']=='ADMITTED'
    assert ra_sha(RA_SOURCE)==RA_SOURCE_SHA,'live shared guardian source changed'
    props=ra_props(RA_GUARD)
    assert props['ActiveState']=='active' and int(props['MainPID'])>0,'shared guardian not active'
    cmd=(RaPath('/proc')/props['MainPID']/'cmdline').read_bytes().split(bytes([0]))
    assert str(RA_SOURCE).encode() in cmd,'shared guardian command identity mismatch'
    state=ra_json.loads(RA_STATE.read_text())
    age=ra_time.time()-state['measured_unix']
    assert 0<=age<5,'shared memory observation stale'
    assert state['schema']=='tunya-global-memory-observer-v4'
    assert state['active_ram_metric_valid'] is True and state['status']=='MONITORING'
    assert state['pause_armed'] is True and state['paused_unit'] is None
    assert not state['existing_frozen_owned_units'] and not state['freeze_state_unknown_units']
    assert (state['pause_bytes'],state['urgent_bytes'],state['user_global_boundary_bytes'])==(825000000000,875000000000,900000000000)
    assert state['started_unix']>=receipt['deployed_unix'],'guardian predates registered epoch'
    active=state['active_ram_estimate_bytes']
    assert isinstance(active,int) and active>=0
    # Conservative whole reservation sum, even though some is already active.
    forecast=active+318000000000
    assert forecast<825000000000,'insufficient shared active-RAM reservation headroom'
    return {'state':state,'age_seconds':age,'forecast_active_plus_all_reservations_bytes':forecast,
            'guardian':props,'source_sha256':RA_SOURCE_SHA}
def ra_slice():
    props=ra_props(RA_SLICE)
    assert props['ActiveState']=='active' and props['FreezerState']=='running'
    assert (props['MemoryMax'],props['TasksMax'],props['MemorySwapMax'],props['CPUQuotaPerSecUSec'])==('40000000000','12','0','8s')
    root=RaPath('/sys/fs/cgroup')/props['ControlGroup'].lstrip('/')
    assert root.name==RA_SLICE and root.is_dir()
    assert (root/'cgroup.freeze').read_text().strip()=='0'
    return props,root
def resource_admit():
    guard=ra_fresh_guard()
    props,root=ra_slice()
    live=[]
    for p in root.iterdir():
        if p.is_dir() and (p/'pids.current').exists() and int((p/'pids.current').read_text())>0:
            live.append(p.name)
    assert not live,('another DSpecM1D service is live',live)
    assert int(props['TasksCurrent'])==0 and int(props['MemoryCurrent'])<40000000000
    return {'guard':guard,'slice':props,'live_dspecm1d_services':live}
def resource_verify_live(unit,kind):
    props=ra_props(unit)
    assert props['ActiveState']=='active' and props['FreezerState']=='running'
    assert (props['MemoryMax'],props['TasksMax'],props['MemorySwapMax'],props['CPUQuotaPerSecUSec'])==(str(RA_CAPS[kind]),'12','0','8s')
    slice_props,root=ra_slice()
    cg=RaPath('/sys/fs/cgroup')/props['ControlGroup'].lstrip('/')
    assert cg.parent==root and cg.name==unit,'service escaped registered kernel slice'
    pid=int(props['MainPID']);assert pid>0
    proc=RaPath('/proc')/str(pid)
    fields=(proc/'stat').read_text().split(') ',1)[1].split()
    start=fields[19]
    assert (cg/'cgroup.freeze').read_text().strip()=='0'
    events=dict(line.split() for line in (cg/'memory.events').read_text().splitlines())
    assert all(int(events.get(k,0))==0 for k in ('max','oom','oom_kill','oom_group_kill'))
    assert int((cg/'pids.current').read_text())<=12
    # Import only; derive its research root from the LIVE guardian's app.slice
    # cgroup, not the SSH caller's different session hierarchy.
    guard=ra_fresh_guard()
    guard_cg=RaPath('/sys/fs/cgroup')/guard['guardian']['ControlGroup'].lstrip('/')
    research=guard_cg.parent.parent/'ccdsmextended20261031.slice'
    assert research.parent==root.parent
    spec=ra_import.spec_from_file_location('registered_guard_readonly',RA_SOURCE)
    module=ra_import.module_from_spec(spec);spec.loader.exec_module(module)
    module.research_cgroup=lambda:research
    inventory=module.owned_units()
    record=next(row for row in inventory if row['name']==unit and RaPath(row['path'])==cg)
    assert record['pids']>0 and record['frozen'] is False and record['freeze_error'] is None
    return {'unit':unit,'properties':props,'slice':slice_props,'main_pid':pid,
            'main_pid_starttime':start,'memory_events':events,'guardian_inventory':record,
            'guard':guard}
def resource_launch(command,cwd,log,kind):
    assert kind in RA_CAPS and all(isinstance(v,str) for v in command)
    unit=RA_PREFIX+kind+'.service'
    admission=resource_admit()
    lock=RaPath('/space/adcm2/research/runs/dspec-yspec-kernels/2026-10-27-tunya/dspec-science-angular-v1-resource-active.lock')
    env={'LD_LIBRARY_PATH':'/opt/gcc-13.2.0/lib64','LD_PRELOAD':'','DSPEC_PROFILE_FTZ':'0','PYTHONDONTWRITEBYTECODE':'1'}
    for name in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
                 'BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','GOTO_NUM_THREADS'):
        env[name]='1'
    args=['systemd-run','--user','--quiet','--unit='+unit,'--slice='+RA_SLICE,
          '--working-directory='+str(cwd),'--property=Type=exec',
          '--property=MemoryMax='+str(RA_CAPS[kind]),'--property=MemorySwapMax=0',
          '--property=TasksMax=12','--property=CPUQuota=800%',
          '--property=MemoryAccounting=yes','--property=TasksAccounting=yes',
          '--property=CPUAccounting=yes','--property=StandardInput=null',
          '--property=StandardOutput=append:'+str(log),'--property=StandardError=append:'+str(log),
          '/usr/bin/env']
    args += [name+'='+value for name,value in sorted(env.items())]
    args += ['/usr/bin/flock','--no-fork','--nonblock',str(lock),*command]
    result=ra_sub.run(args,text=True,capture_output=True,check=True)
    live=resource_verify_live(unit,kind)
    ra_time.sleep(.5)
    later=ra_props(unit)
    progress=int(later['CPUUsageNSec'])-int(live['properties']['CPUUsageNSec'])
    assert progress>0 or (RaPath(log).exists() and RaPath(log).stat().st_size>0),'no observed service progress'
    return {'admission':admission,'systemd_command':args,'stdout':result.stdout,'stderr':result.stderr,
            'live':live,'later_properties':later,'cpu_progress_ns':progress}
'''
