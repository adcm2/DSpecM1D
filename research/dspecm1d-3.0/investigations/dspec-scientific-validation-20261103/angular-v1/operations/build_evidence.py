"""Read-only retrieval and independent byte checks of the angular build."""
from pathlib import Path
import hashlib
import importlib.util
import json
import time

ROOT=Path(__file__).resolve().parents[1]
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
def main():
    operator=load("build_operator",ROOT/"operations/build_stage.py")
    resource=load("resource",ROOT/"operations/resource_admission.py")
    transport=load("transport",ROOT.parents[1]/"tunya-profile-20261027/ssh_transport.py")
    cfg=operator.common_payload()
    cfg["build_config"]=json.loads((ROOT/"configs/angular-build-v1.json").read_text())
    cfg["build_plan_sha256"]=hashlib.sha256((ROOT/"configs/angular-build-v1-plan.json").read_bytes()).hexdigest()
    source="CFG="+repr(cfg)+"\n"+resource.REMOTE_RESOURCE_CODE+r'''
import json,hashlib
from pathlib import Path
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
task=Path(CFG['control'])/'continuation'/CFG['task']
state=json.loads((task/'state.json').read_text())
guard=ra_fresh_guard()
slice_props,slice_root=ra_slice()
unit=ra_props('dspec-science-angular-v1-build.service')
out={'status':state['status'],'supervisor_state':state,'resource_guard':guard,
     'unit_properties':unit,'slice_properties':slice_props,
     'slice_memory_peak_bytes':int((slice_root/'memory.peak').read_text()),
     'slice_memory_events':dict(line.split() for line in (slice_root/'memory.events').read_text().splitlines()),
     'remote_mutations':False}
if state['status']=='completed':
    assert state['plan_sha256']==CFG['build_plan_sha256']
    assert len(state['jobs'])==1
    job=state['jobs'][0]
    assert job['id']==CFG['task'] and job['status']=='completed' and job['returncode']==0
    assert not job['failure_reasons'] and job['inputs_unchanged'] and job['cleanup']['verified_gone']
    for relative in job['outputs']:
        assert sha(task/job['id']/relative)==job['outputs_sha256'][relative]
    assert sha(task/job['id']/job['log'])==job['log_sha256']
    build=Path(CFG['build'])
    campaign=json.loads((build/'build_campaign.json').read_text())
    assert campaign['status']=='completed' and len(campaign['variants'])==1
    variant=CFG['build_config']['variants'][0]
    folder=build/variant['id'];binary=folder/variant['output_name']
    record=json.loads((folder/'command.json').read_text())
    expected=[arg.replace('{OUTPUT}',str(binary)) for arg in variant['command']]
    assert record['command']==expected and record['cwd']==variant['cwd'] and record['returncode']==0
    assert record['tracked_inputs_unchanged'] is True
    assert record['tracked_input_sha256']==record['tracked_input_sha256_after']
    assert record['binary_sha256']==sha(binary)==campaign['variants'][0]['binary_sha256']
    spec=ra_import.spec_from_file_location('frozen_build_tool',CFG['build_tool'])
    tool=ra_import.module_from_spec(spec);spec.loader.exec_module(tool)
    for p,h in record['tracked_input_sha256'].items():assert tool.tree_hash(Path(p))==h,p
    for p,h in CFG['dependencies'].items():assert tool.tree_hash(Path(p))==h,p
    assert sha(folder/'stdout.log')==record['stdout_sha256']
    assert sha(folder/'stderr.log')==record['stderr_sha256']
    record.update(binary_sha256_observed=sha(binary),command_json_sha256=sha(folder/'command.json'),
                  stderr_text=(folder/'stderr.log').read_text())
    out.update(build_campaign=campaign,build_result=json.loads((task/job['id']/'build-result.json').read_text()),
               command_receipts={variant['id']:record},external_dependencies=CFG['dependencies'])
print(json.dumps(out))
'''
    compile(source,"readonly-build-evidence","exec")
    result=json.loads(transport.python(source))
    dest=ROOT/"receipts"/("angular-build-v1-evidence.json" if result["status"]=="completed"
                          else f"angular-build-v1-status-{time.time_ns()}.json")
    with dest.open("x") as f:json.dump(result,f,indent=2);f.write("\n")
    print(json.dumps({"status":result["status"],"receipt":str(dest),
        "sha256":hashlib.sha256(dest.read_bytes()).hexdigest(),
        "slice_peak_bytes":result["slice_memory_peak_bytes"],
        "jobs":[{k:j.get(k) for k in ("id","status","returncode","elapsed_seconds","failure_reasons")}
                for j in result["supervisor_state"]["jobs"]]},indent=2))
if __name__=="__main__":main()
