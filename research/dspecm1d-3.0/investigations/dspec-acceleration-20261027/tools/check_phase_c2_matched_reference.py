#!/usr/bin/env python3
"""Compare one fresh accepted DSpec native L200 run with C2 and saved F8a."""
import argparse, csv, hashlib, json, math, pathlib, struct


def need(ok, text):
    if not ok: raise RuntimeError(text)


def read_json(p): return json.loads(pathlib.Path(p).read_text())
def sha(p): return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def read_csv(p):
    with pathlib.Path(p).open(newline="") as f: return list(csv.DictReader(f))
def rel_l2(a,b):
    need(len(a)==len(b),"array length mismatch")
    err=norm=0.0
    for x,y in zip(a,b):
        need(math.isfinite(x) and math.isfinite(y),"nonfinite value")
        err+=(x-y)**2; norm+=y*y
    return math.sqrt(err/norm) if norm else (0.0 if err==0 else None)
def doubles(p):
    raw=pathlib.Path(p).read_bytes();need(len(raw)%8==0,"invalid kernel binary")
    vals=struct.unpack('<'+'d'*(len(raw)//8),raw)
    need(all(math.isfinite(x) for x in vals),"nonfinite kernel binary")
    return vals
def csv_metrics(a_path,b_path,label):
    a,b=read_csv(a_path),read_csv(b_path)
    need(len(a)==len(b),f"{label} row count mismatch")
    with pathlib.Path(a_path).open(newline="") as f: ah=csv.DictReader(f).fieldnames
    with pathlib.Path(b_path).open(newline="") as f: bh=csv.DictReader(f).fieldnames
    need(ah==bh,f"{label} columns mismatch")
    exact={"bin","time_s","window","frequency_hz","sample"}; metrics={}
    for col in ah:
        av=[r[col] for r in a]; bv=[r[col] for r in b]
        if col in exact:
            need(av==bv,f"{label} exact grid column {col} differs"); continue
        try: x=[float(v) for v in av]; y=[float(v) for v in bv]
        except ValueError:
            need(av==bv,f"{label} metadata column {col} differs"); continue
        e=rel_l2(x,y);need(e is not None and e<=1e-11,f"{label} channel {col} differs: {e}");metrics[col]=e
    return metrics
def verify_outputs(state,root,job_id,names):
    job=next((x for x in state["jobs"] if x["id"]==job_id),None);need(job is not None,f"missing job {job_id}")
    need(job["status"]=="completed" and job["returncode"]==0 and job.get("inputs_unchanged") is True and job["cleanup"]["verified_gone"],f"job receipt failed: {job_id}")
    for name in names: need(sha(root/job_id/name)==job["outputs_sha256"].get(name),f"output hash mismatch: {job_id}/{name}")
    return job
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--task-root',type=pathlib.Path,required=True);ap.add_argument('--candidate',type=pathlib.Path,required=True);ap.add_argument('--candidate-state',type=pathlib.Path,required=True);ap.add_argument('--saved-reference',type=pathlib.Path,required=True);ap.add_argument('--saved-state',type=pathlib.Path,required=True);ap.add_argument('--out',type=pathlib.Path,required=True);a=ap.parse_args()
    new_state=read_json(a.task_root/'state.json');saved_state=read_json(a.saved_state);candidate_state=read_json(a.candidate_state)
    ref_id='accepted-reference-native-l200';candidate_id='phase-c2-l200-native-skip-on';saved_id='F8a-native-L200'
    numeric_names=['data/parameters.txt','data/result.json','data/trace.csv','data/spectrum.csv','data/components.csv','data/nodes.csv','data/along.bin','data/along-angles.csv','data/across.bin','data/across-angles.csv']
    fresh=verify_outputs(new_state,a.task_root,ref_id,numeric_names)
    cand=next((j for j in candidate_state['jobs'] if j['id']==candidate_id),None); need(cand is not None,"missing C2 native job receipt")
    saved=next((j for j in saved_state['jobs'] if j['id']==saved_id),None); need(saved is not None,"missing saved F8a job receipt")
    need(cand['status']=='completed' and cand['returncode']==0 and cand.get('inputs_unchanged') is True and cand['cleanup']['verified_gone'],"C2 native job receipt failed")
    need(saved['status']=='completed' and saved['returncode']==0 and saved.get('inputs_unchanged') is True and saved['cleanup']['verified_gone'],"saved F8a job receipt failed")
    # Saved and C2 artifacts are direct checker inputs; their owning supervisor receipts bind them.
    for name in ('result.json','trace.csv','spectrum.csv','components.csv','nodes.csv','along.bin','along-angles.csv','across.bin','across-angles.csv'):
        need(sha(a.candidate/name)==cand['outputs_sha256'].get('data/'+name),f"C2 output hash mismatch: {name}")
        need(sha(a.saved_reference/name)==saved['outputs_sha256'].get('data/'+name),f"saved F8a output hash mismatch: {name}")
    fresh_dir=a.task_root/ref_id/'data'; c2=a.candidate; old=a.saved_reference
    results={"fresh_accepted_reference":read_json(fresh_dir/'result.json'),"C2_native":read_json(c2/'result.json'),"saved_F8a":read_json(old/'result.json')}
    for label,d in results.items():
        for k,v in {"degree":200,"nq":10,"refinement":0.025,"frequency_count":35,"first_bin":7,"last_bin":41,"native_side_count":293,"along_angles":135,"across_angles":19}.items(): need(d.get(k)==v,f"{label} {k} mismatch")
    need(results['C2_native'].get('phase_c2_preweighted_complex_contraction') is True,"C2 native path flag missing")
    for name in ('nodes.csv','along-angles.csv','across-angles.csv'):
        need(sha(fresh_dir/name)==sha(c2/name)==sha(old/name),f"radial/angular geometry mismatch: {name}")
    waveform={}
    for name in ('trace.csv','spectrum.csv','components.csv'):
        waveform[name]={"C2_native_vs_fresh_reference":csv_metrics(c2/name,fresh_dir/name,'C2/fresh '+name),"fresh_reference_vs_saved_F8a":csv_metrics(fresh_dir/name,old/name,'fresh/saved '+name)}
    kernels={}
    for name in ('along.bin','across.bin'):
        vals=doubles(fresh_dir/name); c2vals=doubles(c2/name); oldvals=doubles(old/name)
        ec=rel_l2(vals,c2vals); eo=rel_l2(vals,oldvals)
        need(ec is not None and ec<=1e-11 and eo is not None and eo<=1e-11,f"fresh reference {name} differs: C2={ec}, F8a={eo}")
        kernels[name]={"C2_native_relative_l2":ec,"saved_F8a_relative_l2":eo,"fresh_sha256":sha(fresh_dir/name)}
    norm={"C2_native_relative_error":abs(results['C2_native']['measurement_rate_norm']-results['fresh_accepted_reference']['measurement_rate_norm'])/max(abs(results['fresh_accepted_reference']['measurement_rate_norm']),1e-300),"saved_F8a_relative_error":abs(results['saved_F8a']['measurement_rate_norm']-results['fresh_accepted_reference']['measurement_rate_norm'])/max(abs(results['fresh_accepted_reference']['measurement_rate_norm']),1e-300)}
    need(max(norm.values())<=1e-11,f"measurement rate norm mismatch: {norm}")
    timing={"C2_native_driver_elapsed_s":results['C2_native'].get('elapsed_s'),"C2_native_supervisor_elapsed_s":cand.get('elapsed_seconds'),"fresh_reference_driver_elapsed_s":results['fresh_accepted_reference'].get('elapsed_s'),"fresh_reference_supervisor_elapsed_s":fresh.get('elapsed_seconds'),"saved_F8a_driver_elapsed_s":results['saved_F8a'].get('elapsed_s'),"saved_F8a_supervisor_elapsed_s":saved.get('elapsed_seconds'),"historical_stage2_native_median_s":16.8176,"interpretation":"single fresh pair; report driver and supervisor boundaries separately; descriptive timing, not a repeated-run distribution"}
    report={"status":"pass","case":{"degree":200,"nq":10,"refinement":0.025,"bins":35,"frequency_workers":8,"inner_workers":1},"waveform_channel_relative_l2":waveform,"kernel_relative_l2":kernels,"measurement_rate_norm_relative_error":norm,"timing":timing,"reference_provenance":{"fresh_binary_sha256":fresh['inputs_sha256_before'].get('/space/adcm2/research/builds/stages12-stage2-v1/doc-native/doc-driver'),"saved_F8a_binary_sha256":saved['inputs_sha256_before'].get('/space/adcm2/research/builds/stages12-stage2-v1/doc-native/doc-driver'),"fresh_state_job_sha256":new_state['jobs'][0]['command_sha256'],"candidate_C2_job_sha256":cand['command_sha256']}}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');print(json.dumps({"status":"pass","out":str(a.out),"timing":timing}))
if __name__=='__main__':main()
