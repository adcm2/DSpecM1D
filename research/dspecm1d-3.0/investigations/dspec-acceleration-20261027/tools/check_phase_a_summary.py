#!/usr/bin/env python3
"""Read saved receipts and preserve complete profiler components without rerunning solvers."""
import argparse,hashlib,json
from pathlib import Path
import check_pair

def summarize(root):
    state=check_pair.load(root/'state.json')
    out={'source_snapshot':str(root.resolve()),'state_sha256':check_pair.sha(root/'state.json'),'status':state['status'],'jobs':{}}
    for job in state['jobs']:
        row={'status':job['status'],'returncode':job.get('returncode'),'process_wall_s':job.get('elapsed_seconds'),'peak_threads':job.get('peak_observed_threads'),'peak_rss_gib_userspace':job.get('peak_observed_rss_gib_userspace'),'inputs_unchanged':job.get('inputs_unchanged'),'cleanup_verified':job.get('cleanup',{}).get('verified_gone'),'command':job['command']}
        folder=root/job['id']
        for rel,digest in {**job.get('outputs_sha256',{}),job['log']:job.get('log_sha256')}.items():
            path=folder/rel
            if path.exists() and digest and check_pair.sha(path)!=digest:raise ValueError('local snapshot hash mismatch: '+str(path))
        result=folder/'data/result.json'
        if result.exists():
            row['result']=check_pair.load(result)
            profile=check_pair.profile(folder/job['log'])
            if profile is not None:row['profile']=profile
        checks=folder/'checks.json'
        if checks.exists():row['checks']=check_pair.load(checks)
        manifest=folder/'data/kernel.bin.json'
        if manifest.exists():
            try:row['global_manifest']=check_pair.load(manifest)
            except json.JSONDecodeError:row['global_manifest_parse_failed']=True
        out['jobs'][job['id']]=row
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('snapshot');p.add_argument('output');a=p.parse_args()
    data=summarize(Path(a.snapshot))
    with Path(a.output).open('x') as stream:json.dump(data,stream,indent=2,allow_nan=False);stream.write('\n')
