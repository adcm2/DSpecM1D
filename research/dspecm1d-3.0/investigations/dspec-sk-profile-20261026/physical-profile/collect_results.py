"""Collect saved physical cases; never interpret a partial output as completed."""
from pathlib import Path
import argparse
import json
import numpy as np

p = argparse.ArgumentParser()
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()
root = Path(__file__).resolve().parent
runs = root / 'runs'

def trace(name):
    return np.genfromtxt(runs / name / 'data/trace.csv', delimiter=',', names=True)

def compare(candidate, reference):
    x, y = trace(candidate), trace(reference)
    m = (y['time_s'] >= 1000) & (y['time_s'] <= 1200)
    t, w = y['time_s'][m], y['window_weight'][m]
    xx, yy = np.interp(t, x['time_s'], x['vertical_m']), y['vertical_m'][m]
    ex = .5 * np.diff(x['time_s'])[0] * np.sum(x['window_weight'] * x['vertical_m']**2)
    ey = .5 * np.diff(y['time_s'])[0] * np.sum(y['window_weight'] * y['vertical_m']**2)
    dw = np.linalg.norm(np.sqrt(w) * (xx-yy)) / np.linalg.norm(np.sqrt(w) * yy)
    de = abs(ex / ey - 1)
    return dict(candidate=candidate, reference=reference, energy_m2_s=float(ex),
                reference_energy_m2_s=float(ey), window_weighted_relative_L2=float(dw),
                energy_relative=float(de), passed_1e_3=bool(max(dw, de) <= 1e-3))

complete, pending = {}, []
for r in sorted(runs.iterdir()):
    if not (r / 'run.json').exists():
        pending.append(r.name)
        continue
    run = json.loads((r / 'run.json').read_text())
    result = None
    if (r / 'data/result.json').exists():
        try:
            result = json.loads((r / 'data/result.json').read_text(), parse_constant=lambda x: None)
        except json.JSONDecodeError:
            pass
    complete[r.name] = dict(returncode=run['returncode'], wall_s=run['wall_s'],
                            peak_rss_bytes=run['peak_rss_bytes'], result=result,
                            profiles=run.get('profiles', []),
                            child_user_cpu_s=run.get('child_user_cpu_s'),
                            child_system_cpu_s=run.get('child_system_cpu_s'))
ys = [(l, f'yspec-degree-l{l}-t7-q1-{v}') for l, v in
      [(64,'v2'),(96,'v2'),(128,'v2'),(192,'v2'),(256,'v2'),(384,'v2'),
       (512,'v3'),(768,'v3'),(1024,'v3'),(1280,'v3'),(1536,'v3'),(1792,'v4'),(1827,'v4')]]
ys = [(l,n) for l,n in ys if n in complete and complete[n]['returncode'] == 0]
adjacent = [dict(L_coarse=l1, L_fine=l2, **compare(n2,n1))
            for (l1,n1),(l2,n2) in zip(ys,ys[1:])]
ref = 'dspec-frozen-v3-l1536-nq6-q1'
refinements = [compare(n,ref) for n in complete if n.startswith('dspec-v4-l1536-')
               and complete[n]['returncode'] == 0]
cross = []
for ds,ysn in [('dspec-frozen-v3-l384-nq6-q1','yspec-degree-l384-t7-q1-v2'),
               ('dspec-frozen-v3-l1024-nq6-q1','yspec-degree-l1024-t7-q1-v3'),
               ('dspec-frozen-v3-l1536-nq6-q1','yspec-degree-l1536-t7-q1-v3'),
               ('dspec-v4-l1827-nq8-q1','yspec-degree-l1827-t7-q1-v4')]:
    if ds in complete and ysn in complete and complete[ds]['returncode'] == complete[ysn]['returncode'] == 0:
        cross.append(compare(ds,ysn))
out = dict(protocol=json.loads((root/'protocol.json').read_text()),
           interpretation='Completed runs only. Accuracy runs used distinct CPUs concurrently and are not controlled speed comparisons. No full high-degree global kernel is yet materialized.',
           pending_or_incomplete=pending, cases=complete,
           degree_adjacent=adjacent, radial_padding_dc_checks=refinements, cross_code=cross)
a.output.write_text(json.dumps(out, indent=2, allow_nan=False) + '\n')
print(json.dumps(dict(completed=len(complete),pending=pending,degree_last=adjacent[-1] if adjacent else None),indent=2))
