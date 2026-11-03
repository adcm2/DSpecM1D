#!/usr/bin/env python3
"""Read-only comparison of two completed, supervisor-receipted DSpec jobs."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

LIMIT = 2e-13
PHYSICAL = (
    'degree', 'sem_degree', 'radial_refinement_maxstep', 'qex', 'dt_s',
    'duration_s', 'fft_samples', 'fft_period_s', 'frequency_bins',
    'epsilon_si_s_inverse', 'df_si_hz', 'source', 'receiver', 'dc_included',
    'nyquist_raw_solved', 'source_contact_enabled', 'nq_gll_points_per_element',
    'polynomial_degree', 'model_length_scale_m', 'model_time_scale_s',
    'sem_elements', 'gll_points_per_element', 'native_element_node_sides',
    'radial_dofs', 'spheroidal_dofs_max_l', 'toroidal_dofs_max_l',
    'source_element', 'source_element_fluid', 'receiver_element',
    'receiver_element_fluid', 'vertical_component', 'spectrum_units',
    'energy_definition', 'pole_only_nonradial_adjoint', 'forward_only_baseline',
)


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def load(path):
    return json.loads(path.read_text())


def finite(value):
    if isinstance(value, dict):
        for v in value.values(): finite(v)
    elif isinstance(value, list):
        for v in value: finite(v)
    elif isinstance(value, (int, float)) and not isinstance(value, bool):
        if not math.isfinite(value): raise ValueError('nonfinite JSON value')


def verified_job(root, job_id):
    if Path(job_id).name != job_id or job_id in ('.', '..'):
        raise ValueError('invalid job id')
    state_bytes = (root / 'state.json').read_bytes()
    state = json.loads(state_bytes)
    matches = [j for j in state['jobs'] if j['id'] == job_id]
    if len(matches) != 1: raise ValueError('missing/duplicate job ' + job_id)
    job = matches[0]
    if not (job['status'] == 'completed' and job['returncode'] == 0
            and job['inputs_unchanged'] and job['cleanup']['verified_gone']
            and not job.get('failure_reasons')):
        raise ValueError('unqualified job ' + job_id)
    if set(job['outputs']) != set(job['outputs_sha256']):
        raise ValueError('incomplete output manifest ' + job_id)
    files = {**job['outputs_sha256'], job['log']: job['log_sha256']}
    before = {}
    for rel, digest in files.items():
        if Path(rel).is_absolute() or '..' in Path(rel).parts:
            raise ValueError('unsafe output path')
        p = root / job_id / rel
        before[rel] = sha(p)
        if before[rel] != digest: raise ValueError('artifact hash mismatch ' + str(p))
        if p.suffix == '.json': finite(load(p))
    return job, before, hashlib.sha256(state_bytes).hexdigest()


def csv_array(path):
    data = np.genfromtxt(path, names=True, delimiter=',', dtype=float)
    if data.size == 0 or not data.dtype.names: raise ValueError('empty CSV ' + str(path))
    data = np.atleast_1d(data)
    a = np.column_stack([data[n] for n in data.dtype.names])
    if not np.isfinite(a).all(): raise ValueError('nonfinite CSV ' + str(path))
    return data.dtype.names, a


def compare_columns(a, b):
    if a.shape != b.shape: raise ValueError('array shape mismatch')
    if not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError('nonfinite array')
    out = []
    for i in range(a.shape[1]):
        error = float(np.max(np.abs(a[:, i] - b[:, i])))
        scale = float(max(np.max(np.abs(a[:, i])), np.max(np.abs(b[:, i]))))
        ratio = error / scale if scale else 0.0
        out.append({'column': i, 'max_abs_difference': error, 'scale': scale,
                    'scaled_max_difference': ratio, 'passed': ratio <= LIMIT})
    return out


def profile(path):
    records = [json.loads(line.split('MATHPROFILE ', 1)[1])
               for line in path.read_text().splitlines() if 'MATHPROFILE ' in line]
    if not records: return None
    if len(records) != 1: raise ValueError('duplicate profiler records')
    r = records[0]; finite(r)
    if r['solver'] != 'dspec' or r['objective'] != -1:
        raise ValueError('unexpected profiler identity')
    total = r['stages']['physical_case_total']['inclusive_s']
    exclusive = 0.0
    for stage in r['stages'].values():
        if any(stage[x] < 0 for x in ('exclusive_s', 'inclusive_s', 'calls')):
            raise ValueError('negative profiler values')
        exclusive += stage['exclusive_s']
    if abs(total - exclusive) > max(1e-8, total * 2e-12):
        raise ValueError('profiler scope accounting fails closure')
    return r


def validate_domain(data, result, domain, hashes):
    """Validate each job independently, rather than merely matching two files."""
    probe = result['mode'] == 'kernel-probe'
    if result['mode'] not in ('kernel-probe', 'kernel'):
        raise ValueError('comparison requires a kernel result')
    if result['source_contact_enabled'] is not False or domain['source_contact_enabled'] is not False:
        raise ValueError('this comparator covers contacts-off volume kernels only')
    nrad = domain['selected_radial_trace_count']; nt = domain['n_theta']; np_ = domain['n_phi']
    angles = domain['sampled_angular_points_per_trace']; points = domain['actual_sampled_point_count']
    for v in (nrad, nt, np_, angles, points):
        if not isinstance(v, int) or isinstance(v, bool) or v <= 0:
            raise ValueError('invalid positive grid count')
    if nt != 2 * result['degree'] + 1 or np_ < 4 * result['degree'] + 1:
        raise ValueError('grid does not resolve complete kernel bandwidth')
    if not (points == nrad * angles and angles <= nt * np_):
        raise ValueError('sampled angle/radial product differs from point count')
    expected = {'full_target_angular_points_per_trace': nt * np_,
                'full_target_point_count': nt * np_ * nrad,
                'full_target_kernel_binary_bytes': nt * np_ * nrad * 48,
                'actual_kernel_binary_bytes': points * 48,
                'sampled_angle_subset_is_complete': angles == nt * np_,
                'diagnostic_angle_subset_only': probe,
                'actual_contact_binary_bytes': 0, 'full_target_contact_binary_bytes': 0,
                'actual_combined_output_bytes': points * 48}
    for k, v in expected.items():
        if domain.get(k) != v: raise ValueError('inconsistent domain field ' + k)
    if not probe and angles != nt * np_: raise ValueError('incomplete non-probe grid')
    for name in ('radial_nodes.csv', 'kernel_radial_nodes.csv', 'angular_grid.csv'):
        if 'data/' + name not in hashes: raise ValueError('unreceipted grid ' + name)
    h, all_sides = csv_array(data / 'radial_nodes.csv')
    if h != ('element', 'node', 'radius_m', 'element_width_m', 'fluid'):
        raise ValueError('unknown native radial table columns')
    h, selected = csv_array(data / 'kernel_radial_nodes.csv')
    if h != ('radial_index', 'element', 'node', 'radius_m', 'element_width_m', 'fluid'):
        raise ValueError('unknown selected radial table columns')
    if len(selected) != nrad or not np.array_equal(selected[:, 0], np.arange(nrad)):
        raise ValueError('selected radial side count/index mismatch')
    if len(all_sides) != domain['available_radial_trace_sides']:
        raise ValueError('available radial side count mismatch')
    if len(all_sides) > result['native_element_node_sides']:
        raise ValueError('too many positive native radial sides')
    if result['native_element_node_sides'] != result['sem_elements'] * result['gll_points_per_element']:
        raise ValueError('native SEM node count mismatch')
    for rows in (all_sides, selected[:, 1:]):
        if not (np.all(rows[:, :2] == np.floor(rows[:, :2]))
                and np.all((rows[:, 0] >= 0) & (rows[:, 0] < result['sem_elements']))
                and np.all((rows[:, 1] >= 0) & (rows[:, 1] < result['gll_points_per_element']))
                and np.all(rows[:, 2:4] > 0) and np.all(np.isin(rows[:, 4], [0, 1]))):
            raise ValueError('invalid native radial/material side')
        if len({tuple(x[:2]) for x in rows}) != len(rows):
            raise ValueError('duplicate element/node side')
    native_by_id = {tuple(row[:2]): row for row in all_sides}
    for row in selected[:, 1:]:
        if tuple(row[:2]) not in native_by_id or not np.array_equal(row, native_by_id[tuple(row[:2])]):
            raise ValueError('selected side differs from native radius/material identity')
    h, grid = csv_array(data / 'angular_grid.csv')
    if h != ('theta_index', 'theta_rad', 'theta_weight') or len(grid) != nt:
        raise ValueError('invalid theta grid shape/header')
    if not (np.array_equal(grid[:, 0], np.arange(nt)) and np.all((grid[:, 1] > 0) & (grid[:, 1] < np.pi))
            and np.all(np.diff(grid[:, 1]) > 0) and np.all(grid[:, 2] > 0)
            and abs(float(np.sum(grid[:, 2])) - 2) <= 1e-12):
        raise ValueError('invalid theta quadrature')
    if probe:
        name = 'kernel_probe_angles.csv'
        if 'data/' + name not in hashes: raise ValueError('unreceipted sampled angles')
        h, g = csv_array(data / name)
        if h != ('sample_index', 'theta_index', 'phi_index', 'theta_rad', 'phi_rad', 'theta_weight', 'phi_weight') or len(g) != angles:
            raise ValueError('sampled angle shape/header mismatch')
        idx = np.arange(angles); ti = idx // np_; pi = idx % np_
        expected_angles = np.column_stack((idx, ti, pi, grid[ti, 1], 2 * np.pi * pi / np_, grid[ti, 2], np.full(angles, 2 * np.pi / np_)))
        if not np.allclose(g, expected_angles, rtol=LIMIT, atol=0):
            raise ValueError('sampled angle coordinates/order/weights differ from native grid')


def assess(reference_root, reference_id, candidate_root, candidate_id, exact=False):
    roots = [reference_root, candidate_root]; ids = [reference_id, candidate_id]
    receipts = [verified_job(r, j) for r, j in zip(roots, ids)]
    dirs = [r / j / 'data' for r, j in zip(roots, ids)]
    results = [load(d / 'result.json') for d in dirs]
    for k in PHYSICAL:
        if k not in results[0] or k not in results[1] or results[0][k] != results[1][k]:
            raise ValueError('physical metadata differs or absent: ' + k)
    for r in results:
        if not (r['energy_m2_s'] > 0 and 0 <= r['energy_adjoint_relative_error_vs_2E'] <= 1e-10):
            raise ValueError('energy/adjoint identity failed')
        recomputed = abs(r['energy_adjoint_qdotp'] - 2 * r['energy_m2_s']) / (2 * r['energy_m2_s'])
        if not math.isfinite(recomputed) or recomputed > 1e-10:
            raise ValueError('recomputed qdot/(2E) identity failed')
    comparison = {'energy': compare_columns(
        np.array([[results[0]['energy_m2_s'], results[0]['energy_adjoint_qdotp']]]),
        np.array([[results[1]['energy_m2_s'], results[1]['energy_adjoint_qdotp']]]))}
    mode = results[0]['mode']
    if results[1]['mode'] != mode: raise ValueError('kernel mode mismatch')
    binary = 'kernel_probe.bin' if mode == 'kernel-probe' else 'kernel.bin'
    domain_name = 'kernel_probe.json' if mode == 'kernel-probe' else 'kernel_domain.json'
    domains = [load(d / domain_name) for d in dirs]
    for d, r, dom, rec in zip(dirs, results, domains, receipts):
        validate_domain(d, r, dom, rec[1])
    for k in ('selected_radial_trace_count', 'n_theta', 'n_phi',
              'actual_sampled_point_count', 'sampled_angular_points_per_trace'):
        if domains[0][k] != domains[1][k]: raise ValueError('domain mismatch ' + k)
    for d, m in zip(dirs, domains):
        if not ('point-major' in m['binary_layout'] and 'doubles[A,C,F,L,N,kappa]' in m['binary_layout']):
            raise ValueError('unrecognized binary layout')
        if (d / binary).stat().st_size != m['actual_sampled_point_count'] * 6 * 8:
            raise ValueError('binary size differs from declared grid')
    names = ['trace.csv', 'spectrum.csv', 'radial_nodes.csv',
             'kernel_radial_nodes.csv', 'angular_grid.csv']
    if mode == 'kernel-probe': names += ['kernel_probe_angles.csv']
    exact_files = {}
    for name in names:
        # The supervisor must bind every compared artifact, including grid CSVs.
        for hashes in (receipts[0][1], receipts[1][1]):
            if 'data/' + name not in hashes: raise ValueError('unreceipted comparison file ' + name)
        ha, a = csv_array(dirs[0] / name); hb, b = csv_array(dirs[1] / name)
        if ha != hb: raise ValueError('CSV headers differ ' + name)
        if name not in ('trace.csv', 'spectrum.csv') and not np.array_equal(a, b):
            raise ValueError('sampling grid/side identity differs ' + name)
        comparison[name] = compare_columns(a, b)
        exact_files[name] = sha(dirs[0] / name) == sha(dirs[1] / name)
    for hashes in (receipts[0][1], receipts[1][1]):
        if 'data/' + binary not in hashes: raise ValueError('unreceipted kernel')
    a = np.memmap(dirs[0] / binary, dtype=np.float64, mode='r').reshape(-1, 6)
    b = np.memmap(dirs[1] / binary, dtype=np.float64, mode='r').reshape(-1, 6)
    comparison[binary] = compare_columns(a, b)
    exact_files[binary] = sha(dirs[0] / binary) == sha(dirs[1] / binary)
    profiles = [profile(r / j / rec[0]['log']) for r, j, rec in zip(roots, ids, receipts)]
    after = [verified_job(r, j) for r, j in zip(roots, ids)]
    if any(a[0:2] != b[0:2] for a, b in zip(receipts, after)):
        raise ValueError('completed-job evidence changed during comparison')
    passed = all(x['passed'] for rows in comparison.values() for x in rows)
    if exact: passed = passed and all(exact_files.values())
    return {'status': 'PASS' if passed else 'FAIL', 'tolerance': LIMIT,
            'require_byte_exact': exact, 'reference_job': reference_id,
            'candidate_job': candidate_id, 'comparison': comparison,
            'byte_exact_files': exact_files, 'profiles': profiles,
            'telemetry': [rec[0] for rec in receipts],
            'state_sha256_at_read': [rec[2] for rec in receipts],
            'verified_artifact_sha256': [rec[1] for rec in receipts]}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('reference-root', 'reference-job', 'candidate-root', 'candidate-job', 'output'):
        p.add_argument('--' + name, required=True)
    p.add_argument('--exact', action='store_true')
    a = p.parse_args()
    try:
        report = assess(Path(a.reference_root), a.reference_job,
                        Path(a.candidate_root), a.candidate_job, a.exact)
    except Exception as e:
        report = {'status': 'FAIL', 'error_type': type(e).__name__, 'error': str(e)}
    with Path(a.output).open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False); stream.write('\n')
    return 0 if report['status'] == 'PASS' else 1


if __name__ == '__main__':
    raise SystemExit(main())
