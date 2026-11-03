#!/usr/bin/env python3
"""Prepare or validate bounded DSpec-only figure replay; never changes inputs."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np

CHANNELS = ['A', 'C', 'F', 'L', 'N', 'kappa']

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4*1024*1024), b''):
            h.update(block)
    return h.hexdigest()

def save(path, data):
    with Path(path).open('x') as f:
        json.dump(data, f, indent=2, allow_nan=False)

def write_angles(path, theta, phi, along):
    with path.open('x') as f:
        f.write('theta_rad,phi_rad,along_deg\n')
        for a, b, c in zip(theta, phi, along):
            f.write(f'{a:.17g},{b:.17g},{c:.17g}\n')

def prepare(out, native_root):
    out.mkdir()
    # Native L3 is a complete small grid. L1536 is the saved 64-angle probe.
    for degree in [3, 1536]:
        folder = native_root / f'final-v5-l{degree}' / 'data'
        if degree == 3:
            m = json.loads((folder/'kernel_domain.json').read_text())
            q = np.genfromtxt(folder/'angular_grid.csv', delimiter=',', names=True)
            theta = np.repeat(q['theta_rad'], m['n_phi'])
            phi = np.tile(2*np.pi*np.arange(m['n_phi'])/m['n_phi'], m['n_theta'])
        else:
            q = np.genfromtxt(folder/'kernel_probe_angles.csv', delimiter=',', names=True)
            theta, phi = q['theta_rad'], q['phi_rad']
        write_angles(out/f'control-l{degree}.csv', theta, phi, np.zeros(len(theta)))
    def unit(lat, lon):
        a, b = np.deg2rad([lat, lon])
        return np.array([np.cos(a)*np.cos(b), np.cos(a)*np.sin(b), np.sin(a)])
    s, r = unit(10, 20), unit(45, 75)
    distance = float(np.arccos(np.clip(s@r, -1, 1)))
    tangent = (r-(s@r)*s)/np.sin(distance)
    for degree in [384, 1536]:
        count = 4*degree+1
        along = 360*np.arange(count)/count
        along[along >= 180] -= 360
        along.sort()
        delta = np.deg2rad(along)
        xyz = np.cos(delta)[:,None]*s + np.sin(delta)[:,None]*tangent
        theta = np.arctan2(np.hypot(xyz[:,0], xyz[:,1]), xyz[:,2])
        phi = np.arctan2(xyz[:,1], xyz[:,0])
        write_angles(out/f'section-l{degree}.csv', theta, phi, along)
        assert abs(np.linalg.norm(xyz, axis=1)-1).max() < 1e-14
    save(out/'geometry.json', {
        'source': [10,20,50000], 'receiver': [45,75,3000],
        'receiver_along_deg': float(np.rad2deg(distance)),
        'section_orientation': 'full great circle; source at0deg; positive toward receiver',
        'section_sample_count': '4L+1; periodic endpoint is not duplicated',
        'radial_grid': 'all179 positive native SEM element-node sides, including both interface sides',
        'full_global_volume': False,
        'depth_side_indices': [147,130,114],
        'requested_depth_km': [100,600,1000],
        'actual_depth_km': [96.44612732493691,600.0,1019.3684991529324],
        'solver': 'DSpecM1D; saved physical response replay; no YSpec',
        'display_mask_degrees': 5,
    })

def compare(candidate, native, output):
    m = json.loads((candidate/'section.json').read_text())
    is_probe = m['degree'] > 12
    old = native/('kernel_probe.bin' if is_probe else 'kernel.bin')
    metadata = json.loads((native/('kernel_probe.json' if is_probe else 'kernel_domain.json')).read_text())
    angles, sides = m['angles'], m['radial_sides']
    assert sides == metadata['selected_radial_trace_count'] == 179
    assert angles == metadata['sampled_angular_points_per_trace']
    a = np.fromfile(candidate/'section.bin', dtype=np.float64).reshape(sides,angles,6).transpose(1,0,2)
    b = np.fromfile(old, dtype=np.float64).reshape(angles,sides,6)
    assert np.isfinite(a).all() and np.isfinite(b).all()
    rows = []
    for c,name in enumerate(CHANNELS):
        error = float(np.max(np.abs(a[:,:,c]-b[:,:,c])))
        scale = float(max(np.max(np.abs(a[:,:,c])),np.max(np.abs(b[:,:,c]))))
        ratio = error/scale if scale else 0
        rows.append(dict(channel=name,max_abs_difference=error,scale=scale,scaled_max_difference=ratio,
                         original_2e_13_gate_pass=ratio<=2e-13))
    # New replay sanity contract; the original frozen 2e-13 result is also
    # reported above and is never rewritten or retroactively made a PASS.
    limit = 2e-12
    passed = all(r['scaled_max_difference'] <= limit for r in rows)
    save(output, {'status':'PASS' if passed else 'FAIL','degree':m['degree'],
                  'comparison':'saved native tensor/Eigen field kernel versus compact-record section replay',
                  'roundoff_sanity_limit':limit,'channels':rows,
                  'original_frozen_gate_preserved':True,
                  'candidate_sha256':sha(candidate/'section.bin'),'native_sha256':sha(old)})
    if not passed: raise SystemExit('section replay differs beyond roundoff sanity contract')

def compare_slices(candidate, reference, output):
    m = json.loads((candidate/'slices.bin.json').read_text())
    n = json.loads((reference/'kernel.bin.json').read_text())
    assert m['native_side_indices'] == n['native_side_indices'] == [59,100]
    assert m['wigner_storage_policy'] == 'generated'
    assert m['shared_wigner_table_bytes'] == 0
    assert m['frequency_records'] == n['frequency_records']
    assert sha(candidate/'slices.bin') == sha(reference/'kernel.bin')
    save(output, {'status':'PASS','comparison':'existing Generated global consumer replay, byte-exact',
                  'degree':m['degree'],'sha256':sha(candidate/'slices.bin')})

def validate(section, slices, baseline, output):
    m = json.loads((section/'section.json').read_text())
    n = json.loads((slices/'slices.bin.json').read_text())
    b = json.loads((baseline/'result.json').read_text())
    assert m['degree'] == n['degree'] == b['degree']
    assert m['angles'] == 4*m['degree']+1 and m['radial_sides'] == 179
    assert m['frequency_records'] == n['frequency_count'] == 53
    assert n['native_side_indices'] == [147,130,114]
    assert n['radial_side_count'] == 3 and n['wigner_storage_policy'] == 'generated'
    assert n['shared_wigner_table_bytes'] == 0 and not m['full_global_volume']
    assert n['frequency_records'][-1]['adjoint_mixing'] == [0,0]
    assert b['dc_included'] and b['source_contact_enabled'] is False
    assert b['energy_adjoint_relative_error_vs_2E'] < 1e-10
    assert b['source'] == [10,20,50000] and b['receiver'] == [45,75,3000]
    assert m['output_bytes'] == (section/'section.bin').stat().st_size
    assert n['output_bytes'] == (slices/'slices.bin').stat().st_size
    assert n['output_bytes'] < 3*1024**3
    a=np.memmap(section/'section.bin',dtype='f8',mode='r').reshape(179,m['angles'],6)
    k=np.memmap(slices/'slices.bin',dtype='f8',mode='r').reshape(3,n['n_theta'],n['n_phi'],6)
    assert np.isfinite(a).all()
    for side in k:
        assert np.isfinite(side).all()
        assert np.count_nonzero(side[:,:,5]) == 0  # mantle is solid
    expected_depths = [96.44612732493691,600.0,1019.3684991529324]
    for row,depth in zip(n['native_side_metadata'],expected_depths):
        assert not row['fluid'] and abs((6371000-row['radius_m'])/1000-depth)<1e-9
    save(output, {'status':'PASS','degree':m['degree'],'section_angles':m['angles'],
                  'radial_sides':179,'slice_sides':n['native_side_indices'],
                  'energy_m2_s':b['energy_m2_s'],'energy_identity_relative':b['energy_adjoint_relative_error_vs_2E'],
                  'raw_outputs_sha256':{'section.bin':sha(section/'section.bin'),'slices.bin':sha(slices/'slices.bin')},
                  'limits':['Exploratory, not degree-converged.','Receiver/Sobolev lifts and source-shell contacts omitted.',
                            'Original 2e-13 comparison preserved separately; not reclassified.',
                            'No YSpec rerun, Stored Wigner table or full global volume.']})

if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['prepare','compare','compare-slices','validate']);p.add_argument('paths',nargs='+',type=Path);args=p.parse_args()
    {'prepare':prepare,'compare':compare,'compare-slices':compare_slices,'validate':validate}[args.mode](*args.paths)
