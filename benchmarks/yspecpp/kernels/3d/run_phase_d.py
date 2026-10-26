#!/usr/bin/env python3
"""Guarded Phase D baseline driver. Perturbed runs are intentionally separate."""
from __future__ import annotations
import argparse
import hashlib
from importlib.metadata import PackageNotFoundError, version
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
YSPEC = Path(os.environ.get('PHASE_D_YSPEC_SOURCE_DIR', str(ROOT.parent / 'YSpecpp'))).expanduser().resolve()
DEPENDENCY_ROOT = os.environ.get('PHASE_D_DEPENDENCY_ROOT')
DSPEC_PLANETMODEL = os.environ.get('PHASE_D_DSPEC_PLANETMODEL_SOURCE_DIR')
YSPEC_PLANETMODEL = os.environ.get('PHASE_D_YSPEC_PLANETMODEL_SOURCE_DIR')
FIXTURE = HERE / 'fixture.json'
SOURCE_DECK = HERE.parent / 'homogeneous_shell.deck'
D_SPEC = ROOT
CONFIG = json.loads(FIXTURE.read_text())
PINS = CONFIG['dependencies']
BREAKS = CONFIG['model']['radial_breaks_m']


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def git(path: Path, *args: str) -> str:
    return subprocess.check_output(['git', '-C', str(path), *args], text=True).strip()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def runtime_versions() -> dict[str, str | None]:
    def package(name: str) -> str | None:
        try:
            return version(name)
        except PackageNotFoundError:
            return None
    return {'python': sys.version, 'python_executable': sys.executable,
            'numpy': package('numpy'), 'scipy': package('scipy'),
            'matplotlib': package('matplotlib')}


def check_empty_first(out: Path) -> None:
    # This is deliberately the first filesystem-sensitive action: a rejected
    # destination is untouched and no configuration or build is created.
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise RuntimeError(f'output must be absent or empty; refusing to touch {out}')


def source_hashes() -> dict[str, str]:
    suffixes = {'.py', '.cpp', '.hpp', '.h', '.json', '.txt', '.md', '.deck', '.cmake'}
    files = sorted(p for p in HERE.rglob('*') if p.is_file() and p.suffix in suffixes)
    return {str(p.relative_to(ROOT)): sha(p) for p in files}


def git_files(path: Path) -> list[str]:
    raw = git(path, 'ls-files', '--modified', '--others', '--exclude-standard')
    return [line for line in raw.splitlines() if line]


def checkout_identity(path: Path, expected: str | None = None) -> dict:
    if not path.is_dir():
        return {'path': str(path), 'available_before_build': False, 'expected_head': expected}
    head = git(path, 'rev-parse', 'HEAD')
    require(expected is None or head == expected, f'dependency pin mismatch at {path}: expected {expected}, got {head}')
    files = git_files(path)
    return {'path': str(path), 'available_before_build': True, 'head': head,
            'dirty_files': files,
            'dirty_sha256': {name: sha(path/name) for name in files if (path/name).is_file()}}


def checks() -> dict:
    require((D_SPEC / 'DSpecM1D/Windowed3DKernels').is_file() and (D_SPEC / 'DSpecM1D/src/Windowed3DKernels.h').is_file(),
            'DSpec opt-in 3D kernel API is unavailable')
    require((YSPEC / 'YSpec/Kernels.hpp').is_file() and (YSPEC / 'src/Spectra.hpp').is_file(),
            'YSpec sampled-kernel/spectrum API is unavailable')
    require((YSPEC / 'CMakeLists.txt').is_file(), f'missing YSpecpp source checkout: {YSPEC}')
    dep_paths = {}
    if DEPENDENCY_ROOT:
        dep_root = Path(DEPENDENCY_ROOT).expanduser().resolve()
        dep_paths = {name: checkout_identity(dep_root/f'{name.lower()}-src', pin)
                     for name, pin in PINS.items() if name != 'planetmodel'}
    if DSPEC_PLANETMODEL or (DEPENDENCY_ROOT and (Path(DEPENDENCY_ROOT).expanduser().resolve()/'planetmodel-src').is_dir()):
        dplanet = Path(DSPEC_PLANETMODEL).expanduser().resolve() if DSPEC_PLANETMODEL else Path(DEPENDENCY_ROOT).expanduser().resolve()/'planetmodel-src'
        dep_paths['DSpec planetmodel'] = checkout_identity(dplanet, CONFIG['dspec_planetmodel'])
    if YSPEC_PLANETMODEL:
        dep_paths['YSpec planetmodel'] = checkout_identity(Path(YSPEC_PLANETMODEL).expanduser().resolve(), PINS['planetmodel'])
    require(SOURCE_DECK.is_file(), f'missing source deck {SOURCE_DECK}')
    return {'dspec_head': git(ROOT, 'rev-parse', 'HEAD'),
            'dspec_dirty_files': git_files(ROOT),
            'dspec_dirty_sha256': {name: sha(ROOT/name) for name in git_files(ROOT) if (ROOT/name).is_file()},
            'yspec_source': checkout_identity(YSPEC),
            'configured_dependency_source_root': str(Path(DEPENDENCY_ROOT).expanduser().resolve()) if DEPENDENCY_ROOT else None,
            'available_dependency_checkouts': dep_paths,
            'dependency_pins': PINS,
            'phase_d_source_sha256': source_hashes(),
            'runtime_versions': runtime_versions(),
            'source_deck_sha256': sha(SOURCE_DECK),
            'origin_preflight_sha256': CONFIG['origin_preflight_sha256'],
            'fixture_sha256': sha(FIXTURE),
            'frozen_planning_protocol_sha256': sha(HERE/'protocol.json')}


def snapshot_sources(out: Path, provenance: dict) -> dict:
    target = out/'source_snapshot'
    target.mkdir()
    snapshots = []
    source_sets = [('dspec', ROOT, provenance['dspec_dirty_sha256'])]
    yspec_identity = provenance['yspec_source']
    source_sets.append(('yspec', Path(yspec_identity['path']), yspec_identity['dirty_sha256']))
    for label, checkout, hashes in source_sets:
        for relative, expected in sorted(hashes.items()):
            source = checkout/relative
            require(source.is_file() and sha(source) == expected,
                    f'source changed before snapshot: {checkout}/{relative}')
            destination = target/'checkouts'/label/relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)
            require(sha(destination) == expected, f'source snapshot hash mismatch: {relative}')
            snapshots.append({'checkout': label, 'path': relative, 'sha256': expected})
    for name, digest in provenance['phase_d_source_sha256'].items():
        require(provenance['dspec_dirty_sha256'].get(name) == digest,
                f'Phase D file is missing or mismatched in the source snapshot manifest: {name}')
    manifest = {'schema': 'dspec-phase-d-source-snapshot-v1',
                'dspec_head': provenance['dspec_head'],
                'dspec_dirty_sha256': provenance['dspec_dirty_sha256'],
                'phase_d_source_sha256': provenance['phase_d_source_sha256'],
                'yspec_source': provenance['yspec_source'],
                'configured_dependency_source_root': provenance['configured_dependency_source_root'],
                'available_dependency_checkouts': provenance['available_dependency_checkouts'],
                'dependency_pins': provenance['dependency_pins'],
                'files': snapshots}
    manifest_path = target/'manifest.json'
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True)+'\n')
    return {'path': str(target), 'manifest_sha256': sha(manifest_path),
            'file_count': len(snapshots)}


def write_deck(path: Path) -> None:
    lines = SOURCE_DECK.read_text().splitlines()
    require(len(lines) >= 4, 'source deck has no material row')
    props = lines[3].split()[1:]
    require(len(props) == 5, 'expected five homogeneous material properties')
    out = lines[:3]
    for i, radius in enumerate(BREAKS):
        row = f'{radius:.1f} ' + ' '.join(props)
        out.append(row)
        if i not in (0, len(BREAKS)-1):
            out.append(row)
    path.write_text('\n'.join(out) + '\n')


def run(cmd: list[str], env: dict[str, str], log: Path) -> dict:
    with log.open('w') as stream:
        print('+', ' '.join(cmd))
        completed = subprocess.run(cmd, check=False, cwd=ROOT, env=env,
                                   stdout=stream, stderr=subprocess.STDOUT)
    reproducibility_keys = {
        'PATH', 'CC', 'CXX', 'CFLAGS', 'CXXFLAGS', 'LDFLAGS', 'CMAKE_GENERATOR',
        'CMAKE_PREFIX_PATH', 'PKG_CONFIG_PATH', 'LD_LIBRARY_PATH', 'PYTHONPATH',
        'PYTHONHOME', 'VIRTUAL_ENV', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
        'MKL_NUM_THREADS', 'BLAS_NUM_THREADS',
    }
    effective_environment = {k: env[k] for k in sorted(env)
                             if k in reproducibility_keys or k.startswith('PHASE_D_')}
    return {'argv': cmd, 'cwd': str(ROOT),
            'environment': {k: env[k] for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLAS_NUM_THREADS')},
            'reproducibility_environment': effective_environment,
            'environment_sha256': hashlib.sha256(json.dumps(sorted(env.items()), separators=(',', ':')).encode()).hexdigest(),
            'log': str(log), 'exit_status': completed.returncode}


def record_command(runfile: Path, state: dict, record: dict) -> None:
    state.setdefault('commands', []).append(record)
    runfile.write_text(json.dumps(state, indent=2) + '\n')
    if record['exit_status'] != 0:
        raise RuntimeError(f"command failed with status {record['exit_status']}: {' '.join(record['argv'])}")



def csvrows(path: Path) -> list[dict[str, str]]:
    import csv
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


def cval(row: dict[str, str], re_key: str, im_key: str) -> complex:
    return complex(float(row[re_key]), float(row[im_key]))


def baseline_metrics(out: Path) -> dict:
    dmeta = json.loads((out/'dspec/meta.json').read_text())
    ymeta = json.loads((out/'yspecpp/meta.json').read_text())
    for key in ('nt', 'dt_s', 'df_hz', 'epsilon_s_inv', 'physical_G', 'radius_m'):
        a, b = float(dmeta[key]), float(ymeta[key])
        require(abs(a-b) <= 1e-12*max(1.0, abs(a), abs(b)), f'baseline metadata mismatch: {key}: {a} vs {b}')
    require(dmeta['solver_bins'] == ymeta['solver_bins'] == [1, 2, 3], 'active solver bins differ from frozen grid')
    require(dmeta['degrees'] == ymeta['degrees'] == [0, 3], 'retained degree range differs from frozen grid')
    ds = csvrows(out/'dspec/spectra.csv'); ys = csvrows(out/'yspecpp/spectra.csv')
    for rows, label in ((ds, 'DSpec'), (ys, 'YSpec')):
        outside = [abs(cval(r, 'raw_real_m_s', 'raw_imag_m_s')) for r in rows if int(r['bin']) not in (1, 2, 3)]
        require(max(outside, default=0.0) == 0.0, f'{label} raw spectrum is nonzero outside active bins')
    spectrum = []
    for obj in ('J_Z', 'J_N'):
        dint = csvrows(out/'dspec/integrals.csv')
        S = float(next(r['S_m_s'] for r in dint if r['objective'] == obj))
        for comp in ('Z', 'N', 'E'):
            da = {(int(r['bin'])): cval(r, 'raw_real_m_s', 'raw_imag_m_s') for r in ds if r['objective']==obj and r['component']==comp}
            ya = {(int(r['bin'])): cval(r, 'raw_real_m_s', 'raw_imag_m_s') for r in ys if r['objective']==obj and r['component']==comp}
            dn = sum(abs(da[k])**2 for k in (1,2,3))**0.5
            yn = sum(abs(ya[k])**2 for k in (1,2,3))**0.5
            diff = sum(abs(da[k]-ya[k])**2 for k in (1,2,3))**0.5
            informative = min(dn,yn) > 1e-6*S
            spectrum.append({'objective':obj,'component':comp,'norm_D_m_s':dn,'norm_Y_m_s':yn,'difference_m_s':diff,'relative_l2':diff/max(dn,yn) if informative else None,'scaled_absolute':diff/S if not informative else None,'informative':informative,'gate_pass':(diff/max(dn,yn) <= 1e-5) if informative else None})
    profile_D=csvrows(out/'dspec/profiles.csv'); profile_Y=csvrows(out/'yspecpp/profiles.csv')
    require(len(profile_D)==len(profile_Y), 'profile sample row count differs')
    profiles=[]
    for obj in ('J_Z','J_N'):
      for mod in ('A','C','F','L','N'):
        aa=[r for r in profile_D if r['objective']==obj and r['modulus']==mod]
        bb=[r for r in profile_Y if r['objective']==obj and r['modulus']==mod]
        require(len(aa)==len(bb) and aa, f'empty or mismatched profile rows for {obj}/{mod}')
        W=0.0; vD=0.0; vY=0.0; delta=0.0; ID=0.0; IY=0.0
        for x,y in zip(aa,bb):
          for key in ('radius_m','theta_rad','phi_rad','dV_m3','deltaX_Pa'):
            require(abs(float(x[key])-float(y[key])) <= 1e-12*max(1.0,abs(float(x[key]))), f'shared profile quadrature mismatch: {key}')
          dv=float(x['dV_m3']); dx=float(x['deltaX_Pa']); kx=float(x['kernel_per_Pa_m3']); ky=float(y['kernel_per_Pa_m3'])
          W += dv*dx*dx; vD += dv*dx*dx*kx*kx; vY += dv*dx*dx*ky*ky; delta += dv*dx*dx*(kx-ky)**2
          ID += dv*dx*kx; IY += dv*dx*ky
        kwD=(vD/W)**0.5; kwY=(vY/W)**0.5; rms=(delta/W)**0.5
        ints=csvrows(out/'dspec/integrals.csv'); rec=next(r for r in ints if r['objective']==obj and r['modulus']==mod)
        S=float(rec['S_m_s']); ks=float(rec['Kscale_per_Pa_m3']); informative=min(kwD,kwY)>1e-6*ks
        ie=abs(ID-IY)/max(abs(ID),abs(IY),1e-14*S); iinformative=min(abs(ID),abs(IY))>1e-6*S
        profiles.append({'objective':obj,'modulus':mod,'samples':len(aa),'Kweighted_D_per_Pa_m3':kwD,'Kweighted_Y_per_Pa_m3':kwY,'Kscale_per_Pa_m3':ks,'profile_informative':informative,'Eprofile':rms/max(kwD,kwY) if informative else None,'RMSdiff_over_Kscale':rms/ks if not informative else None,'profile_gate_pass':(rms/max(kwD,kwY)<=1e-3) if informative else None,'I_D_m_s':ID,'I_Y_m_s':IY,'EJ':ie if iinformative else None,'abs_I_diff_over_S':abs(ID-IY)/S if not iinformative else None,'integrated_informative':iinformative,'integrated_gate_pass':(ie<=1e-3) if iinformative else None})
    require(all(r['gate_pass'] for r in spectrum), 'matched baseline spectrum gate failed or was uninformative')
    require(all(r['profile_gate_pass'] and r['integrated_gate_pass'] for r in profiles), 'matched baseline profile/integrated gate failed or was uninformative')
    mesh=csvrows(out/'dspec/mesh_nodes.csv'); local=[r for r in mesh if float(r['source_local_weight'])>0]
    source=float(dmeta['source_radius_m']); bounds=dmeta['source_element_bounds_m']
    require(bounds[0] <= source <= bounds[1], 'upper source element does not contain source radius')
    lower_nodes=sum(1 for r in local if 6311000<float(r['radius_m'])<source)
    upper_nodes=sum(1 for r in local if source<float(r['radius_m'])<6331000)
    require(lower_nodes>0 and upper_nodes>0, 'source-local envelope lacks interior SEM nodes on one side of source')
    oracle = json.loads((out/'oracle/baseline.json').read_text())
    strain_check = oracle['cartesian_strain_fd_check']
    require(strain_check['passed'] and
            float(strain_check['max_absolute_component_error']) <= float(strain_check['tolerance']),
            'independent spherical-to-Cartesian strain self-check failed')
    oracle_errors = oracle['baseline_spectrum_agreement']['relative_l2_by_quantity_component']
    oracle_gate = float(CONFIG['gates']['oracle_baseline_relative'])
    require(oracle['baseline_spectrum_agreement']['passed'], 'independent oracle baseline spectrum bridge failed')
    require(oracle['baseline_spectrum_agreement']['predeclared_tolerance'] == oracle_gate,
            'oracle baseline gate differs from fixture')
    require(max(oracle_errors.values()) <= oracle_gate, 'oracle baseline spectrum error exceeds frozen gate')
    objectives = oracle['objectives']
    integrals = csvrows(out/'dspec/integrals.csv')
    objective_checks = {}
    for obj in ('J_Z', 'J_N'):
        drow = next(r for r in integrals if r['objective'] == obj)
        scale = float(drow['S_m_s'])
        objective_checks[obj] = {}
        for key in ('J0_m_s', 'S_m_s'):
            a = float(objectives[obj][key]); b = float(drow[key])
            error_over_S = abs(a-b)/scale
            require(error_over_S <= oracle_gate,
                    f'independent oracle baseline mismatch for {obj} {key}: |delta|/S={error_over_S:.3e} exceeds {oracle_gate:.1e}')
            objective_checks[obj][key] = {'oracle': a, 'dspec': b,
                                          'absolute_difference_over_S': error_over_S,
                                          'gate_pass': True}
    born_gate = float(CONFIG['gates']['independent_born_relative'])
    signal_floor = float(CONFIG['gates']['informative_signal_floor_over_S'])
    born_checks = []
    for row in profiles:
        obj, mod = row['objective'], row['modulus']
        prediction = float(oracle['regular_volume_born_m_s_per_h'][obj][mod])
        kernel_integral = float(row['I_D_m_s'])
        scale = float(next(r['S_m_s'] for r in integrals if r['objective'] == obj))
        informative = min(abs(prediction), abs(kernel_integral)) > signal_floor*scale
        difference = abs(prediction-kernel_integral)
        relative = difference/max(abs(prediction), abs(kernel_integral)) if informative else None
        scaled_absolute = difference/scale if not informative else None
        passed = bool(informative and relative <= born_gate)
        born_checks.append({'objective': obj, 'modulus': mod,
                            'regular_born_m_s_per_h': prediction,
                            'dspec_kernel_integral_m_s': kernel_integral,
                            'difference_m_s': difference,
                            'informative': informative,
                            'relative_error': relative,
                            'difference_over_S': scaled_absolute,
                            'gate': born_gate, 'gate_pass': passed})
    born_all_pass = all(row['gate_pass'] for row in born_checks)
    local_rows = csvrows(out/'dspec/source_local_integrals.csv')
    local_profiles = csvrows(out/'dspec/source_local_profiles.csv')
    local_contacts = csvrows(out/'dspec/contacts.csv')
    expected_local = {(obj, mod) for obj in ('J_Z', 'J_N') for mod in ('C', 'F', 'L')}
    require({(r['objective'], r['modulus']) for r in local_rows} == expected_local,
            'source-local production Born table lacks a J/C/F/L row')
    source_local = []
    source_element = int(dmeta['source_element'])
    for row in local_rows:
        obj, mod = row['objective'], row['modulus']
        key_rows = [r for r in local_profiles if r['objective'] == obj and r['modulus'] == mod]
        require(key_rows, f'no source-local production samples for {obj}/{mod}')
        sample_elements = {int(r['element']) for r in key_rows}
        radii = [float(r['radius_m']) for r in key_rows]
        low, high = float(CONFIG['perturbations']['source_local']['support_m'][0]), float(CONFIG['perturbations']['source_local']['support_m'][1])
        radius_source = float(CONFIG['perturbations']['source_local']['source_m'])
        require(source_element in sample_elements and any(r < radius_source for r in radii) and
                any(r > radius_source for r in radii),
                f'{obj}/{mod} source-local samples omit the upper lift or one side of the source')
        endpoint_tolerance = 8.0*sys.float_info.epsilon*max(abs(low), abs(high))
        require(any(abs(float(r['radius_m'])-low) <= endpoint_tolerance and float(r['w']) == 0.0 for r in key_rows) and
                any(abs(float(r['radius_m'])-high) <= endpoint_tolerance and float(r['w']) == 0.0 for r in key_rows),
                f'{obj}/{mod} source-local profile endpoints are missing or nonzero')
        support_volume = float(row['Vsupport_m3'])
        analytic_volume = (4.0*3.141592653589793/3.0)*(high**3-low**3)
        require(abs(support_volume-analytic_volume) <= 1e-12*analytic_volume,
                f'{obj}/{mod} source-local support volume differs from the analytic shell')
        volume = float(row['I_volume_m_s_per_h'])
        contact = float(row['I_contact_m_s_per_h'])
        total = float(row['Kpred_total_m_s_per_h'])
        scale = float(next(r['S_m_s'] for r in integrals if r['objective'] == obj))
        require(abs(total-volume-contact) <= 1e-12*max(abs(total), abs(volume), abs(contact), 1e-14*scale),
                f'{obj}/{mod} source-local volume/contact split does not sum to Kpred')
        contact_sum = sum(float(r['contribution_m_s']) for r in local_contacts
                          if r['objective'] == obj and r['modulus'] == mod)
        require(abs(contact-contact_sum) <= 1e-12*max(abs(contact), abs(contact_sum), 1e-14*scale),
                f'{obj}/{mod} source-local contact differs from the saved contact quadrature')
        volume_sum = sum(float(r['contribution_m_s_per_h']) for r in key_rows)
        require(abs(volume-volume_sum) <= 1e-12*max(abs(volume), abs(volume_sum), 1e-14*scale),
                f'{obj}/{mod} source-local volume differs from its saved point contributions')
        source_local.append({
            'objective': obj, 'modulus': mod, 'samples': len(key_rows),
            'elements': sorted(sample_elements), 'source_element': source_element,
            'source_side_samples': [sum(r < radius_source for r in radii),
                                    sum(r > radius_source for r in radii)],
            'S_m_s': scale, 'Vsupport_m3': support_volume,
            'Kscale_per_Pa_m3': float(row['Kscale_per_Pa_m3']),
            'I_volume_m_s_per_h': volume, 'I_contact_m_s_per_h': contact,
            'Kpred_total_m_s_per_h': total, 'volume_contact_sum_pass': True})
    source_bridge = oracle['source_local_born_vs_dspec_integral']['directions']
    expected_bridge = {f'{obj}_{mod}' for obj in ('J_Z', 'J_N') for mod in ('C', 'F', 'L')}
    require(set(source_bridge) == expected_bridge,
            'independent source-local Born bridge must contain all six objective/modulus rows')
    source_bridge_rows = []
    for key in sorted(source_bridge):
        objective, modulus = key.rsplit('_', 1)
        terms = source_bridge[key]
        for term in ('volume', 'contact', 'total'):
            check = terms[term]
            source_bridge_rows.append({
                'objective': objective, 'modulus': modulus, 'term': term,
                'cartesian_born_m_s_per_h': float(check['cartesian_m_s_per_h']),
                'dspec_integral_m_s_per_h': float(check['dspec_m_s_per_h']),
                'informative': bool(check['informative']),
                'relative_error': check['relative_error'],
                'scaled_absolute_error': float(check['scaled_absolute_error']),
                'gate_pass': bool(check['passed'])})
    source_bridge_all_pass = all(row['gate_pass'] for row in source_bridge_rows)
    require(source_bridge_all_pass,
            'independent source-local volume/contact Born bridge failed or was uninformative')
    return {
        'baseline_metadata_match': True,
        'raw_zero_outside_bins_1_2_3': True,
        'S_J_Z_m_s': float(next(r['S_m_s'] for r in integrals if r['objective'] == 'J_Z')),
        'S_J_N_m_s': float(next(r['S_m_s'] for r in integrals if r['objective'] == 'J_N')),
        'J0_D': {o: float(next(r['J0_m_s'] for r in integrals if r['objective'] == o))
                 for o in ('J_Z', 'J_N')},
        'source_element': dmeta['source_element'],
        'source_element_bounds_m': bounds,
        'source_local_interior_nodes_lower_upper': [lower_nodes, upper_nodes],
        'spectrum': spectrum,
        'regular_profiles': profiles,
        'source_local_production_born': {
            'interpretation': 'DSpec production point-kernel integral over the frozen source-local profile, plus the separately integrated source-contact term.',
            'support_m': [low, high], 'profile': CONFIG['perturbations']['source_local']['profile'],
            'source_element': source_element, 'rows': source_local},
        'source_local_born_comparisons': {
            'definition': 'Independent Cartesian source-local volume and contact Born predictions compared with DSpec production integrals; source-owning upper-element lift is included in volume.',
            'relative_gate': born_gate, 'informative_signal_floor_over_S': signal_floor,
            'rows': source_bridge_rows, 'all_pass': source_bridge_all_pass},
        'cartesian_strain_fd_check': strain_check,
        'oracle_baseline': {
            'status': oracle['status'],
            'relative_l2_by_quantity_component': oracle_errors,
            'max_relative_l2': max(oracle_errors.values()),
            'predeclared_tolerance': oracle_gate,
            'objective_checks_normalized_by_S': objective_checks,
            'total_dofs': oracle['total_dofs'], 'mode_count': oracle['mode_count'],
            'objectives': objectives,
            'regular_volume_born_m_s_per_h': oracle['regular_volume_born_m_s_per_h'],
            'resource_estimate': oracle['resource_estimate']},
        'oracle_born_comparisons': {
            'definition': 'Compare regular-volume oracle Born dJ/dh to the DSpec sampled-kernel integral for the same source-vanishing profile; informative when both magnitudes exceed the configured floor times S.',
            'relative_gate': born_gate, 'informative_signal_floor_over_S': signal_floor,
            'rows': born_checks, 'all_pass': born_all_pass},
        'all_baseline_gates_pass': born_all_pass and source_bridge_all_pass}


def compare_baseline_reference(current_out: Path, reference_out: Path) -> dict:
    reference_out = reference_out.expanduser().resolve()
    require(reference_out.is_dir(), f'baseline reference directory is unavailable: {reference_out}')
    reference_run = json.loads((reference_out/'run.json').read_text())
    require(reference_run.get('state') == 'BASELINES_COMPLETE',
            'baseline reference run is not a completed baseline-only packet')
    reference_fixture = json.loads((reference_out/'fixture.json').read_text())
    current_fixture = json.loads((current_out/'fixture.json').read_text())
    require(current_fixture == reference_fixture, 'full-run fixture differs from the accepted baseline reference')
    current = json.loads((current_out/'baseline_metrics.json').read_text())
    reference = json.loads((reference_out/'baseline_metrics.json').read_text())
    require(reference.get('all_baseline_gates_pass') is True,
            'baseline reference metrics do not record all gates passing')
    for key in ('S_J_Z_m_s', 'S_J_N_m_s'):
        a, b = float(current[key]), float(reference[key])
        require(abs(a-b) <= 1e-12*max(abs(a), abs(b)), f'baseline reference mismatch for {key}: {a} vs {b}')
    def scales(metrics: dict) -> dict[tuple[str, str], float]:
        return {(row['objective'], row['modulus']): float(row['Kscale_per_Pa_m3'])
                for row in metrics['regular_profiles']}
    actual_scales, expected_scales = scales(current), scales(reference)
    require(actual_scales.keys() == expected_scales.keys(), 'baseline reference Kscale rows differ')
    for key, a in actual_scales.items():
        b = expected_scales[key]
        require(abs(a-b) <= 1e-12*max(abs(a), abs(b)), f'baseline reference Kscale mismatch for {key}: {a} vs {b}')
    return {'fixture_match': True, 'S_match_relative_tolerance': 1e-12,
            'Kscale_match_relative_tolerance': 1e-12,
            'S_J_Z_m_s': current['S_J_Z_m_s'], 'S_J_N_m_s': current['S_J_N_m_s'],
            'Kscale_per_Pa_m3': {f'{obj}/{mod}': value for (obj, mod), value in actual_scales.items()},
            'reference_source_sha256': reference_run.get('provenance', {}).get('phase_d_source_sha256'),
            'reference_executable_sha256': {name: record.get('copied_executable_sha256')
                                            for name, record in reference_run.get('builds', {}).items()},
            'reference_run_json_sha256': sha(reference_out/'run.json'),
            'reference_metrics_sha256': sha(reference_out/'baseline_metrics.json')}


def finite_difference_gate_summary(path: Path) -> dict:
    packet = json.loads(path.read_text())
    expected_steps = [float(x) for x in CONFIG['perturbations']['fractional_steps']]
    require([float(x) for x in packet['steps_h']] == expected_steps,
            'finite-difference packet steps differ from the frozen fixture')
    checks = []
    for support, support_record in packet['supports'].items():
        for modulus, objectives in support_record['directions'].items():
            for objective, row in objectives.items():
                if not row.get('informative_kernel', False):
                    checks.append({'support': support, 'objective': objective,
                                   'modulus': modulus, 'informative': False,
                                   'passed': None})
                    continue
                plateau = row.get('plateau') or {}
                comparison = row.get('born_comparison') or {}
                production = row.get('production_comparison') or {}
                born_pass = (not comparison.get('informative') or comparison.get('passed'))
                passed = bool(plateau.get('informative') and plateau.get('passed') and
                              production.get('informative') and production.get('passed') and born_pass)
                checks.append({'support': support, 'objective': objective,
                               'modulus': modulus, 'informative': True,
                               'plateau_pass': bool(plateau.get('passed')),
                               'production_comparison_pass': bool(production.get('passed')),
                               'born_comparison_pass': (bool(comparison.get('passed'))
                                                        if comparison.get('informative') else None),
                               'passed': passed})
    informative = [row for row in checks if row['informative']]
    require(bool(informative), 'finite-difference packet has no informative predictions')
    require(all(row['passed'] for row in informative),
            'finite-difference plateau or independent Born comparison gate failed')
    return {'informative_cases': len(informative),
            'weak_diagnostic_cases': len(checks)-len(informative),
            'all_informative_gates_pass': True,
            'checks': checks}

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', required=True, type=Path)
    ap.add_argument('--baseline-only', action='store_true', help='build and run matched baselines only')
    ap.add_argument('--full', action='store_true', help='run the frozen baseline and finite-difference campaign')
    ap.add_argument('--baseline-reference', type=Path,
                    help='optional completed baseline packet for physics and S/Kscale consistency checks')
    args = ap.parse_args()
    if args.baseline_only == args.full:
        ap.error('select exactly one of --baseline-only or --full')
    if args.baseline_reference and not args.full:
        ap.error('--baseline-reference is only used with --full')
    out = args.output.expanduser().resolve()
    output_guard_passed = False
    try:
        check_empty_first(out)
        output_guard_passed = True
        provenance = checks()
        out.mkdir(parents=True, exist_ok=True)
        source_snapshot = snapshot_sources(out, provenance)
        runfile = out / 'run.json'
        runstate = {'state': 'STARTED',
                    'stage': 'full campaign' if args.full else 'baseline-only; no finite-difference result',
                    'requested_mode': 'full' if args.full else 'baseline-only',
                    'provenance': provenance,
                    'source_snapshot': source_snapshot}
        runfile.write_text(json.dumps(runstate, indent=2) + '\n')
        model = out / 'common_homogeneous.deck'
        shutil.copy2(FIXTURE, out/'fixture.json')
        shutil.copy2(HERE/'protocol.json', out/'frozen-origin-preflight.json')
        write_deck(model)
        env = os.environ.copy()
        default_build_root = ROOT / 'build/phase-d-build-cache'
        build_root = Path(env.get('PHASE_D_BUILD_CACHE_ROOT', str(default_build_root))).expanduser().resolve()
        require(not (build_root == out or out in build_root.parents or build_root in out.parents),
                'mutable build cache and immutable run output must be separate')
        require(not any((parent/'run.json').is_file() for parent in [build_root, *build_root.parents]),
                'refusing to reuse a directory inside an existing run root as mutable build cache')
        build_root.mkdir(parents=True, exist_ok=True)
        env.update({'OMP_NUM_THREADS': str(CONFIG['resources']['omp_threads']),
                    'OPENBLAS_NUM_THREADS': str(CONFIG['resources']['blas_threads']),
                    'MKL_NUM_THREADS': str(CONFIG['resources']['blas_threads']),
                    'BLAS_NUM_THREADS': str(CONFIG['resources']['blas_threads'])})
        build_records = {}
        dependency_root_arg = str(Path(DEPENDENCY_ROOT).expanduser().resolve()) if DEPENDENCY_ROOT else ''
        for solver in ('dspec', 'yspecpp'):
            build = build_root / f'build-{solver}'
            dplanet_arg = (Path(DSPEC_PLANETMODEL).expanduser().resolve() if DSPEC_PLANETMODEL else
                           (Path(DEPENDENCY_ROOT).expanduser().resolve()/'planetmodel-src' if DEPENDENCY_ROOT else None))
            yplanet_arg = Path(YSPEC_PLANETMODEL).expanduser().resolve() if YSPEC_PLANETMODEL else None
            configure_args = ['cmake', '-S', str(HERE), '-B', str(build), f'-DSOLVER={solver}',
                              f'-DDSPEC_SOURCE_DIR={D_SPEC}', f'-DYSPEC_SOURCE_DIR={YSPEC}',
                              f'-DDEPENDENCY_ROOT={dependency_root_arg}']
            if dplanet_arg:
                configure_args.append(f'-DDSPEC_PLANETMODEL_SOURCE_DIR={dplanet_arg}')
            if yplanet_arg:
                configure_args.append(f'-DYSPEC_PLANETMODEL_SOURCE_DIR={yplanet_arg}')
            record_command(runfile, runstate, run(
                configure_args,
                env, out / f'configure-{solver}.log'))
            record_command(runfile, runstate, run(
                ['cmake', '--build', str(build), '--target', 'phase_d_adapter', f'-j{CONFIG["resources"]["jobs"]}'],
                env, out / f'build-{solver}.log'))
            cache = build/'CMakeCache.txt'
            executable = build/'phase_d_adapter'
            require(cache.is_file() and executable.is_file(), f'incomplete {solver} build identity')
            resolved_deps = {}
            for dep_name, pin in PINS.items():
                if dep_name == 'planetmodel':
                    if solver == 'dspec':
                        dep_path = (Path(DSPEC_PLANETMODEL).expanduser().resolve() if DSPEC_PLANETMODEL else
                                    (Path(DEPENDENCY_ROOT).expanduser().resolve()/'planetmodel-src' if DEPENDENCY_ROOT else build/'_deps'/'planetmodel-src'))
                    else:
                        dep_path = Path(YSPEC_PLANETMODEL).expanduser().resolve() if YSPEC_PLANETMODEL else build/'_deps'/'planetmodel-src'
                elif DEPENDENCY_ROOT:
                    dep_path = Path(DEPENDENCY_ROOT).expanduser().resolve()/f'{dep_name.lower()}-src'
                else:
                    dep_path = build/'_deps'/f'{dep_name.lower()}-src'
                if dep_path.is_dir():
                    expected_pin = (CONFIG['dspec_planetmodel'] if dep_name == 'planetmodel' and solver == 'dspec'
                                    else pin)
                    resolved_deps[dep_name] = checkout_identity(dep_path, expected_pin)
                elif dep_name != 'planetmodel':
                    raise RuntimeError(f'{solver} build did not resolve dependency source {dep_path}')
            require(source_hashes() == provenance['phase_d_source_sha256'],
                    'Phase D source/config changed while builds were running; discard this run and review the source snapshot')
            build_records[solver] = {'build_dir': str(build), 'cmake_cache_sha256': sha(cache),
                                     'build_source_hashes': source_hashes(),
                                     'resolved_dependencies': resolved_deps,
                                     'built_executable_sha256': sha(executable)}
        exec_dir = out/'executables'; exec_dir.mkdir()
        for solver in ('dspec', 'yspecpp'):
            shutil.copy2(build_root/f'build-{solver}'/'phase_d_adapter', exec_dir/f'{solver}_adapter')
            build_records[solver]['copied_executable'] = str(exec_dir/f'{solver}_adapter')
            build_records[solver]['copied_executable_sha256'] = sha(exec_dir/f'{solver}_adapter')
        runstate.update({'model': {'path': str(model), 'sha256': sha(model)},
                         'builds': build_records,
                         'build_cache_root': str(build_root),
                         'build_cache_mutable': True,
                         'environment': {k:env[k] for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLAS_NUM_THREADS')}})
        runfile.write_text(json.dumps(runstate, indent=2) + '\n')
        for name in ('dspec', 'yspecpp'):
            (out / name).mkdir()
        record_command(runfile, runstate, run(
            [str(exec_dir/'dspec_adapter'), str(model), str(out/'dspec'), '0.05', '8'], env, out/'run-dspec.log'))
        record_command(runfile, runstate, run(
            [str(exec_dir/'yspecpp_adapter'), str(model), str(out/'yspecpp'),
             str(out/'dspec/radial_nodes.csv'), str(out/'dspec/angles.csv')], env, out/'run-yspecpp.log'))
        oracle_dir = out/'oracle'; oracle_dir.mkdir()
        oracle_cmd = [sys.executable, str(HERE/'cartesian_ti_oracle.py'),
                      '--baseline', str(out/'dspec'), '--output', str(oracle_dir),
                      '--fixture', str(out/'fixture.json')]
        record_command(runfile, runstate, run(oracle_cmd, env, out/'run-oracle-baseline.log'))
        metrics = baseline_metrics(out)
        (out / 'baseline_metrics.json').write_text(json.dumps(metrics, indent=2) + '\n')
        baseline_status = 'PASS' if metrics['all_baseline_gates_pass'] else 'FAIL'
        max_spectrum = max(row['relative_l2'] for row in metrics['spectrum'] if row['informative'])
        max_profile = max(row['Eprofile'] for row in metrics['regular_profiles'] if row['profile_informative'])
        max_integral = max(row['EJ'] for row in metrics['regular_profiles'] if row['integrated_informative'])
        born_failed = sum(not row['gate_pass'] for row in metrics['oracle_born_comparisons']['rows'])
        (out / 'baseline_report.md').write_text(
            f"# Matched Phase D baselines: {baseline_status}\n\n"
            f"DSpec/YSpec maximum errors: spectrum `{max_spectrum:.3g}`, "
            f"weighted profile `{max_profile:.3g}`, integrated derivative `{max_integral:.3g}`. "
            f"Independent oracle bridge maximum: `{metrics['oracle_baseline']['max_relative_l2']:.3g}`. "
            f"Regular Born comparisons failing or uninformative: {born_failed}/{len(metrics['oracle_born_comparisons']['rows'])}.\n\n"
            "See [report.html](report.html) for the readable tables and plots, "
            "[baseline_metrics.json](baseline_metrics.json) for full precision, "
            "[oracle/baseline.json](oracle/baseline.json) for the saved modal packet, "
            "and [source_snapshot/manifest.json](source_snapshot/manifest.json) for archived source hashes.\n")
        baseline_report_record = run([sys.executable, str(HERE/'report_phase_d.py'), '--output', str(out)],
                                     env, out/'report-generation-baseline.log')
        runstate.setdefault('commands', []).append(baseline_report_record)
        require(baseline_report_record['exit_status'] == 0,
                'responsive baseline report or scientific plot generation failed')
        runfile.write_text(json.dumps(runstate, indent=2) + '\n')
        require(metrics['all_baseline_gates_pass'],
                'independent Cartesian regular Born predictions do not pass the frozen DSpec-kernel comparison gate')
        reference_check = compare_baseline_reference(out, args.baseline_reference) if args.baseline_reference else None
        if args.full:
            require(source_hashes() == provenance['phase_d_source_sha256'],
                    'Phase D source/config changed since baseline; stopping before finite differences')
            frozen_packet = {'fixture_sha256': sha(out/'fixture.json'),
                             'baseline_metrics_sha256': sha(out/'baseline_metrics.json'),
                             'oracle_baseline_sha256': sha(out/'oracle/baseline.json'),
                             'source_snapshot_manifest_sha256': source_snapshot['manifest_sha256'],
                             'phase_d_source_sha256': provenance['phase_d_source_sha256'],
                             'dspec_head': provenance['dspec_head'],
                             'dspec_dirty_sha256': provenance['dspec_dirty_sha256'],
                             'yspec_source': provenance['yspec_source'],
                             'dependency_pins': provenance['dependency_pins'],
                             'builds': build_records,
                             'baseline_reference_check': reference_check}
            runstate.update({'state': 'READY_FOR_FINITE_DIFFERENCES',
                             'stage': 'baseline and scale packet frozen before finite differences',
                             'frozen_before_finite_differences': frozen_packet})
            runfile.write_text(json.dumps(runstate, indent=2) + '\n')
            fd_cmd = [sys.executable, str(HERE/'cartesian_ti_oracle.py'),
                      '--finite-differences', '--baseline', str(out/'dspec'),
                      '--output', str(oracle_dir), '--fixture', str(out/'fixture.json')]
            record_command(runfile, runstate, run(fd_cmd, env, out/'run-oracle-finite-differences.log'))
            require((oracle_dir/'finite_differences.json').is_file(),
                    'finite-difference command completed without its result packet')
            require(source_hashes() == provenance['phase_d_source_sha256'],
                    'Phase D source/config changed during finite differences; results are not reproducible from the frozen packet')
            progress_path = oracle_dir/'finite_difference_progress.jsonl'
            progress_identity = ({'path': str(progress_path), 'sha256': sha(progress_path)}
                                 if progress_path.is_file() else None)
            runstate.update({'state': 'FINITE_DIFFERENCES_WRITTEN',
                             'stage': 'finite-difference packet saved; report generated before gate evaluation',
                             'finite_differences_sha256': sha(oracle_dir/'finite_differences.json'),
                             'finite_difference_progress': progress_identity})
            runfile.write_text(json.dumps(runstate, indent=2) + '\n')
        if args.full:
            report_record = run([sys.executable, str(HERE/'report_phase_d.py'), '--output', str(out)],
                                env, out/'report-generation-full.log')
            runstate.setdefault('commands', []).append(report_record)
            require(report_record['exit_status'] == 0,
                    'responsive full-campaign report or scientific plot generation failed')
            runstate['finite_differences_report_sha256'] = sha(out/'report.html')
            runstate.update({'provenance': provenance,
                'model': {'path': str(model), 'sha256': sha(model)},
                'metrics_sha256': sha(out/'baseline_metrics.json'),
                'environment': {k:env[k] for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLAS_NUM_THREADS')},
                'builds': build_records, 'build_cache_root': str(build_root), 'build_cache_mutable': True,
                'report_sha256': sha(out/'report.html'),
                'plot_sha256': {str(p.relative_to(out)): sha(p) for p in sorted((out/'plots').glob('*.png'))},
                'completed_stage_note': 'full campaign requested; report contains every saved baseline and partial/full FD result'})
            runfile.write_text(json.dumps(runstate, indent=2) + '\n')
            fd_gates = finite_difference_gate_summary(oracle_dir/'finite_differences.json')
            runstate.update({'state': 'PHASE_D_COMPLETE',
                             'stage': 'full bounded campaign; all frozen gates passed',
                             'finite_difference_gate_summary': fd_gates})
        else:
            runstate.update({'state': 'BASELINES_COMPLETE', 'stage': 'baseline-only', 'provenance': provenance})
        runstate.update({'provenance': provenance,
            'model': {'path': str(model), 'sha256': sha(model)}, 'metrics_sha256': sha(out/'baseline_metrics.json'),
            'environment': {k:env[k] for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLAS_NUM_THREADS')},
            'builds': build_records, 'build_cache_root': str(build_root), 'build_cache_mutable': True,
            'report_sha256': sha(out/'report.html'),
            'plot_sha256': {str(p.relative_to(out)): sha(p) for p in sorted((out/'plots').glob('*.png'))},
            'completed_stage_note': ('full bounded Phase D campaign and frozen FD packet' if args.full else
                                     'baseline evidence only; no finite-h result is included')})
        runfile.write_text(json.dumps(runstate, indent=2) + '\n')
        return 0
    except Exception as e:
        if not output_guard_passed:
            pass
        elif out.exists() and out.is_dir() and not any(out.iterdir()):
            # Failed preflight left no work products; preserve no empty marker.
            pass
        elif out.exists() and out.is_dir():
            (out / 'failure.txt').write_text(f'{type(e).__name__}: {e}\n')
            runfile = out / 'run.json'
            if runfile.exists():
                try:
                    state = json.loads(runfile.read_text()); state['state'] = 'FAILED'; state['error'] = str(e)
                    runfile.write_text(json.dumps(state, indent=2) + '\n')
                except Exception:
                    pass
        print(f'Phase D baseline stopped: {e}', file=sys.stderr)
        return 2

if __name__ == '__main__':
    raise SystemExit(main())
