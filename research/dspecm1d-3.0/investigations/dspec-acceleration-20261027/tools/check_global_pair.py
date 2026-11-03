#!/usr/bin/env python3
"""Read-only validation of a completed compact-response global-kernel job."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

import check_pair as pair

LIMIT = pair.LIMIT
CHANNELS = ('A', 'C', 'F', 'L', 'N', 'kappa')
PHYSICAL_OUTPUT_BINS = 52
SPOOL_BINS = 53
TOTAL_NATIVE_SIDES = 180
POSITIVE_RADIUS_SIDES = 179
GLOBAL_SELECTED_SIDES_MAX = 8
SUPPORTED_DEGREES = (3, 12)
SUPPORTED_SEM_DEGREE = 6
REFERENCE_FILES = (
    'result.json', 'trace.csv', 'spectrum.csv', 'radial_nodes.csv',
    'kernel_radial_nodes.csv', 'kernel_domain.json', 'kernel.bin',
)
GLOBAL_FILES = (
    'result.json', 'trace.csv', 'spectrum.csv', 'radial_nodes.csv',
    'pole_response_spool.bin', 'kernel.bin', 'kernel.bin.json',
)


def need_receipted(receipt: tuple, names: tuple[str, ...]) -> None:
    hashes = receipt[1]
    for name in names:
        if 'data/' + name not in hashes:
            raise ValueError('unreceipted required artifact ' + name)


def read_native(path: Path) -> tuple[tuple[str, ...], np.ndarray]:
    header, rows = pair.csv_array(path)
    expected = ('element', 'node', 'radius_m', 'element_width_m', 'fluid')
    if header != expected or rows.ndim != 2 or rows.shape[1] != 5:
        raise ValueError('unknown native radial table schema: ' + str(path))
    if (not np.all(rows[:, :2] == np.floor(rows[:, :2]))
            or not np.all(rows[:, :2] >= 0)
            or not np.all(rows[:, 2:4] > 0)
            or not np.all(np.isin(rows[:, 4], (0, 1)))):
        raise ValueError('invalid native radial side table: ' + str(path))
    ids = [tuple(row[:2].astype(np.int64)) for row in rows]
    if len(ids) != len(set(ids)):
        raise ValueError('duplicate native element/node side: ' + str(path))
    return header, rows


def _int(value: Any, name: str, positive: bool = False) -> int:
    if not isinstance(value, int) or isinstance(value, bool) or (positive and value <= 0):
        raise ValueError('invalid integer field ' + name)
    return value


def _finite_pair(value: Any, name: str) -> tuple[float, float]:
    if (not isinstance(value, list) or len(value) != 2
            or any(not isinstance(x, (int, float)) or isinstance(x, bool)
                   or not math.isfinite(x) for x in value)):
        raise ValueError('invalid complex field ' + name)
    return float(value[0]), float(value[1])


def validate_energy(result: dict) -> None:
    energy = result.get('energy_m2_s')
    qdotp = result.get('energy_adjoint_qdotp')
    reported = result.get('energy_adjoint_relative_error_vs_2E')
    if (not isinstance(energy, (int, float)) or not math.isfinite(energy) or energy <= 0
            or not isinstance(qdotp, (int, float)) or not math.isfinite(qdotp)
            or not isinstance(reported, (int, float)) or not math.isfinite(reported)
            or reported < 0 or reported > 1e-10):
        raise ValueError('invalid positive energy or reported qdot/(2E) audit')
    relative = abs(qdotp - 2.0 * energy) / (2.0 * energy)
    if relative > 1e-10:
        raise ValueError('recomputed qdot/(2E) identity failed')


def validate_physical_pair(reference: dict, candidate: dict,
                           candidate_mode: str = 'kernel') -> None:
    for key in pair.PHYSICAL:
        if key not in reference or key not in candidate or reference[key] != candidate[key]:
            raise ValueError('physical metadata differs or absent: ' + key)
    if reference.get('mode') != 'kernel' or candidate.get('mode') != candidate_mode:
        raise ValueError('physical result mode does not match the selected comparison')
    for result in (reference, candidate):
        if result.get('source_contact_enabled') is not False:
            raise ValueError('contacts-on cases are outside this comparison')
        if result.get('forward_only_baseline') is not False:
            raise ValueError('both physical baselines must use forward-only off')
        validate_energy(result)
    degree = reference.get('degree')
    if candidate.get('degree') != degree:
        raise ValueError('reference and candidate degrees differ')
    if candidate_mode == 'kernel-probe':
        if not isinstance(degree, int) or isinstance(degree, bool) or not 1 <= degree <= 1536:
            raise ValueError('probe degree is outside the supported angular-library range')
    elif degree not in SUPPORTED_DEGREES:
        raise ValueError('full global comparisons are restricted to L=3 and L=12')
    if reference.get('sem_degree') != SUPPORTED_SEM_DEGREE or candidate.get('sem_degree') != SUPPORTED_SEM_DEGREE:
        raise ValueError('global comparison cases require SEM degree 6')
    if reference.get('native_element_node_sides') != TOTAL_NATIVE_SIDES:
        raise ValueError('result metadata must report 180 total native element/node sides')
    if candidate.get('native_element_node_sides') != TOTAL_NATIVE_SIDES:
        raise ValueError('result metadata must report 180 total native element/node sides')
    if reference.get('frequency_bins') != PHYSICAL_OUTPUT_BINS or candidate.get('frequency_bins') != PHYSICAL_OUTPUT_BINS:
        raise ValueError('physical result must report 52 output frequency bins')


def validate_global_manifest(manifest: dict, result: dict, data: Path,
                             receipt: tuple, probe64: bool) -> dict:
    need_receipted(receipt, GLOBAL_FILES)
    if manifest.get('mode') != 'compact_response_global_gsh_fft_frequency_accumulator':
        raise ValueError('unexpected compact global sidecar mode')
    degree = _int(manifest.get('degree'), 'degree', positive=True)
    if degree != result.get('degree'):
        raise ValueError('global sidecar degree differs from physical result')
    ntheta = _int(manifest.get('n_theta'), 'n_theta', positive=True)
    nphi = _int(manifest.get('n_phi'), 'n_phi', positive=True)
    full_angles = _int(manifest.get('full_angular_points'), 'full_angular_points', positive=True)
    sampled_angles = _int(manifest.get('sampled_angular_points'), 'sampled_angular_points', positive=True)
    if ntheta != 2 * degree + 1 or nphi < 4 * degree + 1 or full_angles != ntheta * nphi:
        raise ValueError('global angular count does not resolve declared degree bandwidth')
    angular_rule = manifest.get('angular_rule', '')
    if not ('ForBand(L,2,2.0)' in angular_rule and 'FastFFTSize(4L+1)' in angular_rule):
        raise ValueError('global sidecar does not declare the accepted complete-band angular rule')
    complete = manifest.get('angular_subset_is_complete')
    if not isinstance(complete, bool) or complete != (sampled_angles == full_angles):
        raise ValueError('global complete/prefix angular flag disagrees with counts')
    if sampled_angles > full_angles:
        raise ValueError('sampled global angular count exceeds full grid')
    if probe64:
        if result.get('mode') != 'kernel-probe' or sampled_angles != 64 or complete:
            raise ValueError('explicit probe64 mode requires exactly a 64-angle kernel-probe result')
    elif result.get('mode') != 'kernel' or not complete:
        raise ValueError('full global validation requires a complete kernel grid; use --probe64 only for a 64-angle probe')
    nsides = _int(manifest.get('radial_side_count'), 'radial_side_count', positive=True)
    if nsides > GLOBAL_SELECTED_SIDES_MAX:
        raise ValueError('combined-v3 global manifest exceeds the eight-side selection cap')
    indices = manifest.get('native_side_indices')
    metadata = manifest.get('native_side_metadata')
    if (not isinstance(indices, list) or len(indices) != nsides
            or not isinstance(metadata, list) or len(metadata) != nsides):
        raise ValueError('global native-side index/metadata count mismatch')
    if any(not isinstance(i, int) or isinstance(i, bool) for i in indices):
        raise ValueError('global native-side indices must be integers')
    if len(set(indices)) != nsides:
        raise ValueError('duplicate selected native-side index')
    header, native = read_native(data / 'radial_nodes.csv')
    if len(native) != POSITIVE_RADIUS_SIDES:
        raise ValueError('candidate positive-radius native CSV must contain 179 sides')
    selected_count = _int(result.get('global_consumer_radial_side_count'),
                          'global_consumer_radial_side_count', positive=True)
    if selected_count > GLOBAL_SELECTED_SIDES_MAX:
        raise ValueError('combined-v3 global candidate exceeds the eight-side selection cap')
    if nsides != selected_count:
        raise ValueError('result and sidecar selected native-side counts differ')
    if result.get('global_consumer_side_selection_explicit') is not True:
        raise ValueError('global candidate must declare explicit native-side selection')
    if any(i < 0 or i >= len(native) for i in indices):
        raise ValueError('global native-side index is outside candidate all-native CSV')
    scale = result.get('model_length_scale_m')
    if not isinstance(scale, (int, float)) or not math.isfinite(scale) or scale <= 0:
        raise ValueError('invalid model length scale')
    identities = []
    for ordinal, (index, side) in enumerate(zip(indices, metadata)):
        if not isinstance(side, dict):
            raise ValueError('invalid global native-side metadata entry')
        row = native[index]
        for key in ('element', 'node'):
            _int(side.get(key), 'native_side_metadata.' + key)
        for key in ('radius_m', 'radius_model', 'width_m'):
            value = side.get(key)
            if (not isinstance(value, (int, float)) or isinstance(value, bool)
                    or not math.isfinite(value) or value <= 0):
                raise ValueError('invalid global native-side metadata value: ' + key)
        if not isinstance(side.get('fluid'), bool):
            raise ValueError('global native-side fluid metadata must be boolean')
        expected = {
            'element': int(row[0]), 'node': int(row[1]), 'radius_m': float(row[2]),
            'width_m': float(row[3]), 'fluid': bool(int(row[4])),
        }
        for key, value in expected.items():
            if key not in side or side[key] != value:
                raise ValueError('global side metadata disagrees with candidate all-native CSV at output side ' + str(ordinal) + ': ' + key)
        model_radius = side.get('radius_model')
        if (not isinstance(model_radius, (int, float)) or not math.isfinite(model_radius)
                or model_radius <= 0
                or not math.isclose(model_radius * scale, expected['radius_m'], rel_tol=2e-15, abs_tol=0.0)):
            raise ValueError('global model/SI radius identity mismatch at output side ' + str(ordinal))
        identities.append((expected['element'], expected['node']))
    if len(set(identities)) != nsides:
        raise ValueError('duplicate selected native side identity')
    if manifest.get('source_contacts_included') is not False or manifest.get('kernel_is_full_earth_volume_integral') is not False:
        raise ValueError('global output must be contacts-off and explicitly non-volume-integrated')
    if not ('radial-side-major' in manifest.get('output_layout', '')
            and 'theta-major' in manifest.get('output_layout', '')
            and 'A,C,F,L,N,kappa' in manifest.get('output_layout', '')):
        raise ValueError('unrecognized global radial-side-major output layout')
    expected_bytes = nsides * sampled_angles * len(CHANNELS) * 8
    if manifest.get('output_bytes') != expected_bytes or (data / 'kernel.bin').stat().st_size != expected_bytes:
        raise ValueError('global kernel byte count disagrees with side and angle counts')
    if manifest.get('binary_file') != 'kernel.bin':
        raise ValueError('global sidecar names an unexpected binary file')
    if manifest.get('frequency_count') != SPOOL_BINS:
        raise ValueError('global spool must contain ' + str(SPOOL_BINS) + ' frequency records')
    records = manifest.get('frequency_records')
    validate_frequency_records(records, SPOOL_BINS, PHYSICAL_OUTPUT_BINS)
    spool_bytes = (data / 'pole_response_spool.bin').stat().st_size
    if result.get('response_spool_bytes') != spool_bytes:
        raise ValueError('result spool byte count differs from receipted spool file')
    return {
        'native_side_indices': indices,
        'native_side_metadata': metadata,
        'native_identities': identities,
        'native_rows': native,
        'n_theta': ntheta,
        'n_phi': nphi,
        'full_angular_points': full_angles,
        'sampled_angular_points': sampled_angles,
        'radial_side_count': nsides,
        'spool_frequency_count': len(records),
    }


def validate_frequency_records(records: Any, expected_count: int,
                               output_bins: int = PHYSICAL_OUTPUT_BINS) -> None:
    if not isinstance(records, list) or len(records) != expected_count:
        raise ValueError('global spool frequency audit must contain the expected record count')
    for i, record in enumerate(records):
        if not isinstance(record, dict) or record.get('spool_index') != i or record.get('original_bin') != i:
            raise ValueError('global frequency records are not contiguous original bins')
        z = _finite_pair(record.get('z'), 'z')
        mixing = _finite_pair(record.get('adjoint_mixing'), 'adjoint_mixing')
        residual = record.get('receiver_basis_reconstruction_relative')
        if z[1] == 0.0:
            raise ValueError('global frequency record has zero damping')
        q_is_zero = mixing[0] == 0.0 and mixing[1] == 0.0
        if i == output_bins:
            if not q_is_zero:
                raise ValueError('out-of-band guard record must have exactly zero complex q weight')
        elif i < output_bins:
            if q_is_zero:
                raise ValueError('in-band frequency record has zero complex q weight')
        else:
            raise ValueError('global spool record lies beyond the designated guard bin')
        if (not isinstance(residual, (int, float)) or not math.isfinite(residual)
                or residual < 0.0 or residual > 1e-12):
            raise ValueError('global receiver basis reconstruction exceeds 1e-12')


def validate_global_result(result: dict, profile: dict) -> None:
    if result.get('response_record_spool_used') is not True:
        raise ValueError('candidate did not use compact response spool')
    if result.get('response_reuse_used') is not False:
        raise ValueError('candidate response cache must be disabled for spool comparison')
    if result.get('response_cache_entries') != 0 or result.get('response_cache_bytes') != 0:
        raise ValueError('candidate response cache retained entries or bytes')
    expected_factorizations = SPOOL_BINS * (2 * result['degree'] + 1)
    actual = profile['counters'].get('sparse_factorization_calls')
    if actual != expected_factorizations:
        raise ValueError('sparse factorization count differs from one pass: expected '
                         + str(expected_factorizations) + ', got ' + str(actual))


def validate_worker_count(result: dict, manifest: dict, timing: dict,
                          expected_requested: int | None = None) -> dict:
    requested = _int(result.get('global_consumer_workers_requested'),
                     'global_consumer_workers_requested', positive=True)
    if requested > 48:
        raise ValueError('requested global worker count exceeds the producer limit')
    selected = _int(result.get('global_consumer_radial_side_count'),
                    'global_consumer_radial_side_count', positive=True)
    actual = _int(timing.get('actual_workers'), 'actual manifest worker count', positive=True)
    if expected_requested is not None and requested != expected_requested:
        raise ValueError('global worker request differs from the required parallel comparison')
    expected_actual = min(requested, selected, 48)
    if actual != expected_actual or manifest.get('workers') != actual:
        raise ValueError('actual global worker count differs from requested/selected-side count')
    return {'requested': requested, 'actual': actual, 'selected_sides': selected}


def validate_parallel_worker_pair(reference_result: dict, reference_manifest: dict,
                                  reference_timing: dict, candidate_result: dict,
                                  candidate_manifest: dict, candidate_timing: dict) -> dict:
    reference = validate_worker_count(reference_result, reference_manifest,
                                      reference_timing, expected_requested=1)
    candidate = validate_worker_count(candidate_result, candidate_manifest,
                                      candidate_timing, expected_requested=4)
    if (reference['selected_sides'] != GLOBAL_SELECTED_SIDES_MAX
            or candidate['selected_sides'] != GLOBAL_SELECTED_SIDES_MAX):
        raise ValueError('w1-versus-w4 gate requires eight selected native sides')
    return {'reference': reference, 'candidate': candidate}


def validate_global_outputs_equal(reference_data: Path, candidate_data: Path,
                                  reference_result: dict,
                                  candidate_result: dict) -> dict:
    exact = {}
    for name in ('kernel.bin', 'trace.csv', 'spectrum.csv'):
        exact[name] = pair.sha(reference_data / name) == pair.sha(candidate_data / name)
        if not exact[name]:
            raise ValueError('global-to-global ' + name + ' files are not byte-exact')
    scalar_keys = ('energy_m2_s', 'energy_adjoint_qdotp')
    scalar_exact = {key: reference_result.get(key) == candidate_result.get(key)
                    for key in scalar_keys}
    if not all(scalar_exact.values()):
        raise ValueError('global-to-global energy/qdot scalars are not exact')
    return {'byte_exact': exact, 'energy_qdot_exact': scalar_exact}


def validate_direct_domain(domain: dict, result: dict, data: Path, receipt: tuple) -> dict:
    need_receipted(receipt, REFERENCE_FILES + ('angular_grid.csv',))
    if result.get('mode') != 'kernel' or result.get('source_contact_enabled') is not False:
        raise ValueError('reference must be a contacts-off full kernel case')
    if domain.get('source_contact_enabled') is not False:
        raise ValueError('reference domain must be contacts-off')
    if domain.get('actual_contact_binary_bytes') != 0 or domain.get('full_target_contact_binary_bytes') != 0:
        raise ValueError('reference declares source-contact output bytes')
    degree = result['degree']
    nt = _int(domain.get('n_theta'), 'reference n_theta', positive=True)
    np_ = _int(domain.get('n_phi'), 'reference n_phi', positive=True)
    angles = nt * np_
    nsides = _int(domain.get('selected_radial_trace_count'), 'reference side count', positive=True)
    if domain.get('available_radial_trace_sides') != POSITIVE_RADIUS_SIDES:
        raise ValueError('reference domain must report all 179 available positive-radius sides')
    if nt != 2 * degree + 1 or np_ < 4 * degree + 1:
        raise ValueError('reference angular grid does not resolve degree bandwidth')
    if domain.get('full_target_angular_points_per_trace') != angles:
        raise ValueError('reference full angular count disagrees with dimensions')
    if domain.get('sampled_angular_points_per_trace') != angles or domain.get('sampled_angle_subset_is_complete') is not True:
        raise ValueError('reference direct output is not a complete angular grid')
    if domain.get('actual_sampled_point_count') != nsides * angles:
        raise ValueError('reference point count disagrees with radial/angular product')
    if domain.get('actual_kernel_binary_bytes') != nsides * angles * len(CHANNELS) * 8:
        raise ValueError('reference binary byte count disagrees with dimensions')
    layout = domain.get('binary_layout', '')
    if not ('angle_index outer' in layout and 'selected radial side inner' in layout
            and 'A,C,F,L,N,kappa' in layout):
        raise ValueError('unrecognized direct reference angle-major layout')
    if 'ForBand(L,2,2.0)' not in domain.get('angular_rule', ''):
        raise ValueError('reference does not declare the matching ForBand angular rule')
    grid_header, grid = pair.csv_array(data / 'angular_grid.csv')
    if (grid_header != ('theta_index', 'theta_rad', 'theta_weight') or len(grid) != nt
            or not np.array_equal(grid[:, 0], np.arange(nt))
            or not np.all((grid[:, 1] > 0) & (grid[:, 1] < np.pi))
            or not np.all(np.diff(grid[:, 1]) > 0) or not np.all(grid[:, 2] > 0)
            or abs(float(np.sum(grid[:, 2])) - 2.0) > 1e-12):
        raise ValueError('invalid direct reference theta quadrature')
    if (data / 'kernel.bin').stat().st_size != nsides * angles * len(CHANNELS) * 8:
        raise ValueError('reference kernel file size disagrees with declared grid')
    header, selected = pair.csv_array(data / 'kernel_radial_nodes.csv')
    if header != ('radial_index', 'element', 'node', 'radius_m', 'element_width_m', 'fluid'):
        raise ValueError('unknown reference selected-side table schema')
    if len(selected) != nsides or not np.array_equal(selected[:, 0], np.arange(nsides)):
        raise ValueError('reference selected-side table count/order mismatch')
    ids = [tuple(row[1:3].astype(np.int64)) for row in selected]
    if len(ids) != len(set(ids)):
        raise ValueError('duplicate reference selected native side')
    _, allnative = read_native(data / 'radial_nodes.csv')
    if len(allnative) != POSITIVE_RADIUS_SIDES:
        raise ValueError('reference positive-radius native CSV must contain 179 sides')
    if nsides != POSITIVE_RADIUS_SIDES:
        raise ValueError('direct reference must include all 179 positive-radius native sides')
    by_id = {(int(row[0]), int(row[1])): row for row in allnative}
    selected_map = {}
    for ordinal, row in enumerate(selected):
        identity = tuple(row[1:3].astype(np.int64))
        native = by_id.get(identity)
        if native is None or not np.array_equal(row[1:], native):
            raise ValueError('reference selected side differs from state-hash-bound native CSV')
        selected_map[identity] = ordinal
    return {'n_theta': nt, 'n_phi': np_, 'angular_points': angles,
            'radial_side_count': nsides, 'selected_rows': selected,
            'selected_map': selected_map, 'allnative': allnative}


def compare_csv(reference_data: Path, candidate_data: Path, receipts: tuple,
                names: tuple[str, ...]) -> tuple[dict, dict]:
    comparisons, exact = {}, {}
    for name in names:
        need_receipted(receipts[0], (name,))
        need_receipted(receipts[1], (name,))
        ha, a = pair.csv_array(reference_data / name)
        hb, b = pair.csv_array(candidate_data / name)
        if ha != hb:
            raise ValueError('CSV headers differ: ' + name)
        comparisons[name] = pair.compare_columns(a, b)
        exact[name] = pair.sha(reference_data / name) == pair.sha(candidate_data / name)
    return comparisons, exact


def validate_selected_identity_subset(selected: list[tuple[int, int]],
                                      reference: set[tuple[int, int]]) -> None:
    if not 1 <= len(selected) <= GLOBAL_SELECTED_SIDES_MAX:
        raise ValueError('global output side count must be within the one-to-eight selection cap')
    if len(set(selected)) != len(selected):
        raise ValueError('duplicate selected native-side identity')
    if not set(selected).issubset(reference):
        raise ValueError('global native-side selection is not a subset of the direct all-native reference')


def compare_kernel(reference_path: Path, candidate_path: Path,
                   reference_sides: dict, candidate_sides: dict,
                   angles: int) -> list[dict]:
    validate_selected_identity_subset(candidate_sides['native_identities'],
                                      set(reference_sides['selected_map']))
    ref_angles = reference_sides['angular_points']
    if angles > ref_angles:
        raise ValueError('candidate angular prefix exceeds reference grid')
    ref = np.memmap(reference_path, dtype=np.float64, mode='r').reshape(
        ref_angles, reference_sides['radial_side_count'], len(CHANNELS))
    cand = np.memmap(candidate_path, dtype=np.float64, mode='r').reshape(
        candidate_sides['radial_side_count'], angles, len(CHANNELS))
    ref_index = reference_sides['selected_map']
    maximum = np.zeros(6)
    scale = np.zeros(6)
    chunk = 2048
    for side_index, identity in enumerate(candidate_sides['native_identities']):
        rside = ref_index[identity]
        for begin in range(0, angles, chunk):
            end = min(angles, begin + chunk)
            a = np.asarray(ref[begin:end, rside, :])
            b = np.asarray(cand[side_index, begin:end, :])
            if not np.isfinite(a).all() or not np.isfinite(b).all():
                raise ValueError('nonfinite kernel binary value')
            maximum = np.maximum(maximum, np.max(np.abs(a - b), axis=0))
            scale = np.maximum(scale, np.maximum(np.max(np.abs(a), axis=0), np.max(np.abs(b), axis=0)))
    ratio = np.divide(maximum, scale, out=np.zeros_like(maximum), where=scale != 0)
    return [{'column': i, 'channel': CHANNELS[i], 'max_abs_difference': float(maximum[i]),
             'scale': float(scale[i]), 'scaled_max_difference': float(ratio[i]),
             'passed': bool(ratio[i] <= LIMIT)} for i in range(6)]


def validate_worker_timings(manifest: dict) -> dict:
    elapsed = manifest.get('elapsed_s')
    if not isinstance(elapsed, (int, float)) or not math.isfinite(elapsed) or elapsed < 0:
        raise ValueError('invalid global consumer elapsed_s')
    timings = manifest.get('timings')
    if timings is None:
        workers = manifest.get('worker_owned_side_seconds')
        if not isinstance(workers, list) or any(not isinstance(x, (int, float)) or not math.isfinite(x) or x < 0 for x in workers):
            raise ValueError('missing or invalid global worker timing data')
        return {'total_wall_s': float(elapsed), 'worker_execution_wall_s': None,
                'serial_context_setup_wall_s': None,
                'worker_owned_side_seconds_sum_overlapping': float(sum(workers)),
                'actual_workers': len(workers), 'timing_schema': 'global-v2'}
    fields = ('total_wall_s', 'worker_execution_wall_s', 'serial_context_setup_wall_s')
    if any(not isinstance(timings.get(k), (int, float)) or not math.isfinite(timings[k]) or timings[k] < 0 for k in fields):
        raise ValueError('invalid global-v3 wall timing fields')
    if abs(float(timings['total_wall_s']) - float(elapsed)) > max(1e-8, abs(float(elapsed)) * 2e-12):
        raise ValueError('global total wall timing differs from elapsed_s')
    if (float(timings['serial_context_setup_wall_s'])
            + float(timings['worker_execution_wall_s'])
            > float(timings['total_wall_s']) + max(1e-8, abs(float(elapsed)) * 2e-12)):
        raise ValueError('serial setup plus worker wall exceeds total consumer wall')
    stage_sums = timings.get('worker_stage_elapsed_sums')
    stage_keys = ('spool_reads', 'coefficient_rotation_layout', 'gsh_inverse_fft',
                  'point_kernel_contractions', 'frequency_accumulation', 'output_writes')
    if (not isinstance(stage_sums, dict) or any(
            not isinstance(stage_sums.get(k), (int, float)) or isinstance(stage_sums.get(k), bool)
            or not math.isfinite(stage_sums[k]) or stage_sums[k] < 0 for k in stage_keys)):
        raise ValueError('invalid global-v3 worker_stage_elapsed_sums')
    workers = timings.get('workers')
    declared_workers = manifest.get('workers')
    if (not isinstance(workers, list) or not workers
            or declared_workers != len(workers)):
        raise ValueError('global-v3 worker timing list is empty')
    worker_elapsed_sum = timings.get('worker_elapsed_sum_s')
    if (not isinstance(worker_elapsed_sum, (int, float)) or isinstance(worker_elapsed_sum, bool)
            or not math.isfinite(worker_elapsed_sum) or worker_elapsed_sum < 0):
        raise ValueError('missing or invalid global-v3 worker_elapsed_sum_s')
    stage_fields = {
        'spool_reads': 'spool_reads_s',
        'coefficient_rotation_layout': 'coefficient_rotation_layout_s',
        'gsh_inverse_fft': 'gsh_inverse_fft_s',
        'point_kernel_contractions': 'point_kernel_contractions_s',
        'frequency_accumulation': 'frequency_accumulation_s',
        'output_writes': 'output_writes_s',
    }
    for worker in workers:
        if not isinstance(worker, dict) or any(
            not isinstance(worker.get(k), (int, float)) or not math.isfinite(worker[k]) or worker[k] < 0
            for k in ('side_elapsed_s', 'spool_reads_s', 'coefficient_rotation_layout_s',
                      'gsh_inverse_fft_s', 'point_kernel_contractions_s',
                      'frequency_accumulation_s', 'output_writes_s')):
            raise ValueError('invalid per-worker elapsed timing record')
    def rollup_matches(actual: float, expected: float) -> bool:
        return abs(actual - expected) <= max(1e-12, max(abs(actual), abs(expected)) * 2e-12)
    measured_side_sum = sum(float(worker['side_elapsed_s']) for worker in workers)
    if not rollup_matches(float(worker_elapsed_sum), measured_side_sum):
        raise ValueError('worker_elapsed_sum_s disagrees with per-worker side_elapsed_s')
    for aggregate_key, worker_key in stage_fields.items():
        measured = sum(float(worker[worker_key]) for worker in workers)
        if not rollup_matches(float(stage_sums[aggregate_key]), measured):
            raise ValueError('worker_stage_elapsed_sums disagrees with per-worker ' + worker_key)
    return {'total_wall_s': float(timings['total_wall_s']),
            'worker_execution_wall_s': float(timings['worker_execution_wall_s']),
            'serial_context_setup_wall_s': float(timings['serial_context_setup_wall_s']),
            'worker_elapsed_sum_overlapping': float(worker_elapsed_sum),
            'actual_workers': len(workers),
            'timing_schema': 'global-v3'}


def _unchanged(roots: list[Path], jobs: list[str], before: list[tuple]) -> list[tuple]:
    after = [pair.verified_job(root, job) for root, job in zip(roots, jobs)]
    for initial, current in zip(before, after):
        if initial[0] != current[0] or initial[1] != current[1]:
            raise ValueError('completed-job receipt or verified artifact map changed during comparison')
    return after


def assess(reference_root: Path, reference_job: str, candidate_root: Path,
           candidate_job: str, probe64: bool = False,
           global_reference: bool = False) -> dict:
    roots = [reference_root, candidate_root]
    jobs = [reference_job, candidate_job]
    receipts = [pair.verified_job(root, job) for root, job in zip(roots, jobs)]
    dirs = [root / job / 'data' for root, job in zip(roots, jobs)]
    need_receipted(receipts[1], GLOBAL_FILES)
    candidate_result = pair.load(dirs[1] / 'result.json')
    candidate_profile = pair.profile(roots[1] / candidate_job / receipts[1][0]['log'])
    if candidate_profile is None:
        raise ValueError('candidate has no MATHPROFILE record')
    if candidate_profile['solver'] != 'dspec' or candidate_profile['objective'] != -1:
        raise ValueError('unexpected candidate profiler identity')
    validate_global_result(candidate_result, candidate_profile)
    candidate_manifest = pair.load(dirs[1] / 'kernel.bin.json')
    candidate_side = validate_global_manifest(candidate_manifest, candidate_result,
                                                dirs[1], receipts[1], probe64)
    candidate_timing = validate_worker_timings(candidate_manifest)
    candidate_workers = validate_worker_count(candidate_result, candidate_manifest,
                                              candidate_timing)
    reference_result = pair.load(dirs[0] / 'result.json')
    expected_mode = 'kernel-probe' if probe64 else 'kernel'
    validate_physical_pair(reference_result, candidate_result,
                           expected_mode if not global_reference else 'kernel')
    if probe64 and global_reference:
        raise ValueError('--global-reference requires complete kernel outputs')
    reference_profile = pair.profile(roots[0] / reference_job / receipts[0][0]['log'])
    if reference_profile is None:
        raise ValueError('reference has no MATHPROFILE record')
    if global_reference:
        need_receipted(receipts[0], GLOBAL_FILES)
        ref_manifest = pair.load(dirs[0] / 'kernel.bin.json')
        ref_side = validate_global_manifest(ref_manifest, reference_result,
                                             dirs[0], receipts[0], probe64)
        validate_global_result(reference_result, reference_profile)
        reference_timing = validate_worker_timings(ref_manifest)
        worker_pair = validate_parallel_worker_pair(
            reference_result, ref_manifest, reference_timing,
            candidate_result, candidate_manifest, candidate_timing)
        if ref_side['native_identities'] != candidate_side['native_identities']:
            raise ValueError('global-to-global output side ordering differs')
        output_parity = validate_global_outputs_equal(
            dirs[0], dirs[1], reference_result, candidate_result)
        profiles = [reference_profile, candidate_profile]
        after = _unchanged(roots, jobs, receipts)
        return {'status': 'PASS', 'comparison_mode': 'global-to-global-byte-exact',
                'tolerance': LIMIT, 'reference_job': reference_job,
                'candidate_job': candidate_job,
                **output_parity,
                'worker_counts': worker_pair,
                'worker_timings': {'reference': reference_timing,
                                   'candidate': candidate_timing},
                'candidate_global_timing': candidate_timing,
                'profiles': profiles,
                'state_sha256_at_read': [r[2] for r in after],
                'verified_artifact_sha256': [r[1] for r in after]}

    need_receipted(receipts[0], REFERENCE_FILES + ('angular_grid.csv',))
    reference_domain = pair.load(dirs[0] / 'kernel_domain.json')
    reference_side = validate_direct_domain(reference_domain, reference_result,
                                            dirs[0], receipts[0])
    if (candidate_side['n_theta'], candidate_side['n_phi']) != (reference_side['n_theta'], reference_side['n_phi']):
        raise ValueError('global and direct reference angular dimensions differ')
    if reference_side['radial_side_count'] != POSITIVE_RADIUS_SIDES:
        raise ValueError('direct reference must contain all 179 positive-radius native sides')
    validate_selected_identity_subset(candidate_side['native_identities'],
                                      set(reference_side['selected_map']))
    _, candidate_native = read_native(dirs[1] / 'radial_nodes.csv')
    _, reference_native = read_native(dirs[0] / 'radial_nodes.csv')
    if not np.array_equal(candidate_native, reference_native):
        raise ValueError('candidate all-native CSV differs from direct reference native table')
    if candidate_side['sampled_angular_points'] > reference_side['angular_points']:
        raise ValueError('global angular prefix exceeds the direct reference full grid')
    comparisons, exact = compare_csv(dirs[0], dirs[1], receipts,
                                    ('trace.csv', 'spectrum.csv'))
    comparisons['energy'] = pair.compare_columns(
        np.array([[reference_result['energy_m2_s'], reference_result['energy_adjoint_qdotp']]]),
        np.array([[candidate_result['energy_m2_s'], candidate_result['energy_adjoint_qdotp']]]))
    comparisons['kernel.bin'] = compare_kernel(
        dirs[0] / 'kernel.bin', dirs[1] / 'kernel.bin', reference_side,
        candidate_side, candidate_side['sampled_angular_points'])
    exact['kernel.bin'] = pair.sha(dirs[0] / 'kernel.bin') == pair.sha(dirs[1] / 'kernel.bin')
    profiles = [reference_profile, candidate_profile]
    after = _unchanged(roots, jobs, receipts)
    passed = all(row['passed'] for rows in comparisons.values() for row in rows)
    return {'status': 'PASS' if passed else 'FAIL',
            'comparison_mode': 'direct-fullgrid-to-global-spool-prefix' if probe64 else 'direct-fullgrid-to-global-spool-full',
            'tolerance': LIMIT, 'reference_job': reference_job,
            'candidate_job': candidate_job, 'comparison': comparisons,
            'byte_exact_files': exact, 'candidate_global_timing': candidate_timing,
            'candidate_workers': candidate_workers,
            'candidate_native_side_indices': candidate_side['native_side_indices'],
            'angular_points_compared_per_side': candidate_side['sampled_angular_points'],
            'profiles': profiles,
            'state_sha256_at_read': [r[2] for r in after],
            'verified_artifact_sha256': [r[1] for r in after]}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('reference-root', 'reference-job', 'candidate-root', 'candidate-job', 'output'):
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--probe64', action='store_true',
                        help='explicitly allow candidate kernel-probe with exactly 64 prefix angles')
    parser.add_argument('--global-reference', action='store_true',
                        help='compare two compact global-spool jobs byte-exactly instead of using direct reference')
    args = parser.parse_args()
    try:
        report = assess(Path(args.reference_root), args.reference_job,
                        Path(args.candidate_root), args.candidate_job,
                        args.probe64, args.global_reference)
    except Exception as exc:
        report = {'status': 'FAIL', 'error_type': type(exc).__name__, 'error': str(exc)}
    with Path(args.output).open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write('\n')
    return 0 if report['status'] == 'PASS' else 1


if __name__ == '__main__':
    raise SystemExit(main())
