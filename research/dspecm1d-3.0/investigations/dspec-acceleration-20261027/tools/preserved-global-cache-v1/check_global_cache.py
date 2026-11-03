#!/usr/bin/env python3
"""Read-only comparison of completed full-grid global cache/storage jobs."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

import check_pair as pair
import check_global_pair as spool

CHANNELS = spool.CHANNELS
DEGREES = (3, 12, 128, 384, 1536)
RESOURCE_FIELDS = (
    'shared_wigner_table_bytes', 'known_working_bytes_per_worker',
    'estimated_working_bytes_per_worker', 'estimated_total_working_bytes',
    'spool_offset_validation_scratch_bytes', 'output_bytes',
)


def validate_physical_pair(reference: dict, candidate: dict) -> int:
    for key in pair.PHYSICAL:
        if key not in reference or key not in candidate or reference[key] != candidate[key]:
            raise ValueError('physical metadata differs or absent: ' + key)
    degree = reference.get('degree')
    if degree not in DEGREES or candidate.get('degree') != degree:
        raise ValueError('cache comparison degree must be one of L=3, 12, 128, 384, 1536')
    if reference.get('mode') != 'kernel' or candidate.get('mode') != 'kernel':
        raise ValueError('cache comparison requires full kernel mode')
    if reference.get('sem_degree') != 6 or candidate.get('sem_degree') != 6:
        raise ValueError('cache comparison requires SEM degree 6')
    if reference.get('native_element_node_sides') != spool.TOTAL_NATIVE_SIDES or candidate.get('native_element_node_sides') != spool.TOTAL_NATIVE_SIDES:
        raise ValueError('physical results must report 180 total element/node sides')
    if reference.get('frequency_bins') != spool.PHYSICAL_OUTPUT_BINS or candidate.get('frequency_bins') != spool.PHYSICAL_OUTPUT_BINS:
        raise ValueError('physical results must report 52 output bins')
    for result in (reference, candidate):
        if result.get('source_contact_enabled') is not False or result.get('forward_only_baseline') is not False:
            raise ValueError('cache comparison requires contacts off and forward-only off')
        spool.validate_energy(result)
    return degree


def _command_tokens(receipt: tuple) -> list[str]:
    command = receipt[0].get('command')
    if not isinstance(command, list) or not command or any(not isinstance(x, str) for x in command):
        raise ValueError('receipt has no tokenized command provenance')
    return command


def _option(command: list[str], name: str) -> str:
    positions = [i for i, token in enumerate(command) if token == name]
    if len(positions) != 1 or positions[0] + 1 >= len(command):
        raise ValueError('receipt command must contain exactly one ' + name)
    return command[positions[0] + 1]


def validate_invocation(result: dict, manifest: dict, timing: dict, receipt: tuple,
                        expected_workers: int, expected_storage: str,
                        selected_indices: list[int]) -> dict:
    if expected_storage not in ('generated', 'stored'):
        raise ValueError('expected storage must be generated or stored')
    if expected_workers not in (1, 2, 4, 8):
        raise ValueError('expected worker count must be one of 1, 2, 4, 8')
    requested = spool._int(result.get('global_consumer_workers_requested'),
                           'global_consumer_workers_requested', positive=True)
    actual = spool._int(manifest.get('workers'), 'manifest workers', positive=True)
    if requested != expected_workers or actual != expected_workers or timing.get('actual_workers') != expected_workers:
        raise ValueError('requested/result/manifest worker counts differ from explicit expected count')
    if result.get('global_consumer_wigner_storage_requested') != expected_storage:
        raise ValueError('result Wigner storage mode differs from explicit expectation')
    if manifest.get('wigner_storage_policy') != expected_storage:
        raise ValueError('manifest Wigner storage mode differs from explicit expectation')
    command = _command_tokens(receipt)
    if command.count('kernel') != 1:
        raise ValueError('receipted command is not a kernel-mode invocation')
    if _option(command, '--lmax') != str(result.get('degree')):
        raise ValueError('receipted command degree differs from physical result')
    if _option(command, '--nq') != '6' or _option(command, '--radial-traces') != '0':
        raise ValueError('receipted command does not use the required SEM-6/all-native input')
    for option, value in (('--dc', 'include'), ('--contacts', 'off'),
                          ('--pole-only', 'on'), ('--forward-only', 'off'),
                          ('--response-reuse', 'off'), ('--response-record-spool', 'on'),
                          ('--global-workers', str(expected_workers)),
                          ('--global-wigner-storage', expected_storage)):
        if _option(command, option) != value:
            raise ValueError('receipted command option differs from expected: ' + option)
    raw_indices = _option(command, '--global-radial-sides')
    try:
        command_indices = [int(x) for x in raw_indices.split(',')]
    except ValueError as exc:
        raise ValueError('invalid receipted selected-side list') from exc
    if command_indices != selected_indices:
        raise ValueError('receipted selected-side order differs from global manifest')
    if result.get('global_consumer_radial_side_count') != len(selected_indices):
        raise ValueError('result selected-side count differs from global manifest')
    # Bind the producer's result counters to the independently validated sidecar.
    # This catches stale or mixed result/manifest pairs even when their files are
    # individually well-formed.
    mappings = (
        ('global_consumer_full_angular_points', 'full_angular_points'),
        ('global_consumer_frequency_count', 'frequency_count'),
        ('global_consumer_radial_side_count', 'radial_side_count'),
        ('response_spool_kernel_output_bytes', 'output_bytes'),
        ('response_spool_kernel_sampled_angles', 'sampled_angular_points'),
        ('global_consumer_shared_wigner_bytes', 'shared_wigner_table_bytes'),
        ('global_consumer_estimated_working_bytes', 'estimated_total_working_bytes'),
        ('global_consumer_same_spin_batch_size', 'same_spin_channel_batch_size'),
    )
    for result_key, manifest_key in mappings:
        if result.get(result_key) != manifest.get(manifest_key):
            raise ValueError('result/manifest producer field mismatch: ' + result_key)
    if (result.get('global_consumer_angular_subset_complete') is not True
            or result.get('global_consumer_full_transforms_performed') is not True):
        raise ValueError('result does not confirm complete-grid transforms')
    expected_transforms = manifest.get('frequency_count') * len(selected_indices)
    executed = result.get('global_consumer_frequency_side_transforms_executed')
    skipped = result.get('global_consumer_frequency_side_zero_mixing_skipped')
    raw_timing = manifest.get('timings', {})
    timing_executed = raw_timing.get('frequency_side_transforms_executed')
    timing_skipped = raw_timing.get('frequency_side_exact_zero_mixing_skipped')
    if (not isinstance(executed, int) or isinstance(executed, bool)
            or not isinstance(skipped, int) or isinstance(skipped, bool)
            or executed + skipped != expected_transforms
            or executed != timing_executed or skipped != timing_skipped):
        raise ValueError('result/manifest frequency-side transform counters disagree')
    shared = manifest.get('shared_wigner_table_bytes')
    if ((expected_storage == 'generated' and shared != 0)
            or (expected_storage == 'stored' and (not isinstance(shared, int) or shared <= 0))):
        raise ValueError('shared Wigner allocation disagrees with requested storage mode')
    if not 1 <= len(selected_indices) <= spool.GLOBAL_SELECTED_SIDES_MAX:
        raise ValueError('selected native-side count is outside the supported 1..8 range')
    return {'requested_workers': requested, 'actual_workers': actual,
            'storage_mode': expected_storage, 'selected_side_count': len(selected_indices),
            'selected_side_indices': selected_indices}


def _final_timing_adapter(manifest: dict) -> dict:
    """Adapt final-v1 setup key to the already-reviewed v3 timing validator."""
    normalized = dict(manifest)
    timings = dict(manifest.get('timings', {}))
    if ('serial_context_setup_wall_s' not in timings
            and 'serial_grid_and_context_setup_wall_s' in timings):
        timings['serial_context_setup_wall_s'] = timings['serial_grid_and_context_setup_wall_s']
    normalized['timings'] = timings
    return normalized


def validate_pair_structure(reference_manifest: dict, candidate_manifest: dict,
                             reference_result: dict, candidate_result: dict,
                             reference_side: dict, candidate_side: dict) -> None:
    for key in ('degree', 'n_theta', 'n_phi', 'full_angular_points',
                'sampled_angular_points', 'radial_side_count', 'frequency_count',
                'angular_rule', 'output_layout'):
        if reference_manifest.get(key) != candidate_manifest.get(key):
            raise ValueError('global manifests differ in structural field ' + key)
    if (reference_manifest.get('angular_subset_is_complete') is not True
            or candidate_manifest.get('angular_subset_is_complete') is not True
            or reference_side['sampled_angular_points'] != reference_side['full_angular_points']
            or candidate_side['sampled_angular_points'] != candidate_side['full_angular_points']):
        raise ValueError('cache comparison requires complete angular grids')
    if reference_manifest.get('frequency_records') != candidate_manifest.get('frequency_records'):
        raise ValueError('global frequency audits differ')
    if reference_side['native_side_indices'] != candidate_side['native_side_indices']:
        raise ValueError('global native-side indices/order differ')
    if reference_side['native_identities'] != candidate_side['native_identities']:
        raise ValueError('global native-side identities/order differ')
    if reference_side['native_side_metadata'] != candidate_side['native_side_metadata']:
        raise ValueError('global native-side metadata/order differ')
    if not np.array_equal(reference_side['native_rows'], candidate_side['native_rows']):
        raise ValueError('all-native positive-radius tables differ')


def compare_kernel(reference_path: Path, candidate_path: Path,
                   radial_sides: int, angular_points: int,
                   chunk_values: int = 65536) -> list[dict]:
    shape = (radial_sides, angular_points, len(CHANNELS))
    expected_bytes = math.prod(shape) * np.dtype(np.float64).itemsize
    if reference_path.stat().st_size != expected_bytes or candidate_path.stat().st_size != expected_bytes:
        raise ValueError('kernel binary size disagrees with the matched full-grid dimensions')
    ref = np.memmap(reference_path, dtype=np.float64, mode='r').reshape(-1, len(CHANNELS))
    cand = np.memmap(candidate_path, dtype=np.float64, mode='r').reshape(-1, len(CHANNELS))
    maximum = np.zeros(len(CHANNELS))
    scale = np.zeros(len(CHANNELS))
    for begin in range(0, len(ref), chunk_values):
        end = min(len(ref), begin + chunk_values)
        a = np.asarray(ref[begin:end])
        b = np.asarray(cand[begin:end])
        if not np.isfinite(a).all() or not np.isfinite(b).all():
            raise ValueError('nonfinite global kernel channel value')
        maximum = np.maximum(maximum, np.max(np.abs(a-b), axis=0))
        scale = np.maximum(scale, np.maximum(np.max(np.abs(a), axis=0),
                                              np.max(np.abs(b), axis=0)))
    ratio = np.divide(maximum, scale, out=np.zeros_like(maximum), where=scale != 0)
    return [{'channel': CHANNELS[i], 'max_abs_difference': float(maximum[i]),
             'scale': float(scale[i]), 'scaled_max_difference': float(ratio[i]),
             'passed': bool(ratio[i] <= pair.LIMIT)} for i in range(len(CHANNELS))]


def compare_csv(reference_path: Path, candidate_path: Path) -> tuple[list[dict], bool]:
    ha, a = pair.csv_array(reference_path)
    hb, b = pair.csv_array(candidate_path)
    if ha != hb:
        raise ValueError('CSV headers differ: ' + reference_path.name)
    return pair.compare_columns(a, b), pair.sha(reference_path) == pair.sha(candidate_path)


def assess(reference_root: Path, reference_job: str,
           candidate_root: Path, candidate_job: str,
           reference_workers: int, candidate_workers: int,
           reference_storage: str, candidate_storage: str,
           exact: bool = False) -> dict:
    roots, jobs = [reference_root, candidate_root], [reference_job, candidate_job]
    receipts = [pair.verified_job(root, job) for root, job in zip(roots, jobs)]
    dirs = [root / job / 'data' for root, job in zip(roots, jobs)]
    results = [pair.load(directory / 'result.json') for directory in dirs]
    degree = validate_physical_pair(*results)
    profiles = [pair.profile(root / job / receipt[0]['log'])
                for root, job, receipt in zip(roots, jobs, receipts)]
    if any(profile is None for profile in profiles):
        raise ValueError('global cache comparison requires a MATHPROFILE record for each job')
    manifests = [pair.load(directory / 'kernel.bin.json') for directory in dirs]
    for result, profile in zip(results, profiles):
        spool.validate_global_result(result, profile)
    sides = [spool.validate_global_manifest(manifest, result, directory,
                                             receipt, probe64=False)
             for manifest, result, directory, receipt in zip(manifests, results, dirs, receipts)]
    normalized_manifests = [_final_timing_adapter(manifest) for manifest in manifests]
    timings = [spool.validate_worker_timings(manifest) for manifest in normalized_manifests]
    invocation = [validate_invocation(result, manifest, timing, receipt, workers, storage,
                                      side['native_side_indices'])
                 for result, manifest, timing, receipt, workers, storage, side in
                 zip(results, manifests, timings, receipts,
                     (reference_workers, candidate_workers),
                     (reference_storage, candidate_storage), sides)]
    validate_pair_structure(manifests[0], manifests[1], results[0], results[1], sides[0], sides[1])
    comparisons: dict[str, list[dict]] = {}
    exact_files: dict[str, bool] = {}
    for name in ('trace.csv', 'spectrum.csv'):
        comparisons[name], exact_files[name] = compare_csv(dirs[0] / name, dirs[1] / name)
    comparisons['kernel.bin'] = compare_kernel(
        dirs[0] / 'kernel.bin', dirs[1] / 'kernel.bin', sides[0]['radial_side_count'],
        sides[0]['sampled_angular_points'])
    exact_files['kernel.bin'] = pair.sha(dirs[0] / 'kernel.bin') == pair.sha(dirs[1] / 'kernel.bin')
    comparisons['energy'] = pair.compare_columns(
        np.array([[results[0]['energy_m2_s'], results[0]['energy_adjoint_qdotp']]], dtype=float),
        np.array([[results[1]['energy_m2_s'], results[1]['energy_adjoint_qdotp']]], dtype=float))
    scalar_exact = {key: results[0][key] == results[1][key]
                    for key in ('energy_m2_s', 'energy_adjoint_qdotp')}
    if exact and (not all(exact_files.values()) or not all(scalar_exact.values())):
        raise ValueError('--exact comparison requires byte-exact kernel/trace/spectrum and exact energy/qdot')
    after = [pair.verified_job(root, job) for root, job in zip(roots, jobs)]
    for initial, current, root, job in zip(receipts, after, roots, jobs):
        spool._unchanged(root, job, initial, current)
    passed = all(row['passed'] for rows in comparisons.values() for row in rows)
    return {'status': 'PASS' if passed else 'FAIL', 'degree': degree,
            'tolerance': pair.LIMIT, 'exact_requested': exact,
            'comparison': comparisons, 'byte_exact_files': exact_files,
            'energy_qdot_exact': scalar_exact, 'invocations': invocation,
            'timings': timings,
            'resources': [{key: manifest.get(key) for key in RESOURCE_FIELDS}
                          for manifest in manifests],
            'profiles': profiles,
            'state_sha256_at_read': [receipt[2] for receipt in after],
            'verified_artifact_sha256': [receipt[1] for receipt in after]}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('reference-root', 'reference-job', 'candidate-root', 'candidate-job', 'output'):
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--reference-workers', type=int, choices=(1, 2, 4, 8), required=True)
    parser.add_argument('--candidate-workers', type=int, choices=(1, 2, 4, 8), required=True)
    parser.add_argument('--reference-storage-mode', choices=('generated', 'stored'), required=True)
    parser.add_argument('--candidate-storage-mode', choices=('generated', 'stored'), required=True)
    parser.add_argument('--exact', action='store_true',
                        help='also require byte-exact binary/CSV outputs and exact energy/qdot scalars')
    args = parser.parse_args()
    try:
        report = assess(Path(args.reference_root), args.reference_job,
                        Path(args.candidate_root), args.candidate_job,
                        args.reference_workers, args.candidate_workers,
                        args.reference_storage_mode, args.candidate_storage_mode, args.exact)
    except Exception as exc:
        report = {'status': 'FAIL', 'error_type': type(exc).__name__, 'error': str(exc)}
    with Path(args.output).open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write('\n')
    return 0 if report['status'] == 'PASS' else 1


if __name__ == '__main__':
    raise SystemExit(main())
