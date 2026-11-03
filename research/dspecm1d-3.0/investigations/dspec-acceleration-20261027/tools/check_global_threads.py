#!/usr/bin/env python3
"""Compare full-grid global jobs and validate explicit GSH inner-thread settings."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import check_global_cache as cache
import check_pair

CAP = 48
PLANNER = 'FFTWpp::Estimate'
RESULT_THREAD_KEYS = (
    'global_consumer_gsh_threads_requested',
    'global_consumer_gsh_execution_policy',
    'global_consumer_openmp_dynamic_teams_disabled',
    'global_consumer_thread_product_cap',
    'global_consumer_requested_thread_product',
    'global_consumer_worker_count_actual',
    'global_consumer_effective_thread_product_bound',
)
MANIFEST_THREAD_KEYS = (
    'gsh_threads_requested', 'gsh_execution_policy', 'openmp_available',
    'openmp_dynamic_teams_disabled', 'configured_worker_count',
    'requested_thread_product', 'effective_worker_count',
    'effective_thread_product_bound', 'compute_thread_product_cap',
    'generated_wigner_tls_bytes_per_inner_thread',
    'fft_input_output_tls_bytes_per_inner_thread_bound',
    'thread_local_fft_batch_slots_bound',
    'known_tls_bytes_per_inner_thread_bound',
    'inner_thread_tls_bytes_per_worker_bound',
    'inner_thread_tls_bytes_total_bound',
)


def _int_equals(value: Any, expected: int) -> bool:
    return isinstance(value, int) and not isinstance(value, bool) and value == expected


def validate_thread_invocation(result: dict, manifest: dict, receipt: tuple,
                               expected_workers: int, expected_gsh_threads: int,
                               legacy_serial: bool = False) -> dict[str, Any]:
    if (not _int_equals(expected_workers, expected_workers)
            or expected_workers not in (1, 2, 4, 8)
            or not _int_equals(expected_gsh_threads, expected_gsh_threads)
            or expected_gsh_threads < 1 or expected_gsh_threads > CAP
            or expected_workers * expected_gsh_threads > CAP):
        raise ValueError('requested worker/GSH thread product must be within 1..48')
    command = cache._command_tokens(receipt)
    has_thread_flag = '--global-gsh-threads' in command
    if legacy_serial:
        if expected_gsh_threads != 1:
            raise ValueError('legacy serial-reference exception is valid only for one GSH thread')
        if has_thread_flag:
            raise ValueError('legacy serial reference must omit --global-gsh-threads')
        if any(key in result for key in RESULT_THREAD_KEYS) or any(
                key in manifest for key in MANIFEST_THREAD_KEYS):
            raise ValueError('legacy serial reference must have no GSH thread metadata')
        if manifest.get('fftw_planner_policy') != PLANNER:
            raise ValueError('legacy serial reference must declare the Estimate planner policy')
        return {'expected_workers': expected_workers, 'expected_gsh_threads': 1,
                'requested_thread_product': expected_workers,
                'effective_thread_product_bound': manifest['workers'],
                'legacy_serial_reference': True,
                'fftw_planner_policy': PLANNER}

    if not has_thread_flag:
        raise ValueError('new global invocation must explicitly set --global-gsh-threads')
    if cache._option(command, '--global-gsh-threads') != str(expected_gsh_threads):
        raise ValueError('receipted GSH thread count differs from explicit expectation')
    policy = ('sequential' if expected_gsh_threads == 1
              else 'GSHTrans::Execution::Parallel(explicit)')
    requested_product = expected_workers * expected_gsh_threads
    actual_workers = manifest.get('workers')
    if (not _int_equals(result.get('global_consumer_gsh_threads_requested'), expected_gsh_threads)
            or not _int_equals(manifest.get('gsh_threads_requested'), expected_gsh_threads)):
        raise ValueError('CLI/result/manifest GSH thread counts disagree')
    if (result.get('global_consumer_gsh_execution_policy') != policy
            or manifest.get('gsh_execution_policy') != policy):
        raise ValueError('GSH execution policy disagrees with explicit thread count')
    if (not _int_equals(result.get('global_consumer_workers_requested'), expected_workers)
            or not _int_equals(manifest.get('configured_worker_count'), expected_workers)
            or not _int_equals(result.get('global_consumer_worker_count_actual'), actual_workers)
            or not _int_equals(manifest.get('effective_worker_count'), actual_workers)):
        raise ValueError('configured/actual outer worker counts disagree')
    if (not _int_equals(result.get('global_consumer_requested_thread_product'), requested_product)
            or not _int_equals(manifest.get('requested_thread_product'), requested_product)):
        raise ValueError('requested worker-by-GSH thread product disagrees')
    effective_product = actual_workers * expected_gsh_threads
    if (not _int_equals(result.get('global_consumer_effective_thread_product_bound'), effective_product)
            or not _int_equals(manifest.get('effective_thread_product_bound'), effective_product)):
        raise ValueError('effective worker-by-GSH thread product disagrees')
    if (not _int_equals(result.get('global_consumer_thread_product_cap'), CAP)
            or not _int_equals(manifest.get('compute_thread_product_cap'), CAP)
            or requested_product > CAP or effective_product > CAP):
        raise ValueError('thread product exceeds or misreports the producer cap')
    dynamic_disabled = expected_gsh_threads > 1
    if (result.get('global_consumer_openmp_dynamic_teams_disabled') is not dynamic_disabled
            or manifest.get('openmp_dynamic_teams_disabled') is not dynamic_disabled):
        raise ValueError('OpenMP dynamic-team state disagrees with explicit thread count')
    if expected_gsh_threads > 1 and manifest.get('openmp_available') is not True:
        raise ValueError('parallel GSH mode requires OpenMP availability')
    if manifest.get('fftw_planner_policy') != PLANNER:
        raise ValueError('new global job does not declare the Estimate planner policy')
    return {'expected_workers': expected_workers,
            'actual_workers': actual_workers,
            'expected_gsh_threads': expected_gsh_threads,
            'requested_thread_product': requested_product,
            'effective_thread_product_bound': effective_product,
            'legacy_serial_reference': False,
            'gsh_execution_policy': policy,
            'fftw_planner_policy': PLANNER}


def observed_peak(root: Path, job_id: str) -> dict[str, Any]:
    job, _, _ = check_pair.verified_job(root, job_id)
    fields = ('peak_observed_threads', 'peak_observed_rss_gib_userspace',
              'elapsed_seconds')
    return {key: job.get(key) for key in fields}


def assess(reference_root: Path, reference_job: str,
           candidate_root: Path, candidate_job: str,
           reference_workers: int, candidate_workers: int,
           reference_storage: str, candidate_storage: str,
           reference_gsh_threads: int, candidate_gsh_threads: int,
           reference_legacy_serial: bool = False,
           exact: bool = False) -> dict[str, Any]:
    report = cache.assess(reference_root, reference_job, candidate_root, candidate_job,
                          reference_workers, candidate_workers,
                          reference_storage, candidate_storage, exact)
    if report['status'] != 'PASS':
        return report
    receipts = [check_pair.verified_job(root, job)
                for root, job in ((reference_root, reference_job),
                                  (candidate_root, candidate_job))]
    data_dirs = [root / job / 'data' for root, job in
                 ((reference_root, reference_job), (candidate_root, candidate_job))]
    results = [check_pair.load(data / 'result.json') for data in data_dirs]
    manifests = [check_pair.load(data / 'kernel.bin.json') for data in data_dirs]
    thread_reports = [
        validate_thread_invocation(results[0], manifests[0], receipts[0],
                                   reference_workers, reference_gsh_threads,
                                   reference_legacy_serial),
        validate_thread_invocation(results[1], manifests[1], receipts[1],
                                   candidate_workers, candidate_gsh_threads, False),
    ]
    report['thread_configuration'] = {'reference': thread_reports[0],
                                      'candidate': thread_reports[1],
                                      'thread_product_cap': CAP}
    report['observed_process_peaks'] = {
        'reference': observed_peak(reference_root, reference_job),
        'candidate': observed_peak(candidate_root, candidate_job),
        'interpretation': ('Supervisor sampled peaks are reported without a minimum gate; '
                           'short jobs may finish between samples.'),
    }
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('reference-root', 'reference-job', 'candidate-root', 'candidate-job', 'output'):
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--reference-workers', type=int, choices=(1, 2, 4, 8), required=True)
    parser.add_argument('--candidate-workers', type=int, choices=(1, 2, 4, 8), required=True)
    parser.add_argument('--reference-storage-mode', choices=('generated', 'stored'), required=True)
    parser.add_argument('--candidate-storage-mode', choices=('generated', 'stored'), required=True)
    parser.add_argument('--reference-gsh-threads', type=int, required=True)
    parser.add_argument('--candidate-gsh-threads', type=int, required=True)
    parser.add_argument('--reference-legacy-serial', action='store_true',
                        help='allow only a one-thread old reference lacking thread fields/flag and declaring Estimate')
    parser.add_argument('--exact', action='store_true',
                        help='require exact kernel/trace/spectrum bytes and exact energy/qdot scalars')
    args = parser.parse_args()
    try:
        report = assess(Path(args.reference_root), args.reference_job,
                        Path(args.candidate_root), args.candidate_job,
                        args.reference_workers, args.candidate_workers,
                        args.reference_storage_mode, args.candidate_storage_mode,
                        args.reference_gsh_threads, args.candidate_gsh_threads,
                        args.reference_legacy_serial, args.exact)
    except Exception as exc:
        report = {'status': 'FAIL', 'error_type': type(exc).__name__, 'error': str(exc)}
    with Path(args.output).open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write('\n')
    return 0 if report['status'] == 'PASS' else 1


if __name__ == '__main__':
    raise SystemExit(main())
