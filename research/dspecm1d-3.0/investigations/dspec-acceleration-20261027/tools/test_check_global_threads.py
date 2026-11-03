"""Data-only regression and receipt integration checks for GSH thread accounting."""
import json
import struct
from pathlib import Path
import tempfile
import unittest

import check_global_threads as threads
import check_pair
import test_check_global_cache as cache_fixture


def _set_threads(root: Path, job_id: str, gsh_threads: int | None,
                 command_threads: int | None = None, peaks: bool = True):
    data = root / job_id / 'data'
    result = check_pair.load(data / 'result.json')
    manifest = check_pair.load(data / 'kernel.bin.json')
    if gsh_threads is None:
        manifest['fftw_planner_policy'] = threads.PLANNER
        if command_threads is None:
            manifest.pop('gsh_threads_requested', None)
    else:
        workers = 1
        product = workers * gsh_threads
        policy = ('sequential' if gsh_threads == 1
                  else 'GSHTrans::Execution::Parallel(explicit)')
        manifest.update(fftw_planner_policy=threads.PLANNER,
                        gsh_threads_requested=gsh_threads,
                        gsh_execution_policy=policy,
                        openmp_available=True,
                        openmp_dynamic_teams_disabled=gsh_threads > 1,
                        configured_worker_count=workers,
                        requested_thread_product=product,
                        effective_worker_count=workers,
                        effective_thread_product_bound=product,
                        compute_thread_product_cap=48)
        result.update(global_consumer_gsh_threads_requested=gsh_threads,
                      global_consumer_gsh_execution_policy=policy,
                      global_consumer_openmp_dynamic_teams_disabled=gsh_threads > 1,
                      global_consumer_thread_product_cap=48,
                      global_consumer_requested_thread_product=product,
                      global_consumer_worker_count_actual=workers,
                      global_consumer_effective_thread_product_bound=product)
    (data / 'result.json').write_text(json.dumps(result) + '\n')
    (data / 'kernel.bin.json').write_text(json.dumps(manifest) + '\n')
    state_path = root / 'state.json'
    state = json.loads(state_path.read_text())
    job = next(item for item in state['jobs'] if item['id'] == job_id)
    if command_threads is not None:
        job['command'].extend(('--global-gsh-threads', str(command_threads)))
    if peaks:
        job.update(peak_observed_threads=4, peak_observed_rss_gib_userspace=.2,
                   elapsed_seconds=3.)
    for rel in job['outputs']:
        job['outputs_sha256'][rel] = check_pair.sha(root / job_id / rel)
    state_path.write_text(json.dumps(state) + '\n')


class GlobalThreadChecks(unittest.TestCase):
    @staticmethod
    def _completed_jobs(root, specs):
        jobs = [cache_fixture.GlobalCacheChecks._write_completed_job(root, job_id, storage)
                for job_id, storage in specs]
        (root / 'state.json').write_text(json.dumps({'jobs': jobs}) + '\n')

    def test_new_metadata_cli_and_thread_product_match_end_to_end(self):
        with tempfile.TemporaryDirectory(prefix='global-thread-pair-') as tmp:
            root = Path(tmp)
            self._completed_jobs(root, [('ref', 'generated'), ('cand', 'stored')])
            _set_threads(root, 'ref', 1, 1)
            _set_threads(root, 'cand', 32, 32)
            report = threads.assess(root, 'ref', root, 'cand', 1, 1,
                                    'generated', 'stored', 1, 32, exact=True)
            self.assertEqual(report['status'], 'PASS')
            self.assertEqual(report['thread_configuration']['candidate']['requested_thread_product'], 32)
            self.assertEqual(report['observed_process_peaks']['candidate']['peak_observed_threads'], 4)
            self.assertTrue(report['exact_requested'])

    def test_explicit_legacy_serial_reference_only_and_estimate_bound(self):
        with tempfile.TemporaryDirectory(prefix='global-thread-legacy-') as tmp:
            root = Path(tmp)
            self._completed_jobs(root, [('ref', 'generated'), ('cand', 'generated')])
            _set_threads(root, 'ref', None)
            _set_threads(root, 'cand', 32, 32)
            report = threads.assess(root, 'ref', root, 'cand', 1, 1,
                                    'generated', 'generated', 1, 32,
                                    reference_legacy_serial=True)
            self.assertEqual(report['status'], 'PASS')
            self.assertTrue(report['thread_configuration']['reference']['legacy_serial_reference'])
            self.assertEqual(report['thread_configuration']['reference']['fftw_planner_policy'],
                             'FFTWpp::Estimate')

    def test_legacy_exception_rejects_nonserial_or_unbound_reference(self):
        receipt = ({'command': ['driver', 'kernel']}, {}, '')
        result = {'global_consumer_workers_requested': 1}
        manifest = {'workers': 1, 'fftw_planner_policy': threads.PLANNER}
        with self.assertRaisesRegex(ValueError, 'only for one GSH thread'):
            threads.validate_thread_invocation(result, manifest, receipt, 1, 2, True)
        manifest['fftw_planner_policy'] = 'unknown'
        with self.assertRaisesRegex(ValueError, 'Estimate planner'):
            threads.validate_thread_invocation(result, manifest, receipt, 1, 1, True)

    def test_thread_product_and_metadata_mismatches_fail_closed(self):
        with self.assertRaisesRegex(ValueError, 'within 1..48'):
            threads.validate_thread_invocation({}, {}, ({'command': []}, {}, ''), 2, 32)
        receipt = ({'command': ['driver', 'kernel']}, {}, '')
        with self.assertRaisesRegex(ValueError, 'explicitly set'):
            threads.validate_thread_invocation({}, {}, receipt, 1, 1)

    def test_wrong_gsh_cli_and_manifest_product_are_rejected(self):
        workers, gsh = 1, 8
        policy = 'GSHTrans::Execution::Parallel(explicit)'
        result = {'global_consumer_gsh_threads_requested': gsh,
                  'global_consumer_gsh_execution_policy': policy,
                  'global_consumer_workers_requested': workers,
                  'global_consumer_worker_count_actual': workers,
                  'global_consumer_requested_thread_product': 7,
                  'global_consumer_effective_thread_product_bound': 8,
                  'global_consumer_thread_product_cap': 48,
                  'global_consumer_openmp_dynamic_teams_disabled': True}
        manifest = {'workers': 1, 'gsh_threads_requested': gsh,
                    'gsh_execution_policy': policy, 'openmp_available': True,
                    'openmp_dynamic_teams_disabled': True,
                    'configured_worker_count': 1, 'requested_thread_product': 7,
                    'effective_worker_count': 1, 'effective_thread_product_bound': 8,
                    'compute_thread_product_cap': 48,
                    'fftw_planner_policy': threads.PLANNER}
        receipt = ({'command': ['driver', 'kernel', '--global-gsh-threads', '4']}, {}, '')
        with self.assertRaisesRegex(ValueError, 'differs from explicit expectation'):
            threads.validate_thread_invocation(result, manifest, receipt, 1, 8)
        receipt[0]['command'][-1] = '8'
        with self.assertRaisesRegex(ValueError, 'requested worker-by-GSH thread product'):
            threads.validate_thread_invocation(result, manifest, receipt, 1, 8)

    def test_peak_sampling_has_no_minimum_and_reports_missing_values(self):
        with tempfile.TemporaryDirectory(prefix='global-thread-peak-') as tmp:
            root = Path(tmp)
            self._completed_jobs(root, [('short', 'generated')])
            _set_threads(root, 'short', 1, 1, peaks=False)
            peak = threads.observed_peak(root, 'short')
            self.assertIsNone(peak['peak_observed_threads'])
            self.assertIsNone(peak['peak_observed_rss_gib_userspace'])

    def test_exact_gate_rejects_kernel_channel_regression(self):
        with tempfile.TemporaryDirectory(prefix='global-thread-exact-regression-') as tmp:
            root = Path(tmp)
            self._completed_jobs(root, [('ref', 'generated'), ('cand', 'generated')])
            _set_threads(root, 'ref', 1, 1)
            _set_threads(root, 'cand', 8, 8)
            kernel = root / 'cand' / 'data' / 'kernel.bin'
            with kernel.open('r+b') as stream:
                stream.write(struct.pack('<d', 1e-10))
            state_path = root / 'state.json'
            state = json.loads(state_path.read_text())
            candidate_job = next(job for job in state['jobs'] if job['id'] == 'cand')
            candidate_job['outputs_sha256']['data/kernel.bin'] = check_pair.sha(kernel)
            state_path.write_text(json.dumps(state) + '\n')
            with self.assertRaisesRegex(ValueError, '--exact comparison requires'):
                threads.assess(root, 'ref', root, 'cand', 1, 1,
                               'generated', 'generated', 1, 8, exact=True)


if __name__ == '__main__':
    unittest.main()
