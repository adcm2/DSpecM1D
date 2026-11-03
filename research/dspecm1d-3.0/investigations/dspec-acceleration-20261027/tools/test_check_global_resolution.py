import csv
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

import check_global_resolution as resolution
import check_pair


class ResolutionMathTests(unittest.TestCase):
    def test_bilinear_theta_and_periodic_phi(self):
        # At theta halfway between rows and phi halfway across the periodic seam,
        # interpolation blends the 3/0 and 13/10 seam pairs.
        values = np.array([[[0.0], [1.0], [2.0], [3.0]],
                           [[10.0], [11.0], [12.0], [13.0]]])
        result = resolution.bilinear_interpolate(
            values, np.array([0.0, 1.0]), np.array([0.5]),
            np.array([7*np.pi/4]))
        self.assertEqual(result.shape, (1, 1))
        self.assertAlmostEqual(float(result[0, 0]), 6.5)

    def test_fixed_great_circle_excludes_endpoint_neighborhoods(self):
        path = resolution.great_circle((10.0, 20.0), (45.0, 75.0))
        self.assertEqual(path['input_points'], 720)
        self.assertEqual(len(path['theta']), path['included_points'])
        self.assertGreaterEqual(path['first_kept_distance_deg'], 5.0)
        self.assertGreaterEqual(path['last_kept_distance_to_receiver_deg'], 5.0)
        with self.assertRaises(ValueError):
            resolution.great_circle((10.0, 20.0), (10.0, 20.0))

    def test_metrics_mark_only_paired_exact_zeros_structural(self):
        zeros = np.zeros((4, 6))
        metrics = resolution._channel_metrics(zeros, zeros)
        self.assertTrue(all(row['structural_zero'] for row in metrics))
        self.assertTrue(all(row['max_scaled_difference'] == 0.0 for row in metrics))
        changed = zeros.copy()
        changed[0, 0] = 1.0
        metrics = resolution._channel_metrics(zeros, changed)
        self.assertFalse(metrics[0]['structural_zero'])
        self.assertEqual(metrics[0]['max_scaled_difference'], 1.0)

    def test_serial_generated_gate_binds_selected_sides_and_command(self):
        result = {
            'global_consumer_workers_requested': 1,
            'global_consumer_wigner_storage_requested': 'generated',
        }
        manifest = {
            'degree': 128, 'workers': 1, 'wigner_storage_policy': 'generated',
            'fftw_planner_policy': 'FFTWpp::Estimate',
            'sampled_angular_points': 257*520, 'full_angular_points': 257*520,
            'angular_subset_is_complete': True,
        }
        sides = {
            'native_side_indices': [59, 100],
            'native_side_metadata': [{'fluid': True}, {'fluid': False}],
        }
        command = ['driver', 'model', 'data', 'kernel', '--lmax', '128',
                   '--global-workers', '1', '--global-radial-sides', '59,100',
                   '--global-wigner-storage', 'generated', '--forward-only', 'off',
                   '--contacts', 'off', '--pole-only', 'on']
        receipt = ({'command': command}, {}, '')
        resolution._validate_serial_generated(result, manifest, sides, 128, receipt)
        sides['native_side_indices'] = [59, 101]
        with self.assertRaises(ValueError):
            resolution._validate_serial_generated(result, manifest, sides, 128, receipt)


class ResolutionReportFixture(unittest.TestCase):
    @staticmethod
    def _sha(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()

    def _job(self, root, job_id, files, command):
        job_dir = root / job_id
        job_dir.mkdir(parents=True)
        hashes = {}
        for relative, content in files.items():
            path = job_dir / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(content)
            hashes[relative] = self._sha(path)
        log = 'job.log'
        (job_dir / log).write_text(
            'MATHPROFILE ' + json.dumps({'solver': 'dspec', 'objective': -1,
                'stages': {'physical_case_total': {
                    'inclusive_s': 0.0, 'exclusive_s': 0.0, 'calls': 0}}}) + '\n')
        state_path = root / 'state.json'
        if state_path.exists():
            state = json.loads(state_path.read_text())
        else:
            state = {'jobs': []}
        state['jobs'].append({
            'id': job_id, 'status': 'completed', 'returncode': 0,
            'inputs_unchanged': True, 'cleanup': {'verified_gone': True},
            'failure_reasons': [], 'outputs': sorted(hashes),
            'outputs_sha256': hashes, 'log': log,
            'log_sha256': self._sha(job_dir / log), 'command': command,
        })
        state_path.write_text(json.dumps(state))

    @staticmethod
    def _physical_result(degree, mode='kernel'):
        result = {key: 1 for key in check_pair.PHYSICAL}
        result.update({
            'degree': degree, 'mode': mode, 'sem_degree': 6,
            'native_element_node_sides': 180, 'frequency_bins': 52,
            'source_contact_enabled': False, 'forward_only_baseline': False,
            'source': [10.0, 20.0, 50000.0], 'receiver': [45.0, 75.0, 3000.0],
            'energy_m2_s': 1.0, 'energy_adjoint_qdotp': 2.0,
            'energy_adjoint_relative_error_vs_2E': 0.0,
            'global_consumer_workers_requested': 1,
            'global_consumer_wigner_storage_requested': 'generated',
        })
        return result

    def test_actual_schema_compare_and_json_writer(self):
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            global_roots = [base / 'global128', base / 'global384']
            layout_root = base / 'layout'
            layout_root.mkdir()
            for degree, root, job_id in ((128, global_roots[0], 'g128'),
                                         (384, global_roots[1], 'g384')):
                root.mkdir()
                ntheta, nphi = 2*degree+1, 4
                result = self._physical_result(degree)
                manifest = {
                    'degree': degree, 'workers': 1, 'wigner_storage_policy': 'generated',
                    'fftw_planner_policy': 'FFTWpp::Estimate',
                    'n_theta': ntheta, 'n_phi': nphi,
                    'sampled_angular_points': ntheta*nphi,
                    'full_angular_points': ntheta*nphi,
                    'angular_subset_is_complete': True,
                    'radial_side_count': 2,
                    'output_bytes': 2*ntheta*nphi*6*8,
                }
                command = ['driver', 'model', 'data', 'kernel', '--lmax', str(degree),
                    '--global-workers', '1', '--global-radial-sides', '59,100',
                    '--global-wigner-storage', 'generated', '--forward-only', 'off',
                    '--contacts', 'off', '--pole-only', 'on']
                kernel = np.ones((2, ntheta, nphi, 6), dtype=np.float64).tobytes()
                self._job(root, job_id, {
                    'data/result.json': json.dumps(result).encode(),
                    'data/kernel.bin.json': json.dumps(manifest).encode(),
                    'data/kernel.bin': kernel,
                    'data/radial_nodes.csv': b'element,node,radius_m,element_width_m,fluid\n',
                    'data/pole_response_spool.bin': b'fixture',
                }, command)

            theta_jobs = []
            for degree, job_id in ((128, 't128'), (384, 't384')):
                ntheta = 2*degree+1
                theta = np.linspace(0.001, np.pi-0.001, ntheta)
                weights = np.full(ntheta, 2.0/ntheta)
                rows = ['theta_index,theta_rad,theta_weight'] + [
                    f'{i},{t:.17g},{w:.17g}' for i, (t, w) in enumerate(zip(theta, weights))]
                result = self._physical_result(degree, 'kernel-probe')
                probe = {'n_theta': ntheta, 'sampled_angular_points_per_trace': 64,
                         'sampled_angle_subset_is_complete': False}
                command = ['driver', 'model', 'data', 'kernel-probe', '--lmax', str(degree),
                           '--probe-angle-count', '64']
                self._job(layout_root, job_id, {
                    'data/result.json': json.dumps(result).encode(),
                    'data/kernel_probe.json': json.dumps(probe).encode(),
                    'data/angular_grid.csv': ('\n'.join(rows)+'\n').encode(),
                }, command)
                theta_jobs.append(job_id)

            def fake_manifest(manifest, result, data, receipt, probe64):
                return {
                    'native_side_indices': [59, 100],
                    'native_side_metadata': [
                        {'fluid': True, 'radius_m': 1.0},
                        {'fluid': False, 'radius_m': 0.5}],
                    'native_identities': [(10, 1), (20, 2)],
                    'radial_side_count': 2, 'n_theta': manifest['n_theta'],
                    'n_phi': manifest['n_phi'],
                }

            with (mock.patch.object(resolution.global_check, 'validate_global_result'),
                  mock.patch.object(resolution.global_check, 'validate_global_manifest', fake_manifest),
                  mock.patch.object(resolution.global_check, '_unchanged',
                                    side_effect=lambda roots, jobs, before: before)):
                report = resolution.compare(global_roots[0], 'g128', global_roots[1],
                    'g384', layout_root, theta_jobs[0], theta_jobs[1])
            output = base / 'report.json'
            resolution.write_report(report, output)
            loaded = json.loads(output.read_text())
            self.assertEqual(loaded['status'], 'DESCRIPTIVE_COMPARISON')
            self.assertIsNone(loaded['acceptance_gate'])
            self.assertIsInstance(loaded['path']['theta'], list)
            self.assertEqual(len(loaded['path']['theta']), loaded['path']['included_points'])

            # The helper must bind the result hash before it reads probe metadata.
            state_path = layout_root / 'state.json'
            state = json.loads(state_path.read_text())
            theta_job = next(job for job in state['jobs'] if job['id'] == theta_jobs[0])
            theta_job['outputs_sha256']['data/result.json'] = '0'*64
            state_path.write_text(json.dumps(state))
            with self.assertRaises(ValueError):
                resolution.read_theta_grid(layout_root, theta_jobs[0], 128)

            # A result file alone cannot make a receipted full-kernel or wrong-L
            # invocation pass as the requested kernel-probe reference.
            state = json.loads(state_path.read_text())
            theta_job = next(job for job in state['jobs'] if job['id'] == theta_jobs[1])
            theta_job['command'] = ['driver', 'model', 'data', 'kernel', '--lmax', '128',
                                    '--probe-angle-count', '64']
            state_path.write_text(json.dumps(state))
            with self.assertRaises(ValueError):
                resolution.read_theta_grid(layout_root, theta_jobs[1], 384)


if __name__ == '__main__':
    unittest.main()
