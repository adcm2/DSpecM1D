"""Synthetic receipt and scientific-domain regression checks; no solver calls."""
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import check_pair as check


class CandidateChecks(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='dspec-accel-check-')
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name)
        self.a = self.fixture('reference')
        self.b = self.fixture('candidate')

    def fixture(self, label):
        root = self.base / label; d = root / 'case/data'; d.mkdir(parents=True)
        r = {k: 1 for k in check.PHYSICAL}
        r.update(mode='kernel-probe', degree=3, sem_elements=1, gll_points_per_element=6,
                 native_element_node_sides=6, source=[10, 20, 50000], receiver=[45, 75, 3000],
                 source_contact_enabled=False, energy_m2_s=2., energy_adjoint_qdotp=4.,
                 energy_adjoint_relative_error_vs_2E=0.)
        (d / 'result.json').write_text(json.dumps(r))
        dom = dict(selected_radial_trace_count=1, available_radial_trace_sides=1, n_theta=7,
                   n_phi=14, actual_sampled_point_count=2, sampled_angular_points_per_trace=2,
                   binary_layout='point-major doubles[A,C,F,L,N,kappa]', source_contact_enabled=False,
                   full_target_angular_points_per_trace=98, full_target_point_count=98,
                   full_target_kernel_binary_bytes=98*48, actual_kernel_binary_bytes=2*48,
                   sampled_angle_subset_is_complete=False, diagnostic_angle_subset_only=True,
                   actual_contact_binary_bytes=0, full_target_contact_binary_bytes=0,
                   actual_combined_output_bytes=2*48)
        (d / 'kernel_probe.json').write_text(json.dumps(dom))
        np.arange(12, dtype=np.float64).tofile(d / 'kernel_probe.bin')
        for name in ('trace.csv', 'spectrum.csv'):
            (d / name).write_text('a,b\n0,1\n1,2\n')
        (d / 'radial_nodes.csv').write_text('element,node,radius_m,element_width_m,fluid\n0,1,100,10,0\n')
        (d / 'kernel_radial_nodes.csv').write_text('radial_index,element,node,radius_m,element_width_m,fluid\n0,0,1,100,10,0\n')
        x, weights = np.polynomial.legendre.leggauss(7); theta = np.arccos(-x)
        grid = np.column_stack((np.arange(7), theta, weights))
        np.savetxt(d / 'angular_grid.csv', grid, delimiter=',', header='theta_index,theta_rad,theta_weight', comments='')
        angles = np.array([[i, 0, i, theta[0], 2*np.pi*i/14, weights[0], 2*np.pi/14] for i in range(2)])
        np.savetxt(d / 'kernel_probe_angles.csv', angles, delimiter=',',
                   header='sample_index,theta_index,phi_index,theta_rad,phi_rad,theta_weight,phi_weight', comments='')
        p = dict(solver='dspec', objective=-1, stages={'physical_case_total': dict(inclusive_s=1., exclusive_s=1., calls=1)}, counters={})
        (root / 'case/job.log').write_text('MATHPROFILE ' + json.dumps(p) + '\n')
        self.receipt(root)
        return root

    def receipt(self, root):
        d = root / 'case/data'; outputs = {'data/' + p.name: check.sha(p) for p in d.iterdir()}
        job = dict(id='case', status='completed', returncode=0, inputs_unchanged=True,
                   cleanup={'verified_gone': True}, failure_reasons=[], outputs=list(outputs),
                   outputs_sha256=outputs, log='job.log', log_sha256=check.sha(root / 'case/job.log'))
        (root / 'state.json').write_text(json.dumps({'jobs': [job]}))

    def assess(self, exact=False):
        return check.assess(self.a, 'case', self.b, 'case', exact)

    def test_exact_pair_passes(self):
        self.assertEqual(self.assess(True)['status'], 'PASS')

    def test_unreceipted_artifact_tamper_rejected(self):
        (self.b / 'case/data/kernel_probe.bin').write_bytes(b'bad')
        with self.assertRaisesRegex(ValueError, 'hash mismatch'): self.assess()

    def test_receipted_numerical_regression_fails(self):
        p = self.b / 'case/data/kernel_probe.bin'; a = np.fromfile(p, dtype=np.float64); a[1] += .01; a.tofile(p)
        self.receipt(self.b)
        self.assertEqual(self.assess()['status'], 'FAIL')

    def test_forged_reported_adjoint_error_rejected(self):
        p = self.b / 'case/data/result.json'; r = check.load(p); r['energy_adjoint_qdotp'] = 3.; p.write_text(json.dumps(r)); self.receipt(self.b)
        with self.assertRaisesRegex(ValueError, 'recomputed'): self.assess()

    def test_inconsistent_point_product_rejected(self):
        p = self.b / 'case/data/kernel_probe.json'; r = check.load(p); r['selected_radial_trace_count'] = 2; p.write_text(json.dumps(r)); self.receipt(self.b)
        with self.assertRaisesRegex(ValueError, 'product'): self.assess()

    def test_receipted_wrong_radial_side_rejected(self):
        p = self.b / 'case/data/kernel_radial_nodes.csv'; p.write_text(p.read_text().replace('100,10,0', '101,10,0')); self.receipt(self.b)
        with self.assertRaisesRegex(ValueError, 'native radius/material'): self.assess()

    def test_unverified_cleanup_rejected(self):
        p = self.b / 'state.json'; r = check.load(p); r['jobs'][0]['cleanup']['verified_gone'] = False; p.write_text(json.dumps(r))
        with self.assertRaisesRegex(ValueError, 'unqualified'): self.assess()

    def test_changed_solve_flags_rejected(self):
        p = self.b / 'case/data/result.json'
        original = check.load(p)
        for flag in ('pole_only_nonradial_adjoint', 'forward_only_baseline'):
            r = dict(original); r[flag] = 0; p.write_text(json.dumps(r)); self.receipt(self.b)
            with self.assertRaisesRegex(ValueError, flag): self.assess(True)


if __name__ == '__main__':
    unittest.main()
