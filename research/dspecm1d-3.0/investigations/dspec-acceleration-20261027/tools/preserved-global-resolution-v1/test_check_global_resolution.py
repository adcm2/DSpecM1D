import unittest

import numpy as np

import check_global_resolution as resolution


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
        resolution._validate_serial_generated(result, manifest, sides, 128,
                                              {'command': command})
        sides['native_side_indices'] = [59, 101]
        with self.assertRaises(ValueError):
            resolution._validate_serial_generated(result, manifest, sides, 128,
                                                  {'command': command})


if __name__ == '__main__':
    unittest.main()
