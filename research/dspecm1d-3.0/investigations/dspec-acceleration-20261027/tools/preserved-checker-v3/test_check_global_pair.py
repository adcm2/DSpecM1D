"""Small data-only checks for the global spool comparison helpers."""
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

import check_global_pair as check
import check_pair as pair


class GlobalPairChecks(unittest.TestCase):
    @staticmethod
    def _frequency_records():
        records = [dict(spool_index=i, original_bin=i, z=[float(i), -0.1],
                        adjoint_mixing=[0.25, -0.5],
                        receiver_basis_reconstruction_relative=1e-13)
                   for i in range(53)]
        records[52]['adjoint_mixing'] = [0.0, 0.0]
        return records

    @staticmethod
    def _native_table(path):
        with path.open('w') as stream:
            stream.write('element,node,radius_m,element_width_m,fluid\n')
            for i in range(179):
                stream.write(f'{i},1,{1000+i},10,0\n')
        return [dict(element=i, node=1, radius_m=float(1000+i),
                     radius_model=float(1000+i)/2.0, width_m=10.0, fluid=False)
                for i in range(8)]

    def test_transpose_and_native_side_permutation_align_exactly(self):
        with tempfile.TemporaryDirectory(prefix='global-pair-layout-') as tmp:
            root = Path(tmp)
            angles, sides = 7, 8
            direct = np.arange(angles * sides * 6, dtype=np.float64).reshape(angles, sides, 6)
            identities = [(4, i) for i in range(sides)]
            permutation = [7, 0, 3, 1, 6, 2, 5, 4]
            global_layout = direct[:, permutation, :].transpose(1, 0, 2).copy()
            direct_path, global_path = root / 'direct.bin', root / 'global.bin'
            direct.tofile(direct_path)
            global_layout.tofile(global_path)
            ref = {'angular_points': angles, 'radial_side_count': sides,
                   'selected_map': {identity: i for i, identity in enumerate(identities)}}
            candidate = {'radial_side_count': sides,
                         'native_identities': [identities[i] for i in permutation]}
            result = check.compare_kernel(direct_path, global_path, ref, candidate, angles)
            self.assertTrue(all(row['passed'] for row in result))
            self.assertEqual([row['max_abs_difference'] for row in result], [0.] * 6)

    def test_wrong_native_side_set_is_rejected(self):
        ref = {'angular_points': 1, 'radial_side_count': 8,
               'selected_map': {(1, i): i for i in range(8)}}
        candidate = {'radial_side_count': 8,
                     'native_identities': [(1, i) for i in range(7)] + [(1, 99)]}
        with tempfile.TemporaryDirectory(prefix='global-pair-side-') as tmp:
            path = Path(tmp) / 'empty.bin'
            np.zeros(8 * 6, dtype=np.float64).tofile(path)
            with self.assertRaisesRegex(ValueError, 'subset'):
                check.compare_kernel(path, path, ref, candidate, 1)

    def test_frequency_audit_requires_contiguous_bins_complex_q_and_small_residual(self):
        records = self._frequency_records()
        check.validate_frequency_records(records, 53)
        records[11]['original_bin'] = 12
        with self.assertRaisesRegex(ValueError, 'contiguous'):
            check.validate_frequency_records(records, 53)
        records[11]['original_bin'] = 11
        records[11]['adjoint_mixing'] = [0., 0.]
        with self.assertRaisesRegex(ValueError, 'in-band'):
            check.validate_frequency_records(records, 53)
        records[11]['adjoint_mixing'] = [0.25, -0.5]
        records[11]['receiver_basis_reconstruction_relative'] = 1.1e-12
        with self.assertRaisesRegex(ValueError, '1e-12'):
            check.validate_frequency_records(records, 53)

    def test_combined_v3_global_manifest_accepts_180_total_179_positive_and_eight_selected(self):
        with tempfile.TemporaryDirectory(prefix='global-pair-v3-schema-') as tmp:
            data = Path(tmp)
            sides = self._native_table(data / 'radial_nodes.csv')
            (data / 'kernel.bin').write_bytes(bytes(8 * 91 * 6 * 8))
            (data / 'pole_response_spool.bin').write_bytes(b'spool')
            result = dict(mode='kernel', degree=3, sem_degree=6,
                          native_element_node_sides=180, frequency_bins=52,
                          model_length_scale_m=2.0,
                          global_consumer_radial_side_count=8,
                          global_consumer_side_selection_explicit=True,
                          response_spool_bytes=5)
            metadata = [dict(side, radius_m=float(side['radius_m'])) for side in sides]
            manifest = dict(
                mode='compact_response_global_gsh_fft_frequency_accumulator',
                degree=3, n_theta=7, n_phi=13, full_angular_points=91,
                sampled_angular_points=91, angular_subset_is_complete=True,
                angular_rule='ForBand(L,2,2.0) nodes and FastFFTSize(4L+1) longitudes',
                radial_side_count=8, native_side_indices=list(range(8)),
                native_side_metadata=metadata, frequency_count=53,
                frequency_records=self._frequency_records(),
                source_contacts_included=False,
                kernel_is_full_earth_volume_integral=False,
                output_layout='radial-side-major; theta-major angular prefix; point-major channels A,C,F,L,N,kappa',
                output_bytes=8 * 91 * 6 * 8, binary_file='kernel.bin')
            receipted = {f'data/{name}': 'unused' for name in check.GLOBAL_FILES}
            validated = check.validate_global_manifest(
                manifest, result, data, ({}, receipted, 'statehash'), False)
            self.assertEqual(validated['radial_side_count'], 8)
            self.assertEqual(len(validated['native_rows']), 179)
            self.assertEqual(validated['spool_frequency_count'], 53)
            manifest['radial_side_count'] = 2
            manifest['native_side_indices'] = manifest['native_side_indices'][:2]
            manifest['native_side_metadata'] = manifest['native_side_metadata'][:2]
            manifest['output_bytes'] = 2 * 91 * 6 * 8
            (data / 'kernel.bin').write_bytes(bytes(manifest['output_bytes']))
            result['global_consumer_radial_side_count'] = 2
            self.assertEqual(check.validate_global_manifest(
                manifest, result, data, ({}, receipted, 'statehash'), False)['radial_side_count'], 2)
            manifest['radial_side_count'] = 9
            with self.assertRaisesRegex(ValueError, 'eight-side selection cap'):
                check.validate_global_manifest(
                    manifest, result, data, ({}, receipted, 'statehash'), False)
            manifest['radial_side_count'] = 2
            manifest['native_side_metadata'][0]['element'] = 99
            with self.assertRaisesRegex(ValueError, 'disagrees with candidate'):
                check.validate_global_manifest(
                    manifest, result, data, ({}, receipted, 'statehash'), False)

    def test_combined_v3_selected_side_counts_within_cap_are_valid(self):
        selected = [(i, 1) for i in range(2)]
        check.validate_selected_identity_subset(selected, {(i, 1) for i in range(179)})

    def test_combined_v3_result_counter_matches_53_by_2l_plus_1_solves(self):
        result = dict(degree=3, response_record_spool_used=True,
                      response_reuse_used=False, response_cache_entries=0,
                      response_cache_bytes=0)
        profile = {'counters': {'sparse_factorization_calls': 53 * 7}}
        check.validate_global_result(result, profile)
        profile['counters']['sparse_factorization_calls'] -= 1
        with self.assertRaisesRegex(ValueError, 'one pass'):
            check.validate_global_result(result, profile)

    def test_parallel_gate_requires_actual_w1_and_w4_for_eight_sides(self):
        reference = dict(global_consumer_workers_requested=1,
                         global_consumer_radial_side_count=8)
        candidate = dict(global_consumer_workers_requested=4,
                         global_consumer_radial_side_count=8)
        pair_result = check.validate_parallel_worker_pair(
            reference, {'workers': 1}, {'actual_workers': 1},
            candidate, {'workers': 4}, {'actual_workers': 4})
        self.assertEqual(pair_result['reference']['actual'], 1)
        self.assertEqual(pair_result['candidate']['actual'], 4)
        with self.assertRaisesRegex(ValueError, 'actual global worker count'):
            check.validate_parallel_worker_pair(
                reference, {'workers': 1}, {'actual_workers': 1},
                candidate, {'workers': 3}, {'actual_workers': 3})

    def test_global_to_global_checks_kernel_trace_spectrum_and_energy_qdot_exact(self):
        with tempfile.TemporaryDirectory(prefix='global-pair-exact-') as tmp:
            root = Path(tmp)
            left, right = root / 'left', root / 'right'
            left.mkdir(); right.mkdir()
            for name in ('kernel.bin', 'trace.csv', 'spectrum.csv'):
                (left / name).write_bytes(b'producer bytes')
                (right / name).write_bytes(b'producer bytes')
            result = dict(energy_m2_s=2.0, energy_adjoint_qdotp=4.0)
            parity = check.validate_global_outputs_equal(left, right, result, dict(result))
            self.assertTrue(all(parity['byte_exact'].values()))
            self.assertTrue(all(parity['energy_qdot_exact'].values()))
            (right / 'spectrum.csv').write_bytes(b'different spectrum')
            with self.assertRaisesRegex(ValueError, 'spectrum.csv'):
                check.validate_global_outputs_equal(left, right, result, dict(result))
            (right / 'spectrum.csv').write_bytes(b'producer bytes')
            other = dict(result, energy_adjoint_qdotp=4.000000000000001)
            with self.assertRaisesRegex(ValueError, 'energy/qdot'):
                check.validate_global_outputs_equal(left, right, result, other)

    def test_combined_v3_physical_matrix_is_l3_l12_sem6_52_output_bins(self):
        for degree in (3, 12):
            reference = {key: 1 for key in pair.PHYSICAL}
            reference.update(mode='kernel', degree=degree, sem_degree=6,
                             native_element_node_sides=180, frequency_bins=52,
                             source_contact_enabled=False,
                             forward_only_baseline=False, energy_m2_s=2.0,
                             energy_adjoint_qdotp=4.0,
                             energy_adjoint_relative_error_vs_2E=0.0)
            candidate = dict(reference)
            check.validate_physical_pair(reference, candidate)
            reference['frequency_bins'] = 53
            candidate['frequency_bins'] = 53
            with self.assertRaisesRegex(ValueError, '52 output frequency'):
                check.validate_physical_pair(reference, candidate)

    def test_higher_degree_is_allowed_only_on_explicit_probe_route(self):
        reference = {key: 1 for key in pair.PHYSICAL}
        reference.update(mode='kernel', degree=32, sem_degree=6,
                         native_element_node_sides=180, frequency_bins=52,
                         source_contact_enabled=False, forward_only_baseline=False,
                         energy_m2_s=2.0, energy_adjoint_qdotp=4.0,
                         energy_adjoint_relative_error_vs_2E=0.0)
        candidate = dict(reference, mode='kernel-probe')
        check.validate_physical_pair(reference, candidate, 'kernel-probe')
        candidate['mode'] = 'kernel'
        with self.assertRaisesRegex(ValueError, 'L=3 and L=12'):
            check.validate_physical_pair(reference, candidate)

    def test_guard_q_must_be_zero_and_inband_q_must_be_nonzero(self):
        records = self._frequency_records()
        check.validate_frequency_records(records, 53)
        records[52]['adjoint_mixing'] = [1.0, 0.0]
        with self.assertRaisesRegex(ValueError, 'guard record'):
            check.validate_frequency_records(records, 53)
        records = self._frequency_records()
        records[17]['adjoint_mixing'] = [0.0, 0.0]
        with self.assertRaisesRegex(ValueError, 'in-band'):
            check.validate_frequency_records(records, 53)

    def test_eight_global_sides_are_a_subset_of_full_reference_not_equal_to_it(self):
        full = {(i, 1) for i in range(179)}
        selected = [(i, 1) for i in range(8)]
        check.validate_selected_identity_subset(selected, full)
        with self.assertRaisesRegex(ValueError, 'subset'):
            check.validate_selected_identity_subset(selected[:-1] + [(999, 1)], full)
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            check.validate_selected_identity_subset(selected[:-1] + [selected[0]], full)

    def test_csv_receipt_lookup_uses_data_prefix_exactly_once(self):
        with tempfile.TemporaryDirectory(prefix='global-pair-csv-receipt-') as tmp:
            root = Path(tmp)
            candidate_root = root / 'candidate'
            candidate_root.mkdir()
            a, b = root / 'trace.csv', candidate_root / 'trace.csv'
            a.write_text('x,y\n1,2\n')
            b.write_text('x,y\n1,2\n')
            receipts = (({}, {'data/trace.csv': 'a'}, ''),
                        ({}, {'data/trace.csv': 'b'}, ''))
            comparisons, exact = check.compare_csv(root, candidate_root, receipts, ('trace.csv',))
            self.assertTrue(all(row['passed'] for row in comparisons['trace.csv']))
            self.assertTrue(exact['trace.csv'])

    def test_completed_receipt_mutation_is_rejected_by_recheck(self):
        with tempfile.TemporaryDirectory(prefix='global-pair-receipt-') as tmp:
            root = Path(tmp)
            data = root / 'job' / 'data'
            data.mkdir(parents=True)
            artifact = data / 'native.csv'
            artifact.write_text('x\n1\n')
            log = root / 'job' / 'job.log'
            log.write_text('receipt log\n')
            job = dict(id='job', status='completed', returncode=0,
                       inputs_unchanged=True, cleanup={'verified_gone': True},
                       failure_reasons=[], outputs=['data/native.csv'],
                       outputs_sha256={'data/native.csv': pair.sha(artifact)},
                       log='job.log', log_sha256=pair.sha(log))
            (root / 'state.json').write_text(json.dumps({'jobs': [job]}))
            before = pair.verified_job(root, 'job')
            artifact.write_text('x\n2\n')
            with self.assertRaisesRegex(ValueError, 'hash mismatch'):
                check._unchanged([root], ['job'], [before])

    def test_qdot_energy_identity_is_recomputed(self):
        valid = dict(energy_m2_s=2., energy_adjoint_qdotp=4.,
                     energy_adjoint_relative_error_vs_2E=0.)
        check.validate_energy(valid)
        valid['energy_adjoint_qdotp'] = 3.9
        with self.assertRaisesRegex(ValueError, 'recomputed'):
            check.validate_energy(valid)

    def test_worker_wall_is_reported_separately_from_worker_elapsed_sums(self):
        stage = dict(side_elapsed_s=.8, spool_reads_s=.2,
                     coefficient_rotation_layout_s=.2, gsh_inverse_fft_s=.2,
                     point_kernel_contractions_s=.1, frequency_accumulation_s=.1,
                     output_writes_s=.1)
        manifest = dict(elapsed_s=3., workers=2, timings=dict(
            total_wall_s=3., worker_execution_wall_s=1.5,
            serial_context_setup_wall_s=.5, worker_elapsed_sum_s=1.6,
            worker_stage_elapsed_sums=dict(
                spool_reads=.4, coefficient_rotation_layout=.4,
                gsh_inverse_fft=.4, point_kernel_contractions=.2,
                frequency_accumulation=.2, output_writes=.2),
            workers=[stage, stage]))
        result = check.validate_worker_timings(manifest)
        self.assertEqual(result['worker_execution_wall_s'], 1.5)
        self.assertEqual(result['serial_context_setup_wall_s'], .5)
        self.assertEqual(result['worker_elapsed_sum_overlapping'], 1.6)
        self.assertEqual(result['actual_workers'], 2)
        manifest['timings']['worker_elapsed_sum_s'] = 2.1
        with self.assertRaisesRegex(ValueError, 'worker_elapsed_sum_s disagrees'):
            check.validate_worker_timings(manifest)
        manifest['timings']['worker_elapsed_sum_s'] = 1.6
        manifest['timings']['worker_stage_elapsed_sums']['spool_reads'] = .5
        with self.assertRaisesRegex(ValueError, 'worker_stage_elapsed_sums disagrees'):
            check.validate_worker_timings(manifest)
        manifest['timings']['worker_stage_elapsed_sums']['spool_reads'] = .4
        del manifest['timings']['worker_elapsed_sum_s']
        with self.assertRaisesRegex(ValueError, 'missing or invalid'):
            check.validate_worker_timings(manifest)
        manifest['timings']['worker_elapsed_sum_s'] = 1.6
        manifest['timings']['serial_context_setup_wall_s'] = 2.
        with self.assertRaisesRegex(ValueError, 'exceeds total'):
            check.validate_worker_timings(manifest)


if __name__ == '__main__':
    unittest.main()
