"""Synthetic data-only checks for generated/stored full-grid pair validation."""
import tempfile
import json
from pathlib import Path
import unittest

import numpy as np

import check_global_cache as check
import check_global_pair as spool
import check_pair as pair


class GlobalCacheChecks(unittest.TestCase):
    @staticmethod
    def _records():
        rows = [dict(spool_index=i, original_bin=i, z=[float(i), -0.1],
                     adjoint_mixing=[0.25, -0.5],
                     receiver_basis_reconstruction_relative=1e-13)
                for i in range(53)]
        rows[52]['adjoint_mixing'] = [0.0, 0.0]
        return rows

    @staticmethod
    def _result(degree=3, workers=1, storage='generated', sides=2):
        result = {key: 1 for key in pair.PHYSICAL}
        result.update(mode='kernel', degree=degree, sem_degree=6,
                      native_element_node_sides=180, frequency_bins=52,
                      source_contact_enabled=False, forward_only_baseline=False,
                      energy_m2_s=2., energy_adjoint_qdotp=4.,
                      energy_adjoint_relative_error_vs_2E=0.,
                      model_length_scale_m=2., global_consumer_radial_side_count=sides,
                      global_consumer_side_selection_explicit=True,
                      global_consumer_workers_requested=workers,
                      global_consumer_wigner_storage_requested=storage,
                      response_record_spool_used=True, response_reuse_used=False,
                      response_cache_entries=0, response_cache_bytes=0,
                      response_spool_bytes=5,
                      response_spool_kernel_output_bytes=0,
                      response_spool_kernel_sampled_angles=0,
                      global_consumer_full_angular_points=0,
                      global_consumer_angular_subset_complete=True,
                      global_consumer_full_transforms_performed=True,
                      global_consumer_frequency_count=53,
                      global_consumer_frequency_side_transforms_executed=104,
                      global_consumer_frequency_side_zero_mixing_skipped=2,
                      global_consumer_shared_wigner_bytes=0,
                      global_consumer_estimated_working_bytes=1000,
                      global_consumer_same_spin_batch_size=1)
        return result

    @staticmethod
    def _command(degree=3, workers=1, storage='generated', indices=(0, 2)):
        args = ['/usr/bin/env', 'DSPEC_SYMBOLIC_REUSE=0', 'physical-driver',
                '/model/PREM500', '/job/data', 'kernel']
        values = (
            ('--lmax', str(degree)), ('--nq', '6'), ('--radial-traces', '0'),
            ('--dc', 'include'), ('--contacts', 'off'), ('--pole-only', 'on'),
            ('--forward-only', 'off'), ('--response-reuse', 'off'),
            ('--response-record-spool', 'on'), ('--global-workers', str(workers)),
            ('--global-wigner-storage', storage),
            ('--global-radial-sides', ','.join(str(i) for i in indices)),
        )
        for key, value in values:
            args.extend((key, value))
        return ({'command': args}, {}, 'state-hash')

    @classmethod
    def _fixture(cls, data: Path, degree=3, workers=1, storage='generated',
                 indices=(0, 2)):
        data.mkdir(parents=True, exist_ok=True)
        with (data / 'radial_nodes.csv').open('w') as stream:
            stream.write('element,node,radius_m,element_width_m,fluid\n')
            for i in range(179):
                stream.write(f'{i},1,{1000+i},10,0\n')
        identities = [(i, 1) for i in indices]
        metadata = [dict(element=i, node=1, radius_m=float(1000+i),
                         radius_model=float(1000+i)/2., width_m=10., fluid=False)
                    for i in indices]
        ntheta, nphi = 2*degree+1, 4*degree+1
        angles = ntheta*nphi
        output_bytes = len(indices)*angles*6*8
        (data / 'kernel.bin').write_bytes(bytes(output_bytes))
        (data / 'pole_response_spool.bin').write_bytes(b'spool')
        result = cls._result(degree, workers, storage, len(indices))
        shared = 0 if storage == 'generated' else 640
        result.update(response_spool_kernel_output_bytes=output_bytes,
                      response_spool_kernel_sampled_angles=angles,
                      global_consumer_full_angular_points=angles,
                      global_consumer_shared_wigner_bytes=shared)
        manifest = dict(
            mode='compact_response_global_gsh_fft_frequency_accumulator',
            degree=degree,
            angular_rule='ForBand(L,2,2.0) nodes and FastFFTSize(4L+1) longitudes',
            n_theta=ntheta, n_phi=nphi, full_angular_points=angles,
            sampled_angular_points=angles, angular_subset_is_complete=True,
            radial_side_count=len(indices), native_side_indices=list(indices),
            native_side_metadata=metadata, frequency_count=53,
            frequency_records=cls._records(),
            source_contacts_included=False,
            kernel_is_full_earth_volume_integral=False,
            output_layout='radial-side-major; theta-major; point-major channels A,C,F,L,N,kappa',
            output_bytes=output_bytes, binary_file='kernel.bin',
            wigner_storage_policy=storage, shared_wigner_table_bytes=shared,
            same_spin_channel_batch_size=1,
            known_working_bytes_per_worker=100,
            estimated_working_bytes_per_worker=200,
            estimated_total_working_bytes=1000,
            spool_offset_validation_scratch_bytes=50,
            workers=workers, elapsed_s=3., timings=dict(
                serial_grid_and_context_setup_wall_s=.5,
                frequency_side_transforms_executed=104,
                frequency_side_exact_zero_mixing_skipped=2,
                worker_execution_wall_s=2.5, total_wall_s=3.,
                worker_elapsed_sum_s=.8*workers,
                worker_stage_elapsed_sums=dict(
                    spool_reads=.2*workers, coefficient_rotation_layout=.2*workers,
                    gsh_inverse_fft=.2*workers, point_kernel_contractions=.05*workers,
                    frequency_accumulation=.1*workers, output_writes=.05*workers),
                workers=[dict(side_elapsed_s=.8, spool_reads_s=.2,
                              coefficient_rotation_layout_s=.2,
                              gsh_inverse_fft_s=.2, point_kernel_contractions_s=.05,
                              frequency_accumulation_s=.1, output_writes_s=.05)
                         for _ in range(workers)]))
        receipt_hashes = {f'data/{name}': 'synthetic' for name in spool.GLOBAL_FILES}
        return result, manifest, ({'command': cls._command(degree, workers, storage, indices)[0]['command']},
                                  receipt_hashes, 'state-hash'), identities

    @classmethod
    def _write_completed_job(cls, root: Path, job_id: str, storage: str):
        job_dir = root / job_id
        data = job_dir / 'data'
        result, manifest, _, _ = cls._fixture(data, storage=storage)
        (data / 'trace.csv').write_text('sample,value\n0,1.25\n1,-0.5\n')
        (data / 'spectrum.csv').write_text('frequency,real,imag\n0,2,0\n1,1,-0.25\n')
        (data / 'result.json').write_text(json.dumps(result, allow_nan=False) + '\n')
        (data / 'kernel.bin.json').write_text(json.dumps(manifest, allow_nan=False) + '\n')
        factorisations = 53 * (2 * result['degree'] + 1)
        profile = {'solver': 'dspec', 'objective': -1,
                   'stages': {'physical_case_total': {
                       'exclusive_s': 1.0, 'inclusive_s': 1.0, 'calls': 1}},
                   'counters': {'sparse_factorization_calls': factorisations}}
        (job_dir / 'job.log').write_text('MATHPROFILE ' + json.dumps(profile) + '\n')
        outputs = ['data/' + name for name in spool.GLOBAL_FILES]
        output_hashes = {name: pair.sha(root / job_id / name) for name in outputs}
        command = cls._command(storage=storage)[0]['command']
        job = dict(id=job_id, status='completed', returncode=0,
                   inputs_unchanged=True, cleanup={'verified_gone': True},
                   failure_reasons=[], outputs=outputs, outputs_sha256=output_hashes,
                   log='job.log', log_sha256=pair.sha(job_dir / 'job.log'),
                   command=command)
        return job

    def test_supported_degree_matrix_and_physical_count_schema(self):
        for degree in check.DEGREES:
            check.validate_physical_pair(self._result(degree), self._result(degree))
        bad = self._result(64)
        with self.assertRaisesRegex(ValueError, 'degree'):
            check.validate_physical_pair(bad, bad)

    def test_generated_stored_and_worker_pair_expectations_are_explicit(self):
        for workers, storage_a, storage_b in ((1, 'generated', 'stored'),
                                               (2, 'generated', 'stored'),
                                               (1, 'generated', 'generated'),
                                               (2, 'generated', 'generated'),
                                               (4, 'generated', 'generated')):
            with tempfile.TemporaryDirectory(prefix='global-cache-mode-') as tmp:
                data = Path(tmp)
                result, manifest, receipt, _ = self._fixture(
                    data, workers=workers, storage=storage_a)
                timing = spool.validate_worker_timings(check._final_timing_adapter(manifest))
                got = check.validate_invocation(result, manifest, timing, receipt,
                                                workers, storage_a, [0, 2])
                self.assertEqual(got['actual_workers'], workers)
                other, other_manifest, other_receipt, _ = self._fixture(
                    data / 'second', workers=workers, storage=storage_b)
                other_timing = spool.validate_worker_timings(
                    check._final_timing_adapter(other_manifest))
                got2 = check.validate_invocation(other, other_manifest, other_timing,
                                                 other_receipt, workers, storage_b, [0, 2])
                self.assertEqual(got2['storage_mode'], storage_b)

    def test_parallel_generated_jobs_can_be_compared_across_worker_counts(self):
        for candidate_workers in (2, 4):
            with tempfile.TemporaryDirectory(prefix='global-cache-workers-') as tmp:
                data = Path(tmp)
                for workers in (1, candidate_workers):
                    result, manifest, receipt, _ = self._fixture(
                        data / str(workers), workers=workers)
                    timing = spool.validate_worker_timings(
                        check._final_timing_adapter(manifest))
                    checked = check.validate_invocation(
                        result, manifest, timing, receipt, workers, 'generated', [0, 2])
                    self.assertEqual(checked['actual_workers'], workers)

    def test_result_counters_must_match_sidecar_and_storage_bytes(self):
        with tempfile.TemporaryDirectory(prefix='global-cache-counters-') as tmp:
            result, manifest, receipt, _ = self._fixture(Path(tmp))
            timing = spool.validate_worker_timings(check._final_timing_adapter(manifest))
            result['global_consumer_frequency_side_transforms_executed'] -= 1
            with self.assertRaisesRegex(ValueError, 'frequency-side transform counters'):
                check.validate_invocation(result, manifest, timing, receipt, 1,
                                          'generated', [0, 2])

    def test_mutated_native_side_metadata_and_order_are_rejected(self):
        with tempfile.TemporaryDirectory(prefix='global-cache-native-') as tmp:
            result, manifest, receipt, _ = self._fixture(Path(tmp))
            manifest['native_side_metadata'][0]['element'] = 99
            with self.assertRaisesRegex(ValueError, 'disagrees with candidate'):
                spool.validate_global_manifest(manifest, result, Path(tmp), receipt, False)

    def test_nonfinite_q_bad_basis_and_bad_energy_are_rejected(self):
        records = self._records()
        records[5]['adjoint_mixing'] = [float('nan'), 0.]
        with self.assertRaisesRegex(ValueError, 'complex field'):
            spool.validate_frequency_records(records, 53)
        records = self._records()
        records[8]['receiver_basis_reconstruction_relative'] = 1.1e-12
        with self.assertRaisesRegex(ValueError, '1e-12'):
            spool.validate_frequency_records(records, 53)
        energy = dict(energy_m2_s=2., energy_adjoint_qdotp=3.,
                      energy_adjoint_relative_error_vs_2E=0.)
        with self.assertRaisesRegex(ValueError, 'recomputed'):
            spool.validate_energy(energy)

    def test_storage_and_workers_must_match_command_result_and_manifest(self):
        with tempfile.TemporaryDirectory(prefix='global-cache-invocation-') as tmp:
            result, manifest, receipt, _ = self._fixture(Path(tmp))
            timing = spool.validate_worker_timings(check._final_timing_adapter(manifest))
            result['global_consumer_wigner_storage_requested'] = 'stored'
            with self.assertRaisesRegex(ValueError, 'result Wigner storage'):
                check.validate_invocation(result, manifest, timing, receipt, 1,
                                          'generated', [0, 2])
            result['global_consumer_wigner_storage_requested'] = 'generated'
            manifest['workers'] = 2
            with self.assertRaisesRegex(ValueError, 'worker counts'):
                check.validate_invocation(result, manifest, timing, receipt, 1,
                                          'generated', [0, 2])

    def test_kernel_channel_regression_exceeds_unchanged_numeric_gate(self):
        with tempfile.TemporaryDirectory(prefix='global-cache-channel-') as tmp:
            a, b = Path(tmp)/'a.bin', Path(tmp)/'b.bin'
            reference = np.ones((2, 5, 6), dtype=np.float64)
            candidate = reference.copy()
            candidate[1, 2, 4] += 1e-10
            reference.tofile(a); candidate.tofile(b)
            comparisons = check.compare_kernel(a, b, 2, 5)
            by_channel = {row['channel']: row for row in comparisons}
            self.assertTrue(by_channel['A']['passed'])
            self.assertFalse(by_channel['N']['passed'])

    def test_pair_requires_exact_audits_and_matching_grid_dimensions(self):
        with tempfile.TemporaryDirectory(prefix='global-cache-structure-') as tmp:
            root = Path(tmp)
            a, ma, _, ids_a = self._fixture(root/'a')
            b, mb, _, ids_b = self._fixture(root/'b')
            side_a = spool.validate_global_manifest(ma, a, root/'a',
                ({}, {f'data/{x}':'ok' for x in spool.GLOBAL_FILES}, ''), False)
            side_b = spool.validate_global_manifest(mb, b, root/'b',
                ({}, {f'data/{x}':'ok' for x in spool.GLOBAL_FILES}, ''), False)
            check.validate_pair_structure(ma, mb, a, b, side_a, side_b)
            mb['frequency_records'][0]['z'][0] += 1.
            with self.assertRaisesRegex(ValueError, 'frequency audits differ'):
                check.validate_pair_structure(ma, mb, a, b, side_a, side_b)

    def test_assess_end_to_end_with_completed_receipts_and_exact_outputs(self):
        with tempfile.TemporaryDirectory(prefix='global-cache-assess-') as tmp:
            root = Path(tmp)
            jobs = [self._write_completed_job(root, 'generated', 'generated'),
                    self._write_completed_job(root, 'stored', 'stored')]
            (root / 'state.json').write_text(json.dumps({'jobs': jobs}) + '\n')
            report = check.assess(root, 'generated', root, 'stored',
                                  1, 1, 'generated', 'stored', exact=True)
            self.assertEqual(report['status'], 'PASS')
            self.assertTrue(all(report['byte_exact_files'].values()))
            self.assertEqual(report['comparison']['kernel.bin'][0]['scaled_max_difference'], 0.)
            self.assertEqual(report['invocations'][1]['storage_mode'], 'stored')


if __name__ == '__main__':
    unittest.main()
