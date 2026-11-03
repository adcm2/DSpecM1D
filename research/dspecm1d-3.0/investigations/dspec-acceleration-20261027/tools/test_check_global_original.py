"""End-to-end synthetic receipt check for the direct-to-global wrapper."""
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

import check_global_original as wrapper
import check_global_pair as reviewed
import check_pair


def _profile():
    return {'solver': 'dspec', 'objective': -1,
            'stages': {'physical_case_total': {
                'exclusive_s': 1.0, 'inclusive_s': 1.0, 'calls': 1}},
            'counters': {'sparse_factorization_calls': 53 * 7}}


def _physical(mode='kernel'):
    result = {key: 1 for key in check_pair.PHYSICAL}
    result.update(mode=mode, degree=3, sem_degree=6,
                  native_element_node_sides=180, frequency_bins=52,
                  source_contact_enabled=False, forward_only_baseline=False,
                  energy_m2_s=2.0, energy_adjoint_qdotp=4.0,
                  energy_adjoint_relative_error_vs_2E=0.0,
                  model_length_scale_m=2.0)
    return result


def _records():
    rows = [dict(spool_index=i, original_bin=i, z=[float(i), -0.1],
                 adjoint_mixing=[0.25, -0.5],
                 receiver_basis_reconstruction_relative=1e-13)
            for i in range(53)]
    rows[-1]['adjoint_mixing'] = [0.0, 0.0]
    return rows


def _write_common(data: Path, result: dict):
    data.mkdir(parents=True, exist_ok=True)
    (data / 'result.json').write_text(json.dumps(result, allow_nan=False) + '\n')
    (data / 'trace.csv').write_text('sample,value\n0,1.25\n1,-0.5\n')
    (data / 'spectrum.csv').write_text('frequency,real,imag\n0,2,0\n1,1,-0.25\n')


def _write_direct(root: Path, job_id: str, selected_tensor: np.ndarray):
    job = root / job_id
    data = job / 'data'
    result = _physical()
    _write_common(data, result)
    rows = [(i, 1, float(1000+i), 10.0, 0) for i in range(179)]
    with (data / 'radial_nodes.csv').open('w') as stream:
        stream.write('element,node,radius_m,element_width_m,fluid\n')
        stream.writelines(f'{e},{n},{r},{w},{f}\n' for e,n,r,w,f in rows)
    with (data / 'kernel_radial_nodes.csv').open('w') as stream:
        stream.write('radial_index,element,node,radius_m,element_width_m,fluid\n')
        stream.writelines(f'{i},{e},{n},{r},{w},{f}\n' for i,(e,n,r,w,f) in enumerate(rows))
    theta = np.arange(1, 8, dtype=float) * np.pi / 8
    with (data / 'angular_grid.csv').open('w') as stream:
        stream.write('theta_index,theta_rad,theta_weight\n')
        stream.writelines(f'{i},{x},{2/7}\n' for i,x in enumerate(theta))
    direct = selected_tensor
    direct.tofile(data / 'kernel.bin')
    ntheta, nphi, angles, sides = 7, 13, 91, 179
    domain = dict(
        source_contact_enabled=False, actual_contact_binary_bytes=0,
        full_target_contact_binary_bytes=0, selected_radial_trace_count=sides,
        available_radial_trace_sides=179, n_theta=ntheta, n_phi=nphi,
        full_target_angular_points_per_trace=angles,
        sampled_angular_points_per_trace=angles,
        sampled_angle_subset_is_complete=True,
        actual_sampled_point_count=sides*angles,
        full_target_point_count=sides*angles,
        actual_kernel_binary_bytes=sides*angles*6*8,
        full_target_kernel_binary_bytes=sides*angles*6*8,
        actual_combined_output_bytes=sides*angles*6*8,
        binary_layout='angle_index outer; selected radial side inner; doubles[A,C,F,L,N,kappa]',
        angular_rule='ForBand(L,2,2.0) complete grid')
    (data / 'kernel_domain.json').write_text(json.dumps(domain) + '\n')
    _write_job_receipt(root, job_id, reviewed.REFERENCE_FILES + ('angular_grid.csv',),
                       ['/usr/bin/env', 'direct-driver', 'kernel'])


def _write_global(root: Path, job_id: str, selected_tensor: np.ndarray):
    job = root / job_id
    data = job / 'data'
    result = _physical()
    result.update(global_consumer_radial_side_count=8,
                  global_consumer_side_selection_explicit=True,
                  response_spool_bytes=5, response_record_spool_used=True,
                  response_reuse_used=False, response_cache_entries=0,
                  response_cache_bytes=0, global_consumer_workers_requested=1)
    _write_common(data, result)
    rows = [(i, 1, float(1000+i), 10.0, 0) for i in range(179)]
    with (data / 'radial_nodes.csv').open('w') as stream:
        stream.write('element,node,radius_m,element_width_m,fluid\n')
        stream.writelines(f'{e},{n},{r},{w},{f}\n' for e,n,r,w,f in rows)
    (data / 'pole_response_spool.bin').write_bytes(b'spool')
    selected_tensor.tofile(data / 'kernel.bin')
    metadata = [dict(element=e, node=n, radius_m=r, radius_model=r/2.,
                     width_m=w, fluid=bool(f)) for e,n,r,w,f in rows[:8]]
    manifest = dict(
        mode='compact_response_global_gsh_fft_frequency_accumulator', degree=3,
        angular_rule='ForBand(L,2,2.0) nodes and FastFFTSize(4L+1) longitudes',
        n_theta=7, n_phi=13, full_angular_points=91, sampled_angular_points=91,
        angular_subset_is_complete=True, radial_side_count=8,
        native_side_indices=list(range(8)), native_side_metadata=metadata,
        frequency_count=53, frequency_records=_records(),
        source_contacts_included=False, kernel_is_full_earth_volume_integral=False,
        output_layout='radial-side-major; theta-major; point-major channels A,C,F,L,N,kappa',
        output_bytes=8*91*6*8, binary_file='kernel.bin', workers=1, elapsed_s=3.,
        timings=dict(serial_grid_and_context_setup_wall_s=.5,
                     worker_execution_wall_s=2.5, total_wall_s=3.,
                     worker_elapsed_sum_s=.8,
                     worker_stage_elapsed_sums=dict(
                         spool_reads=.2, coefficient_rotation_layout=.2,
                         gsh_inverse_fft=.2, point_kernel_contractions=.05,
                         frequency_accumulation=.1, output_writes=.05),
                     workers=[dict(side_elapsed_s=.8, spool_reads_s=.2,
                                   coefficient_rotation_layout_s=.2, gsh_inverse_fft_s=.2,
                                   point_kernel_contractions_s=.05,
                                   frequency_accumulation_s=.1, output_writes_s=.05)]))
    (data / 'kernel.bin.json').write_text(json.dumps(manifest) + '\n')
    cmd = ['/usr/bin/env', 'global-driver', '/model', '/job/data', 'kernel',
           '--lmax', '3', '--nq', '6', '--radial-traces', '0', '--dc', 'include',
           '--contacts', 'off', '--pole-only', 'on', '--forward-only', 'off',
           '--response-reuse', 'off', '--response-record-spool', 'on',
           '--global-workers', '1', '--global-wigner-storage', 'generated',
           '--global-radial-sides', '0,1,2,3,4,5,6,7']
    _write_job_receipt(root, job_id, reviewed.GLOBAL_FILES, cmd)


def _write_job_receipt(root: Path, job_id: str, output_names, command):
    job = root / job_id
    (job / 'job.log').write_text('MATHPROFILE ' + json.dumps(_profile()) + '\n')
    outputs = ['data/' + name for name in output_names]
    job_record = dict(id=job_id, status='completed', returncode=0,
                      inputs_unchanged=True, cleanup={'verified_gone': True},
                      failure_reasons=[], outputs=outputs,
                      outputs_sha256={name: check_pair.sha(root / job_id / name) for name in outputs},
                      log='job.log', log_sha256=check_pair.sha(job / 'job.log'),
                      command=command)
    state_path = root / 'state.json'
    prior = json.loads(state_path.read_text())['jobs'] if state_path.exists() else []
    state_path.write_text(json.dumps({'jobs': prior + [job_record]}) + '\n')


class OriginalGlobalWrapperIntegration(unittest.TestCase):
    def test_full_direct_179_to_global_8_through_reviewed_assess(self):
        with tempfile.TemporaryDirectory(prefix='global-original-wrapper-') as tmp:
            root = Path(tmp)
            full = np.arange(91*179*6, dtype=np.float64).reshape(91, 179, 6)
            selected = full[:, :8, :].transpose(1, 0, 2).copy()
            _write_direct(root, 'direct', full)
            _write_global(root, 'global', selected)
            report = wrapper.assess(root, 'direct', root, 'global')
            self.assertEqual(report['status'], 'PASS')
            self.assertEqual(report['comparison_mode'], 'direct-fullgrid-to-global-spool-full')
            self.assertEqual(report['candidate_native_side_indices'], list(range(8)))
            self.assertTrue(all(row['passed'] for rows in report['comparison'].values() for row in rows))
            self.assertEqual(len(report['timing_adaptation']), 1)
            adaptation = report['timing_adaptation'][0]
            self.assertTrue(adaptation['applied'])
            self.assertTrue(adaptation['original_fields_preserved'])
            self.assertEqual(adaptation['source_value_s'], .5)
            self.assertEqual(adaptation['compatibility_value_s'], .5)
            self.assertEqual(report['tolerance'], 2e-13)
            self.assertEqual(len(report['reviewed_helper_sha256']), 64)


if __name__ == '__main__':
    unittest.main()
