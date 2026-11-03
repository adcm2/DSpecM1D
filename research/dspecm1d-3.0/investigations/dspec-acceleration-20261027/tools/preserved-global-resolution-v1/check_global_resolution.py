#!/usr/bin/env python3
"""Describe L128/L384 global-kernel traces on a masked source-receiver great circle."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

import check_global_pair as global_check
import check_pair

CHANNELS = ('A', 'C', 'F', 'L', 'N', 'kappa')
SIDE_INDICES = (59, 100)
POINT_COUNT = 720
ENDPOINT_MASK_DEG = 5.0
EXPECTED_SOURCE = (10.0, 20.0)
EXPECTED_RECEIVER = (45.0, 75.0)


def great_circle(source_latlon: tuple[float, float], receiver_latlon: tuple[float, float],
                 count: int = POINT_COUNT,
                 endpoint_mask_deg: float = ENDPOINT_MASK_DEG) -> dict[str, Any]:
    if count != POINT_COUNT or endpoint_mask_deg != ENDPOINT_MASK_DEG:
        raise ValueError('resolution comparison uses its fixed 720-point, 5-degree mask')
    def unit(latlon):
        lat, lon = (math.radians(float(x)) for x in latlon)
        if not (math.isfinite(lat) and math.isfinite(lon) and -math.pi/2 <= lat <= math.pi/2):
            raise ValueError('invalid great-circle endpoint')
        return np.array([math.cos(lat)*math.cos(lon), math.cos(lat)*math.sin(lon), math.sin(lat)])
    start, end = unit(source_latlon), unit(receiver_latlon)
    angle = math.acos(float(np.clip(np.dot(start, end), -1.0, 1.0)))
    if angle <= 1e-14 or abs(math.sin(angle)) <= 1e-14:
        raise ValueError('great-circle endpoints are coincident or antipodal')
    fraction = np.linspace(0.0, 1.0, count)
    vectors = (np.sin((1-fraction)*angle)[:, None] / math.sin(angle) * start
               + np.sin(fraction*angle)[:, None] / math.sin(angle) * end)
    vectors /= np.linalg.norm(vectors, axis=1)[:, None]
    lat = np.arcsin(np.clip(vectors[:, 2], -1.0, 1.0))
    lon = np.mod(np.arctan2(vectors[:, 1], vectors[:, 0]), 2*np.pi)
    theta = np.pi/2 - lat
    endpoint_mask = math.radians(endpoint_mask_deg)
    arc_distance = fraction * angle
    keep = (arc_distance >= endpoint_mask) & ((angle - arc_distance) >= endpoint_mask)
    if np.count_nonzero(keep) == 0:
        raise ValueError('endpoint mask removes the entire great-circle arc')
    return {'theta': theta[keep], 'phi': lon[keep], 'fraction': fraction[keep],
            'arc_distance_deg': math.degrees(angle), 'input_points': count,
            'included_points': int(np.count_nonzero(keep)),
            'endpoint_mask_deg': endpoint_mask_deg,
            'first_kept_distance_deg': float(math.degrees(arc_distance[keep][0])),
            'last_kept_distance_to_receiver_deg': float(math.degrees((angle-arc_distance[keep])[-1]))}


def bilinear_interpolate(values: np.ndarray, theta_nodes: np.ndarray,
                         theta: np.ndarray, phi: np.ndarray) -> np.ndarray:
    if values.ndim != 3 or values.shape[0] != len(theta_nodes):
        raise ValueError('interpolation grid must have shape (theta, phi, channel)')
    ntheta, nphi, _ = values.shape
    if ntheta < 2 or nphi < 2 or np.any(np.diff(theta_nodes) <= 0):
        raise ValueError('interpolation requires increasing theta nodes and at least two phi nodes')
    if theta.shape != phi.shape or not np.isfinite(theta).all() or not np.isfinite(phi).all():
        raise ValueError('invalid interpolation query coordinates')
    if np.any(theta < theta_nodes[0]) or np.any(theta > theta_nodes[-1]):
        raise ValueError('great-circle theta lies outside the exact reference grid')
    upper = np.searchsorted(theta_nodes, theta, side='right')
    lower = np.clip(upper - 1, 0, ntheta - 2)
    upper = lower + 1
    ft = (theta - theta_nodes[lower]) / (theta_nodes[upper] - theta_nodes[lower])
    phi_coordinate = np.mod(phi, 2*np.pi) * (nphi / (2*np.pi))
    p0 = np.floor(phi_coordinate).astype(np.int64) % nphi
    p1 = (p0 + 1) % nphi
    fp = phi_coordinate - np.floor(phi_coordinate)
    v00, v01 = values[lower, p0, :], values[lower, p1, :]
    v10, v11 = values[upper, p0, :], values[upper, p1, :]
    low = v00 * (1-fp[:, None]) + v01 * fp[:, None]
    high = v10 * (1-fp[:, None]) + v11 * fp[:, None]
    out = low * (1-ft[:, None]) + high * ft[:, None]
    if not np.isfinite(out).all():
        raise ValueError('nonfinite bilinear interpolation result')
    return out


def read_theta_grid(root: Path, job_id: str, expected_degree: int,
                    receipt: tuple | None = None) -> tuple[np.ndarray, tuple]:
    receipt = receipt or check_pair.verified_job(root, job_id)
    if 'data/angular_grid.csv' not in receipt[1]:
        raise ValueError('reference angular grid is not receipt-bound')
    data = root / job_id / 'data'
    result = check_pair.load(data / 'result.json')
    if result.get('mode') != 'kernel' or result.get('degree') != expected_degree:
        raise ValueError('theta-grid reference must be a full kernel job at the matching degree')
    header, rows = check_pair.csv_array(data / 'angular_grid.csv')
    ntheta = 2*expected_degree+1
    if (header != ('theta_index', 'theta_rad', 'theta_weight') or rows.shape != (ntheta, 3)
            or not np.array_equal(rows[:, 0], np.arange(ntheta))
            or np.any(rows[:, 1] <= 0) or np.any(rows[:, 1] >= np.pi)
            or np.any(np.diff(rows[:, 1]) <= 0) or np.any(rows[:, 2] <= 0)
            or abs(float(np.sum(rows[:, 2])) - 2.0) > 1e-12):
        raise ValueError('invalid receipt-bound exact theta grid')
    return rows[:, 1].copy(), receipt


def _physical_pair(result128: dict, result384: dict) -> None:
    for key in check_pair.PHYSICAL:
        if key == 'degree':
            continue
        if key not in result128 or key not in result384 or result128[key] != result384[key]:
            raise ValueError('L128/L384 physical metadata differs or is missing: ' + key)
    if (result128.get('degree') != 128 or result384.get('degree') != 384
            or result128.get('mode') != 'kernel' or result384.get('mode') != 'kernel'
            or result128.get('sem_degree') != 6 or result384.get('sem_degree') != 6
            or result128.get('native_element_node_sides') != 180
            or result384.get('native_element_node_sides') != 180
            or result128.get('frequency_bins') != 52 or result384.get('frequency_bins') != 52
            or result128.get('source_contact_enabled') is not False
            or result384.get('source_contact_enabled') is not False):
        raise ValueError('resolution pair is not the required contacts-off L128/L384 SEM-6 case')
    if tuple(result128.get('source', [])[:2]) != EXPECTED_SOURCE or tuple(result384.get('source', [])[:2]) != EXPECTED_SOURCE:
        raise ValueError('source must be 10 N, 20 E in both results')
    if tuple(result128.get('receiver', [])[:2]) != EXPECTED_RECEIVER or tuple(result384.get('receiver', [])[:2]) != EXPECTED_RECEIVER:
        raise ValueError('receiver must be 45 N, 75 E in both results')
    global_check.validate_energy(result128)
    global_check.validate_energy(result384)


def _argv_options(command: list[str]) -> dict[str, str]:
    options = {}
    for index, token in enumerate(command[:-1]):
        if token.startswith('--'):
            options[token] = command[index + 1]
    return options


def _validate_serial_generated(result: dict, manifest: dict,
                               sides: dict, degree: int,
                               receipt_job: dict) -> None:
    if (manifest.get('degree') != degree or manifest.get('workers') != 1
            or result.get('global_consumer_workers_requested') != 1
            or result.get('global_consumer_wigner_storage_requested') != 'generated'
            or manifest.get('wigner_storage_policy') != 'generated'
            or manifest.get('fftw_planner_policy') != 'FFTWpp::Estimate'):
        raise ValueError('resolution comparison requires serial Generated/Estimate global jobs')
    command = receipt_job.get('command')
    if not isinstance(command, list):
        raise ValueError('completed-job receipt lacks its argv command')
    options = _argv_options(command)
    if (options.get('--lmax') != str(degree)
            or options.get('--global-workers') != '1'
            or options.get('--global-radial-sides') != ','.join(map(str, SIDE_INDICES))
            or options.get('--global-wigner-storage') != 'generated'
            or options.get('--forward-only') != 'off'
            or options.get('--contacts') != 'off'
            or options.get('--pole-only') != 'on'):
        raise ValueError('receipted command differs from the required full-grid physical case')
    if (sides['native_side_indices'] != list(SIDE_INDICES)
            or sides['native_side_metadata'][0].get('fluid') is not True
            or sides['native_side_metadata'][1].get('fluid') is not False):
        raise ValueError('global output must contain native sides 59-fluid and 100-solid in order')
    if (manifest.get('sampled_angular_points') != manifest.get('full_angular_points')
            or not manifest.get('angular_subset_is_complete')):
        raise ValueError('resolution comparison requires complete global angular grids')


def _channel_metrics(a: np.ndarray, b: np.ndarray) -> list[dict]:
    if a.shape != b.shape or a.ndim != 2 or a.shape[1] != len(CHANNELS):
        raise ValueError('interpolated native trace shapes differ')
    if not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError('nonfinite interpolated native trace')
    metrics = []
    for channel, name in enumerate(CHANNELS):
        x, y = a[:, channel], b[:, channel]
        diff = x-y
        max_abs = float(np.max(np.abs(diff)))
        scale = float(max(np.max(np.abs(x)), np.max(np.abs(y))))
        norm_scale = float(max(np.linalg.norm(x), np.linalg.norm(y)))
        structural_zero = bool(np.all(x == 0.0) and np.all(y == 0.0))
        if scale == 0.0 or norm_scale == 0.0:
            if not structural_zero:
                raise ValueError('zero comparison scale without a paired structural zero')
            max_scaled = scaled_l2 = 0.0
        else:
            max_scaled = max_abs/scale
            scaled_l2 = float(np.linalg.norm(diff)/norm_scale)
        metrics.append({'channel': name, 'max_abs_difference': max_abs,
                        'max_scaled_difference': max_scaled,
                        'norm_scaled_l2_difference': scaled_l2,
                        'scale_max_abs': scale, 'scale_max_l2': norm_scale,
                        'structural_zero': structural_zero})
    return metrics


def _sample_job(root: Path, job_id: str, expected_degree: int,
                theta_nodes: np.ndarray) -> tuple[dict, dict, dict, np.ndarray, tuple]:
    receipt = check_pair.verified_job(root, job_id)
    data = root / job_id / 'data'
    result = check_pair.load(data / 'result.json')
    manifest = check_pair.load(data / 'kernel.bin.json')
    profile = check_pair.profile(root / job_id / receipt[0]['log'])
    if profile is None:
        raise ValueError('global resolution job lacks profiler evidence')
    global_check.validate_global_result(result, profile)
    sides = global_check.validate_global_manifest(manifest, result, data, receipt, probe64=False)
    _validate_serial_generated(result, manifest, sides, expected_degree, receipt[0])
    ntheta, nphi = manifest['n_theta'], manifest['n_phi']
    if len(theta_nodes) != ntheta:
        raise ValueError('receipt-bound layout theta grid differs from global grid dimension')
    expected = sides['radial_side_count']*ntheta*nphi*len(CHANNELS)*8
    if (manifest.get('output_bytes') != expected
            or (data / 'kernel.bin').stat().st_size != expected):
        raise ValueError('kernel byte count does not match declared side/theta/phi/channel dimensions')
    kernel = np.memmap(data / 'kernel.bin', dtype=np.float64, mode='r',
                       shape=(sides['radial_side_count'], ntheta, nphi, len(CHANNELS)))
    return result, manifest, sides, kernel, receipt


def compare(l128_root: Path, l128_job: str, l384_root: Path, l384_job: str,
            layout_root: Path, theta128_job: str, theta384_job: str) -> dict[str, Any]:
    l128_receipt = check_pair.verified_job(l128_root, l128_job)
    l384_receipt = check_pair.verified_job(l384_root, l384_job)
    theta128, grid128_receipt = read_theta_grid(layout_root, theta128_job, 128)
    theta384, grid384_receipt = read_theta_grid(layout_root, theta384_job, 384)
    result128, manifest128, sides128, kernel128, l128_receipt = _sample_job(
        l128_root, l128_job, 128, theta128)
    result384, manifest384, sides384, kernel384, l384_receipt = _sample_job(
        l384_root, l384_job, 384, theta384)
    _physical_pair(result128, result384)
    if sides128['native_identities'] != sides384['native_identities']:
        raise ValueError('L128/L384 selected native-side identities/order differ')
    if sides128['native_side_metadata'] != sides384['native_side_metadata']:
        raise ValueError('L128/L384 selected native-side radii/materials differ')
    path = great_circle(EXPECTED_SOURCE, EXPECTED_RECEIVER)
    traces: dict[str, dict[str, Any]] = {}
    for side_ordinal, side_metadata in enumerate(sides128['native_side_metadata']):
        a = bilinear_interpolate(kernel128[side_ordinal], theta128, path['theta'], path['phi'])
        b = bilinear_interpolate(kernel384[side_ordinal], theta384, path['theta'], path['phi'])
        if not np.any(a) and not np.any(b):
            raise ValueError('sampled native trace is all zero at both resolutions')
        traces[str(SIDE_INDICES[side_ordinal])] = {
            'native_identity': sides128['native_identities'][side_ordinal],
            'radius_m': side_metadata['radius_m'], 'fluid': side_metadata['fluid'],
            'channels': _channel_metrics(a, b),
        }
    receipts = [(l128_root, l128_job, l128_receipt),
                (l384_root, l384_job, l384_receipt),
                (layout_root, theta128_job, grid128_receipt),
                (layout_root, theta384_job, grid384_receipt)]
    after = global_check._unchanged([x[0] for x in receipts], [x[1] for x in receipts],
                                    [x[2] for x in receipts])
    return {
        'status': 'DESCRIPTIVE_COMPARISON',
        'acceptance_gate': None,
        'interpretation': ('Approximate bilinear interpolation along the specified great-circle section; '
                           'descriptive only and not a convergence certificate.'),
        'degrees': {'coarse': 128, 'fine': 384},
        'theta_grid_sources': {
            'l128': {'job': theta128_job, 'sha256': grid128_receipt[1]['data/angular_grid.csv']},
            'l384': {'job': theta384_job, 'sha256': grid384_receipt[1]['data/angular_grid.csv']},
        },
        'kernel_sources': {
            'l128': {'job': l128_job, 'sha256': l128_receipt[1]['data/kernel.bin']},
            'l384': {'job': l384_job, 'sha256': l384_receipt[1]['data/kernel.bin']},
        },
        'path': {'source_lat_lon_deg': list(EXPECTED_SOURCE),
                 'receiver_lat_lon_deg': list(EXPECTED_RECEIVER), **path},
        'angular_grid_rule': 'longitude phi=2*pi*k/n_phi; theta from receipt-bound reference angular_grid.csv',
        'output_layout': 'radial-side-major, theta-major, phi-major, channels A,C,F,L,N,kappa',
        'traces': traces,
        'state_sha256_at_read': [receipt[2] for receipt in after],
        'verified_artifact_sha256': [receipt[1] for receipt in after],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('l128-root', 'l128-job', 'l384-root', 'l384-job',
                 'layout-root', 'theta128-job', 'theta384-job', 'output'):
        parser.add_argument('--' + name, required=True)
    args = parser.parse_args()
    try:
        report = compare(Path(args.l128_root), args.l128_job,
                         Path(args.l384_root), args.l384_job,
                         Path(args.layout_root), args.theta128_job,
                         args.theta384_job)
    except Exception as exc:
        report = {'status': 'FAIL', 'error_type': type(exc).__name__, 'error': str(exc),
                  'acceptance_gate': None,
                  'interpretation': 'No scientific acceptance or convergence conclusion is produced.'}
    with Path(args.output).open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write('\n')
    return 0 if report['status'] == 'DESCRIPTIVE_COMPARISON' else 1


if __name__ == '__main__':
    raise SystemExit(main())
