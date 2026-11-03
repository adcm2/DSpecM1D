#!/usr/bin/env python3
"""Apply the frozen field-parity gate, then verify measured reuse was exercised."""
import argparse
import hashlib
import json
from pathlib import Path

import check_pair


def reuse_checks(result, profile, symbolic):
    if result.get('response_reuse_used') is not True or result.get('response_record_spool_used') is not False:
        raise ValueError('direct response reuse was not exercised')
    if (result.get('response_reuse_capacity_fallback') is not False
            or result.get('response_reuse_preflight_fallback') is not False):
        raise ValueError('response reuse fell back to the reference route')
    if result.get('forward_only_baseline') is not False:
        raise ValueError('baseline must report its unit-adjoint capture')
    if result.get('symbolic_reuse_enabled') is not symbolic:
        raise ValueError('symbolic mode differs from the requested experiment')
    # This fixed benchmark includes radial L0 and both nonradial families.
    if result.get('frequency_bins') != 52:
        raise ValueError('reuse attribution requires the fixed 52-bin physical case')
    families = 2 * result['degree'] + 1
    misses = (result['frequency_bins'] + 1) * families
    hits = result['frequency_bins'] * families
    if not (result['response_cache_entries'] == misses == result['response_cache_misses']
            == result['response_cache_preflight_entries'] and result['response_cache_hits'] == hits):
        raise ValueError('cache entry/hit/miss counts differ from the fixed physical case')
    if not (0 < result['response_cache_bytes'] <= result['response_cache_preflight_estimated_bytes']
            <= result['response_cache_budget_bytes']):
        raise ValueError('cache storage exceeds its preflight/budget')
    if not 0 <= result['response_basis_reconstruction_max_relative'] <= 1e-12:
        raise ValueError('unit-receiver basis reconstruction gate failed')
    if profile is None: raise ValueError('reuse attribution requires an instrumented job')
    counters = profile['counters']
    factors = counters['sparse_factorization_calls']
    if factors != misses: raise ValueError('factorisation count does not show one capture pass')
    analyses = counters['sparse_symbolic_analysis_calls']
    reuses = counters.get('sparse_symbolic_analysis_reuses', 0)
    if analyses + reuses != factors or (symbolic and reuses <= 0) or (not symbolic and reuses != 0):
        raise ValueError('symbolic analysis/reuse counts do not match the factorisations')
    return {'expected_capture_factorisations': misses, 'expected_cache_hits': hits,
            'observed_factorisations': factors, 'symbolic_analyses': analyses,
            'symbolic_reuses': reuses, 'receiver_basis_max_relative': result['response_basis_reconstruction_max_relative'],
            'cache_bytes': result['response_cache_bytes']}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('reference-root', 'reference-job', 'candidate-root', 'candidate-job', 'output'):
        p.add_argument('--' + key, required=True)
    p.add_argument('--symbolic', choices=('on', 'off'), required=True)
    a = p.parse_args()
    try:
        # Both runs use --forward-only off, so captured-unit methodology is
        # explicit and the reviewed metadata gate remains unchanged.
        report = check_pair.assess(Path(a.reference_root), a.reference_job,
                                   Path(a.candidate_root), a.candidate_job)
        if report['status'] == 'PASS':
            raw = (Path(a.candidate_root) / a.candidate_job / 'data/result.json').read_bytes()
            if hashlib.sha256(raw).hexdigest() != report['verified_artifact_sha256'][1]['data/result.json']:
                raise ValueError('candidate result changed before reuse attribution')
            result = json.loads(raw)
            report['reuse_attribution'] = reuse_checks(result, report['profiles'][1], a.symbolic == 'on')
    except Exception as e:
        report = {'status': 'FAIL', 'error_type': type(e).__name__, 'error': str(e)}
    with Path(a.output).open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False); stream.write('\n')
    return 0 if report['status'] == 'PASS' else 1


if __name__ == '__main__':
    raise SystemExit(main())
