#!/usr/bin/env python3
"""Apply strict field parity, then attribute nonradial reuse with radial-L0 fallback."""
import argparse
import hashlib
import json
from pathlib import Path

import check_pair


def reuse_checks(result, profile, symbolic):
    if result.get('response_reuse_requested') is not True:
        raise ValueError('direct response reuse was not explicitly requested')
    if result.get('response_reuse_used') is not True or result.get('response_record_spool_used') is not False:
        raise ValueError('direct nonradial response reuse was not exercised')
    if (result.get('response_reuse_capacity_fallback') is not False
            or result.get('response_reuse_preflight_fallback') is not False):
        raise ValueError('response reuse fell back to the reference route')
    if result.get('forward_only_baseline') is not False:
        raise ValueError('baseline must report its unit-adjoint capture')
    if result.get('symbolic_reuse_enabled') is not symbolic:
        raise ValueError('symbolic mode differs from the requested experiment')
    if result.get('frequency_bins') != 52:
        raise ValueError('reuse attribution requires the fixed 52-bin physical case')
    if result.get('response_cache_radial_l0_bypassed') is not True:
        raise ValueError('radial L0 must use the original weighted-solve path')
    l0_calls = result.get('response_cache_radial_l0_bypass_calls')
    if l0_calls != 105:
        raise ValueError('radial L0 bypass must account for 53 baseline and 52 kernel calls')

    degree = result.get('degree')
    if not isinstance(degree, int) or isinstance(degree, bool) or degree < 1:
        raise ValueError('invalid maximum degree')
    nonradial_families = 2 * degree
    expected_misses = 53 * nonradial_families
    expected_hits = 52 * nonradial_families
    if not (result.get('response_cache_entries') == expected_misses
            and result.get('response_cache_misses') == expected_misses
            and result.get('response_cache_preflight_entries') == expected_misses
            and result.get('response_cache_hits') == expected_hits):
        raise ValueError('nonradial cache entries/hits/misses differ from the fixed physical case')
    if not (0 < result.get('response_cache_bytes', 0)
            <= result.get('response_cache_preflight_estimated_bytes', 0)
            <= result.get('response_cache_budget_bytes', 0)):
        raise ValueError('cache storage exceeds its preflight/budget')
    residual = result.get('response_basis_reconstruction_max_relative')
    if not isinstance(residual, (int, float)) or not 0 <= residual <= 1e-12:
        raise ValueError('unit-receiver basis reconstruction gate failed')
    if profile is None:
        raise ValueError('reuse attribution requires an instrumented job')
    counters = profile['counters']
    factors = counters.get('sparse_factorization_calls')
    expected_factors = 53 * (2 * degree + 1) + 52
    if factors != expected_factors:
        raise ValueError('factorisation count differs from nonradial capture plus original radial-L0 solves')
    profile_l0 = counters.get('response_cache_radial_l0_bypass_calls')
    if profile_l0 != l0_calls:
        raise ValueError('radial L0 result/profile bypass counts disagree')
    analyses = counters.get('sparse_symbolic_analysis_calls')
    reuses = counters.get('sparse_symbolic_analysis_reuses', 0)
    if (not isinstance(analyses, int) or not isinstance(reuses, int)
            or analyses + reuses != factors
            or (symbolic and reuses <= 0) or (not symbolic and reuses != 0)):
        raise ValueError('symbolic analysis/reuse counts do not match factorisations')
    return {
        'nonradial_capture_misses': expected_misses,
        'nonradial_kernel_hits': expected_hits,
        'radial_l0_bypass_calls': l0_calls,
        'expected_factorisations': expected_factors,
        'observed_factorisations': factors,
        'symbolic_analyses': analyses,
        'symbolic_reuses': reuses,
        'receiver_basis_max_relative': residual,
        'cache_bytes': result['response_cache_bytes'],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('reference-root', 'reference-job', 'candidate-root', 'candidate-job', 'output'):
        parser.add_argument('--' + key, required=True)
    parser.add_argument('--symbolic', choices=('on', 'off'), required=True)
    args = parser.parse_args()
    try:
        report = check_pair.assess(Path(args.reference_root), args.reference_job,
                                   Path(args.candidate_root), args.candidate_job)
        if report['status'] == 'PASS':
            raw = (Path(args.candidate_root) / args.candidate_job / 'data/result.json').read_bytes()
            digest = hashlib.sha256(raw).hexdigest()
            if digest != report['verified_artifact_sha256'][1]['data/result.json']:
                raise ValueError('candidate result changed before reuse attribution')
            result = json.loads(raw)
            report['reuse_attribution'] = reuse_checks(
                result, report['profiles'][1], args.symbolic == 'on')
    except Exception as exc:
        report = {'status': 'FAIL', 'error_type': type(exc).__name__, 'error': str(exc)}
    with Path(args.output).open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write('\n')
    return 0 if report['status'] == 'PASS' else 1


if __name__ == '__main__':
    raise SystemExit(main())
