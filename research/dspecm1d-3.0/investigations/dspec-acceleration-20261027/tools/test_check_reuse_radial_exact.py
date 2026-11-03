"""Pure-data regression checks for the radial-L0 exact reuse gate."""
import unittest

import check_reuse_radial_exact as checker


def fixture(symbolic=True):
    degree = 3
    misses = 53 * 2 * degree
    hits = 52 * 2 * degree
    factors = 53 * (2 * degree + 1) + 52
    result = dict(
        response_reuse_requested=True, response_reuse_used=True,
        response_record_spool_used=False,
        response_reuse_capacity_fallback=False,
        response_reuse_preflight_fallback=False,
        forward_only_baseline=False, symbolic_reuse_enabled=symbolic,
        frequency_bins=52, degree=degree,
        response_cache_radial_l0_bypassed=True,
        response_cache_radial_l0_bypass_calls=105,
        response_cache_entries=misses, response_cache_misses=misses,
        response_cache_preflight_entries=misses, response_cache_hits=hits,
        response_cache_bytes=500, response_cache_preflight_estimated_bytes=600,
        response_cache_budget_bytes=700,
        response_basis_reconstruction_max_relative=5e-13,
    )
    profile = {'counters': {
        'sparse_factorization_calls': factors,
        'response_cache_radial_l0_bypass_calls': 105,
        'sparse_symbolic_analysis_calls': factors - 40 if symbolic else factors,
        'sparse_symbolic_analysis_reuses': 40 if symbolic else 0,
    }}
    return result, profile


class RadialExactReuseChecks(unittest.TestCase):
    def test_valid_nonradial_cache_and_original_radial_l0_route(self):
        result, profile = fixture()
        report = checker.reuse_checks(result, profile, True)
        self.assertEqual(report['nonradial_capture_misses'], 318)
        self.assertEqual(report['nonradial_kernel_hits'], 312)
        self.assertEqual(report['radial_l0_bypass_calls'], 105)
        self.assertEqual(report['expected_factorisations'], 423)

    def test_exact_nonradial_counter_formula_rejects_single_count_drift(self):
        result, profile = fixture()
        result['response_cache_hits'] -= 1
        with self.assertRaisesRegex(ValueError, 'entries/hits/misses'):
            checker.reuse_checks(result, profile, True)
        result, profile = fixture()
        profile['counters']['sparse_factorization_calls'] += 1
        with self.assertRaisesRegex(ValueError, 'factorisation count'):
            checker.reuse_checks(result, profile, True)

    def test_radial_bypass_count_must_include_both_passes(self):
        result, profile = fixture()
        result['response_cache_radial_l0_bypass_calls'] = 53
        with self.assertRaisesRegex(ValueError, '53 baseline and 52 kernel'):
            checker.reuse_checks(result, profile, True)

    def test_missing_or_true_fallback_and_wrong_route_fail_closed(self):
        for key, value in (('response_reuse_capacity_fallback', None),
                           ('response_reuse_preflight_fallback', True)):
            result, profile = fixture()
            result[key] = value
            with self.assertRaisesRegex(ValueError, 'fell back'):
                checker.reuse_checks(result, profile, True)
        result, profile = fixture()
        result['response_record_spool_used'] = True
        with self.assertRaisesRegex(ValueError, 'not exercised'):
            checker.reuse_checks(result, profile, True)

    def test_basis_storage_and_symbolic_closure_remain_gated(self):
        result, profile = fixture()
        result['response_basis_reconstruction_max_relative'] = 1.1e-12
        with self.assertRaisesRegex(ValueError, 'basis reconstruction'):
            checker.reuse_checks(result, profile, True)
        result, profile = fixture(symbolic=False)
        with self.assertRaisesRegex(ValueError, 'symbolic mode'):
            checker.reuse_checks(result, profile, True)
        result, profile = fixture()
        profile['counters']['sparse_symbolic_analysis_reuses'] = 0
        profile['counters']['sparse_symbolic_analysis_calls'] = 423
        with self.assertRaisesRegex(ValueError, 'symbolic analysis/reuse'):
            checker.reuse_checks(result, profile, True)


if __name__ == '__main__':
    unittest.main()
