from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import binary32_covariance_coercivity as C


class Binary32CoercivityTests(unittest.TestCase):
    def test_uniform_format_bounds_are_strict_and_exact(self):
        b = C.bounds()
        self.assertGreater(b['covariance_lambda_min_lower'], 0)
        self.assertGreater(b['joint_storage_lambda_min_lower'], 0)
        self.assertEqual(b['joint_storage_lambda_min_lower'], 1/(21*C.MAX_FINITE))
        self.assertEqual(b['covariance_lambda_min_lower']*(21*C.MAX_FINITE)**20,
                         C.QUANTUM**21)
        self.assertFalse(b['shipping_SPD_preservation_proved'])
        self.assertFalse(b['source_uniform_rho_certified'])
        self.assertFalse(b['storage_search_allowed'])

    def test_smallest_subnormal_is_valid_and_coercive(self):
        r = C.exact_spd([[C.QUANTUM, 0], [0, C.MAX_FINITE]])
        self.assertTrue(r['exact_SPD'])
        self.assertEqual(r['determinant'], C.QUANTUM*C.MAX_FINITE)

    def test_positive_diagonal_and_symmetry_do_not_imply_SPD(self):
        with self.assertRaisesRegex(ValueError, 'not SPD'):
            C.exact_spd([[1, 2], [2, 1]])
        with self.assertRaisesRegex(ValueError, 'not SPD'):
            C.exact_spd([[1, 1], [1, 1]])

    def test_near_singular_representable_matrix_uses_exact_pivots(self):
        eps = F(1, 2**23)
        r = C.exact_spd([[1, 1], [1, 1+eps]])
        self.assertEqual(r['determinant'], eps)
        self.assertEqual(r['minimum_exact_LDL_pivot'], eps)

    def test_grid_membership_is_not_covariance_rounding(self):
        for x in (C.QUANTUM/2, F(1, 3), C.MAX_FINITE+1):
            with self.assertRaises(ValueError):
                C.binary32_value(x)
        with self.assertRaisesRegex(ValueError, 'symmetry'):
            C.exact_spd([[1, 0], [C.QUANTUM, 1]])


    def test_committed_exact_native_audit_preserves_infinite_boundary(self):
        import json
        from pathlib import Path
        path = Path(__file__).resolve().parents[2]/'reports/results/ou3_alt_storage/binary32-coercivity.json'
        report = json.loads(path.read_text())
        self.assertEqual(sum(len(h['endpoints']) for h in report['histories'].values()), 36)
        traces = report['retained_service_finite_trace_audit']
        self.assertEqual(sum(h['covariance_records_checked'] for h in traces.values()), 18586)
        for h in traces.values():
            self.assertGreater(F(h['minimum_exact_LDL_pivot']), 0)
            self.assertEqual(h['sample_boundaries'], 1200)
            self.assertFalse(h['infinite_SPD_preservation_proved'])
        self.assertFalse(report['shipping_infinite_domain_membership_proved'])


if __name__ == '__main__':
    unittest.main()
