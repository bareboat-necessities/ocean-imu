from __future__ import annotations

from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'tools'/'stability'))
import ou3_brmm_infinite_continuation as C


class InfiniteContinuationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report = C.build()

    def test_full_radial_source_column_not_an_independent_S_ball(self):
        h = C.symbol('h')
        for n in (18, 21, 24):
            root = C.ambiguity_map(n, {})
            after = C.ambiguity_map(n, h)
            for j in range(3):
                self.assertEqual(root[9+j][j], C.constant(1))
                self.assertEqual(root[12+j][j], {})
                self.assertEqual(after[12+j][j], h)
            for i in list(range(9))+list(range(15, n)):
                self.assertTrue(all(not x for x in after[i]))
        with self.assertRaises(ValueError):
            C.ambiguity_map(3, h)

    def test_all_event_identities_are_coefficientwise_exact(self):
        records = self.report['literal_event_certificates']
        self.assertEqual([r['dimension'] for r in records], [18, 21, 24])
        self.assertEqual(sum(r['polynomial_identities_checked'] for r in records), 945)
        for record in records:
            self.assertEqual(record['nonzero_polynomial_residuals'], [])
            for rows in record['Joseph_residual_root_coefficients'].values():
                self.assertTrue(all(not x for row in rows for x in row))

    def test_dropping_or_reversing_physical_S_destroys_cancellation(self):
        for n in (18, 21, 24):
            for bad_sign in (0, -1):
                bad = C.event_certificate(n, physical_S_sign=bad_sign)
                self.assertTrue(bad['nonzero_polynomial_residuals'])
                self.assertNotEqual(bad['Joseph_residual_root_coefficients']['S_zero'][0][0], {})

    def test_one_common_estimate_cannot_bound_both_histories(self):
        for h in (F(0), F(3), F(2400), F(480001, 200)):
            for estimate in (F(-10000), F(-1, 7), F(0), F(19, 3), F(10000)):
                p, m = C.paired_errors(h, F(1, 8), estimate)
                self.assertEqual(p-m, h/4)
                self.assertGreaterEqual(max(p*p, m*m), h*h/64)
                self.assertEqual(p*p+m*m, 2*estimate*estimate+h*h/32)
        # This is a retention radius, never a replacement entry radius.
        h = F(480001, 200)
        p, m = C.paired_errors(h, F(1, 8), F(7, 9))
        self.assertGreater(max(abs(p), abs(m)), 300)

    def test_short_window_increment_does_not_bound_absolute_supply(self):
        d = F(1, 8)
        self.assertEqual(3*d, F(3, 8))
        for h in (F(0), F(3), F(3000)):
            self.assertEqual(C.source_energy(h, 3, d), (3*h*h+9*h+9)/64)
        self.assertGreater(C.source_energy(3000, 3, d), 100000)
        with self.assertRaises(ValueError):
            C.source_energy(-1, 3, d)
        with self.assertRaises(ValueError):
            C.paired_errors(-1, d, 0)

    def test_each_bias_family_was_separately_admitted(self):
        families = self.report['separate_bias_witnesses']
        self.assertEqual(set(families), {'BIAS0', 'BIAS1', 'BIAS2'})
        self.assertEqual(len({x['qualification'] for x in families.values()}), 3)
        for value in families.values():
            self.assertTrue(value['own_family_admission'])
            self.assertEqual(value['recurrence_residual'], {})
            self.assertEqual(value['true_bias'], '0')
            self.assertEqual(value['driver'], '0')

    def test_scope_does_not_claim_an_internal_instability_or_P4(self):
        self.assertEqual(C.validate(self.report), [])
        self.assertEqual(self.report['classification'], 'B')
        self.assertTrue(self.report['bounded_all18_indefinite_target_refuted_under_finite_window_definition'])
        for key in ('nominal_filter_instability_established',
                    'conditional_ISS_with_unbounded_S_input_refuted',
                    'full_Normal_Live_source_admission_claimed',
                    'full_source_transition_cover_materialized',
                    'finite_startup_runtime_handoff_refuted', 'P4_PASS',
                    'P5_MAY_START', 'deployment_arithmetic_qualification_claimed'):
            self.assertFalse(self.report[key], key)
        self.assertEqual(self.report['P3_delta'], 1e-18)

    def test_mutated_proof_objects_and_promotion_are_rejected(self):
        for key, value in (('P4_PASS', True), ('P3_delta', 1e-17),
                           ('parallelogram_identity_residual', {'1': '1'})):
            bad = deepcopy(self.report)
            bad[key] = value
            self.assertTrue(C.validate(bad), key)
        bad = deepcopy(self.report)
        bad['literal_event_certificates'][0]['coefficient_map_B_h'][12][0] = {}
        self.assertTrue(C.validate(bad))


if __name__ == '__main__':
    unittest.main()
