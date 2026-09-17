import unittest

from tools.stability.ou3_alt_contraction.informative_service_theorem import InformativeServiceAssumption as C


class InformativeServiceTests(unittest.TestCase):
    def setUp(self):
        self.contract = C(1.0, 3.0, .01, 1.0)
        self.row = {'start_sample': 0, 'end_sample': 600,
                    'accepted_magnetic_events': 3,
                    'max_observed_magnetic_gap_s': 1.0,
                    'minimum_transported_heading_response': .02,
                    'transported_heading_bias_min_eigenvalue': 2.0}

    def test_membership_does_not_prove_stability_or_recurrence(self):
        result = self.contract.audit_native_window(self.row, .005)
        self.assertTrue(result['finite_window_membership_point_pass'])
        self.assertFalse(result['infinite_recurrence_certified'])
        self.assertFalse(result['stability_implied_by_membership'])
        declared = self.contract.declaration()
        self.assertTrue(declared['assumption_explicitly_selected'])
        self.assertFalse(declared['unconditional_source_class_restricted'])
        self.assertFalse(declared['storage_search_allowed'])

    def test_rejected_degenerate_sparse_and_wrong_duration_fail(self):
        for field, value in (('accepted_magnetic_events', 0),
                             ('max_observed_magnetic_gap_s', 1.01),
                             ('minimum_transported_heading_response', 0),
                             ('transported_heading_bias_min_eigenvalue', 0),
                             ('end_sample', 1200)):
            with self.subTest(field=field):
                row = dict(self.row, **{field: value})
                self.assertFalse(self.contract.audit_native_window(row, .005)
                                 ['finite_window_membership_point_pass'])

    def test_invalid_parameter_and_measurement_rejected(self):
        for value in (0, -1, float('nan'), float('inf')):
            with self.assertRaises(ValueError):
                C(value, 3, .01, 1)
        with self.assertRaises(ValueError):
            C(2, 3, .01, 1)
        with self.assertRaises(ValueError):
            self.contract.audit_native_window(dict(self.row,
                transported_heading_bias_min_eigenvalue=float('nan')), .005)


class QuotientObstructionTests(unittest.TestCase):
    def test_quotient_storage_minimizes_over_gauge_fibre(self):
        import numpy as np
        from tools.stability.ou3_alt_contraction.service_regime_diagnostic import quotient_metric
        M = np.array([[2., 1.], [1., 2.]])
        reduced = quotient_metric(M, np.eye(2)[:, :1], np.eye(2)[:, 1:])
        self.assertAlmostEqual(reduced[0, 0], 1.5)
        self.assertNotEqual(reduced[0, 0], M[1, 1])


    def test_yaw_removal_retains_neutral_bias_and_transverse_split_checks_leakage(self):
        import numpy as np
        from tools.stability.ou3_alt_contraction.service_regime_diagnostic import quotient_diagnostic
        A = .5*np.eye(24)
        A[2, 2] = A[5, 5] = 1
        A[2, 5] = 3
        result = quotient_diagnostic(A)
        self.assertEqual(result['heading_only']['motion_spectral_radius'], 1)
        self.assertEqual(result['heading_and_axial_bias']['motion_spectral_radius'], .5)
        # A moving source can couple centre bias back into translation.
        A[6, 5] = .01
        result = quotient_diagnostic(A)
        self.assertFalse(result['heading_and_axial_bias']['point_quotient_exists'])
        self.assertFalse(result['heading_and_axial_bias']['axial_bias_is_heading_gauge'])


class RecordedServiceEvidenceTests(unittest.TestCase):
    def test_recorded_native_family_and_scope_guards(self):
        import json
        from pathlib import Path
        from tools.stability.ou3_alt_contraction.service_regime_diagnostic import validate
        path = Path(__file__).resolve().parents[2]/'reports/results/ou3_alt_storage/service-regimes.json'
        report = json.loads(path.read_text())
        self.assertEqual(validate(report), [])
        report['histories']['H:stride200']['retained_MAG_CALL_SCHEDULE_v1_satisfied'] = True
        self.assertIn('H:stride200 callback profile misclassified', validate(report))
        report['source_uniform_rho_certified'] = True
        self.assertIn('source_uniform_rho_certified incorrectly promoted', validate(report))


if __name__ == '__main__':
    unittest.main()
