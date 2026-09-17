import os
from pathlib import Path
import tempfile
import unittest

import numpy as np

from tools.stability.ou3_alt_contraction import carried_storage_rho_diagnostic as D
from tools.stability.ou3_alt_contraction import magnetic_service_formulation as F


class CarriedStorageTests(unittest.TestCase):
    def test_terminal_high_precision_keeps_active_bias_output_energy(self):
        A = np.eye(24)
        A[:18, :18] *= .5
        A[18, 0] = 2
        P = np.eye(21)
        # Dropping active-bias output would incorrectly report 0.25.
        rho = D.high_precision_ratio(A, P, P)
        self.assertAlmostEqual(float(rho), 4.25)

    def test_storage_keeps_held_covariance_cross_terms(self):
        P = np.eye(21)
        P[0, 18] = P[18, 0] = .5
        M = D.covariance_metric(P)
        self.assertAlmostEqual(M[0, 0], 4 / 3)
        self.assertAlmostEqual(M[0, 18], -2 / 3)
        self.assertFalse(np.array_equal(M[:18, :18], np.linalg.inv(P[:18, :18])))

    def test_indefinite_covariance_has_no_floor_fallback(self):
        P = np.eye(21)
        P[0, 0] = -1e-9
        with self.assertRaises(np.linalg.LinAlgError):
            D.covariance_metric(P)

    def test_bias_prediction_has_physical_bias_column(self):
        A = np.eye(21)
        A[18:21, 18:21] *= .75
        J = D.lift(A, 2)
        np.testing.assert_array_equal(J[18:21, 21:24], .25 * np.eye(3))
        np.testing.assert_array_equal(J[21:24, 21:24], np.eye(3))

    def test_covariance_growth_can_fake_decay_without_coercivity(self):
        # An unchanged error can shrink in an unbounded-covariance metric.
        A = np.eye(24)
        M0 = D.covariance_metric(np.eye(21))
        M1 = D.covariance_metric(2 * np.eye(21))
        result = F.projected_storage_ratio(A, M0, M1, F.joint24_motion_injection())
        self.assertAlmostEqual(result['rho_point'], .5)
        self.assertFalse(result['source_uniform_rho_certified'])

    def test_detached_event_rejected(self):
        P = ' '.join(map(str, np.eye(21).ravel()))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'broken.txt'
            path.write_text(f'0 0 0 {P}\n1 1011 0 {P}\n')
            with self.assertRaisesRegex(ValueError, 'does not complete'):
                D.analyze(path, {'samples': 600, 'dt_s': .005})

    @unittest.skipUnless(os.environ.get('OU3_ALT_REQUIRE_NATIVE') == '1',
                         'explicit native feasibility requested')
    def test_native_carried_histories_and_nonpromotion(self):
        with tempfile.TemporaryDirectory() as directory:
            report = D.run(directory, samples=1200)
        self.assertTrue(report['passive_instrumentation_bit_identity'])
        self.assertEqual(D.validate(report), [])
        self.assertEqual(set(report['histories']), {'H', 'A', 'HA'})
        for history in report['histories'].values():
            self.assertEqual(len(history['windows']), 2)
            for row in history['windows']:
                self.assertGreater(row['accepted_magnetic_events'], 0)
                self.assertGreater(row['transported_heading_bias_min_eigenvalue'], 0)
                self.assertLess(row['energy_accounting_residual'], 1e-8)
        self.assertEqual(report['histories']['HA']['windows'][1]['release_edges'], 1)
        for key in ('source_uniform_rho_certified', 'uniform_metric_coercivity_certified',
                    'storage_search_allowed', 'ALT_STARTUP_PASS', 'ALT_LIVE_PASS',
                    'ALT_END_TO_END_PASS', 'interval_enclosure_authorized'):
            self.assertFalse(report[key])
        report['source_uniform_rho_certified'] = True
        self.assertIn('source_uniform_rho_certified not false', D.validate(report))


if __name__ == '__main__':
    unittest.main()
