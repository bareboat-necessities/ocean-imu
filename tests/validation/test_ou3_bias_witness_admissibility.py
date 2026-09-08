from pathlib import Path
import json
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tests/kalman_ou_iii"))

import ou3_p4_bias_witness_admissibility as AUDIT  # noqa: E402


class BiasWitnessAdmissibilityTests(unittest.TestCase):
    def test_perfect_cancellation_has_no_positive_BIAS2_separation(self):
        d = AUDIT.separation_ratios(2.0, 2.0, -2.0, 0.0, 1.0)
        self.assertEqual(d["kappa_point"], 1.0)
        self.assertEqual(d["mu_sufficient_point"], 0.0)
        self.assertFalse(d["positive_sufficient_separation_at_this_point"])

    def test_dense_cross_term_is_retained_with_the_correct_sign(self):
        d = AUDIT.separation_ratios(4.0, 1.0, -1.5, 2.0, 2.0)
        self.assertEqual(d["energy_identity_residual"], 0.0)
        self.assertAlmostEqual(d["kappa_point"], 0.6)
        self.assertEqual(d["mu_actual_point"], 1.0)
        self.assertEqual(d["mu_sufficient_point"], 1.0)
        self.assertFalse(d["source_uniform_mu_certified"])

    def test_positive_cross_term_does_not_become_a_false_lower_bound(self):
        d = AUDIT.separation_ratios(4.0, 1.0, 1.5, 8.0, 2.0)
        self.assertEqual(d["mu_actual_point"], 4.0)
        self.assertEqual(d["mu_sufficient_point"], 1.0)
        with self.assertRaises(ValueError):
            AUDIT.separation_ratios(1.0, 1.0, 0.0, 2.0, 0.0)

    def test_bias_compatibility_never_self_certifies_source_membership(self):
        self.assertEqual(AUDIT.case_verdict(True, True, True),
                         "BIAS_COMPATIBLE_FULL_BRMM_ADMISSIBILITY_UNRESOLVED")
        self.assertEqual(AUDIT.case_verdict(False, True, True),
                         "OUTSIDE_CHECKED_CONDITIONAL_POINT_DOMAIN")
        self.assertEqual(AUDIT.case_verdict(True, True, False),
                         "BIAS1_POINT_RECURRENCE_FAILED")

    def _one_event(self, initial_bias):
        # A correction moves an initially out-of-domain estimate back inside.
        # Checking only event outputs would incorrectly admit the initial state.
        n = 21
        x = np.zeros(n)
        x[18] = initial_bias
        H = np.zeros((3, n))
        H[:, 15:18] = np.eye(3)
        H[:, 18:21] = np.eye(3)
        K = np.zeros((n, 3))
        K[18:21] = 0.5 * np.eye(3)
        event = {"type": AUDIT.BASE.EV_ACC, "name": "accelerometer", "time": 0.0,
                 "H": H, "R": np.eye(3)}
        linear = {"direction": x, "Q0": np.eye(n), "QN": np.eye(n),
                  "path": [{"Pafter": np.eye(n)}]}
        domain = json.loads(AUDIT.BIAS.DEFAULT_DOMAIN.read_text())
        return AUDIT.audit_case(
            {"mode_dim": n, "t0": 0.0, "events": [event]}, linear,
            [{"K": K, "Qafter": np.eye(n)}], domain, 1.0,
            Path("wave_data_pmstokes_H1.500_L50.710_A-30.00_P120.00.csv"),
        )

    def test_initial_estimate_interior_is_checked_even_if_correction_repairs_it(self):
        d = self._one_event(0.38)
        self.assertFalse(d["estimated_bias_interior_retained"])
        self.assertEqual(d["first_estimated_bias_interior_failure"]["event_index"], -1)
        self.assertTrue(d["legacy_state_domain_retained_including_initial"])

    def test_measurement_changes_error_without_becoming_external_bias_forcing(self):
        d = self._one_event(0.2)
        self.assertTrue(d["bias_recurrence_point_check_pass"])
        self.assertTrue(d["conditional_point_checks_pass"])
        self.assertAlmostEqual(d["max_measurement_bias_correction"], 0.1)
        self.assertAlmostEqual(d["max_corrected_error_minus_free_GM_path"], 0.1)
        self.assertTrue(d["external_bias_forcing_zero"])
        self.assertFalse(d["canonical_P4_falsified_here"])


if __name__ == "__main__":
    unittest.main()
