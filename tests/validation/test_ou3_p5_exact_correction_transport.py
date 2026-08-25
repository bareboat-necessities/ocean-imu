import math
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_exact_correction_transport as T


class Ou3P5ExactCorrectionTransportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = T.build()

    def test_source_bound_exact_transport_algebra_passes(self):
        d = self.d
        self.assertEqual(T.validate(d), [])
        self.assertEqual(d["P5_EXACT_CORRECTION_TRANSPORT_ALGEBRA_CERTIFICATE"], "PASS")
        self.assertTrue(d["source_generated_not_trajectory_fit"])
        self.assertFalse(d["source_replay_used"])
        self.assertFalse(d["filter_changed"])
        self.assertTrue(d["full_state_cross_terms_retained"])

    def test_reset_congruence_has_no_condition_number_penalty(self):
        d = self.d
        self.assertFalse(d["condition_number_multiplier_used_for_reset_transport"])
        self.assertIn("Ge^-1 rho", d["exact_reset_congruence_identity"])
        for row in d["nodes"].values():
            self.assertEqual(row["reset_inverse_operator_norm_upper"], 1.0)
            self.assertEqual(row["reset_min_singular_value_lower"], 1.0)
            self.assertGreaterEqual(row["reset_determinant_lower"], 1.0)

    def test_current_first_S_correction_keeps_both_gauged_nodes_off_antipode(self):
        d = self.d
        self.assertGreater(d["first_due_S_correction_norm_upper_rad"], 1.0)
        normal = d["nodes"]["normal_gauged"]
        timeout = d["nodes"]["timeout_gauged"]
        self.assertTrue(normal["chart_safe"])
        self.assertTrue(timeout["chart_safe"])
        self.assertGreater(normal["cayley_composition_denominator_lower"], 0.0)
        self.assertGreater(timeout["cayley_composition_denominator_lower"], 0.0)
        self.assertGreater(normal["cayley_composition_denominator_lower"], timeout["cayley_composition_denominator_lower"])
        self.assertLess(timeout["injected_cayley_norm_upper"], 3.0)

    def test_deployed_axis_angle_cayley_bound_is_finite_and_not_linearized(self):
        row = T.correction_cayley_norm_bounds(1.0)
        self.assertIn("axis_angle", row["source_branch_family"])
        self.assertGreater(row["injected_cayley_norm_upper"], 1.0)
        self.assertGreater(row["injected_cayley_minus_delta_norm_upper"], 0.0)
        self.assertTrue(math.isfinite(row["injected_cayley_norm_upper"]))

    def test_small_source_series_branch_is_enclosed_separately(self):
        row = T.correction_cayley_norm_bounds(1.0e-3)
        self.assertEqual(row["source_branch_family"], "source_polynomial_series")
        self.assertLess(row["injected_cayley_minus_delta_norm_upper"], 1.0e-9)
        self.assertGreater(row["series_cayley_coefficient_upper"], 0.0)

    def test_algebra_does_not_promote_complete_word(self):
        self.assertFalse(self.d["complete_word_numerical_budget_closed_here"])
        self.assertIn("complete source-correlated 1 s word", self.d["next_obligation"])


if __name__ == "__main__":
    unittest.main()
