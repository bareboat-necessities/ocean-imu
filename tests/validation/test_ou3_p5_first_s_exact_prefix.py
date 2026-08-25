import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_first_s_exact_prefix as P


class Ou3P5FirstSExactPrefixTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = P.build()

    def test_first_S_exact_prefix_passes_on_widened_finite_chart(self):
        d = self.d
        self.assertEqual(P.validate(d), [])
        self.assertEqual(d["P5_FIRST_DUE_S_EXACT_CAYLEY_PREFIX_CERTIFICATE"], "PASS_WIDENED_CHART")
        self.assertTrue(d["current_first_S_prefix_inside_widened_chart"])
        self.assertEqual(d["first_failure"], "NONE_AT_FIRST_S_ATTITUDE_CHART")
        self.assertFalse(d["diagnostic_q_lt_1_is_promotion_gate"])

    def test_exact_source_quaternion_composition_keeps_antipodal_denominator_positive(self):
        for row in self.d["nodes"].values():
            self.assertTrue(row["chart_safe"])
            self.assertGreater(row["cayley_composition_denominator_lower"], 0.0)
            self.assertGreater(row["post_injection_cayley_norm_upper"], 1.0)
            self.assertFalse(row["inside_diagnostic_cayley_lt_1"])
            self.assertTrue(row["inside_widened_prefix_chart"])

    def test_widening_is_exact_chart_not_filter_change_or_old_local_gate(self):
        d = self.d
        self.assertGreater(d["widened_prefix_cayley_norm_upper"], d["required_first_S_post_cayley_norm_upper"])
        self.assertLessEqual(d["widened_prefix_cayley_norm_upper"], 16.0)
        self.assertGreater(d["widened_prefix_antipodal_one_plus_cosine_margin_lower"], 0.0)
        self.assertGreater(d["widened_prefix_exact_vector_residual_factor_lower"], 0.0)
        self.assertGreater(d["widened_prefix_pair_information_vs_goLive_attitude_metric_lower"], 0.0)
        self.assertFalse(d["filter_changed"])
        self.assertTrue(d["full_S_to_attitude_gain_retained"])

    def test_q_lt_1_tightening_target_is_retained_as_diagnostic_only(self):
        n = self.d["nodes"]["normal_gauged"]
        t = self.d["nodes"]["timeout_gauged"]
        self.assertGreater(
            n["certified_correction_radius_for_diagnostic_cayley_lt_1_rad"],
            t["certified_correction_radius_for_diagnostic_cayley_lt_1_rad"],
        )
        self.assertLess(t["certified_correction_radius_for_diagnostic_cayley_lt_1_rad"], 0.30)
        self.assertGreater(t["certified_correction_radius_for_diagnostic_cayley_lt_1_rad"], 0.20)
        self.assertGreater(t["diagnostic_Ptheta_tightening_factor"], 1.0)
        self.assertLess(
            t["directional_Ptheta_upper_target_for_diagnostic_cayley_lt_1_if_other_factors_unchanged"],
            t["current_directional_Ptheta_upper"],
        )

    def test_no_replay_or_S_cross_gain_removal(self):
        self.assertFalse(self.d["source_replay_used"])
        self.assertTrue(self.d["full_S_to_attitude_gain_retained"])


if __name__ == "__main__":
    unittest.main()
