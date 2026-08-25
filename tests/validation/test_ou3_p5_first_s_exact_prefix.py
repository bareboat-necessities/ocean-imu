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

    def test_fail_closed_first_S_prefix_obstruction_is_machine_named(self):
        d = self.d
        self.assertEqual(P.validate(d), [])
        self.assertEqual(d["P5_FIRST_DUE_S_EXACT_CAYLEY_PREFIX_CERTIFICATE"], "NOT_ESTABLISHED")
        self.assertFalse(d["current_first_S_prefix_inside_outer_bootstrap"])
        self.assertEqual(
            d["first_failure"],
            "FIRST_DUE_S_EXACT_CAYLEY_PREFIX_NOT_CERTIFIED_WITH_CURRENT_STAGED_BOUND",
        )

    def test_exact_source_quaternion_composition_keeps_antipodal_denominator_positive(self):
        for row in self.d["nodes"].values():
            self.assertTrue(row["chart_safe"])
            self.assertGreater(row["cayley_composition_denominator_lower"], 0.0)
            self.assertGreater(row["post_injection_cayley_norm_upper"], 1.0)
            self.assertFalse(row["inside_common_outer_bootstrap"])

    def test_timeout_node_sets_the_limiting_tightening_target(self):
        n = self.d["nodes"]["normal_gauged"]
        t = self.d["nodes"]["timeout_gauged"]
        self.assertGreater(
            n["certified_correction_radius_for_cayley_lt_1_rad"],
            t["certified_correction_radius_for_cayley_lt_1_rad"],
        )
        self.assertLess(t["certified_correction_radius_for_cayley_lt_1_rad"], 0.30)
        self.assertGreater(t["certified_correction_radius_for_cayley_lt_1_rad"], 0.20)
        self.assertGreater(t["Ptheta_tightening_factor_required"], 1.0)
        self.assertLess(
            t["directional_Ptheta_upper_target_if_other_first_S_factors_unchanged"],
            t["current_directional_Ptheta_upper"],
        )

    def test_no_filter_or_S_cross_gain_change_is_used(self):
        self.assertFalse(self.d["filter_changed"])
        self.assertTrue(self.d["full_S_to_attitude_gain_retained"])
        self.assertFalse(self.d["source_replay_used"])


if __name__ == "__main__":
    unittest.main()
