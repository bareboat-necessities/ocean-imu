import math
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_45deg_first_accel_q8_bridge as B


class Ou3P545DegFirstAccelQ8BridgeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = B.build(source_pieces=2)

    def test_bridge_validates_and_is_source_bound(self):
        d = self.d
        self.assertEqual(B.validate(d), [])
        self.assertTrue(d["source_generated_not_trajectory_fit"])
        self.assertFalse(d["source_replay_used"])
        self.assertFalse(d["filter_changed"])
        self.assertEqual(d["P5_45DEG_FIRST_ACCEL_Q8_BRIDGE_CERTIFICATE"], "PASS")

    def test_45deg_first_packet_is_inside_q8_without_favorable_sign(self):
        d = self.d
        p = d["P5_45deg_entrance_first_accel"]
        self.assertFalse(d["accepted_correction_sign_assumed_favorable"])
        self.assertTrue(d["all_correction_directions_covered_by_dot_product_extremum"])
        self.assertTrue(p["accepted_or_rejected_branch_family_inside_q8"])
        self.assertGreater(p["product"]["homogeneous_product_scalar_lower"], 0.0)
        self.assertGreater(
            p["product"]["q8_test_lhs_lower_17W2"],
            p["product"]["q8_test_rhs_upper_4plusq2"],
        )
        self.assertLess(p["product"]["post_update_q_upper_from_scalar"], 8.0)
        self.assertLess(p["max_first_accelerometer_correction_norm_upper_rad"], 6.0)

    def test_candidate_ladder_is_also_inside_q8(self):
        rows = self.d["P4_candidate_first_accel_rows"]
        self.assertEqual([r["angle_deg"] for r in rows], [30.0, 25.0, 20.0, 15.0])
        self.assertTrue(self.d["all_P4_candidate_first_accel_branches_inside_q8"])
        for r in rows:
            self.assertTrue(r["product"]["inside_q8_for_every_correction_direction"])
            self.assertLess(r["product"]["post_update_q_upper_from_scalar"], 8.0)

    def test_reset_covariance_bound_keeps_joseph_and_identity_branches(self):
        p = self.d["P5_45deg_entrance_first_accel"]
        r = p["reset"]
        self.assertTrue(r["accepted_Joseph_posterior_Loewner_below_prior"])
        self.assertTrue(r["rejected_branch_is_identity"])
        self.assertGreaterEqual(r["accepted_or_rejected_covariance_multiplier_upper"], 1.0)
        expected = 1.0 + 0.25 * p["max_first_accelerometer_correction_norm_upper_rad"] ** 2
        self.assertGreaterEqual(r["reset_covariance_spectral_multiplier_upper"], expected)

    def test_no_transcendental_or_promotion_shortcut(self):
        d = self.d
        self.assertTrue(d["validated_sin_cos_used"])
        self.assertFalse(d["atan_tan_used_in_promoted_q8_test"])
        self.assertFalse(d["deployed_correction_limit_increased"])
        self.assertFalse(d["detailed_post_reset_cross_covariance_propagated_here"])
        self.assertFalse(d["returned_to_30deg_P4_sector_here"])
        self.assertFalse(d["P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE"])
        self.assertFalse(d["P5_FINITE_CAPTURE_TO_P4_ESTABLISHED_HERE"])


if __name__ == "__main__":
    unittest.main()
