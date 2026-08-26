import math
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_first_accel_post_reset as G


class Ou3P5FirstAccelPostResetTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = G.build(source_pieces=2)

    def test_post_reset_stage_validates(self):
        d = self.d
        self.assertEqual(G.validate(d), [])
        self.assertEqual(d["P5_FIRST_ACCEL_POST_RESET_PREFIX_CERTIFICATE"], "PASS")
        self.assertTrue(d["source_generated_not_trajectory_fit"])
        self.assertFalse(d["source_replay_used"])
        self.assertFalse(d["filter_changed"])

    def test_sign_complete_injection_stays_in_declared_q8_chart(self):
        d = self.d
        self.assertTrue(d["SO3_triangle_inequality_used_for_sign_complete_composition"])
        self.assertTrue(d["rejected_accel_correction_is_identity"])
        self.assertTrue(d["post_accel_complete_branch_family_inside_q8"])
        self.assertGreater(d["post_accel_cayley_norm_upper"], d["pre_accel_cayley_norm_upper"])
        self.assertLess(d["post_accel_cayley_norm_upper"], d["q_chart_target"])
        self.assertLess(d["post_accel_geodesic_angle_upper_rad"], math.pi)

    def test_reset_bound_matches_deployed_left_error_reset_identity(self):
        d = self.d
        dmax = float(d["accepted_accel_correction_norm_upper_rad"])
        expected_op = math.sqrt(1.0 + 0.25 * dmax * dmax)
        self.assertTrue(d["reset_is_nonsingular_for_all_first_accel_children"])
        self.assertGreaterEqual(d["reset_operator_norm_upper"], expected_op)
        self.assertGreaterEqual(
            d["reset_covariance_quadratic_multiplier_upper"],
            d["reset_operator_norm_upper"] ** 2,
        )

    def test_no_theorem_or_deployment_gate_is_relaxed(self):
        d = self.d
        self.assertEqual(d["deployed_correction_limit_rad"], 6.0)
        self.assertFalse(d["deployed_correction_limit_increased"])
        self.assertFalse(d["full_matrix_Joseph_reset_children_propagated_here"])
        self.assertFalse(d["whole_word_promoted_here"])
        self.assertFalse(d["N_H_words_set_here"])
        self.assertEqual(
            d["next_obligation"],
            "PROPAGATE_EXACT_FIRST_ACCEL_JOSEPH_RESET_COVARIANCE_AND_STATE_CHILDREN_TO_SAMPLE_1",
        )


if __name__ == "__main__":
    unittest.main()
