from __future__ import annotations
import math
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import ou3_p5_sample1_structured_full_gain_v12c as V12C


class StructuredFullGainV12CTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Fast fixture.  The authoritative focused CI uses the full 24^3 V10
        # subdivision.  At 2^3 V10 may intentionally remain fail closed.
        cls.d = V12C.build(source_pieces=4, source_cell_index=0,
                           p_pieces=2, tangent_pieces=2, axial_pieces=2)

    def test_one_plus_two_structure_is_active(self):
        for k in (
            "nominal_one_plus_two_innovation_structure_used",
            "scalar_positive_innovation_identity_used",
            "two_by_two_positive_determinant_identity_used",
            "two_by_two_lambda_floor_det_over_trace_used",
            "attitude_gain_rows_only_bounded_for_q_gate",
            "actual_shipping_covariance_PSD_used",
            "actual_innovation_noise_floor_fallback_retained",
        ):
            self.assertIs(self.d[k], True)
        self.assertGreater(self.d["minimum_nominal_block_lambda_lower"], 0.0)

    def test_no_generic_ldlt_or_broad_inverse_regression(self):
        self.assertIs(self.d["generic_interval_innovation_LDLT_floor_used"], False)
        self.assertIs(self.d["broad_sample1_3x3_interval_inverse_reintroduced"], False)
        self.assertTrue(math.isfinite(self.d["max_actual_innovation_inverse_operator_upper"]))

    def test_result_is_fail_closed_or_closed(self):
        self.assertIn(self.d["P5_SAMPLE1_ONE_PLUS_TWO_ATTITUDE_RESOLVENT_V12C"],
                      ("PASS", "NOT_ESTABLISHED"))
        vf = V12C.validate(self.d)
        v10_missing = any("V10 prerequisite did not pass" in x for x in vf)
        if v10_missing:
            self.assertEqual(self.d["P5_SAMPLE1_ONE_PLUS_TWO_ATTITUDE_RESOLVENT_V12C"],
                             "NOT_ESTABLISHED")
            self.assertTrue(any("V10 prerequisite did not pass" in x
                                for x in self.d["failures"]))
        else:
            self.assertEqual(vf, [])
            if self.d["P5_SAMPLE1_ONE_PLUS_TWO_ATTITUDE_RESOLVENT_V12C"] == "PASS":
                self.assertEqual(self.d["unclosed_joint_cells"], 0)
                self.assertIsNone(self.d["first_unclosed_joint_cell"])
            else:
                self.assertIsNotNone(self.d["first_unclosed_joint_cell"])

    def test_no_promotion_or_shipping_limit_change(self):
        self.assertEqual(float(self.d["deployed_correction_limit_rad"]), 6.0)
        for k in (
            "source_replay_used", "filter_changed",
            "deployed_correction_limit_increased",
            "complete_sample1_branch_closed_here",
            "signed_cayley_q8_composed_here", "q8_word_promoted_here",
            "whole_word_promoted_here", "N_H_words_set_here",
        ):
            self.assertIs(self.d[k], False)


if __name__ == "__main__":
    unittest.main()
