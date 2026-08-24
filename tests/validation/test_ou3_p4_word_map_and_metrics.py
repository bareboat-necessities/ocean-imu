from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p4_exact_word_map as WORD
import ou3_p4_node_metrics as METRIC


class Ou3P4WordMapAndMetricTests(unittest.TestCase):
    def test_exact_word_map_uses_shipping_order_and_immediate_resets(self):
        d = WORD.build()
        self.assertEqual(WORD.validate(d), [])
        self.assertEqual(d["shipping_operation_order"], WORD.EXPECTED_ORDER)
        for mode in ("H", "A"):
            p = d[mode]["correction_policy"]
            self.assertTrue(p["deployed_normalized_quaternion_map"])
            self.assertTrue(p["left_error_reset_after_each_accepted_correction"])
            self.assertTrue(p["full_S_to_attitude_cross_gain"])
            self.assertFalse(p["linearized_attitude_injection_allowed"])
            self.assertFalse(p["one_shared_reset_after_all_measurements"])

    def test_S_pseudo_precedes_accelerometer_in_same_imu_update(self):
        d = WORD.build()
        order = d["shipping_operation_order"]
        self.assertLess(
            order.index("periodic_S_zero_when_due_then_immediate_quaternion_injection_and_left_error_reset"),
            order.index("accelerometer_correction_or_rejection_then_immediate_quaternion_injection_and_left_error_reset_if_accepted"),
        )

    def test_A_mode_retains_exact_bias_projection_semantic(self):
        d = WORD.build()
        self.assertIsNone(d["H"]["active_bias_projection"])
        p = d["A"]["active_bias_projection"]
        self.assertEqual(p["kind"], "exact_Euclidean_projection_onto_closed_ball")
        self.assertTrue(p["nonexpansive"])

    def test_metric_family_is_group_compatible_and_ou_scaled(self):
        d = METRIC.build()
        self.assertEqual(METRIC.validate(d), [])
        self.assertTrue(d["selection_requires_validated_endpoint_word_bound"])
        self.assertIsNone(d["selected_candidate"])
        for candidate in d["candidates"].values():
            self.assertTrue(candidate["same_multiplier_on_every_graph_node"])
            for mode, dim in (("H",18),("A",21)):
                for node in candidate["nodes"][mode]:
                    m = node["metric"]
                    self.assertEqual(m["kind"], "GROUP_COMPATIBLE_NODE_METRIC")
                    self.assertFalse(m["attitude_linear_cross_terms"])
                    self.assertFalse(m["equals_Kalman_inverse_covariance"])
                    self.assertEqual(len(m["Pbar_diagonal"]), dim)
                    labels = [b["label"] for b in m["P_xi_blocks"]]
                    self.assertEqual(labels[:5], ["b_g","v","p","S","a_w"])
                    if mode == "A":
                        self.assertEqual(labels[-1], "b_a")

    def test_no_old_metric_fallback_names_are_present(self):
        text = (ROOT / "tools" / "ou3_p4_node_metrics.py").read_text(encoding="utf-8").casefold()
        self.assertNotIn("common euclidean fallback", text)
        self.assertNotIn("schmidt", text)


if __name__ == "__main__":
    unittest.main()
