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

    def test_metric_is_mode_normalized_exact_cayley_information_geometry(self):
        d = METRIC.build()
        self.assertEqual(METRIC.validate(d), [])
        self.assertTrue(d["single_quantitative_metric_route"])
        self.assertFalse(d["retired_block_diagonal_route_available"])
        self.assertTrue(d["global_scaling_does_not_change_physical_level_sets"])
        for mode, dim in (("H",18),("A",21)):
            m = d["modes"][mode]
            self.assertEqual(m["dimension"], dim)
            self.assertEqual(m["kind"], "CAYLEY_LIFTED_SOURCE_INFORMATION_METRIC")
            self.assertTrue(m["source_covariance_inverse"])
            self.assertGreater(m["mode_global_positive_scale"], 0.0)
            self.assertTrue(m["same_scale_on_every_source_node_in_mode"])
            self.assertTrue(m["full_attitude_linear_cross_terms_retained"])
            self.assertFalse(m["block_diagonal_metric_used"])
            self.assertFalse(m["common_Euclidean_metric_used"])
            self.assertTrue(m["local_coordinate_matches_P3_delta_theta"])
            self.assertTrue(m["local_quadratic_is_positive_scalar_multiple_of_P3_information_metric"])
            self.assertTrue(m["endpoint_metric_must_match_endpoint_source_covariance"])
            self.assertGreater(m["metric_lambda_min_lower"], 0.0)
            self.assertLessEqual(m["metric_lambda_min_lower"], 1.0)
            self.assertGreaterEqual(m["metric_lambda_max_upper"], m["metric_lambda_min_lower"])
            self.assertGreater(m["P3_word_endpoint_margin_lower"], 0.0)

    def test_no_retired_block_metric_or_schmidt_route_remains(self):
        text = (ROOT / "tools" / "ou3_p4_node_metrics.py").read_text(encoding="utf-8").casefold()
        self.assertNotIn("schmidt", text)
        self.assertIn("retired_block_diagonal_route_available", text)
        self.assertIn("full_attitude_linear_cross_terms_retained", text)
        self.assertIn("same_scale_on_every_source_node_in_mode", text)


if __name__ == "__main__":
    unittest.main()
