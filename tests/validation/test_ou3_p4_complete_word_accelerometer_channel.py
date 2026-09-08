from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "tools" / "stability"))

from ou3_interval import matrix_point, symmetric_positive_definite_ldlt  # noqa: E402
import ou3_p4_complete_word_accelerometer_channel as mod  # noqa: E402


class CompleteWordAccelerometerChannelTest(unittest.TestCase):
    def test_contract_is_joint_complete_word_and_nonpromoting(self):
        d = mod.build()
        self.assertEqual(mod.validate(d), [])
        self.assertEqual(d["canonical_source"], "COMPLETE_BRMM_NORMAL_LIVE_WORD")
        self.assertTrue(d["accelerometer_measurement_noise_is_PSD_final_covariance_component"])
        self.assertTrue(d["all_later_due_S_updates_remain_inside_suffix"])
        self.assertTrue(d["actual_applied_RS_required_for_every_later_S_suffix_event"])
        self.assertTrue(d["actual_RS_regularization_not_removed_by_channel_reduction"])
        self.assertFalse(d["packetwise_norm_sum_used"])
        self.assertFalse(d["packet_count_multiplier_used"])
        self.assertFalse(d["inverse_metric_floor_claim_used"])
        self.assertFalse(d["correction_radius_claim_used"])
        self.assertFalse(d["stacked_nonlinear_residual_graph_closed_here"])
        self.assertFalse(d["P4_promoted_here"])

    def test_covariance_component_and_precision_channel_are_full_matrix(self):
        P = matrix_point([[1.0, 0.0], [0.0, 1.0]])
        A = [
            matrix_point([[0.2], [0.0]]),
            matrix_point([[0.0], [0.3]]),
        ]
        R = [matrix_point([[1.0]]), matrix_point([[1.0]])]
        residual = mod.covariance_component_margin(P, A, R)
        ok_residual, pivots = symmetric_positive_definite_ldlt(residual)
        self.assertTrue(ok_residual)
        self.assertGreater(min(x.lo for x in pivots), 0.0)
        ok_precision, precision_pivots = mod.certify_strict_fixture_precision_domination(P, A, R)
        self.assertTrue(ok_precision)
        self.assertGreater(min(precision_pivots), 0.0)

    def test_channel_rejects_dimension_or_block_count_mismatch(self):
        P = matrix_point([[1.0]])
        A = [matrix_point([[0.2]])]
        R = [matrix_point([[1.0]])]
        self.assertEqual(len(mod.covariance_component_margin(P, A, R)), 1)
        with self.assertRaises(ValueError):
            mod.channel_covariance(A, [])
        with self.assertRaises(ValueError):
            mod.channel_covariance([matrix_point([[0.2, 0.1]])], R)
        with self.assertRaises(ValueError):
            mod.covariance_component_margin(matrix_point([[1.0, 0.0], [0.0, 1.0]]), A, R)


if __name__ == "__main__":
    unittest.main()
