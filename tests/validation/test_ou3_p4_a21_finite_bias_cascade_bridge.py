from __future__ import annotations

import math
import unittest

import ou3_p4_a21_finite_bias_cascade_bridge as BRIDGE


class TestA21FiniteTauBCascadeBridge(unittest.TestCase):
    def test_exact_schur_budget_is_strict(self):
        dh = 1.0e-18
        beta = 1.0e-3
        da = 5.0e-19
        lim = BRIDGE.cascade_g2_limit_conservative(dh, beta, da)
        self.assertGreater(lim, 0.0)
        g2 = 0.25 * lim
        self.assertGreater(
            BRIDGE.cascade_determinant_slack_lower(dh, beta, da, g2), 0.0
        )

    def test_target_must_be_below_both_diagonal_gaps(self):
        with self.assertRaises(ValueError):
            BRIDGE.cascade_g2_limit_conservative(1.0e-18, 1.0e-3, 1.0e-18)
        with self.assertRaises(ValueError):
            BRIDGE.cascade_g2_limit_conservative(1.0e-3, 1.0e-4, 1.0e-4)

    def test_mu_logarithm_does_not_overflow(self):
        x = BRIDGE.required_mu_log10(1.0e300, 1.0e-22)
        self.assertTrue(math.isfinite(x))
        self.assertGreater(x, 300.0)
        self.assertEqual(BRIDGE.required_mu_log10(0.0, 1.0e-22), -math.inf)

    def test_repository_build_is_fail_closed_at_joint_coupling(self):
        d = BRIDGE.build()
        self.assertEqual(d["canonical_source"], "COMPLETE_SEA3_NORMAL_LIVE_WORD")
        self.assertEqual(d["P3_frozen_gate"], 1.0e-18)
        self.assertTrue(d["P3_frozen_not_modified"])
        self.assertTrue(d["actual_applied_R_S_retained"])
        self.assertTrue(d["full_H18_state_retained_including_aw"])
        self.assertTrue(d["full_ba_state_retained"])
        self.assertFalse(d["packet_count_coupling_bound_consumed"])
        self.assertFalse(d["joint_same_history_C_Hb_metric_bound_closed"])
        self.assertIsNone(d["joint_C_Hb_c2_upper"])
        self.assertFalse(d["full_rank_A21_cascade_metric_closed"])
        self.assertFalse(d["P4_promoted_here"])
        self.assertFalse(d["P5_may_start"])
        self.assertEqual(BRIDGE.validate(d), [])

    def test_metric_lower_is_full_matrix_not_marginal_inverse_floor(self):
        d = BRIDGE.build()
        m = d["H18_metric_equivalence"]
        self.assertGreater(m["posterior_covariance_lambda_min_lower"], 0.0)
        self.assertGreater(m["M_H_lambda_max_upper"], 0.0)
        self.assertEqual(m["reset_sigma_min_lower"], 1.0)
        self.assertFalse(m["reset_lower_bound_requires_correction_radius"])
        self.assertFalse(m["marginal_inverse_metric_floor_inference_used"])


if __name__ == "__main__":
    unittest.main()
