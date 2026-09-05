from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_sea3_spread_s_regularizer as REG  # noqa: E402


class Sea3SpreadSRegularizerTest(unittest.TestCase):
    def test_complete_sea3_actual_rs_spread_regularizer_is_strict(self):
        d = REG.build()
        self.assertEqual(REG.validate(d), [])
        self.assertEqual(d["canonical_source"], "COMPLETE_SEA3_NORMAL_LIVE_WORD")
        self.assertTrue(d["complete_SEA3_source_consumed"])
        self.assertFalse(d["source_family_replaced"])
        self.assertFalse(d["trajectory_replay_used"])
        self.assertFalse(d["independent_tau_sigma_RS_TS_extrema_product_used"])
        self.assertFalse(d["selected_S_events_replace_full_scheduler_word"])
        self.assertTrue(d["all_due_S_updates_remain_in_literal_word"])
        self.assertTrue(d["actual_applied_SpectralMSE_R_S_consumed"])
        self.assertTrue(d["selected_firings_are_guaranteed_members_of_full_word"])
        self.assertLessEqual(d["scheduler_uniform_gap_s_upper"], 0.151)
        self.assertEqual(len(d["selected_backward_lag_windows_s"]), 3)
        self.assertGreater(d["endpoint_integrator_vandermonde_det_abs_lower"], 0.0)
        self.assertGreater(d["endpoint_integrator_sigma_min_lower"], 0.0)
        info = d["integrator_information_lambda_min_lower"]
        self.assertTrue(math.isfinite(info))
        self.assertGreater(info, 0.0)
        self.assertTrue(d["spread_S_regularizer_pass"])
        self.assertFalse(d["P3_promoted"])

    def test_rs_is_applied_std_not_target_or_variance(self):
        d = REG.build()
        eff = d["effective_observation_covariance"]
        rs_std = eff["actual_applied_R_S_axis_std_upper"]
        rs_var = eff["actual_applied_R_S_variance_upper"]
        self.assertGreater(rs_std, 0.0)
        self.assertGreaterEqual(rs_var, rs_std * rs_std)
        self.assertEqual(d["R_S_axis_std_factors"], [0.72, 0.72, 1.0])


if __name__ == "__main__":
    unittest.main()
