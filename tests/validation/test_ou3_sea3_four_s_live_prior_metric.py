from pathlib import Path
import math
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_sea3_four_s_live_prior_metric as mod  # noqa: E402


class Sea3FourSLivePriorMetricTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = mod.build()

    def test_complete_source_and_shipping_live_seed(self):
        self.assertEqual(mod.validate(self.d), [])
        self.assertEqual(self.d["canonical_source"], "COMPLETE_SEA3_NORMAL_LIVE_WORD")
        self.assertTrue(self.d["component_of_complete_SEA3_full_word"])
        self.assertTrue(self.d["shipping_Live_covariance_seed_consumed"])
        self.assertTrue(self.d["Live_aw_reset_to_same_committed_stationary_covariance_consumed"])
        self.assertFalse(self.d["selected_four_S_events_replace_complete_word"])
        self.assertFalse(self.d["source_family_replaced"])

    def test_translation_prior_metric_clears_exact_batch_requirement(self):
        info = self.d["translation_Live_prior_metric_information_lambda_min_lower"]
        req = self.d["batch_required_information_ratio"]
        self.assertTrue(math.isfinite(info) and info > 0.0)
        self.assertTrue(math.isfinite(req) and req > 0.0)
        self.assertGreaterEqual(info, req)
        self.assertTrue(self.d["translation_batch_prior_condition_pass"])
        self.assertEqual(self.d["useful_gate"], 1.0e-18)

    def test_actual_applied_rs_and_all_due_updates_remain(self):
        self.assertTrue(self.d["actual_applied_SpectralMSE_R_S_consumed"])
        self.assertTrue(self.d["same_tight_four_S_record_covariance_consumed"])
        self.assertTrue(self.d["all_due_S_updates_remain_in_literal_word"])

    def test_no_retired_promotion_shortcuts(self):
        self.assertFalse(self.d["determinant_used_for_metric_bound"])
        self.assertFalse(self.d["blockwise_contraction_ratio_used"])
        self.assertFalse(self.d["D_W_L_W_product_used"])
        self.assertFalse(self.d["trajectory_replay_used"])
        self.assertFalse(self.d["independent_tau_sigma_RS_source_created"])
        self.assertFalse(self.d["P3_promoted"])


if __name__ == "__main__":
    unittest.main()
