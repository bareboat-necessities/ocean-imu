#!/usr/bin/env python3
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import ou3_sea3_a21_prior_free_completion as A21


class Sea3A21MixedCompletionStatusTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = A21.build()
        cls.failures = A21.validate(cls.d)

    def test_complete_sea3_and_paper_detectability_route_are_retained(self):
        d = self.d
        self.assertEqual([], self.failures)
        self.assertEqual("COMPLETE_SEA3_NORMAL_LIVE_WORD", d["canonical_source"])
        self.assertEqual(3.0, d["complete_word_horizon_s"])
        self.assertEqual("ETA6_PLUS_FINITE_RESIDUAL_BIAS_CORRELATION", d["paper_active_bias_route"])
        self.assertFalse(d["eta9_point_packet_shortcut_used"])
        self.assertTrue(d["A21_finite_bias_detectability_consumed"])
        self.assertTrue(d["H18_prior_free_completion_consumed"])

    def test_invalid_full_D_inverse_promotion_is_rejected(self):
        d = self.d
        self.assertTrue(d["finite_tau_detectability_does_not_imply_full_A21_information_inverse"])
        self.assertFalse(d["full_A21_D_inverse_available"])
        self.assertFalse(d["full_A21_prior_free_D_inverse_identity_used"])
        self.assertTrue(d["invalid_append_Qba_to_H18_full_D_completion_rejected"])
        self.assertFalse(d["full_21x21_Omega_minus_delta_P_LDLT_closed"])
        self.assertFalse(d["A21_prior_free_completion_closed"])
        self.assertFalse(d["full_21x21_interval_LDLT_used"])
        self.assertFalse(d["P3_promoted"])

    def test_finite_bias_prior_is_quantitatively_tied_to_shipping_process(self):
        d = self.d
        rel = d["source_generated_bias_release"]
        self.assertTrue(rel["H18_ba_cross_covariances_zero_at_release"])
        self.assertGreater(rel["bias_prior_variance"], 0.0)
        self.assertGreater(rel["first_active_prediction_ba_process_lambda_min_lower"], 0.0)
        dom = d["stable_bias_prior_process_domination"]
        self.assertGreater(dom["ratio_c_b_upper"], 0.0)
        self.assertTrue(math.isfinite(dom["ratio_c_b_upper"]))
        self.assertGreater(dom["delta_times_ratio_upper"], 0.0)
        self.assertLess(dom["delta_times_ratio_upper"], 1.0)
        self.assertTrue(dom["delta_times_ratio_is_small"])
        self.assertTrue(dom["same_congruence_preserves_tag_domination_after_first_prediction"])
        self.assertTrue(dom["later_total_Omega_contains_tag_plus_additional_PSD_noise"])
        self.assertFalse(dom["this_lemma_alone_closes_full_A21_Riccati_inequality"])

    def test_actual_rs_and_no_replacement_source(self):
        d = self.d
        self.assertTrue(d["actual_applied_SpectralMSE_R_S_retained_through_H18_component"])
        self.assertTrue(d["event_algebra_preserves_margin_after_closure"])
        self.assertEqual(1.0e-18, d["useful_gate"])
        for key in (
            "source_family_replaced",
            "trajectory_replay_used",
            "independent_tau_sigma_RS_source_created",
        ):
            self.assertFalse(d[key], key)


if __name__ == "__main__":
    unittest.main()