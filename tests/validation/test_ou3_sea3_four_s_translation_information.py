from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_sea3_four_s_translation_information as FOUR  # noqa: E402


class Sea3FourSTranslationInformationTest(unittest.TestCase):
    def test_full_v_p_s_aw_information_uses_complete_sea3_and_applied_rs(self):
        d = FOUR.build()
        self.assertEqual(FOUR.validate(d), [])
        self.assertEqual(d["canonical_source"], "COMPLETE_SEA3_NORMAL_LIVE_WORD")
        self.assertTrue(d["component_of_complete_SEA3_full_word"])
        self.assertFalse(d["P3_architecture_replaced"])
        self.assertFalse(d["source_family_replaced"])
        self.assertFalse(d["trajectory_replay_used"])
        self.assertFalse(d["selected_four_S_events_replace_full_scheduler_word"])
        self.assertTrue(d["all_due_S_updates_remain_in_literal_word"])
        self.assertTrue(d["complete_SEA3_same_realization_requirement_retained"])
        self.assertTrue(d["actual_applied_SpectralMSE_R_S_consumed"])
        self.assertFalse(d["instantaneous_R_S_target_substituted_for_applied_R_S"])
        self.assertTrue(d["time_varying_tau_allowed_inside_selected_subword"])
        self.assertTrue(d["time_varying_sigma_allowed_inside_selected_subword"])
        self.assertEqual(
            d["dimensionless_translation_state"],
            ["S", "g*p", "g^2*v", "g^3*a_w"],
        )
        self.assertEqual(len(d["four_S_windows_s"]), 4)
        self.assertLessEqual(d["uniform_S_gap_s_upper"], 0.151)
        self.assertGreater(
            d["scaled_observation_determinant_abs_lower_rank_witness_only"], 0.0
        )

    def test_newton_information_is_genuine_4x4_not_determinant_scalarization(self):
        d = FOUR.build()
        ni = d["newton_coordinate_information"]
        self.assertTrue(ni["full_4x4_matrix_inequality_closed"])
        self.assertFalse(ni["determinant_used_for_information_lower"])
        self.assertFalse(ni["frobenius_singular_value_conversion_used"])
        self.assertFalse(ni["scalar_information_beta_used"])
        self.assertTrue(ni["exact_rational_arithmetic_used_before_outward_float_conversion"])
        info = ni["D_S_newton_lambda_min_lower"]
        self.assertTrue(math.isfinite(info))
        self.assertGreater(info, 0.0)
        self.assertTrue(d["P3_RS_TRANSLATION_INFORMATION_MATRIX_CLOSED"])
        self.assertFalse(d["P3_promoted"])

    def test_shipping_rs_cap_and_axis_factors_are_consumed(self):
        d = FOUR.build()
        self.assertEqual(d["R_S_axis_std_factors"], [0.72, 0.72, 1.0])
        self.assertLessEqual(d["R_S_applied_base_std"][1], 100.0001)
        vars_ = d["selected_S_record_noise"]["measurement_variance_axis_upper"]
        self.assertEqual(len(vars_), 3)
        self.assertGreater(max(vars_), 0.0)


if __name__ == "__main__":
    unittest.main()
