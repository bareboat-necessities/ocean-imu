from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import ou3_sea3_rs_word_information as mod  # noqa: E402


class Sea3RsWordInformationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = mod.build()

    def test_source_uniform_four_s_geometry_and_matrix_information_close(self):
        self.assertEqual(mod.validate(self.d), [])
        self.assertTrue(self.d["P3_RS_TRANSLATION_OBSERVATION_GEOMETRY_CLOSED"])
        self.assertTrue(self.d["P3_RS_BATCH_NOISE_UPPER_CLOSED"])
        self.assertTrue(self.d["P3_RS_TRANSLATION_INFORMATION_MATRIX_CLOSED"])
        self.assertTrue(self.d["four_S_translation_observation_operator_full_rank"])
        self.assertGreater(self.d["aw_scaled_third_divided_difference_lower"], 0.0)
        self.assertGreater(
            self.d["scaled_observation_determinant_abs_lower_rank_witness_only"], 0.0
        )
        self.assertEqual(len(self.d["four_S_windows"]), 4)

    def test_newton_information_avoids_determinant_frobenius_scalarization(self):
        ni = self.d["newton_coordinate_information"]
        self.assertTrue(ni["full_4x4_matrix_inequality_closed"])
        self.assertFalse(ni["determinant_used_for_information_lower"])
        self.assertFalse(ni["frobenius_singular_value_conversion_used"])
        self.assertFalse(ni["scalar_information_beta_used"])
        self.assertLessEqual(ni["L_inverse_one_norm_upper"], 2.400000000000001)
        self.assertLessEqual(ni["L_inverse_infinity_norm_upper"], 2.000000000000001)
        self.assertGreater(ni["L_transpose_L_lambda_min_lower"], 0.17)
        self.assertGreater(ni["D_S_newton_lambda_min_lower"], 0.0)

    def test_no_retired_source_graph_or_final_gate_reenters(self):
        self.assertFalse(self.d["old_P2_800_state_graph_consumed"])
        self.assertFalse(self.d["source_history_graph_consumed"])
        self.assertFalse(self.d["predecessor_path_enumeration_consumed"])
        self.assertTrue(
            self.d["determinant_used_only_as_rank_witness_not_eigenvalue_scalarization"]
        )
        self.assertFalse(self.d["P3_PROMOTED"])

    def test_applied_rs_is_not_replaced_by_instantaneous_target(self):
        lag = self.d["SEA3_target_vs_applied_RS_contract"]
        self.assertTrue(lag["SpectralMSE_target_same_tau_sigma_TS_point"])
        self.assertTrue(lag["applied_RS_has_separate_EMA"])
        self.assertFalse(lag["instantaneous_target_formula_substituted_for_applied_RS"])
        self.assertTrue(lag["active_applied_RS_safety_ceiling_used_here"])
        self.assertTrue(lag["complete_word_still_requires_same_joint_source_path"])
        self.assertGreater(
            self.d["selected_S_record_noise"]["Sigma_S_inverse_scalar_lower"], 0.0
        )


if __name__ == "__main__":
    unittest.main()
