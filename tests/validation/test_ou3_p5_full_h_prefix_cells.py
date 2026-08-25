import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_full_h_prefix_cells as F


class Ou3P5FullHPrefixCellsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = F.build()

    def test_full_matrix_stage_is_source_bound_and_validates(self):
        d = self.d
        self.assertEqual(F.validate(d), [])
        self.assertTrue(d["source_generated_not_trajectory_fit"])
        self.assertFalse(d["source_replay_used"])
        self.assertFalse(d["filter_changed"])
        self.assertEqual(d["mode"], "H")
        self.assertEqual(d["dimension"], 18)
        self.assertTrue(d["full_18x18_covariance_propagated"])

    def test_same_cell_recomputes_full_kalman_payload(self):
        d = self.d
        self.assertTrue(d["H_R_S_K_r_d_eff_recomputed_in_same_prefix_cell"])
        self.assertTrue(d["shipping_Joseph_update_used"])
        self.assertTrue(d["immediate_left_error_reset_congruence_used"])
        self.assertTrue(d["physical_attitude_correction_is_minus_Etheta_Kr"])
        self.assertTrue(d["signed_cayley_primitive_consumes_actual_interval_d"])
        self.assertFalse(d["signed_a_dot_c_replaced_by_independent_abs_product"])

    def test_vector_reductions_are_active_not_standalone_eta_budget(self):
        d = self.d
        self.assertTrue(d["magnetometer_radial_K_action_exact_zero"])
        self.assertTrue(d["accelerometer_effective_aw_input_used"])
        self.assertFalse(d["standalone_vector_eta_penalty_used"])

    def test_numerical_nonclosure_is_fail_closed_with_witness(self):
        d = self.d
        if d["complete_q_le_8_prefix_family_closed"]:
            self.assertEqual(d["P5_FULL_H_PREFIX_MATRIX_CERTIFICATE"], "PASS")
            self.assertIsNone(d["first_failure"])
            self.assertIsNotNone(d["smaller_source_reachable_chart_upper"])
            self.assertLess(d["smaller_source_reachable_chart_upper"], 8.0)
        else:
            self.assertEqual(d["P5_FULL_H_PREFIX_MATRIX_CERTIFICATE"], "NOT_ESTABLISHED")
            self.assertIsNotNone(d["first_failure"])
            self.assertIn("sample", d["first_failure"])
            self.assertIn("reason", d["first_failure"])

    def test_source_cell_is_joint_complete_word_invariant(self):
        s = self.d["source_cell"]
        self.assertTrue(s["joint_invariant_cell_over_complete_word"])
        for key in ("tau_s", "sigma_aw_mps2", "R_S_filter_std", "pseudo_period_s"):
            self.assertEqual(len(s[key]), 2)
            self.assertLessEqual(s[key][0], s[key][1])


if __name__ == "__main__":
    unittest.main()
