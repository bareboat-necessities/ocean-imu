import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_joint_prefix_cells as J


class Ou3P5JointPrefixCellsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = J.build()

    def test_active_payload_is_full_matrix_not_directional_fallback(self):
        d = self.d
        self.assertEqual(J.validate(d), [])
        self.assertEqual(d["active_P_payload"], "OUTWARD_FULL_18X18_H_COVARIANCE_CELL")
        self.assertTrue(d["full_signed_matrix_covariance_cells_available"])
        self.assertFalse(d["directional_P_payload_retained_as_active_backend"])
        self.assertFalse(d["old_directional_scalar_route_used_for_promotion"])
        self.assertEqual(d["P5_JOINT_PREFIX_SCALAR_CELL_CERTIFICATE"], "RETIRED_AS_ACTIVE_ROUTE")

    def test_exact_signed_measurement_reset_calculus_is_mandatory(self):
        d = self.d
        self.assertTrue(d["P_H_R_K_S_r_d_eff_recomputed_in_same_prefix_cell"])
        self.assertTrue(d["shipping_Joseph_update_used"])
        self.assertTrue(d["immediate_left_error_reset_congruence_used"])
        self.assertTrue(d["physical_attitude_correction_is_minus_Etheta_Kr"])
        self.assertTrue(d["signed_cayley_primitive_consumes_actual_interval_d"])
        self.assertFalse(d["signed_a_dot_c_replaced_by_independent_abs_product"])

    def test_tangent_mag_and_effective_acc_routes_remain_active(self):
        d = self.d
        self.assertTrue(d["magnetometer_radial_K_action_exact_zero"])
        self.assertTrue(d["magnetometer_radial_Joseph_information_exact_zero"])
        self.assertTrue(d["accelerometer_effective_aw_input_used"])
        self.assertFalse(d["standalone_vector_eta_penalty_used"])

    def test_no_NH_words_before_complete_matrix_word_closes(self):
        d = self.d
        self.assertFalse(d["N_H_words_set_here"])
        if d["P5_JOINT_PREFIX_FULL_MATRIX_CERTIFICATE"] == "PASS":
            self.assertTrue(d["signed_cayley_prefix_composition_closed"])
            self.assertTrue(d["P5_numerical_status_can_promote_from_this_stage"])
        else:
            self.assertFalse(d["P5_numerical_status_can_promote_from_this_stage"])
            self.assertNotEqual(d["first_unclosed_numerical_obligation"], "NONE_AT_COMPLETE_GAUGED_H_PREFIX")


if __name__ == "__main__":
    unittest.main()
