import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_effective_vector_input as E


class Ou3P5EffectiveVectorInputTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = E.build()

    def test_source_bound_reduction_passes_without_filter_change(self):
        d = self.d
        self.assertEqual(E.validate(d), [])
        self.assertEqual(d["P5_EFFECTIVE_VECTOR_INPUT_CERTIFICATE"], "PASS")
        self.assertTrue(d["source_generated_not_trajectory_fit"])
        self.assertFalse(d["source_replay_used"])
        self.assertFalse(d["filter_changed"])
        self.assertTrue(d["standalone_vector_eta_penalty_retired_from_P5_numerical_route"])
        self.assertFalse(d["eta_declared_identically_zero"])
        self.assertTrue(d["joseph_information_identity_remains_valid"])

    def test_configured_magnetometer_radial_residual_is_exactly_killed(self):
        d = self.d
        src = d["source_semantics"]
        self.assertTrue(src["configured_magnetometer_covariance_isotropic"])
        m = d["magnetometer"]
        self.assertTrue(m["H_theta_transpose_v_exact_zero"])
        self.assertTrue(m["cross_covariance_action_on_v_exact_zero"])
        self.assertTrue(m["kalman_gain_radial_action_exact_zero"])
        self.assertIn("K_m y_m=K_m H_m d_m", m["exact_state_correction_identity"])
        self.assertTrue(m["effective_coordinate_nonexpansive"])
        self.assertFalse(m["standalone_radial_eta_changes_state"])
        self.assertFalse(m["standalone_eta_information_penalty_required_for_state_correction"])

    def test_accelerometer_eta_is_exact_effective_aw_input(self):
        d = self.d
        src = d["source_semantics"]
        self.assertTrue(src["normal_live_accelerometer_Jaw_is_Rwb"])
        self.assertTrue(src["normal_live_accelerometer_Jaw_orthogonal_full_rank"])
        self.assertTrue(src["configured_imu_lever_arm_disabled"])
        a = d["accelerometer"]
        self.assertEqual(a["shipping_J_aw"], "R_wb")
        self.assertTrue(a["J_aw_orthogonal_full_row_rank"])
        self.assertTrue(a["effective_aw_defect_norm_equals_eta_norm"])
        self.assertIn("H_a E_aw e_eta=eta_a", a["exact_measurement_range_identity"])
        self.assertIn("K_a(H_a z+eta_a)=K_a H_a(z+E_aw e_eta)", a["exact_state_correction_identity"])
        self.assertFalse(a["standalone_eta_information_penalty_required_for_state_correction"])
        self.assertGreaterEqual(a["widened_effective_aw_defect_norm_upper_mps2"], 0.0)

    def test_cellwise_effective_coordinates_have_strict_finite_bounds(self):
        cells = self.d["annular_effective_input_cells"]
        self.assertEqual(len(cells), self.d["subdivision_cell_count"])
        self.assertGreater(len(cells), 1)
        prev = -1.0
        for row in cells:
            self.assertGreater(row["mag_effective_tangent_coordinate_gain_upper"], 0.0)
            self.assertLessEqual(row["mag_effective_tangent_coordinate_gain_upper"], 1.0)
            defect = row["mag_effective_vs_tangent_defect_ratio_upper"]
            self.assertGreaterEqual(defect, 0.0)
            self.assertLess(defect, 1.0)
            self.assertGreaterEqual(defect, prev)
            prev = defect
            self.assertGreaterEqual(row["acc_effective_aw_attitude_eta_per_vector_norm_upper"], 0.0)
            self.assertGreaterEqual(row["acc_effective_aw_latent_cross_gain_upper"], 0.0)
            self.assertLess(row["acc_effective_aw_latent_cross_gain_upper"], 2.0)

    def test_gravity_quotient_uses_same_accelerometer_effective_input(self):
        q = self.d["gravity_quotient"]
        self.assertTrue(q["accelerometer_effective_aw_input_descends_to_quotient"])
        self.assertEqual(q["axial_gyro_bias_role_unchanged"], "NEUTRAL_BOUNDED_INPUT")
        self.assertFalse(q["standalone_accelerometer_eta_penalty_required"])


if __name__ == "__main__":
    unittest.main()
