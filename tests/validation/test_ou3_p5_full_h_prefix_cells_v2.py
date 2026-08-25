import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_full_h_prefix_cells_v2 as F


class Ou3P5FullHPrefixCellsV2Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = F.build()

    def test_dependency_preserving_backend_validates(self):
        d = self.d
        self.assertEqual(F.validate(d), [])
        self.assertEqual(d["active_full_matrix_backend"], "DEPENDENCY_PRESERVING_OU_KERNEL_BOUNDS")
        self.assertFalse(d["broad_tau_natural_interval_product_used"])
        self.assertTrue(d["integrated_ou_transition_monotone_endpoint_hull_used"])
        self.assertTrue(d["integrated_ou_process_positive_kernel_moment_bounds_used"])
        self.assertTrue(d["goLive_gyro_bias_covariance_includes_full_startup_RW_upper"])
        self.assertFalse(d["old_v1_natural_interval_prediction_is_promotion_route"])

    def test_full_matrix_signed_transport_semantics_survive_tightening(self):
        d = self.d
        self.assertTrue(d["full_18x18_covariance_propagated"])
        self.assertTrue(d["H_R_S_K_r_d_eff_recomputed_in_same_prefix_cell"])
        self.assertTrue(d["shipping_Joseph_update_used"])
        self.assertTrue(d["immediate_left_error_reset_congruence_used"])
        self.assertTrue(d["physical_attitude_correction_is_minus_Etheta_Kr"])
        self.assertTrue(d["signed_cayley_primitive_consumes_actual_interval_d"])
        self.assertFalse(d["signed_a_dot_c_replaced_by_independent_abs_product"])

    def test_nonclosure_still_has_concrete_matrix_witness(self):
        d = self.d
        if d["complete_q_le_8_prefix_family_closed"]:
            self.assertEqual(d["P5_FULL_H_PREFIX_MATRIX_CERTIFICATE"], "PASS")
            self.assertIsNone(d["first_failure"])
        else:
            self.assertEqual(d["P5_FULL_H_PREFIX_MATRIX_CERTIFICATE"], "NOT_ESTABLISHED")
            self.assertIsNotNone(d["first_failure"])
            self.assertIn("sample", d["first_failure"])
            self.assertIn("reason", d["first_failure"])


if __name__ == "__main__":
    unittest.main()
