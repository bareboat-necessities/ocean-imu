import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_full_h_prefix_cells_v3 as F


class Ou3P5FullHPrefixCellsV3Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = F.build()

    def test_v3_backend_validates_and_keeps_full_matrix_transport(self):
        self.assertEqual(F.validate(self.d), [])
        self.assertTrue(self.d["full_18x18_covariance_propagated"])
        self.assertTrue(self.d["H_R_S_K_r_d_eff_recomputed_in_same_prefix_cell"])
        self.assertTrue(self.d["shipping_Joseph_update_used"])
        self.assertTrue(self.d["immediate_left_error_reset_congruence_used"])

    def test_retired_three_rad_correction_gate_is_absent(self):
        self.assertFalse(self.d["correction_norm_three_rad_is_promotion_gate"])
        self.assertFalse(self.d["correction_cayley_coordinate_formed_before_group_composition"])
        self.assertTrue(self.d["deployed_quaternion_composed_before_result_cayley"])
        self.assertTrue(self.d["only_resulting_error_antipode_is_group_chart_gate"])
        self.assertGreaterEqual(self.d["maximum_validated_deployed_correction_norm_rad"], 6.0)

    def test_numerical_nonclosure_if_any_is_not_hidden(self):
        if self.d["P5_FULL_H_PREFIX_MATRIX_CERTIFICATE"] != "PASS":
            self.assertFalse(self.d["complete_q_le_8_prefix_family_closed"])
            self.assertIsNotNone(self.d["first_failure"])
        else:
            self.assertTrue(self.d["complete_q_le_8_prefix_family_closed"])
            self.assertIsNone(self.d["first_failure"])


if __name__ == "__main__":
    unittest.main()
