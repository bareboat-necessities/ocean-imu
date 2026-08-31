from pathlib import Path
import math
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import ou3_p4_operation_matched_sector_certificate as S


class P4OperationMatchedSectorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = S.build()

    def test_validates(self):
        self.assertEqual(S.validate(self.d), [])

    def test_sector_is_physical_not_microscopic(self):
        self.assertGreaterEqual(self.d["design_full_attitude_angle_rad"], 0.80)
        self.assertGreater(self.d["design_full_attitude_angle_deg"], 45.0)
        self.assertLess(self.d["design_cayley_norm_upper"], 1.0)

    def test_exact_vector_sector_has_material_margin(self):
        self.assertGreater(self.d["exact_vector_strong_monotonicity_factor_lower"], 0.80)
        self.assertLess(self.d["exact_eta_to_rotational_residual_information_ratio_upper"], 0.25)
        self.assertLess(self.d["exact_effective_tangent_defect_ratio_upper"], 0.45)

    def test_p1_gauged_handoffs_are_inside(self):
        o = self.d["P1_overlap"]
        self.assertTrue(o["normal_gauged_inside_sector"])
        self.assertTrue(o["timeout_gauged_inside_sector"])
        self.assertLess(o["normal_gauged_angle_upper_rad"], 0.80)
        self.assertLess(o["timeout_gauged_angle_upper_rad"], 0.80)

    def test_old_global_defect_route_is_not_reintroduced(self):
        self.assertFalse(self.d["global_packet_count_times_lipschitz_defect_used"])
        self.assertFalse(self.d["whole_word_weakest_P3_delta_used_as_attitude_sector_margin"])
        self.assertFalse(self.d["P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE"])


if __name__ == "__main__":
    unittest.main()
