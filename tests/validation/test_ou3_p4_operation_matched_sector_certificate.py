from copy import deepcopy
from pathlib import Path
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
        self.assertTrue(self.d["validated_transcendentals_used_for_design_boundary"])
        self.assertTrue(self.d["validated_half_angle_interval_used_for_design_boundary"])

    def test_design_boundary_has_two_sided_cosine_enclosure(self):
        self.assertLessEqual(
            self.d["design_full_attitude_cosine_lower"],
            self.d["design_full_attitude_cosine_upper"],
        )
        h = self.d["design_geometry"]["validated_half_angle_interval_rad"]
        self.assertLess(h[0], h[1])
        self.assertLessEqual(h[0], 0.4)
        self.assertGreaterEqual(h[1], 0.4)

    def test_exact_vector_sector_has_material_margin(self):
        self.assertGreater(self.d["exact_vector_strong_monotonicity_factor_lower"], 0.80)
        self.assertLess(self.d["exact_eta_to_rotational_residual_information_ratio_upper"], 0.25)
        self.assertLess(self.d["exact_effective_tangent_defect_ratio_upper"], 0.45)

    def test_p1_gauged_handoffs_are_inside(self):
        o = self.d["P1_overlap"]
        self.assertTrue(o["normal_gauged_inside_sector"])
        self.assertTrue(o["timeout_gauged_inside_sector"])
        self.assertLess(o["normal_gauged_cayley_norm_upper"], self.d["design_cayley_norm_upper"])
        self.assertLess(o["timeout_gauged_cayley_norm_upper"], self.d["design_cayley_norm_upper"])

    def test_old_global_defect_route_is_not_reintroduced(self):
        self.assertFalse(self.d["global_packet_count_times_lipschitz_defect_used"])
        self.assertFalse(self.d["whole_word_weakest_P3_delta_used_as_attitude_sector_margin"])
        self.assertFalse(self.d["P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE"])

    def test_boundary_or_dead_route_mutations_fail_closed(self):
        for key, value in (
            ("validated_half_angle_interval_used_for_design_boundary", False),
            ("global_packet_count_times_lipschitz_defect_used", True),
            ("whole_word_weakest_P3_delta_used_as_attitude_sector_margin", True),
        ):
            with self.subTest(key=key):
                d = deepcopy(self.d)
                d[key] = value
                self.assertNotEqual(S.validate(d), [])
        d = deepcopy(self.d)
        d["design_full_attitude_cosine_lower"] = 0.8
        d["design_full_attitude_cosine_upper"] = 0.7
        self.assertNotEqual(S.validate(d), [])


if __name__ == "__main__":
    unittest.main()
