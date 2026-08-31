from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import ou3_p5_outer_sector_capture_certificate as C


class P5OuterSectorCaptureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = C.build()

    def test_validates(self):
        self.assertEqual(C.validate(self.d), [])

    def test_capture_is_immediate_and_non_microscopic(self):
        self.assertEqual(self.d["P5_OUTER_SECTOR_CAPTURE_CERTIFICATE"], "PASS")
        self.assertEqual(self.d["N_outer_words"], 0)
        self.assertGreaterEqual(self.d["outer_sector_angle_rad"], 0.80)
        self.assertTrue(self.d["all_source_handoff_branches_enter_outer_sector"])

    def test_all_branches_are_source_faithful(self):
        b = self.d["branches"]
        self.assertTrue(b["normal_gauged"]["inside_outer_sector"])
        self.assertTrue(b["timeout_gauged"]["inside_outer_sector"])
        self.assertTrue(b["timeout_ungauged"]["inside_outer_gravity_sector"])
        self.assertFalse(b["timeout_ungauged"]["full_heading_radius_assigned"])

    def test_retired_microscopic_route_is_not_capture_definition(self):
        self.assertFalse(self.d["legacy_microscopic_inner_seed_used_as_outer_capture_target"])
        self.assertFalse(self.d["legacy_uniform_transport_route_used"])
        self.assertFalse(self.d["P5_INNER_FUNNEL_FINITE_CAPTURE_ESTABLISHED_HERE"])


if __name__ == "__main__":
    unittest.main()
