import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
STABILITY = ROOT / "tools" / "stability"
if str(STABILITY) not in sys.path:
    sys.path.insert(0, str(STABILITY))

import ou3_p5_late_north_attitude_landing as L


class LateNorthAttitudeLandingTest(unittest.TestCase):
    def test_lands_inside_widest_p4_attitude_sector(self):
        d = L.build()
        self.assertEqual([], L.validate(d))
        self.assertTrue(d["LATE_NORTH_ATTITUDE_LANDING_INTO_WIDEST_P4_SECTOR_CLOSED"])
        self.assertEqual(30.0, d["P4_widest_full_attitude_candidate_deg"])
        self.assertEqual(19.0, d["P5_tilt_capture_core_deg"])
        self.assertLess(d["post_reset_full_attitude_strict_upper_deg"], 30.0)
        self.assertGreater(d["strict_margin_to_widest_P4_attitude_sector_deg"], 0.0)

    def test_does_not_shrink_p4_or_promote_full_capture(self):
        d = L.build()
        self.assertFalse(d["P4_attitude_basin_shrunk"])
        self.assertTrue(d["P5_tilt_core_is_not_P4_basin_radius"])
        self.assertFalse(d["other_P4_coordinates_landing_closed_here"])
        self.assertFalse(d["finite_time_reach_of_19deg_tilt_core_closed_here"])
        self.assertFalse(d["P4_PASS"])
        self.assertFalse(d["P5_PASS"])


if __name__ == "__main__":
    unittest.main()
