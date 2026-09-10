import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
STABILITY = ROOT / "tools" / "stability"
if str(STABILITY) not in sys.path:
    sys.path.insert(0, str(STABILITY))

import ou3_p5_late_north_yaw_reset as R


class LateNorthYawResetTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = R.build()

    def test_reset_structure_is_closed(self):
        self.assertEqual([], R.validate(self.d))
        self.assertTrue(self.d["LATE_NORTH_HYBRID_RESET_MAP_STRUCTURE_CLOSED"])
        self.assertTrue(self.d["gravity_quotient_tilt_error_exactly_invariant"])
        self.assertTrue(self.d["covariance_exactly_unchanged_by_setter"])

    def test_nonattitude_states_are_unchanged(self):
        self.assertTrue(self.d["linear_navigation_state_exactly_unchanged"])
        self.assertTrue(self.d["gyro_bias_state_exactly_unchanged"])
        self.assertTrue(self.d["accelerometer_bias_state_exactly_unchanged"])
        self.assertTrue(self.d["local_attitude_error_state_zeroed_after_rewrite"])

    def test_yaw_landing_error_remains_separate(self):
        self.assertFalse(self.d["magnetic_reference_yaw_error_bound_closed_here"])
        self.assertFalse(self.d["full_attitude_P4_landing_closed_here"])
        self.assertFalse(self.d["P5_PASS"])


if __name__ == "__main__":
    unittest.main()
