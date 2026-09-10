import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
STABILITY = ROOT / "tools" / "stability"
if str(STABILITY) not in sys.path:
    sys.path.insert(0, str(STABILITY))

import ou3_p5_magnetic_progress_admission as M


class P5MagneticProgressAdmissionTest(unittest.TestCase):
    def test_default_shipping_progress_closes_finite_north_event(self):
        d = M.build()
        self.assertEqual([], M.validate(d))
        self.assertTrue(d["DEFAULT_SHIPPING_PATH_PROGRESS_IMPLIES_MAGAUTOTUNER_READY"])
        self.assertTrue(d["FINITE_NORTH_EVENT_UNDER_ADMISSION_CLOSED"])
        self.assertEqual(150.0, d["finite_north_event_time_upper_from_progress_origin_s"])

    def test_progress_is_observability_only(self):
        d = M.build()
        self.assertFalse(d["BRMM_motion_caps_changed"])
        self.assertFalse(d["BIAS_family_changed"])
        self.assertFalse(d["admission"]["with_mag_false_requires_progress"])
        self.assertTrue(d["admission"]["packet_recurrence_alone_is_not_substitute"])

    def test_horizon_can_be_widened_without_changing_shipping_semantics(self):
        d = M.build(progress_horizon_s=300.0)
        self.assertEqual([], M.validate(d))
        self.assertEqual(300.0, d["finite_north_event_time_upper_from_progress_origin_s"])
        self.assertTrue(d["FINITE_NORTH_EVENT_UNDER_ADMISSION_CLOSED"])
        self.assertFalse(d["P5_PASS"])


if __name__ == "__main__":
    unittest.main()
