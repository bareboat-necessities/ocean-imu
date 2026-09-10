import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
STABILITY = ROOT / "tools" / "stability"
if str(STABILITY) not in sys.path:
    sys.path.insert(0, str(STABILITY))

import ou3_startup_magnetic_observability as OBS


class StartupMagneticObservabilityTest(unittest.TestCase):
    def test_current_source_has_structural_yaw_obstruction(self):
        d = OBS.build()
        self.assertEqual([], OBS.validate(d))
        self.assertTrue(d["structural_observability_obstruction_proved"])
        self.assertFalse(d["UNCONDITIONAL_FULL_ATTITUDE_FINITE_CAPTURE_FROM_CURRENT_SOURCE"])
        self.assertEqual(45.0, d["yaw_indistinguishability_witness"]["minimax_full_attitude_error_lower_deg"])

    def test_shipping_acquisition_thresholds_are_exact(self):
        d = OBS.build()
        m = d["shipping_magnetic_acquisition"]
        self.assertEqual(128, m["min_accepted_samples"])
        self.assertEqual(15.0, m["min_accepted_window_s"])
        self.assertEqual(150.0, m["proxy_startup_timeout_s"])
        self.assertFalse(m["timeout_requires_north"])
        self.assertTrue(m["timeout_allows_ungauged_live"])

    def test_theorem_split_does_not_shrink_motion_domain(self):
        d = OBS.build()
        ext = d["weakest_required_source_extension"]
        self.assertTrue(ext["same_history_required"])
        self.assertTrue(ext["does_not_shrink_BRMM_motion_caps"])
        split = d["theorem_split_required"]
        self.assertTrue(split["gravity_quotient_capture_before_north"])
        self.assertTrue(split["full_attitude_P4_capture_only_after_north_gauge_event"])
        self.assertTrue(split["late_north_event_is_hybrid_yaw_reset"])


if __name__ == "__main__":
    unittest.main()
