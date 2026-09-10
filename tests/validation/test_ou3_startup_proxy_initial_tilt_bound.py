from pathlib import Path
import math
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools" / "stability"))

import ou3_startup_proxy_initial_tilt_bound as mod  # noqa: E402


class StartupProxyInitialTiltBoundTests(unittest.TestCase):
    def test_declared_physical_seed_is_inside_mahony_chart(self):
        d = mod.build()
        self.assertEqual(mod.validate(d), [])
        self.assertAlmostEqual(d["physical_initial_tilt_error_upper_deg"], 24.0721360794, places=8)
        self.assertTrue(d["physical_initialization_strictly_inside_chart"])
        self.assertGreater(d["physical_chart_margin_rad_before_sensor_fp_charges"], 0.0)

    def test_initial_lemma_does_not_overclaim_capture(self):
        d = mod.build()
        self.assertFalse(d["sensor_direction_error_charge_closed_here"])
        self.assertFalse(d["binary32_direction_error_charge_closed_here"])
        self.assertFalse(d["whole_startup_trajectory_capture_closed_here"])
        self.assertFalse(d["P5_PASS"])


if __name__ == "__main__":
    unittest.main()
