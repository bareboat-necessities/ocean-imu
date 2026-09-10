import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
STABILITY = ROOT / "tools" / "stability"
if str(STABILITY) not in sys.path:
    sys.path.insert(0, str(STABILITY))

import ou3_startup_handoff_tilt_hemisphere as H


class StartupHandoffTiltHemisphereTest(unittest.TestCase):
    def test_shipping_branch_closes_wide_tilt_bound(self):
        d = H.build()
        self.assertEqual([], H.validate(d))
        self.assertTrue(d["HANDOFF_TILT_HEMISPHERE_BOUND_CLOSED"])
        self.assertAlmostEqual(math.degrees(math.pi / 2 + 0.02), d["true_gravity_quotient_tilt_strict_upper_deg"], places=12)

    def test_bound_does_not_assume_old_chart_or_yaw(self):
        d = H.build()
        self.assertFalse(d["historical_60_deg_mahony_chart_consumed"])
        self.assertFalse(d["yaw_gauge_required"])
        self.assertTrue(d["applies_to_quality_and_timeout_handoff_samples"])
        self.assertTrue(d["early_live_must_contract_to_P4"])


if __name__ == "__main__":
    unittest.main()
