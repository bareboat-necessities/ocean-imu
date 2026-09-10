import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
STABILITY = ROOT / "tools" / "stability"
if str(STABILITY) not in sys.path:
    sys.path.insert(0, str(STABILITY))

import ou3_startup_postlive_regauge_contract as C


class StartupPostLiveRegaugeContractTest(unittest.TestCase):
    def test_contract(self):
        d = C.build()
        self.assertEqual([], C.validate(d))
        self.assertTrue(d["FINITE_POSTLIVE_NORTH_EVENT_UNDER_ADMISSION_CLOSED"])
        self.assertFalse(d["normal_live_PE_alone_implies_first_north"])
        self.assertFalse(d["startup_magnetic_PE_required_for_timeout_branch"])
        self.assertTrue(d["late_north_yaw_landing_bound_still_required"])
        self.assertFalse(d["P5_PASS"])


if __name__ == "__main__":
    unittest.main()
