import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
STABILITY = ROOT / "tools" / "stability"
if str(STABILITY) not in sys.path:
    sys.path.insert(0, str(STABILITY))

import ou3_p5_wide_handoff_h18_backbone as W


class WideHandoffH18BackboneTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = W.build()

    def test_wide_sector_contains_shipping_handoff(self):
        self.assertEqual([], W.validate(self.d))
        self.assertTrue(self.d["strictly_contains_complete_handoff_tilt_section"])
        self.assertLess(self.d["shipping_handoff_tilt_strict_upper_deg"], self.d["transient_outer_angle_deg"])
        self.assertGreater(self.d["shipping_handoff_tilt_strict_upper_deg"], 91.0)

    def test_prior_free_h18_ldlt_survives_wide_sector(self):
        self.assertTrue(self.d["WIDE_HANDOFF_H18_PRIOR_FREE_LDLT_CLOSED"])
        self.assertGreater(self.d["wide_vector_information_retention_lower"], 0.0)
        self.assertGreater(self.d["wide_finite_angle_H18_information_lower"], 1e-18)
        self.assertGreater(self.d["worst_full_H18_LDLT_pivot_lower"], 0.0)
        self.assertFalse(self.d["finite_angle_information_loss_above_45deg_is_P5_blocker"])

    def test_backbone_does_not_claim_finite_map_capture(self):
        self.assertFalse(self.d["finite_map_contraction_closed_here"])
        self.assertFalse(self.d["finite_capture_time_closed_here"])
        self.assertFalse(self.d["P4_PASS"])
        self.assertFalse(self.d["P5_PASS"])


if __name__ == "__main__":
    unittest.main()
