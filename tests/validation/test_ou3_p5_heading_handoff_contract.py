import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_heading_handoff_contract as HEADING


class Ou3P5HeadingHandoffContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = HEADING.build()

    def test_gravity_cosine_is_not_reused_as_full_attitude_bound(self):
        d = self.d
        self.assertEqual(HEADING.validate(d), [])
        self.assertTrue(d["P1_gravity_cosines_are_tilt_only"])
        q = d["gauged_quality_handoff"]
        self.assertLess(q["full_attitude_cosine_lower"], q["tilt_cosine_lower"])
        self.assertGreater(q["full_attitude_cayley_norm_upper"], 0.20)
        self.assertLess(q["full_attitude_cayley_norm_upper"], 0.35)

    def test_gauged_timeout_is_finite_but_wider_than_quality(self):
        d = self.d
        qn = d["gauged_quality_handoff"]["full_attitude_cayley_norm_upper"]
        qt = d["gauged_timeout_subbranch"]["full_attitude_cayley_norm_upper"]
        self.assertGreater(qt, qn)
        self.assertGreater(qt, 0.5)
        self.assertLess(qt, 0.7)

    def test_ungauged_timeout_routes_to_yaw_quotient(self):
        u = self.d["ungauged_timeout_subbranch"]
        self.assertFalse(u["source_timeout_requires_north_ready"])
        self.assertTrue(u["source_can_handoff_with_pending_yaw_abs_nan"])
        self.assertFalse(u["full_heading_cayley_bound_available"])
        self.assertIn("YAW_QUOTIENT", u["required_route"])
        self.assertFalse(self.d["timeout_full_heading_node_covers_complete_timeout_family"])


if __name__ == "__main__":
    unittest.main()
