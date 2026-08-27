from __future__ import annotations
import math
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

from ou3_interval import Interval
import ou3_p5_sample1_signed_cayley_q8_v18b as V18B


class Sample1SignedCayleyQ8V18BTests(unittest.TestCase):
    def test_signed_negative_interval_is_accepted_and_bounded(self):
        s, c, broad = V18B._signed_full_angle_trig_interval(
            Interval.outward_bounds(-1.0, -0.5))
        self.assertFalse(broad)
        self.assertGreaterEqual(s.lo, -1.0)
        self.assertLessEqual(s.hi, 1.0)
        self.assertGreaterEqual(c.lo, -1.0)
        self.assertLessEqual(c.hi, 1.0)

    def test_zero_angle_encloses_identity_rotation(self):
        s, c, broad = V18B._signed_full_angle_trig_interval(Interval.point(0.0))
        self.assertFalse(broad)
        self.assertLessEqual(s.lo, 0.0)
        self.assertGreaterEqual(s.hi, 0.0)
        self.assertLessEqual(c.lo, 1.0)
        self.assertGreaterEqual(c.hi, 1.0)

    def test_quarter_turn_point_is_tight_without_half_angle_assumption(self):
        x = Interval.point(math.pi / 2.0)
        s, c, broad = V18B._signed_full_angle_trig_interval(x)
        self.assertFalse(broad)
        self.assertGreater(s.lo, 0.999999999999)
        self.assertLess(s.hi, 1.000000000001)
        self.assertLess(abs(c.lo), 1.0e-12)
        self.assertLess(abs(c.hi), 1.0e-12)

    def test_outside_audited_midpoint_fails_wide_not_unsafe(self):
        s, c, broad = V18B._signed_full_angle_trig_interval(
            Interval.outward_bounds(5.0, 5.2))
        self.assertTrue(broad)
        self.assertEqual((s.lo, s.hi), (-1.0, 1.0))
        self.assertEqual((c.lo, c.hi), (-1.0, 1.0))

    def test_rotation_helper_accepts_signed_angle(self):
        cy = Interval.outward_bounds(-0.3, 0.4)
        cz = Interval.outward_bounds(0.2, 0.5)
        y, z = V18B._rotate_yz_rx_transpose(
            cy, cz, Interval.outward_bounds(-1.0, -0.5))
        self.assertLessEqual(y.lo, y.hi)
        self.assertLessEqual(z.lo, z.hi)

    def test_v18b_keeps_shipping_limit_and_q_target(self):
        self.assertEqual(V18B.Q_TARGET, 8.0)
        self.assertEqual(V18B.AUDITED_POINT_ABS_MAX, 4.5)


if __name__ == "__main__":
    unittest.main()
