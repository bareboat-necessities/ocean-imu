from __future__ import annotations
import math
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

from ou3_interval import Interval
import ou3_p5_sample1_signed_cayley_q8_v14 as V14
import ou3_p5_sample1_signed_cayley_q8_v17 as V17


class Sample1SignedCayleyQ8V17Tests(unittest.TestCase):
    def test_product_scalar_radius_matches_v14_identity(self):
        q = 0.6593778441001633
        W = Interval.outward_bounds(0.75, 1.25)
        _wmin, parent = V14._qplus_from_product_scalar(q, W)
        child = V17._q_upper_from_product_scalar(q, W)
        self.assertEqual(parent, child)
        self.assertTrue(math.isfinite(child))

    def test_zero_crossing_product_scalar_fails_closed(self):
        q = V17._q_upper_from_product_scalar(
            0.5, Interval.outward_bounds(-0.1, 0.2))
        self.assertTrue(math.isinf(q))

    def test_component_is_intersected_with_tighter_radius(self):
        x = Interval.outward_bounds(-2.0, 3.0)
        y = V17._clip_component_to_radius(x, 0.75)
        self.assertGreaterEqual(y.lo, -0.75 - 1e-15)
        self.assertLessEqual(y.hi, 0.75 + 1e-15)
        self.assertGreaterEqual(y.lo, x.lo)
        self.assertLessEqual(y.hi, x.hi)

    def test_product_radius_can_strictly_improve_triangle_parent(self):
        qpre = 2.0
        W = Interval.outward_bounds(2.4, 2.5)
        q_product = V17._q_upper_from_product_scalar(qpre, W)
        q_triangle = V14.PREFIX2._post_correction_q_upper(qpre, 1.0)
        self.assertLess(q_product, q_triangle)

    def test_v17_keeps_shipping_limit_and_q_target(self):
        self.assertEqual(V17.Q_TARGET, 8.0)
        self.assertEqual(V17.V16.Q_TARGET, 8.0)


if __name__ == "__main__":
    unittest.main()
