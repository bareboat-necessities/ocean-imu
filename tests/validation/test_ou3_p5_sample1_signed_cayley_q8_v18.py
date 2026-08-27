from __future__ import annotations
import math
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

from ou3_interval import Interval
import ou3_p5_sample1_signed_cayley_q8_v14 as V14
import ou3_p5_sample1_signed_cayley_q8_v18 as V18


class Sample1SignedCayleyQ8V18Tests(unittest.TestCase):
    def test_zero_proof_gauge_rotation_preserves_yz(self):
        cy = Interval.outward_bounds(-0.3, 0.4)
        cz = Interval.outward_bounds(0.2, 0.5)
        y, z = V18._rotate_yz_rx_transpose(cy, cz, Interval.point(0.0))
        self.assertLessEqual(y.lo, cy.lo)
        self.assertGreaterEqual(y.hi, cy.hi)
        self.assertLessEqual(z.lo, cz.lo)
        self.assertGreaterEqual(z.hi, cz.hi)

    def test_signed_negative_proof_gauge_rotation_is_supported(self):
        cy = Interval.outward_bounds(-0.3, 0.4)
        cz = Interval.outward_bounds(0.2, 0.5)
        y, z = V18._rotate_yz_rx_transpose(
            cy, cz, Interval.outward_bounds(-1.0, -0.5))
        self.assertLessEqual(y.lo, y.hi)
        self.assertLessEqual(z.lo, z.hi)

    def test_quarter_turn_proof_gauge_rotation_has_expected_orientation(self):
        one = Interval.point(1.0)
        zero = Interval.point(0.0)
        y, z = V18._rotate_yz_rx_transpose(
            one, zero, Interval.point(math.pi / 2.0))
        self.assertLessEqual(abs(y.lo), 2e-16)
        self.assertLessEqual(abs(y.hi), 2e-16)
        self.assertLessEqual(z.lo, -1.0)
        self.assertGreaterEqual(z.hi, -1.0)

    def test_componentwise_yz_support_strictly_tightens_cauchy_parent(self):
        parent_W = Interval.outward_bounds(1.0, 3.0)
        wd = Interval.point(1.0)
        vd = [Interval.point(0.0), Interval.point(0.5), Interval.point(0.0)]
        chart = {
            "cx": Interval.point(0.0),
            "cy": Interval.outward_bounds(0.2, 0.3),
            "cz": Interval.outward_bounds(-0.1, 0.1),
            "cyz_norm_upper": 1.0,
        }
        W, yz_box, yz_joint = V18._support_product_scalar(
            parent_W, wd, vd, chart)
        self.assertGreater(W.lo, parent_W.lo)
        self.assertLess(W.hi, parent_W.hi)
        self.assertGreater(yz_joint.lo, -0.5)
        self.assertLess(yz_joint.hi, 0.5)
        self.assertGreaterEqual(yz_joint.lo, yz_box.lo)
        self.assertLessEqual(yz_joint.hi, yz_box.hi)

    def test_tighter_product_scalar_cannot_worsen_q_bound(self):
        parent_W = Interval.outward_bounds(0.2, 3.0)
        support_W = Interval.outward_bounds(1.0, 2.0)
        wp, qp = V14._qplus_from_product_scalar(2.0, parent_W)
        ws, qs = V14._qplus_from_product_scalar(2.0, support_W)
        self.assertGreaterEqual(ws, wp)
        self.assertLessEqual(qs, qp)

    def test_v18_keeps_shipping_limit_and_q_target(self):
        self.assertEqual(V18.Q_TARGET, 8.0)
        self.assertEqual(V14.Q_TARGET, 8.0)


if __name__ == "__main__":
    unittest.main()
