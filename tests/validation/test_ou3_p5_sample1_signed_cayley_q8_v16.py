from __future__ import annotations
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

from ou3_interval import Interval
import ou3_p5_sample1_signed_cayley_q8_v14d as V14D
import ou3_p5_sample1_signed_cayley_q8_v16 as V16


class Sample1SignedCayleyQ8V16Tests(unittest.TestCase):
    def test_completed_416_worst_witness_retains_positive_dominant_x_axis(self):
        dbox = [
            Interval.outward_bounds(3.2214254998487846, 5.268556327935069),
            Interval.outward_bounds(-0.2384221965655676, 0.161345668536564),
            Interval.outward_bounds(-0.14602278716528508, 0.3028520932014494),
        ]
        v = V16._axis_vector_cone(
            dbox, 3.2214254998487837, 5.269730943810957)
        self.assertGreater(v[0].lo, 0.45)
        self.assertLess(v[1].abs_upper(), 0.09)
        self.assertLess(v[2].abs_upper(), 0.11)

    def test_negative_dominant_axis_keeps_sign(self):
        dbox = [
            Interval.outward_bounds(-1.5974490399960035, -0.7508816110412719),
            Interval.outward_bounds(-0.6086080220702467, 0.18364306824297555),
            Interval.outward_bounds(-0.287281129385631, 2.7678068724031064),
        ]
        u = V16._axis_component_interval(
            dbox, 0, 0.7508816110412716, 2.025610674675054)
        self.assertLess(u.hi, 0.0)
        self.assertLess(u.lo, u.hi)
        self.assertLessEqual(abs(u.lo), 1.0 + 1e-15)

    def test_axis_cone_intersects_and_tightens_v14d_parent(self):
        dbox = [
            Interval.outward_bounds(3.2214254998487846, 5.268556327935069),
            Interval.outward_bounds(-0.2384221965655676, 0.161345668536564),
            Interval.outward_bounds(-0.14602278716528508, 0.3028520932014494),
        ]
        w0, v0, b0 = V14D.radial_sinc_normalized_shipping_quaternion(
            dbox, radial_lower=3.2214254998487837,
            radial_upper=5.269730943810957)
        w1, v1, b1, narrowed = V16.axis_cone_normalized_shipping_quaternion(
            dbox, radial_lower=3.2214254998487837,
            radial_upper=5.269730943810957,
            parent=V14D.radial_sinc_normalized_shipping_quaternion)
        self.assertTrue(narrowed)
        self.assertEqual(w0.as_list(), w1.as_list())
        self.assertEqual(b0, b1)
        self.assertGreater(v1[0].lo, v0[0].lo)
        self.assertLessEqual(v1[0].hi, v0[0].hi)
        for old, new in zip(v0, v1):
            self.assertGreaterEqual(new.lo, old.lo)
            self.assertLessEqual(new.hi, old.hi)

    def test_series_branch_is_left_unchanged(self):
        dbox = [
            Interval.outward_bounds(-0.003, 0.003),
            Interval.outward_bounds(-0.003, 0.003),
            Interval.outward_bounds(-0.003, 0.003),
        ]
        w0, v0, b0 = V14D.radial_sinc_normalized_shipping_quaternion(
            dbox, radial_lower=0.0, radial_upper=0.005)
        w1, v1, b1, narrowed = V16.axis_cone_normalized_shipping_quaternion(
            dbox, radial_lower=0.0, radial_upper=0.005,
            parent=V14D.radial_sinc_normalized_shipping_quaternion)
        self.assertFalse(narrowed)
        self.assertEqual(w0.as_list(), w1.as_list())
        self.assertEqual([x.as_list() for x in v0], [x.as_list() for x in v1])
        self.assertEqual(b0, b1)


if __name__ == "__main__":
    unittest.main()
