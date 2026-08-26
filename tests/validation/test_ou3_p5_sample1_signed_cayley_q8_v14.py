from __future__ import annotations
import math
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
from ou3_interval import Interval
import ou3_p5_sample1_signed_cayley_q8_v14 as V14


class Sample1SignedCayleyQ8V14Tests(unittest.TestCase):
    def test_normalized_axis_quaternion_keeps_signed_direction(self):
        d = [Interval.outward_bounds(-7.2, -6.8),
             Interval.outward_bounds(-0.1, 0.1),
             Interval.outward_bounds(-0.1, 0.1)]
        w, v, branches = V14._normalized_shipping_quaternion(
            d, radial_lower=6.8, radial_upper=7.21)
        self.assertIn("AXIS_ANGLE_UNIT", branches)
        for x in (w, *v):
            self.assertTrue(math.isfinite(x.lo))
            self.assertTrue(math.isfinite(x.hi))
        # Beyond 2*pi, sin(theta/2)<0.  With negative d_x the normalized
        # quaternion vector x component is therefore positive.
        self.assertGreater(v[0].lo, 0.0)

    def test_series_branch_is_source_normalized(self):
        d = [Interval.outward_bounds(-0.005, 0.006),
             Interval.outward_bounds(-0.001, 0.001),
             Interval.outward_bounds(-0.001, 0.001)]
        w, v, branches = V14._normalized_shipping_quaternion(
            d, radial_lower=0.0, radial_upper=0.0062)
        self.assertIn("SERIES_NORMALIZED", branches)
        self.assertGreater(w.lo, 0.0)
        self.assertLessEqual(max(x.abs_upper() for x in v), 0.01)

    def test_qplus_identity_is_strict_when_product_scalar_is_separated(self):
        W = Interval.outward_bounds(-2.0, -1.8)
        wmin, q = V14._qplus_from_product_scalar(3.2, W)
        self.assertGreater(wmin, 0.0)
        self.assertTrue(math.isfinite(q))
        self.assertLess(q, 8.0)

    def test_unknown_small_rotation_widens_without_losing_finiteness(self):
        cx = Interval.outward_bounds(0.8, 1.2)
        y = V14._widen_cx_by_unknown_rotation(cx, 3.0, 1e-3)
        self.assertLess(y.lo, cx.lo)
        self.assertGreater(y.hi, cx.hi)
        self.assertTrue(math.isfinite(y.lo) and math.isfinite(y.hi))

    def test_coarse_family_fails_closed_or_closes(self):
        d = V14.build(source_pieces=4, source_cell_index=0,
                      p_pieces=2, tangent_pieces=2, axial_pieces=2,
                      residual_x_pieces=2, parallel_pieces=2)
        st = d["P5_SAMPLE1_SIGNED_CAYLEY_Q8_V14"]
        self.assertIn(st, ("PASS", "NOT_ESTABLISHED"))
        self.assertEqual(float(d["deployed_correction_limit_rad"]), 6.0)
        self.assertIs(d["deployed_correction_limit_increased"], False)
        self.assertIs(d["q8_word_promoted_here"], False)
        self.assertIs(d["whole_word_promoted_here"], False)
        self.assertIs(d["N_H_words_set_here"], False)
        if st == "PASS":
            self.assertEqual(V14.validate(d), [])
            self.assertLess(d["max_post_sample1_cayley_norm_upper"], 8.0)
        else:
            # Coarse V13E/V12D fixtures are intentionally allowed to stop at
            # their prerequisite or to return a concrete V14 q8 witness.
            vf = V14.validate(d)
            self.assertTrue(vf or d["first_unclosed_q8_cell"] is not None)


if __name__ == "__main__":
    unittest.main()
