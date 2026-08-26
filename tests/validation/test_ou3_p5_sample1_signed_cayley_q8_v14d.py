from __future__ import annotations
import math
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
from ou3_interval import Interval
import ou3_p5_sample1_signed_cayley_q8_v14d as V14D


class Sample1SignedCayleyQ8V14DTests(unittest.TestCase):
    def test_broad_sub6_axis_cell_keeps_quaternion_components_physical(self):
        d = [Interval.outward_bounds(-0.35, 0.02),
             Interval.outward_bounds(-0.02, 0.35),
             Interval.outward_bounds(-0.42, 0.73)]
        w, v, branches = V14D.radial_sinc_normalized_shipping_quaternion(
            d, radial_lower=0.0, radial_upper=0.75)
        self.assertIn("SERIES_NORMALIZED_RADIAL_CLIP", branches)
        self.assertIn("AXIS_ANGLE_MONOTONE_SINC_RADIAL", branches)
        self.assertLessEqual(w.abs_upper(), 1.001)
        for x in v:
            self.assertLessEqual(x.abs_upper(), 1.001)

    def test_winding_branch_is_finite_and_uses_v2_trig(self):
        d = [Interval.outward_bounds(-7.4, -6.4),
             Interval.outward_bounds(-0.1, 0.1),
             Interval.outward_bounds(-0.1, 0.1)]
        w, v, branches = V14D.radial_sinc_normalized_shipping_quaternion(
            d, radial_lower=6.4, radial_upper=7.4)
        self.assertEqual(branches, ["AXIS_ANGLE_V2_WINDING_RADIAL"])
        for x in (w, *v):
            self.assertTrue(math.isfinite(x.lo) and math.isfinite(x.hi))
        self.assertLessEqual(w.abs_upper(), 1.001)
        for x in v:
            self.assertLessEqual(x.abs_upper(), 1.3)

    def test_no_filter_limit_or_promotion_change(self):
        d = V14D.build(source_pieces=4, source_cell_index=0,
                       p_pieces=2, tangent_pieces=2, axial_pieces=2,
                       residual_x_pieces=2, parallel_pieces=2)
        self.assertIn(d["P5_SAMPLE1_SIGNED_CAYLEY_Q8_V14D"],
                      ("PASS", "NOT_ESTABLISHED"))
        self.assertEqual(float(d["deployed_correction_limit_rad"]), 6.0)
        self.assertIs(d["deployed_correction_limit_increased"], False)
        self.assertIs(d["source_replay_used"], False)
        self.assertIs(d["filter_changed"], False)
        self.assertIs(d["q8_word_promoted_here"], False)
        self.assertIs(d["whole_word_promoted_here"], False)
        self.assertIs(d["N_H_words_set_here"], False)
        if d["P5_SAMPLE1_SIGNED_CAYLEY_Q8_V14D"] == "PASS":
            self.assertEqual(V14D.validate(d), [])
            self.assertLess(d["max_post_sample1_cayley_norm_upper"], 8.0)
        else:
            self.assertTrue(V14D.validate(d) or d["first_unclosed_q8_cell"] is not None)


if __name__ == "__main__":
    unittest.main()
