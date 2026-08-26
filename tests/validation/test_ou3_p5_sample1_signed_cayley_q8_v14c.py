from __future__ import annotations
import math
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
from ou3_interval import Interval
import ou3_p5_sample1_signed_cayley_q8_v14c as V14C


class Sample1SignedCayleyQ8V14CTests(unittest.TestCase):
    def test_component_box_is_clipped_by_radial_upper(self):
        d = [Interval.outward_bounds(-0.7, 0.7),
             Interval.outward_bounds(-0.5, 0.6),
             Interval.outward_bounds(-0.9, 0.8)]
        c = V14C._clip_component_box(d, 0.01)
        self.assertIsNotNone(c)
        for x in c:
            self.assertLessEqual(x.abs_upper(), 0.01000000000000001)

    def test_series_branch_cannot_use_large_cartesian_components(self):
        d = [Interval.outward_bounds(-0.7, 0.7),
             Interval.outward_bounds(-0.5, 0.6),
             Interval.outward_bounds(-0.9, 0.8)]
        w, v, branches = V14C.branch_local_normalized_shipping_quaternion(
            d, radial_lower=0.0, radial_upper=0.005)
        self.assertEqual(branches, ["SERIES_NORMALIZED_RADIAL_CLIP"])
        self.assertGreater(w.lo, 0.99)
        for x in v:
            self.assertLess(x.abs_upper(), 0.01)

    def test_mixed_series_axis_quaternion_records_superseded_sinc_dependency_loss(self):
        """V14C is historical: broad axis division loses radius/sine dependency.

        The source quaternion is unit, but V14C independently divides an
        interval sine by an interval half-angle.  This synthetic cell therefore
        must remain an explicit fail-closed witness rather than being relabeled
        as a valid unit enclosure.  V14D replaces this operation by the validated
        radial sinc backend.
        """
        d = [Interval.outward_bounds(-0.35, 0.02),
             Interval.outward_bounds(-0.02, 0.35),
             Interval.outward_bounds(-0.42, 0.73)]
        w, v, branches = V14C.branch_local_normalized_shipping_quaternion(
            d, radial_lower=0.0, radial_upper=0.75)
        self.assertIn("SERIES_NORMALIZED_RADIAL_CLIP", branches)
        self.assertIn("AXIS_ANGLE_UNIT_RADIAL_CLIP", branches)
        self.assertTrue(all(math.isfinite(x.lo) and math.isfinite(x.hi)
                            for x in (w, *v)))
        self.assertGreater(max(x.abs_upper() for x in v), 1.0)

    def test_coarse_family_is_fail_closed_or_closed(self):
        d = V14C.build(source_pieces=4, source_cell_index=0,
                       p_pieces=2, tangent_pieces=2, axial_pieces=2,
                       residual_x_pieces=2, parallel_pieces=2)
        self.assertIn(d["P5_SAMPLE1_SIGNED_CAYLEY_Q8_V14C"],
                      ("PASS", "NOT_ESTABLISHED"))
        self.assertEqual(float(d["deployed_correction_limit_rad"]), 6.0)
        self.assertIs(d["deployed_correction_limit_increased"], False)
        self.assertIs(d["q8_word_promoted_here"], False)
        self.assertIs(d["whole_word_promoted_here"], False)
        self.assertIs(d["N_H_words_set_here"], False)
        if d["P5_SAMPLE1_SIGNED_CAYLEY_Q8_V14C"] == "PASS":
            self.assertEqual(V14C.validate(d), [])
            self.assertLess(d["max_post_sample1_cayley_norm_upper"], 8.0)
        else:
            vf = V14C.validate(d)
            self.assertTrue(vf or d["first_unclosed_q8_cell"] is not None)


if __name__ == "__main__":
    unittest.main()
