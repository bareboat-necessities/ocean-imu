from __future__ import annotations
import math
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
from ou3_interval import Interval
import ou3_p5_sample1_signed_cayley_q8_v14b as V14B


class Sample1SignedCayleyQ8V14BTests(unittest.TestCase):
    def test_two_component_interval_norm_is_outward(self):
        v = [Interval.outward_bounds(-3.0, 3.0),
             Interval.outward_bounds(-4.0, 4.0)]
        q = V14B.interval_euclidean_norm_upper(v)
        self.assertTrue(math.isfinite(q))
        self.assertGreaterEqual(q, 5.0)

    def test_three_component_parent_use_remains_supported(self):
        v = [Interval.outward_bounds(-1.0, 1.0),
             Interval.outward_bounds(-2.0, 2.0),
             Interval.outward_bounds(-2.0, 2.0)]
        self.assertGreaterEqual(V14B.interval_euclidean_norm_upper(v), 3.0)

    def test_empty_interval_vector_is_rejected(self):
        with self.assertRaises(ValueError):
            V14B.interval_euclidean_norm_upper([])

    def test_coarse_family_is_fail_closed_or_closed(self):
        d = V14B.build(source_pieces=4, source_cell_index=0,
                       p_pieces=2, tangent_pieces=2, axial_pieces=2,
                       residual_x_pieces=2, parallel_pieces=2)
        self.assertIn(d["P5_SAMPLE1_SIGNED_CAYLEY_Q8_V14B"],
                      ("PASS", "NOT_ESTABLISHED"))
        self.assertEqual(float(d["deployed_correction_limit_rad"]), 6.0)
        self.assertIs(d["deployed_correction_limit_increased"], False)
        self.assertIs(d["q8_word_promoted_here"], False)
        self.assertIs(d["whole_word_promoted_here"], False)
        self.assertIs(d["N_H_words_set_here"], False)
        if d["P5_SAMPLE1_SIGNED_CAYLEY_Q8_V14B"] == "PASS":
            self.assertEqual(V14B.validate(d), [])
            self.assertLess(d["max_post_sample1_cayley_norm_upper"], 8.0)
        else:
            vf = V14B.validate(d)
            self.assertTrue(vf or d["first_unclosed_q8_cell"] is not None)


if __name__ == "__main__":
    unittest.main()
