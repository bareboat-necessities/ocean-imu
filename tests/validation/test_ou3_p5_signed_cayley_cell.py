import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

from ou3_interval import Interval
import ou3_p5_signed_cayley_cell as C


class Ou3P5SignedCayleyCellTests(unittest.TestCase):
    def test_audit_primitive_passes_without_promoting_word(self):
        d = C.build()
        self.assertEqual(C.validate(d), [])
        self.assertEqual(d["P5_SIGNED_CAYLEY_CELL_PRIMITIVE"], "PASS")
        self.assertTrue(d["signed_a_dot_c_retained"])
        self.assertFalse(d["independent_abs_a_abs_c_denominator_used"])
        self.assertFalse(d["complete_word_promoted_here"])
        self.assertFalse(d["filter_changed"])

    def test_corrective_signed_cell_improves_exact_denominator(self):
        c = [Interval.outward_bounds(0.55, 0.65), Interval.outward_bounds(-0.02, 0.02), Interval.outward_bounds(-0.02, 0.02)]
        d = [Interval.outward_bounds(-0.45, -0.35), Interval.outward_bounds(-0.005, 0.005), Interval.outward_bounds(-0.005, 0.005)]
        row = C.compose_cell(c, d)
        self.assertLess(row["a_dot_c"].hi, 0.0)
        self.assertGreater(row["denominator"].lo, 1.0)
        self.assertGreater(row["correction_scale"].lo, 0.0)
        self.assertGreater(row["c_plus_norm_upper"], 0.0)

    def test_polynomial_and_axis_angle_scale_are_both_positive(self):
        small = C.correction_cayley_scale_interval(0.005)
        large = C.correction_cayley_scale_interval(1.5)
        self.assertGreater(small.lo, 0.0)
        self.assertGreater(large.lo, 0.0)
        self.assertGreater(large.hi, small.lo)

    def test_cell_refuses_antipodal_denominator_crossing(self):
        c = [Interval.outward_bounds(3.9, 4.1), Interval.point(0.0), Interval.point(0.0)]
        d = [Interval.outward_bounds(0.95, 1.05), Interval.point(0.0), Interval.point(0.0)]
        with self.assertRaises(RuntimeError):
            C.compose_cell(c, d)


if __name__ == "__main__":
    unittest.main()
