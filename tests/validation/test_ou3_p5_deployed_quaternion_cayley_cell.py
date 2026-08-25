import math
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

from ou3_interval import Interval
import ou3_p5_deployed_quaternion_cayley_cell as Q


class Ou3P5DeployedQuaternionCayleyCellTests(unittest.TestCase):
    def test_source_bound_primitive_passes_without_promoting_word(self):
        d = Q.build()
        self.assertEqual(Q.validate(d), [])
        self.assertEqual(d["P5_DEPLOYED_QUATERNION_CAYLEY_CELL_PRIMITIVE"], "PASS")
        self.assertTrue(d["only_resulting_error_antipode_is_gate"])
        self.assertFalse(d["correction_cayley_singularity_at_pi_used"])
        self.assertFalse(d["filter_changed"])
        self.assertFalse(d["complete_word_promoted_here"])

    def test_correction_above_retired_three_rad_range_composes(self):
        row = Q.compose_cell(
            [Interval.outward_bounds(-0.62, -0.58), Interval.outward_bounds(-0.005, 0.005), Interval.outward_bounds(-0.005, 0.005)],
            [Interval.outward_bounds(3.18, 3.22), Interval.outward_bounds(-0.002, 0.002), Interval.outward_bounds(-0.002, 0.002)],
        )
        self.assertGreater(row["correction_norm_upper"], 3.0)
        self.assertFalse(row["correction_cayley_coordinate_formed"])
        self.assertFalse(row["product_scalar"].lo <= 0.0 <= row["product_scalar"].hi)
        self.assertGreater(row["c_plus_norm_upper"], 0.0)

    def test_genuine_resulting_antipode_is_refused(self):
        with self.assertRaises(RuntimeError):
            Q.compose_cell(
                [Interval.point(0.0), Interval.point(0.0), Interval.point(0.0)],
                [Interval.point(math.pi), Interval.point(0.0), Interval.point(0.0)],
            )

    def test_polynomial_branch_still_uses_shipping_homogeneous_quaternion(self):
        row = Q.compose_cell(
            [Interval.outward_bounds(0.19, 0.21), Interval.point(0.0), Interval.point(0.0)],
            [Interval.outward_bounds(-0.006, -0.004), Interval.point(0.0), Interval.point(0.0)],
        )
        self.assertLess(row["correction_norm_upper"], 1.0e-2)
        self.assertTrue(row["source_quaternion_normalization_cancels_exactly"])
        self.assertGreater(row["product_scalar"].lo, 0.0)


if __name__ == "__main__":
    unittest.main()
