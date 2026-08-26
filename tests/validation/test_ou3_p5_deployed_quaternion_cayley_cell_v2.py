import math
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

from ou3_interval import Interval
import ou3_p5_deployed_quaternion_cayley_cell_v2 as Q


class Ou3P5DeployedQuaternionCayleyCellV2Tests(unittest.TestCase):
    def test_extended_primitive_passes_without_promoting_word(self):
        d = Q.build()
        self.assertEqual(Q.validate(d), [])
        self.assertEqual(d["P5_DEPLOYED_QUATERNION_CAYLEY_CELL_V2_PRIMITIVE"], "PASS")
        self.assertGreaterEqual(d["maximum_validated_correction_norm_rad"], 9.0)
        self.assertTrue(d["radial_subdivision_required_above_6_rad"])
        self.assertTrue(d["nonmonotone_half_angle_trig_enclosed_without_monotonicity"])
        self.assertFalse(d["filter_changed"])
        self.assertFalse(d["complete_word_promoted_here"])

    def test_cell_above_two_pi_composes(self):
        row = Q.compose_cell(
            [Interval.point(0.0), Interval.point(0.0), Interval.point(0.0)],
            [Interval.outward_bounds(8.0, 8.2), Interval.point(0.0), Interval.point(0.0)],
            d_norm_lower=math.nextafter(8.0, -math.inf),
            d_norm_upper=math.nextafter(8.2, math.inf),
        )
        self.assertGreater(row["correction_norm_lower"], 2.0 * math.pi)
        self.assertGreater(row["correction_norm_upper"], 6.0)
        self.assertEqual(row["quaternion_enclosure_backend"], "RADIAL_SUBCELL_NONMONOTONE_TRIG")
        self.assertFalse(row["product_scalar"].lo <= 0.0 <= row["product_scalar"].hi)
        self.assertTrue(math.isfinite(row["c_plus_norm_upper"]))

    def test_unsubdivided_winding_box_is_refused(self):
        with self.assertRaises(ValueError):
            Q.compose_cell(
                [Interval.point(0.0), Interval.point(0.0), Interval.point(0.0)],
                [Interval.outward_bounds(-8.2, 8.2), Interval.point(0.0), Interval.point(0.0)],
                d_norm_upper=math.nextafter(8.2, math.inf),
            )

    def test_old_range_still_delegates_to_v1(self):
        row = Q.compose_cell(
            [Interval.outward_bounds(-0.62, -0.58), Interval.point(0.0), Interval.point(0.0)],
            [Interval.outward_bounds(3.18, 3.22), Interval.point(0.0), Interval.point(0.0)],
        )
        self.assertEqual(row["quaternion_enclosure_backend"], "V1_MONOTONE_THROUGH_6_RAD")
        self.assertFalse(row["product_scalar"].lo <= 0.0 <= row["product_scalar"].hi)


if __name__ == "__main__":
    unittest.main()
