from __future__ import annotations
import math
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_sample1_signed_cayley_q8_v15 as V15


class Sample1SignedCayleyQ8V15Tests(unittest.TestCase):
    def test_principal_axis_angle_bound_tracks_winding(self):
        self.assertLessEqual(V15._principal_axis_angle_upper(0.6, 1.7), 1.700000000000001)
        self.assertGreaterEqual(V15._principal_axis_angle_upper(3.0, 3.4), math.pi)
        self.assertLess(V15._principal_axis_angle_upper(5.0, 5.5), 1.29)
        self.assertLess(V15._principal_axis_angle_upper(6.2, 6.4), 0.117)
        self.assertLess(V15._principal_axis_angle_upper(7.0, 8.0), 1.718)
        self.assertIsNone(V15._principal_axis_angle_upper(0.0, 0.005))

    def test_completed_416_first_v14d_witness_closes_by_geodesic_metric(self):
        # Exact values emitted by the merged #416 fine V14D job.  The old
        # componentwise product bound gave q+=8.404756..., just outside q<8.
        geo = V15._geodesic_q_and_scalar_lower(
            0.6593778441001633,
            0.641119274810348,
            1.7248624819426428,
        )
        self.assertIsNotNone(geo)
        q_plus, w_lower = geo
        self.assertLess(q_plus, 8.0)
        self.assertGreater(w_lower, 0.0)
        self.assertLess(q_plus, 4.87)

    def test_completed_416_worst_high_angle_witness_still_requires_direction(self):
        # This cell can span a principal correction angle near pi, so V15 must
        # not pretend the metric triangle inequality alone closes it.
        geo = V15._geodesic_q_and_scalar_lower(
            1.9302592522019741,
            3.2214254998487837,
            5.269730943810957,
        )
        self.assertIsNone(geo)

    def test_q8_target_angle_is_stricter_than_antipode_only(self):
        target_angle = 2.0 * math.atan(V15.Q_TARGET / 2.0)
        self.assertLess(target_angle, math.pi)
        geo = V15._geodesic_q_and_scalar_lower(0.5, 0.5, 1.0)
        self.assertIsNotNone(geo)
        self.assertLess(geo[0], V15.Q_TARGET)


if __name__ == "__main__":
    unittest.main()
