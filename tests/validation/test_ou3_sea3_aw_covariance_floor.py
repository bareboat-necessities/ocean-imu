#!/usr/bin/env python3
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_sea3_aw_covariance_floor as FLOOR
from ou3_interval import Interval


class Sea3AwCovarianceFloorTest(unittest.TestCase):
    def diag(self, values):
        z = FLOOR.I(0.0)
        M = [[z for _ in range(3)] for _ in range(3)]
        for i, value in enumerate(values):
            M[i][i] = value if isinstance(value, Interval) else FLOOR.I(value)
        return M

    def test_status_and_shipping_parity(self):
        d = FLOOR.build()
        self.assertEqual(FLOOR.validate(d), [])
        self.assertTrue(d["shipping_source_parity_pass"])
        self.assertTrue(d["floor_request_is_source_event"])
        self.assertFalse(d["floor_increment_is_source_coordinate"])
        self.assertFalse(d["P3_promoted"])

    def test_strict_positive_argument_is_exact(self):
        target = self.diag([1.0, 1.0, 1.0])
        current = self.diag([0.25, 0.5, 0.75])
        delta, case = FLOOR.positive_part_enclosure(target, current)
        self.assertEqual(case, "strict_positive_definite_exact")
        for i, expected in enumerate((0.75, 0.5, 0.25)):
            self.assertTrue(delta[i][i].contains(expected))
        for i in range(3):
            for j in range(3):
                if i != j:
                    self.assertTrue(delta[i][j].contains(0.0))

    def test_strict_negative_argument_is_exact_zero(self):
        target = self.diag([1.0, 1.0, 1.0])
        current = self.diag([1.5, 2.0, 3.0])
        delta, case = FLOOR.positive_part_enclosure(target, current)
        self.assertEqual(case, "strict_negative_definite_exact_zero")
        self.assertTrue(all(x.lo == 0.0 and x.hi == 0.0 for row in delta for x in row))

    def test_mixed_argument_uses_rigorous_outer_box(self):
        target = self.diag([1.0, 1.0, 1.0])
        current = self.diag([0.5, 1.5, Interval(0.8, 1.2)])
        delta, case = FLOOR.positive_part_enclosure(target, current)
        self.assertEqual(case, "mixed_spectrum_frobenius_outer")
        # Point D=diag(.5,-.5,0) has Pi_+(D)=diag(.5,0,0).
        self.assertTrue(delta[0][0].contains(0.5))
        self.assertTrue(delta[1][1].contains(0.0))
        self.assertTrue(delta[2][2].contains(0.0))
        self.assertTrue(all(x.lo <= x.hi for row in delta for x in row))

    def test_aw_block_is_mode_independent_offset(self):
        z = FLOOR.I(0.0)
        for n in (18, 21):
            P = [[z for _ in range(n)] for _ in range(n)]
            for i, value in enumerate((1.0, 2.0, 3.0)):
                P[15 + i][15 + i] = FLOOR.I(value)
            b = FLOOR.aw_block(P)
            self.assertTrue(b[0][0].contains(1.0))
            self.assertTrue(b[1][1].contains(2.0))
            self.assertTrue(b[2][2].contains(3.0))


if __name__ == "__main__":
    unittest.main()
