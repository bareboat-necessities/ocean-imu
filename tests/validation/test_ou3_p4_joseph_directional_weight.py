#!/usr/bin/env python3
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from ou3_interval import Interval
import ou3_p4_joseph_directional_weight as J


class JosephDirectionalWeightTests(unittest.TestCase):
    def I(self, x):
        return Interval.point(float(x))

    def test_scalar_fixture_matches_exact_attenuation(self):
        # P=2, H=1, R=3 => S=5 and H' S^-1 H = (3/5) H' R^-1 H.
        P = [[self.I(2.0)]]
        H = [[self.I(1.0)]]
        R = [[self.I(3.0)]]
        d = J.same_cell_bridge(P, H, R)
        c = d["joseph_vs_R_inverse_directional_attenuation_lower"]
        self.assertGreater(c, 0.59)
        self.assertLessEqual(c, 0.6)
        self.assertTrue(d["rank_and_nullspaces_preserved"])
        self.assertFalse(d["K_interval_matrix_materialized"])
        self.assertFalse(d["condition_number_conversion_used"])

    def test_rank_deficient_direction_is_not_scalarized(self):
        # A zero second column remains an exact null direction.  The bridge is a
        # scalar on H'R^-1H, not a replacement by c*I.
        P = [
            [self.I(1.0), self.I(0.0)],
            [self.I(0.0), self.I(1.0)],
        ]
        H = [[self.I(1.0), self.I(0.0)]]
        R = [[self.I(1.0)]]
        d = J.same_cell_bridge(P, H, R)
        self.assertGreater(d["joseph_vs_R_inverse_directional_attenuation_lower"], 0.0)
        self.assertIn("H^T S^-1 H", d["directional_inequality"])
        self.assertTrue(d["rank_and_nullspaces_preserved"])

    def test_source_bound_initial_HA_bridge_is_positive_and_fail_closed(self):
        d = J.build()
        self.assertEqual(J.validate(d), [])
        self.assertTrue(d["P4_JOSEPH_DIRECTIONAL_WEIGHT_PRIMITIVE_ESTABLISHED"])
        self.assertFalse(d["P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED"])
        self.assertFalse(d["P5_FINITE_CAPTURE_ESTABLISHED"])
        self.assertFalse(d["per_packet_scalar_full_state_margin_claimed"])
        self.assertFalse(d["interval_K_materialized"])
        self.assertFalse(d["covariance_condition_number_conversion_used"])
        for mode, n in (("H", 18), ("A", 21)):
            row = d["modes"][mode]
            self.assertEqual(row["dimension"], n)
            self.assertGreater(row["minimum_initial_directional_attenuation_lower"], 0.0)
            self.assertFalse(row["complete_word_prefix_covariances_propagated_here"])
            self.assertFalse(row["complete_word_directional_credit_accumulated_here"])
            self.assertFalse(row["P4_PROMOTED"])
            for op in ("accelerometer", "magnetometer", "S_zero"):
                c = row[op]["joseph_vs_R_inverse_directional_attenuation_lower"]
                self.assertTrue(math.isfinite(c))
                self.assertGreater(c, 0.0)
                self.assertLessEqual(c, 1.00000000000001)

    def test_diagonal_R_is_required(self):
        P = [
            [self.I(1.0), self.I(0.0)],
            [self.I(0.0), self.I(1.0)],
        ]
        H = [
            [self.I(1.0), self.I(0.0)],
            [self.I(0.0), self.I(1.0)],
        ]
        R = [
            [self.I(1.0), self.I(0.1)],
            [self.I(0.1), self.I(1.0)],
        ]
        with self.assertRaises(RuntimeError):
            J.same_cell_bridge(P, H, R)


if __name__ == "__main__":
    unittest.main()
