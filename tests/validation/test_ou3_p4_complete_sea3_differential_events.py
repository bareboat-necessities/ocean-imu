#!/usr/bin/env python3
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
STABILITY = ROOT / "tools" / "stability"
if str(STABILITY) not in sys.path:
    sys.path.insert(0, str(STABILITY))

from ou3_interval import Interval, matrix_identity, matrix_mul, matrix_sub
import ou3_p4_complete_sea3_differential_events as EVENTS
import ou3_sea3_full_normal_live_word as WORD


def point_matrix(rows):
    return [[Interval.point(float(x)) for x in row] for row in rows]


def zero_state(n):
    return [Interval.point(0.0) for _ in range(n)]


def point_K(n):
    K = [[Interval.point(0.0) for _ in range(3)] for _ in range(n)]
    vals = [0.02, -0.01, 0.015]
    for i in range(min(n, 21)):
        K[i][i % 3] = Interval.point(vals[i % 3] / (1.0 + i))
    return K


def expected_I_minus_KH(n, K, H):
    return matrix_sub(matrix_identity(n), matrix_mul(K, H))


def assert_contains(test, actual, expected):
    test.assertEqual((len(actual), len(actual[0])), (len(expected), len(expected[0])))
    for i in range(len(expected)):
        for j in range(len(expected[0])):
            e = expected[i][j]
            test.assertLessEqual(actual[i][j].lo, e.lo)
            test.assertGreaterEqual(actual[i][j].hi, e.hi)


class CompleteSea3DifferentialEventTests(unittest.TestCase):
    def test_S_zero_zero_error_jacobian_matches_literal_tangent(self):
        for mode, n in (("H", 18), ("A", 21)):
            K = point_K(n)
            J = EVENTS.joseph_event_jacobian(
                mode, zero_state(n), K, "S_zero", actual_applied_RS=True
            )
            expected = expected_I_minus_KH(n, K, WORD.H_S_zero(mode))
            assert_contains(self, J, expected)

    def test_magnetometer_zero_error_jacobian_matches_literal_tangent(self):
        m = [Interval.point(20.0), Interval.point(-5.0), Interval.point(40.0)]
        for mode, n in (("H", 18), ("A", 21)):
            K = point_K(n)
            J = EVENTS.joseph_event_jacobian(
                mode, zero_state(n), K, "magnetometer", m_body=m
            )
            expected = expected_I_minus_KH(n, K, WORD.H_magnetometer(mode, m))
            assert_contains(self, J, expected)

    def test_accelerometer_zero_error_jacobian_matches_literal_tangent(self):
        f = [Interval.point(1.5), Interval.point(-0.7), Interval.point(-9.2)]
        R = point_matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        for mode, n in (("H", 18), ("A", 21)):
            K = point_K(n)
            J = EVENTS.joseph_event_jacobian(
                mode, zero_state(n), K, "accelerometer", f_hat=f, R_hat=R
            )
            expected = expected_I_minus_KH(n, K, WORD.H_accelerometer(mode, f, R))
            assert_contains(self, J, expected)

    def test_finite_cell_uses_outward_exact_cayley_not_linear_reset(self):
        n = 18
        state = zero_state(n)
        state[0] = Interval(-0.08, 0.08)
        state[1] = Interval(-0.05, 0.05)
        state[15] = Interval(-0.2, 0.2)
        K = point_K(n)
        f = [Interval.point(2.0), Interval.point(0.5), Interval.point(-9.0)]
        R = point_matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        J = EVENTS.joseph_event_jacobian(
            "H", state, K, "accelerometer", f_hat=f, R_hat=R
        )
        self.assertEqual((len(J), len(J[0])), (18, 18))
        self.assertTrue(any(x.lo != x.hi for row in J for x in row))

    def test_A21_projection_hybrid_fails_closed(self):
        with self.assertRaisesRegex(RuntimeError, "projection hybrid"):
            EVENTS.joseph_event_jacobian(
                "A", zero_state(21), point_K(21), "S_zero",
                actual_applied_RS=True, bias_projection_inactive=False,
            )

    def test_status_keeps_projection_and_source_obligations_open(self):
        d = EVENTS.build()
        self.assertEqual(EVENTS.validate(d), [])
        self.assertEqual(d["canonical_source"], "COMPLETE_SEA3_NORMAL_LIVE_WORD")
        self.assertTrue(d["actual_applied_RS_required_for_S_event"])
        self.assertTrue(d["exact_Cayley_attitude_state"])
        self.assertGreater(d["A21_bias_projection_nominal_margin_mps2"], 0.0)
        self.assertFalse(d["A21_bias_projection_inactive_source_uniformly_proved_here"])
        self.assertFalse(d["projection_hybrid_silently_ignored"])
        self.assertFalse(d["source_uniform_finite_angle_event_Jacobians_closed"])
        self.assertFalse(d["P4_promoted_here"])


if __name__ == "__main__":
    unittest.main()
