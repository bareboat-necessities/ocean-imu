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


def point_covariance(n):
    P = [[Interval.point(0.0) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        P[i][i] = Interval.point(1.0 + 0.03 * i)
    # Retain a few benign cross-covariances so K is genuinely full-state.
    for i, j, x in ((0, 12, 0.02), (1, 15, -0.015), (6, 12, 0.01)):
        if i < n and j < n:
            P[i][j] = P[j][i] = Interval.point(x)
    return P


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
    def test_S_zero_zero_error_jacobian_matches_same_cell_literal_tangent(self):
        for mode, n in (("H", 18), ("A", 21)):
            P = point_covariance(n)
            R = WORD.R_S_zero([Interval.point(0.72), Interval.point(0.72), Interval.point(1.0)])
            event = EVENTS.source_joseph_event(
                mode, zero_state(n), P, R, "S_zero",
                R_provenance=EVENTS.ACTUAL_RS_PROVENANCE,
            )
            expected = expected_I_minus_KH(n, event["K"], event["H"])
            assert_contains(self, event["J_state"], expected)
            self.assertTrue(event["same_P_H_R_cell"])
            self.assertEqual(event["R_provenance"], EVENTS.ACTUAL_RS_PROVENANCE)

    def test_S_zero_rejects_target_or_unprovenanced_R(self):
        P = point_covariance(18)
        R = WORD.R_S_zero([Interval.point(0.72), Interval.point(0.72), Interval.point(1.0)])
        with self.assertRaisesRegex(ValueError, "actual applied SpectralMSE R_S"):
            EVENTS.source_joseph_event("H", zero_state(18), P, R, "S_zero")
        with self.assertRaisesRegex(ValueError, "actual applied SpectralMSE R_S"):
            EVENTS.source_joseph_event(
                "H", zero_state(18), P, R, "S_zero", R_provenance="TARGET_RS"
            )

    def test_magnetometer_zero_error_jacobian_matches_same_cell_literal_tangent(self):
        m = [Interval.point(20.0), Interval.point(-5.0), Interval.point(40.0)]
        Rm = WORD.diagonal_R([0.3, 0.3, 0.3])
        for mode, n in (("H", 18), ("A", 21)):
            event = EVENTS.source_joseph_event(
                mode, zero_state(n), point_covariance(n), Rm,
                "magnetometer", m_body=m,
            )
            expected = expected_I_minus_KH(n, event["K"], event["H"])
            assert_contains(self, event["J_state"], expected)

    def test_accelerometer_zero_error_jacobian_matches_same_cell_literal_tangent(self):
        f = [Interval.point(1.5), Interval.point(-0.7), Interval.point(-9.2)]
        Rhat = point_matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        Ra = WORD.diagonal_R([0.2, 0.2, 0.2])
        for mode, n in (("H", 18), ("A", 21)):
            event = EVENTS.source_joseph_event(
                mode, zero_state(n), point_covariance(n), Ra,
                "accelerometer", f_hat=f, R_hat=Rhat,
            )
            expected = expected_I_minus_KH(n, event["K"], event["H"])
            assert_contains(self, event["J_state"], expected)

    def test_finite_cell_uses_outward_exact_cayley_not_linear_reset(self):
        n = 18
        state = zero_state(n)
        state[0] = Interval(-0.08, 0.08)
        state[1] = Interval(-0.05, 0.05)
        state[15] = Interval(-0.2, 0.2)
        f = [Interval.point(2.0), Interval.point(0.5), Interval.point(-9.0)]
        Rhat = point_matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        event = EVENTS.source_joseph_event(
            "H", state, point_covariance(n), WORD.diagonal_R([0.2, 0.2, 0.2]),
            "accelerometer", f_hat=f, R_hat=Rhat,
        )
        J = event["J_state"]
        self.assertEqual((len(J), len(J[0])), (18, 18))
        self.assertTrue(any(x.lo != x.hi for row in J for x in row))

    def test_A21_projection_hybrid_fails_closed(self):
        P = point_covariance(21)
        R = WORD.R_S_zero([Interval.point(0.72), Interval.point(0.72), Interval.point(1.0)])
        with self.assertRaisesRegex(RuntimeError, "projection hybrid"):
            EVENTS.source_joseph_event(
                "A", zero_state(21), P, R, "S_zero",
                R_provenance=EVENTS.ACTUAL_RS_PROVENANCE,
                bias_projection_inactive=False,
            )

    def test_status_keeps_projection_and_source_obligations_open(self):
        d = EVENTS.build()
        self.assertEqual(EVENTS.validate(d), [])
        self.assertEqual(d["canonical_source"], "COMPLETE_SEA3_NORMAL_LIVE_WORD")
        self.assertTrue(d["same_P_H_R_cell_derives_S_and_K"])
        self.assertFalse(d["independent_K_input_allowed_for_theorem"])
        self.assertTrue(d["actual_applied_RS_required_for_S_event"])
        self.assertEqual(d["actual_applied_RS_provenance_token"], EVENTS.ACTUAL_RS_PROVENANCE)
        self.assertTrue(d["exact_Cayley_attitude_state"])
        self.assertGreater(d["A21_bias_projection_nominal_margin_mps2"], 0.0)
        self.assertFalse(d["A21_bias_projection_inactive_source_uniformly_proved_here"])
        self.assertFalse(d["projection_hybrid_silently_ignored"])
        self.assertFalse(d["source_uniform_finite_angle_event_Jacobians_closed"])
        self.assertFalse(d["P4_promoted_here"])


if __name__ == "__main__":
    unittest.main()
