#!/usr/bin/env python3
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
STABILITY = ROOT / "tools" / "stability"
if str(STABILITY) not in sys.path:
    sys.path.insert(0, str(STABILITY))

from ou3_interval import Interval, matrix_identity
import ou3_p4_complete_brmm_differential_prediction as PRED
import ou3_brmm_shipping_prediction_primitives as SHIPPING


def I(x): return Interval.point(float(x))

def zero_state(n): return [I(0.0) for _ in range(n)]


def expected_F(mode, omega, h, tau_aw, tau_ba):
    Faa, _ = SHIPPING.attitude_gyro_bias_F_Q(
        omega, h, [I(1e-6), I(1e-6), I(1e-6)], I(1e-10)
    )
    Fll, _ = SHIPPING.translation_F_Q(tau_aw, h, [I(1.0), I(1.0), I(1.0)])
    n = 18 if mode == "H" else 21
    F = matrix_identity(n)
    for i in range(6):
        for j in range(6):
            F[i][j] = Faa[i][j]
    for i in range(12):
        for j in range(12):
            F[6+i][6+j] = Fll[i][j]
    if mode == "A":
        phi, _ = SHIPPING.active_accel_bias_F_Q(h, tau_ba, I(1e-8))
        for i in range(3):
            F[18+i][18+i] = phi
    return F


def assert_contains(test, actual, expected, tol=2e-10):
    test.assertEqual((len(actual), len(actual[0])), (len(expected), len(expected[0])))
    for i in range(len(expected)):
        for j in range(len(expected[0])):
            e = expected[i][j]
            test.assertLessEqual(actual[i][j].lo, e.hi + tol, (i,j,actual[i][j],e))
            test.assertGreaterEqual(actual[i][j].hi, e.lo - tol, (i,j,actual[i][j],e))


class CompleteBrmmDifferentialPredictionTests(unittest.TestCase):
    def test_zero_error_jacobian_contains_literal_shipping_prediction_H18_A21(self):
        omega = [I(0.2), I(-0.1), I(0.05)]
        h = I(0.005)
        tau = I(1.7)
        tau_ba = I(60.0)
        for mode, n in (("H",18),("A",21)):
            event = PRED.prediction_event(
                mode, zero_state(n), omega, h, tau, tau_ba=tau_ba if mode == "A" else None
            )
            expected = expected_F(mode, omega, h, tau, tau_ba)
            assert_contains(self, event["J_state"], expected)
            self.assertTrue(event["same_source_omega_h_tau"])
            self.assertTrue(event["full_F_Eaw_translation_rows_retained"])
            self.assertEqual(len(event["state_out"]), n)
            self.assertTrue(all(x.contains(0.0) for x in event["state_out"]))

    def test_finite_cayley_and_gyro_bias_cell_returns_state_and_jacobian_from_same_map(self):
        state = zero_state(18)
        state[0] = Interval(-0.15,0.15)
        state[1] = Interval(-0.10,0.10)
        state[3] = Interval(-0.01,0.01)
        event = PRED.prediction_event(
            "H", state, [I(0.3),I(-0.2),I(0.1)], I(0.005), I(2.0)
        )
        J = event["J_state"]
        out = event["state_out"]
        self.assertEqual((len(J),len(J[0])),(18,18))
        self.assertEqual(len(out),18)
        self.assertTrue(any(x.lo != x.hi for row in J for x in row))
        self.assertTrue(any(x.lo != x.hi for x in out))

    def test_A21_requires_configured_bias_GM_factor(self):
        with self.assertRaisesRegex(ValueError, "tau_ba"):
            PRED.prediction_event(
                "A", zero_state(21), [I(0),I(0),I(0)], I(0.005), I(2.0)
            )

    def test_status_retains_complete_brmm_and_full_F_Eaw(self):
        d=PRED.build()
        self.assertEqual(PRED.validate(d),[])
        self.assertEqual(d["canonical_source"],"COMPLETE_BRMM_NORMAL_LIVE_WORD")
        self.assertTrue(d["same_complete_BRMM_omega_h_tau_required"])
        self.assertFalse(d["independent_F_input_allowed_for_theorem"])
        self.assertTrue(d["exact_relative_attitude_prediction_differentiated"])
        self.assertTrue(d["finite_physical_state_output_available"])
        self.assertTrue(d["finite_state_and_Jacobian_share_one_AD_composition"])
        self.assertTrue(d["full_F_Eaw_v_p_S_aw_rows_retained"])
        self.assertTrue(d["process_Q_remains_in_same_source_Riccati_metric"])
        self.assertFalse(d["packetwise_norm_sum_used"])
        self.assertFalse(d["source_uniform_prediction_Jacobian_closed"])
        self.assertFalse(d["P4_promoted_here"])


if __name__ == "__main__": unittest.main()
