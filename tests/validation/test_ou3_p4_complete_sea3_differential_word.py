#!/usr/bin/env python3
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
STABILITY = ROOT / "tools" / "stability"
if str(STABILITY) not in sys.path:
    sys.path.insert(0, str(STABILITY))

from ou3_interval import Interval, matrix_identity
import ou3_p4_complete_sea3_differential_events as EVENTS
import ou3_p4_complete_sea3_differential_word as WORD
import ou3_sea3_full_normal_live_word as LITERAL


def point_covariance(n):
    P = [[Interval.point(0.0) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        P[i][i] = Interval.point(1.0 + 0.02 * i)
    return P


def zero_state(n):
    return [Interval.point(0.0) for _ in range(n)]


class CompleteSea3DifferentialWordTests(unittest.TestCase):
    def test_same_history_source_token_is_mandatory(self):
        w = WORD.initialize("H", "sea3-history-A")
        with self.assertRaisesRegex(ValueError, "same complete SEA3 source"):
            WORD.apply_event(w, WORD.DifferentialEvent(
                "prediction", matrix_identity(18), "sea3-history-B"
            ))

    def test_due_s_event_requires_same_cell_and_actual_rs_provenance(self):
        w = WORD.initialize("H", "same")
        with self.assertRaisesRegex(ValueError, "same P/H/R cell"):
            WORD.apply_event(w, WORD.DifferentialEvent(
                "S_zero", matrix_identity(18), "same",
                R_provenance=EVENTS.ACTUAL_RS_PROVENANCE,
                same_P_H_R_cell=False,
            ))
        with self.assertRaisesRegex(ValueError, "actual applied SpectralMSE R_S"):
            WORD.apply_event(w, WORD.DifferentialEvent(
                "S_zero", matrix_identity(18), "same",
                R_provenance="TARGET_RS", same_P_H_R_cell=True,
            ))

        R = LITERAL.R_S_zero([
            Interval.point(0.72), Interval.point(0.72), Interval.point(1.0)
        ])
        source = EVENTS.source_joseph_event(
            "H", zero_state(18), point_covariance(18), R, "S_zero",
            R_provenance=EVENTS.ACTUAL_RS_PROVENANCE,
        )
        WORD.apply_event(w, WORD.event_from_source_joseph("S_zero", "same", source))
        self.assertEqual(w.S_updates, 1)

    def test_covariance_floor_has_identity_state_jacobian(self):
        w = WORD.initialize("A", "same")
        bad = matrix_identity(21)
        bad[15][15] = Interval.point(1.001)
        with self.assertRaisesRegex(ValueError, "identity state Jacobian"):
            WORD.apply_event(w, WORD.DifferentialEvent("aw_floor", bad, "same"))
        WORD.apply_event(w, WORD.DifferentialEvent("aw_floor", matrix_identity(21), "same"))
        self.assertEqual(w.floors, 1)

    def test_unique_rectangular_h_to_a_lift_composes_full_state(self):
        w = WORD.initialize("H", "same")
        WORD.apply_event(w, WORD.DifferentialEvent("prediction", matrix_identity(18), "same"))
        L = [[Interval.point(0.0) for _ in range(18)] for _ in range(21)]
        for i in range(18):
            L[i][i] = Interval.point(1.0)
        WORD.apply_event(w, WORD.DifferentialEvent("H_to_A", L, "same"))
        self.assertEqual((len(w.J_word), len(w.J_word[0])), (21, 18))
        self.assertEqual(w.current_dim, 21)
        self.assertEqual(w.hybrid_lifts, 1)
        with self.assertRaisesRegex(ValueError, "unique 21x18"):
            WORD.apply_event(w, WORD.DifferentialEvent("H_to_A", L, "same"))

    def test_full_matrix_cocycle_feeds_ldlt_gate(self):
        w = WORD.initialize("A", "same")
        half = Interval.point(0.5)
        J = [[Interval.point(0.0) for _ in range(21)] for _ in range(21)]
        for i in range(21):
            J[i][i] = half
        WORD.apply_event(w, WORD.DifferentialEvent("prediction", J, "same"))

        f = [Interval.point(0.0), Interval.point(0.0), Interval.point(-9.80665)]
        Rhat = matrix_identity(3)
        source = EVENTS.source_joseph_event(
            "A", zero_state(21), point_covariance(21), LITERAL.diagonal_R([0.2, 0.2, 0.2]),
            "accelerometer", f_hat=f, R_hat=Rhat,
        )
        # For this algebra-only LDLT smoke test use identity event after proving
        # theorem measurement events themselves require same-cell provenance.
        WORD.apply_event(w, WORD.DifferentialEvent(
            "accelerometer", matrix_identity(21), "same", same_P_H_R_cell=True
        ))
        self.assertTrue(source["same_P_H_R_cell"])
        ok, pivots = WORD.certify_word(
            w,
            matrix_identity(21), matrix_identity(21),
            matrix_identity(21), matrix_identity(21),
            0.5,
        )
        self.assertTrue(ok)
        self.assertEqual(len(pivots), 21)

    def test_event_from_source_joseph_rejects_unbound_event(self):
        with self.assertRaisesRegex(ValueError, "same P/H/R cell"):
            WORD.event_from_source_joseph(
                "accelerometer", "same", {"J_state": matrix_identity(18)}
            )

    def test_status_is_complete_sea3_and_nonpromoting(self):
        d = WORD.build()
        self.assertEqual(WORD.validate(d), [])
        self.assertEqual(d["canonical_source"], "COMPLETE_SEA3_NORMAL_LIVE_WORD")
        self.assertEqual(d["P3_delta_consumed"], 1.0e-18)
        self.assertTrue(d["same_complete_SEA3_source_token_required_for_every_event"])
        self.assertTrue(d["same_P_H_R_cell_required_for_every_Joseph_event"])
        self.assertFalse(d["independent_K_input_allowed_for_theorem"])
        self.assertTrue(d["actual_applied_per_axis_RS_required_on_every_S_event"])
        self.assertEqual(d["actual_RS_provenance_token"], EVENTS.ACTUAL_RS_PROVENANCE)
        self.assertTrue(d["H_to_A_unique_rectangular_event_required"])
        self.assertFalse(d["A21_bias_projection_hybrid_closed"])
        self.assertFalse(d["packetwise_norm_sum_used"])
        self.assertFalse(d["packet_count_multiplier_used"])
        self.assertFalse(d["independent_RS_schedule_used"])
        self.assertFalse(d["state_elimination_used"])
        self.assertFalse(d["source_uniform_complete_word_Jacobian_enclosed"])
        self.assertFalse(d["P4_promoted_here"])


if __name__ == "__main__":
    unittest.main()
