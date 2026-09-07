#!/usr/bin/env python3
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
STABILITY = ROOT / "tools" / "stability"
if str(STABILITY) not in sys.path:
    sys.path.insert(0, str(STABILITY))

from ou3_interval import Interval, matrix_identity
import ou3_p4_complete_sea3_differential_word as WORD


class CompleteSea3DifferentialWordTests(unittest.TestCase):
    def test_same_history_source_token_is_mandatory(self):
        w = WORD.initialize("H", "sea3-history-A")
        with self.assertRaisesRegex(ValueError, "same complete SEA3 source"):
            WORD.apply_event(w, WORD.DifferentialEvent(
                "prediction", matrix_identity(18), "sea3-history-B"
            ))

    def test_due_s_event_requires_actual_applied_rs_provenance(self):
        w = WORD.initialize("H", "same")
        with self.assertRaisesRegex(ValueError, "actual-applied R_S"):
            WORD.apply_event(w, WORD.DifferentialEvent(
                "S_zero", matrix_identity(18), "same", actual_applied_RS=False
            ))
        WORD.apply_event(w, WORD.DifferentialEvent(
            "S_zero", matrix_identity(18), "same", actual_applied_RS=True
        ))
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
        WORD.apply_event(w, WORD.DifferentialEvent("accelerometer", matrix_identity(21), "same"))
        ok, pivots = WORD.certify_word(
            w,
            matrix_identity(21), matrix_identity(21),
            matrix_identity(21), matrix_identity(21),
            0.5,
        )
        self.assertTrue(ok)
        self.assertEqual(len(pivots), 21)

    def test_status_is_complete_sea3_and_nonpromoting(self):
        d = WORD.build()
        self.assertEqual(WORD.validate(d), [])
        self.assertEqual(d["canonical_source"], "COMPLETE_SEA3_NORMAL_LIVE_WORD")
        self.assertEqual(d["P3_delta_consumed"], 1.0e-18)
        self.assertTrue(d["same_complete_SEA3_source_token_required_for_every_event"])
        self.assertTrue(d["actual_applied_per_axis_RS_required_on_every_S_event"])
        self.assertTrue(d["H_to_A_unique_rectangular_event_required"])
        self.assertFalse(d["packetwise_norm_sum_used"])
        self.assertFalse(d["packet_count_multiplier_used"])
        self.assertFalse(d["independent_RS_schedule_used"])
        self.assertFalse(d["state_elimination_used"])
        self.assertFalse(d["source_uniform_complete_word_Jacobian_enclosed"])
        self.assertFalse(d["P4_promoted_here"])


if __name__ == "__main__":
    unittest.main()
