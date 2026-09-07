#!/usr/bin/env python3
import sys
from pathlib import Path
import unittest

REPO = Path(__file__).resolve().parents[2]
STABILITY = REPO / "tools" / "stability"
if str(STABILITY) not in sys.path:
    sys.path.insert(0, str(STABILITY))

from ou3_interval import Interval, matrix_identity
import ou3_p4_complete_sea3_phi_differential_metric as DIFF


class PhiDifferentialMetricTests(unittest.TestCase):
    def test_zero_error_phi_jacobian_contains_identity_in_H18_and_A21(self):
        z = [Interval.point(0.0) for _ in range(3)]
        f = [Interval.point(1.0), Interval.point(-2.0), Interval.point(3.0)]
        R = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        for n in (18, 21):
            T = DIFF.full_phi_jacobian_interval(n, z, z, f, R)
            self.assertEqual((len(T), len(T[0])), (n, n))
            for i in range(n):
                for j in range(n):
                    target = 1.0 if i == j else 0.0
                    self.assertLessEqual(T[i][j].lo, target)
                    self.assertGreaterEqual(T[i][j].hi, target)

    def test_finite_cell_keeps_aw_diagonal_rotation_bounded(self):
        c = [Interval(-0.10, 0.10), Interval(-0.08, 0.08), Interval(-0.06, 0.06)]
        da = [Interval(-0.2, 0.2) for _ in range(3)]
        f = [Interval(-5.0, 5.0), Interval(-4.0, 4.0), Interval(-6.0, 6.0)]
        R = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        T = DIFF.full_phi_jacobian_interval(21, c, da, f, R)
        for i in range(15, 18):
            for j in range(15, 18):
                self.assertTrue(T[i][j].lo <= T[i][j].hi)
                self.assertTrue(abs(T[i][j].lo) < 2.0)
                self.assertTrue(abs(T[i][j].hi) < 2.0)

    def test_interval_ldlt_certifies_full_differential_matrix_not_scalar_ratio(self):
        I = matrix_identity(2)
        half = Interval.point(0.5)
        Jword = [[half if i == j else Interval.point(0.0) for j in range(2)] for i in range(2)]
        ok, pivots = DIFF.certify_differential_contraction(Jword, I, I, I, I, 0.5)
        self.assertTrue(ok)
        self.assertEqual(len(pivots), 2)
        self.assertTrue(all(p.lo > 0.0 for p in pivots))
        bad, _ = DIFF.certify_differential_contraction(Jword, I, I, I, I, 0.2)
        self.assertFalse(bad)

    def test_status_is_complete_sea3_full_rank_and_nonpromoting(self):
        d = DIFF.build()
        self.assertEqual(DIFF.validate(d), [])
        self.assertEqual(d["canonical_source"], "COMPLETE_SEA3_NORMAL_LIVE_WORD")
        self.assertEqual(d["P3_delta_consumed"], 1.0e-18)
        self.assertTrue(d["H18_full_rank"])
        self.assertTrue(d["A21_full_rank"])
        self.assertEqual(d["Phi_jacobian_determinant_exact"], 1.0)
        self.assertTrue(d["all_due_S_updates_with_actual_applied_RS_required"])
        self.assertFalse(d["finite_Phi_storage_used_as_Lyapunov_function"])
        self.assertFalse(d["state_elimination_used"])
        self.assertFalse(d["packet_count_remainder_budget_used"])
        self.assertFalse(d["correction_radius_claim_used"])
        self.assertFalse(d["inverse_metric_floor_claim_used"])
        self.assertFalse(d["P4_promoted_here"])


if __name__ == "__main__":
    unittest.main()
