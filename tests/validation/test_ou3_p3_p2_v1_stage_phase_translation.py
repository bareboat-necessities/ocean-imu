#!/usr/bin/env python3
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from ou3_interval import symmetric_positive_definite_ldlt
import ou3_p3_p2_v1_stage_phase_translation as P
import ou3_source_reachable_matrix_p3 as BASE


class P2V1StagePhaseTranslationTests(unittest.TestCase):
    def test_finite_phase_cover_and_frozen_tail_cut_are_exact(self):
        self.assertEqual(P.PHASES, tuple(range(26)))
        self.assertEqual(P.FROZEN_FRESH_HISTORY_SAMPLES, 13)
        self.assertEqual(P.FROZEN_TRANSIENT_SAMPLES, tuple(range(1, 13)))
        self.assertEqual(BASE.MIN_USEFUL_DELTA, 1.0e-18)

    def test_identity_floor_is_strict_spd(self):
        M = P._identity_floor(0.25)
        ok, _ = symmetric_positive_definite_ldlt(M)
        self.assertTrue(ok)
        self.assertEqual(M[0][0].lo, 0.25)
        self.assertEqual(M[0][1].lo, 0.0)

    def test_certified_delta_matches_simple_diagonal_case(self):
        M = P._identity_floor(1.0)
        delta = P._certified_delta(M, 1.0, [2.0, 4.0, 5.0, 10.0])
        self.assertGreater(delta, 0.099999999999)
        self.assertLess(delta, 0.100000000001)

    def test_gate_fails_above_simple_diagonal_margin(self):
        M = P._identity_floor(1.0)
        A = P._physical_gate_matrix(M, 1.0, [2.0, 4.0, 5.0, 10.0], 0.11)
        self.assertFalse(symmetric_positive_definite_ldlt(A)[0])

    def test_fixed_physical_scaling_penalizes_S_by_h6(self):
        M = P._identity_floor(1.0)
        h = 0.005
        A = P._physical_gate_matrix(M, h, [1.0, 1.0, 1.0, 1.0], 0.0)
        self.assertGreater(A[0][0].lo, 0.0)
        self.assertGreater(A[1][1].lo, 0.0)
        self.assertGreater(A[2][2].lo, 0.0)
        self.assertGreater(A[3][3].lo, 0.0)
        self.assertLess(A[2][2].hi, A[1][1].hi)
        self.assertLess(A[1][1].hi, A[0][0].hi)

    def test_thirteen_sample_tail_argument_is_a_fresh_window(self):
        # The frozen-clock reduction must not accidentally depend on the 26
        # finite-stage phase bound.  After 13 held samples, the last 13 samples
        # form a complete fresh constant-source segment on their own.
        self.assertLess(P.FROZEN_FRESH_HISTORY_SAMPLES, max(P.PHASES) + 1)
        self.assertEqual(max(P.FROZEN_TRANSIENT_SAMPLES) + 1, P.FROZEN_FRESH_HISTORY_SAMPLES)


if __name__ == "__main__":
    unittest.main()
