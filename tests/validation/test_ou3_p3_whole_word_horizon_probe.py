#!/usr/bin/env python3
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import ou3_p3_whole_word_horizon_probe as P
import ou3_p3_p2_v1_stage_phase_translation as STAGE


class WholeWordHorizonProbeTests(unittest.TestCase):
    def test_635_sample_word_guarantees_23_complete_segments(self):
        self.assertEqual(P.guaranteed_complete_segments(635), 23)
        # 22 maximum-length segments plus two maximum partial tails cannot
        # cover the target: 22*26 + 2*25 = 622 < 635.
        self.assertLess(22 * 26 + 2 * P.PHASE_MAX, 635)
        self.assertGreaterEqual(23 * 26 + 2 * P.PHASE_MAX, 635)

    def test_metric_projection_maps_back_below_certified_physical_diagonal(self):
        h = 0.005
        upper = [11.0, 17.0, 23.0, 5.0]
        delta = 2.0e-7
        L = P._metric_lower(delta, h, upper)
        physical = (h, h*h, h*h*h, 1.0)
        for i in range(4):
            self.assertGreater(L[i][i].lo, 0.0)
            mapped = physical[i] * physical[i] * L[i][i].lo
            self.assertLessEqual(mapped, delta * upper[i])
            for j in range(4):
                if i != j:
                    self.assertEqual(L[i][j].lo, 0.0)
                    self.assertEqual(L[i][j].hi, 0.0)

    def test_metric_projection_is_accepted_by_the_canonical_gate_at_smaller_delta(self):
        h = 0.005
        upper = [11.0, 17.0, 23.0, 5.0]
        delta = 2.0e-7
        L = P._metric_lower(delta, h, upper)
        # Outward rounding in _metric_lower deliberately shaves the exact
        # boundary, so test a slightly smaller value rather than equality.
        certified = STAGE._certified_delta(L, h, upper)
        self.assertTrue(math.isfinite(certified))
        self.assertGreater(certified, 0.99 * delta)
        self.assertLessEqual(certified, delta)

    def test_projection_constants_match_frozen_p2_partition(self):
        self.assertEqual(P.TAU_CELLS, 10)
        self.assertEqual(P.SIGMA_CELLS, 8)
        self.assertEqual(P.RS_CELLS, 10)
        self.assertEqual(P.TAU_STRIDE, 80)
        self.assertEqual([i * P.TAU_STRIDE for i in range(P.TAU_CELLS)],
                         [0, 80, 160, 240, 320, 400, 480, 560, 640, 720])


if __name__ == "__main__":
    unittest.main()
