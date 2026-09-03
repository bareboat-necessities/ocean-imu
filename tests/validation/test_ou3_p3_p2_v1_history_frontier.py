#!/usr/bin/env python3
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import ou3_p3_p2_v1_history_frontier as H


class P2V1HistoryFrontierTests(unittest.TestCase):
    def test_componentwise_adverse_dominance(self):
        self.assertTrue(H.dominates((2, 3, 4, 5), (1, 3, 0, 5)))
        self.assertFalse(H.dominates((2, 2, 4, 5), (1, 3, 0, 5)))

    def test_pareto_insert_removes_only_dominated_labels(self):
        front = {(2, 1, 4, 0), (1, 3, 2, 2)}
        self.assertTrue(H.pareto_insert(front, (3, 2, 4, 1)))
        self.assertIn((3, 2, 4, 1), front)
        self.assertNotIn((2, 1, 4, 0), front)
        self.assertIn((1, 3, 2, 2), front)
        before = set(front)
        self.assertFalse(H.pareto_insert(front, (2, 1, 3, 1)))
        self.assertEqual(front, before)

    def test_max_update_is_continuation_monotone(self):
        a = (1, 3, 2, 0)
        b = (2, 3, 4, 1)
        self.assertTrue(H.dominates(b, a))
        node = (5, 0, 3, 7)
        aa = H.update_label(a, node)
        bb = H.update_label(b, node)
        self.assertTrue(H.dominates(bb, aa))
        self.assertEqual(aa, (5, 3, 3, 7))
        self.assertEqual(bb, (5, 3, 4, 7))

    def test_empty_label_becomes_complete_after_one_segment(self):
        self.assertEqual(H.update_label(H.EMPTY_LABEL, (0, 1, 2, 3)), (0, 1, 2, 3))

    def test_label_summary_preserves_no_cartesian_extrema_contract(self):
        fr = {
            "stats": {
                "tables": [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0], [7.0, 8.0]],
                "cadence_lower_global_safe": 0.5,
                "node_ranks": [(0, 0, 0, 0), (1, 1, 1, 1)],
            },
            "target": {
                "history_duration_lower_s": 4.0,
                "terminal_history_duration_upper_s": 4.2,
            },
        }
        s = H.label_summary((1, 0, 1, 0), fr)
        self.assertTrue(s["all_statistics_from_one_legal_P2_history"])
        self.assertFalse(s["independent_global_source_extrema_used"])
        self.assertEqual(s["pseudo_update_cadence_s"], [0.5, 2.0])
        self.assertEqual(s["sigma_squared_upper"], 3.0)
        self.assertEqual(s["q_c_upper"], 6.0)
        self.assertEqual(s["S_measurement_variance_upper"], 7.0)
        self.assertTrue(math.isfinite(s["history_duration_s"][1]))


if __name__ == "__main__":
    unittest.main()
