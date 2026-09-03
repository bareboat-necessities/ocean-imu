#!/usr/bin/env python3
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))

import ou3_p3_matched_history_cost_frontier as M


class MatchedHistoryCostFrontierTests(unittest.TestCase):
    def test_cost_dominance_requires_more_adverse_and_no_more_cost(self):
        a = (3, 2, 4, 1)
        b = (2, 2, 3, 1)
        self.assertTrue(M.cost_dominates(a, 20, b, 20))
        self.assertTrue(M.cost_dominates(a, 19, b, 20))
        self.assertFalse(M.cost_dominates(a, 21, b, 20))
        self.assertFalse(M.cost_dominates(b, 19, a, 20))

    def test_insert_keeps_cheapest_same_label(self):
        q = (1, 2, 3, 4)
        front = {q: 20}
        self.assertFalse(M.insert_cost_frontier(front, q, 21))
        self.assertTrue(M.insert_cost_frontier(front, q, 19))
        self.assertEqual(front, {q: 19})

    def test_insert_removes_only_cost_adverse_dominated_state(self):
        weak = (1, 1, 1, 1)
        strong = (2, 2, 2, 2)
        front = {weak: 30}
        self.assertTrue(M.insert_cost_frontier(front, strong, 25))
        self.assertEqual(front, {strong: 25})

        # A more adverse state that costs more cannot replace the cheaper weak
        # state because it leaves less finite-word continuation budget.
        front = {weak: 20}
        self.assertTrue(M.insert_cost_frontier(front, strong, 25))
        self.assertEqual(front, {weak: 20, strong: 25})

    def test_min_gap_successors_uses_cheapest_exact_support(self):
        rt = {
            "nodes": [{}, {}, {}],
            "gaps": [2, 3],
            "labelled_successors": [
                [{1}, {1, 2}],
                [{1}, {2}],
                [{2}, {0}],
            ],
        }
        edges = M.min_gap_successors(rt)
        self.assertEqual(edges[0], ((1, 2), (2, 3)))
        self.assertEqual(edges[1], ((1, 2), (2, 3)))
        self.assertEqual(edges[2], ((0, 3), (2, 2)))


if __name__ == "__main__":
    unittest.main()
