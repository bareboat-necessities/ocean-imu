#!/usr/bin/env python3
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))

import ou3_p3_matched_history_label_frontier as M


class MatchedHistoryLabelFrontierTests(unittest.TestCase):
    def test_maximum_complete_segments_matches_635_over_13(self):
        self.assertEqual(M.maximum_complete_segments(635, 13), 49)
        self.assertEqual(M.maximum_complete_segments(26, 13), 2)
        with self.assertRaises(ValueError):
            M.maximum_complete_segments(0, 13)

    def test_reduce_keeps_only_componentwise_adverse_labels(self):
        labels = {
            (1, 1, 1, 1),
            (2, 1, 1, 1),
            (1, 2, 1, 1),
            (2, 2, 1, 1),
        }
        self.assertEqual(M._reduce(labels), {(2, 2, 1, 1)})

    def test_allowed_node_set_expands_with_more_adverse_label(self):
        ranks = [
            (0, 0, 0, 0),
            (1, 0, 1, 0),
            (1, 1, 2, 1),
        ]
        a = M._allowed_nodes((1, 0, 1, 0), ranks)
        b = M._allowed_nodes((1, 1, 2, 1), ranks)
        self.assertEqual(a, [0, 1])
        self.assertEqual(b, [0, 1, 2])
        self.assertTrue(set(a).issubset(b))

    def test_union_successor_step_updates_source_before_transition(self):
        rt = {"union_successors": [[1], [1]]}
        ranks = [(1, 2, 3, 4), (4, 3, 2, 1)]
        front = {0: {(-1, -1, -1, -1)}}
        out = M._step(front, rt, ranks)
        self.assertEqual(out, {1: {(1, 2, 3, 4)}})

    def test_global_reduction_includes_terminal_source_before_forgetting_endpoint(self):
        ranks = [(1, 2, 3, 4), (4, 3, 2, 1)]
        front = {1: {(1, 2, 3, 4)}}
        out = M._global_reduce_with_terminal_source(front, ranks)
        self.assertEqual(out, {(4, 3, 3, 4)})


if __name__ == "__main__":
    unittest.main()
