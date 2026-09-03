#!/usr/bin/env python3
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))

import ou3_p3_four_max_global_label_witness as M


class FourMaxGlobalLabelWitnessTests(unittest.TestCase):
    def test_global_rank_and_source_masks(self):
        ranks = [(0, 2, 1, 4), (3, 1, 5, 2), (2, 2, 4, 4)]
        g = M.global_rank_tuple(ranks)
        self.assertEqual(g, (3, 2, 5, 4))
        self.assertEqual(M.source_max_mask(ranks[0], g), 0b1010)
        self.assertEqual(M.source_max_mask(ranks[1], g), 0b0101)
        self.assertEqual(M.source_max_mask(ranks[2], g), 0b1010)

    def test_global_mask_is_all_four_coordinates(self):
        self.assertEqual(M.FULL_MASK, 0b1111)

    def test_min_gap_successors_keeps_cheapest_exact_support(self):
        rt = {
            "nodes": [{}, {}],
            "gaps": [13, 14],
            "labelled_successors": [
                [{1}, {0, 1}],
                [{1}, {0}],
            ],
        }
        self.assertEqual(
            M.min_gap_successors(rt),
            [((0, 14), (1, 13)), ((0, 14), (1, 13))],
        )


if __name__ == "__main__":
    unittest.main()
