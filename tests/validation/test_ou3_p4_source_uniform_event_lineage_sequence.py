#!/usr/bin/env python3
from __future__ import annotations

import sys
import unittest
from pathlib import Path

TOOLS = Path(__file__).resolve().parents[2] / "tools" / "stability"
sys.path.insert(0, str(TOOLS))

import ou3_p4_source_uniform_event_lineage_sequence as LINEAGE


class SourceUniformEventLineageSequenceTest(unittest.TestCase):
    def test_cross_sample_literal_ancestry_is_contiguous(self) -> None:
        d = LINEAGE.build()
        self.assertEqual(LINEAGE.validate(d), [])
        self.assertTrue(d["global_literal_event_lineage_constructor_closed"])
        self.assertTrue(d["first_event_of_next_sample_points_to_previous_literal_event"])
        self.assertTrue(d["estimator_predecessor_remains_previous_selector_child"])
        self.assertTrue(d["sample_selector_and_literal_event_ancestry_kept_distinct"])
        self.assertFalse(d["finite_source_enumeration_used"])
        self.assertFalse(d["favorable_branch_selection_part_of_constructor"])
        self.assertFalse(d["P4_MOTION_PASS"])
        self.assertFalse(d["P4_PASS"])
        self.assertFalse(d["P5_MAY_START"])

    def test_two_sample_smoke_preserves_estimator_parent(self) -> None:
        s = LINEAGE._two_sample_smoke()
        self.assertGreater(s["H_literal_events"], 0)
        self.assertGreater(s["A_literal_events"], 0)
        self.assertTrue(s["H_cross_sample_predecessor_is_prior_literal_event"])
        self.assertTrue(s["A_cross_sample_predecessor_is_prior_literal_event"])
        self.assertTrue(s["H_estimator_parent_remains_prior_selector"])
        self.assertTrue(s["A_estimator_parent_remains_prior_selector"])


if __name__ == "__main__":
    unittest.main()
