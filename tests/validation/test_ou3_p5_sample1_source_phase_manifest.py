from __future__ import annotations
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_sample1_source_phase_manifest as M


class Sample1SourcePhaseManifestTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = M.build(source_pieces=4)

    def test_manifest_validates_and_is_source_complete_semantically(self):
        self.assertEqual(M.validate(self.d), [])
        self.assertEqual(self.d["P5_SAMPLE1_SOURCE_PHASE_MANIFEST"], "PASS")
        self.assertTrue(self.d["same_partition_as_first_accel_rotation_gauge"])
        self.assertGreater(self.d["source_phase_child_count"], 0)
        self.assertGreater(self.d["due_child_count"], 0)
        self.assertGreater(self.d["not_due_child_count"], 0)
        self.assertEqual(
            self.d["due_child_count"] + self.d["not_due_child_count"],
            self.d["source_phase_child_count"],
        )

    def test_due_and_not_due_indices_partition_every_child_exactly_once(self):
        due = self.d["due_source_cell_indices"]
        not_due = self.d["not_due_source_cell_indices"]
        self.assertFalse(set(due) & set(not_due))
        self.assertEqual(
            sorted(due + not_due),
            list(range(self.d["source_phase_child_count"])),
        )
        for row in self.d["source_phase_children"]:
            self.assertEqual(
                row["current_V14D_due_route_eligible"],
                row["sample0_pseudo_phase"] == "due",
            )

    def test_focused_v14d_cell_zero_is_due_but_does_not_promote_family(self):
        self.assertIn(0, self.d["due_source_cell_indices"])
        self.assertTrue(self.d["current_V14D_requires_sample0_due"])
        self.assertTrue(self.d["not_due_identity_branch_must_be_certified_separately"])
        self.assertFalse(self.d["all_source_phase_children_numerically_closed_here"])
        self.assertFalse(self.d["complete_sample1_source_family_promoted_here"])
        self.assertFalse(self.d["whole_word_promoted_here"])
        self.assertFalse(self.d["N_H_words_set_here"])


if __name__ == "__main__":
    unittest.main()
