from __future__ import annotations
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import ou3_p5_sample1_signed_radial_subcells_v13d as V13D


class Sample1SignedRadialSubcellsV13DTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = V13D.build(source_pieces=4, source_cell_index=0,
                           p_pieces=2, tangent_pieces=2, axial_pieces=2,
                           residual_x_pieces=4, parallel_pieces=4)

    def test_adapter_changes_no_filter_or_shipping_limit(self):
        self.assertIs(self.d["source_generated_not_trajectory_fit"], True)
        self.assertIs(self.d["source_replay_used"], False)
        self.assertIs(self.d["filter_changed"], False)
        self.assertEqual(float(self.d["deployed_correction_limit_rad"]), 6.0)
        self.assertIs(self.d["deployed_correction_limit_increased"], False)

    def test_v12d_adapter_is_explicit(self):
        self.assertIs(self.d["V12D_tangent_channel_prerequisite_used"], True)
        self.assertIs(self.d["V12D_attitude_gain_perturbation_used"], True)
        self.assertIs(self.d["V12D_radial_parent_upper_used"], True)

    def test_coarse_fixture_fails_closed_until_v12d(self):
        st = self.d["P5_SAMPLE1_SIGNED_RADIAL_SUBCELLS_V13D"]
        self.assertIn(st, ("PASS", "NOT_ESTABLISHED"))
        if not self.d["V12D_prerequisite_passed"]:
            self.assertEqual(st, "NOT_ESTABLISHED")
            self.assertEqual(self.d["evaluated_signed_subcells"], 0)
            self.assertTrue(any("V12D prerequisite did not pass" in x
                                for x in V13D.validate(self.d)))
        else:
            self.assertEqual(V13D.validate(self.d), [])

    def test_never_promotes_q8_or_word(self):
        for k in (
            "signed_cayley_q8_composed_here", "complete_sample1_branch_closed_here",
            "q8_word_promoted_here", "whole_word_promoted_here", "N_H_words_set_here",
        ):
            self.assertIs(self.d[k], False)


if __name__ == "__main__":
    unittest.main()
