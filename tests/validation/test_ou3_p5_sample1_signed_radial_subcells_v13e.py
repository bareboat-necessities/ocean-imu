from __future__ import annotations
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import ou3_p5_sample1_signed_radial_subcells_v13e as V13E


class Sample1SignedRadialSubcellsV13ETests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = V13E.build(source_pieces=4, source_cell_index=0,
                           p_pieces=2, tangent_pieces=2, axial_pieces=2,
                           residual_x_pieces=4, parallel_pieces=4)

    def test_component_box_and_correlated_norm_have_separate_roles(self):
        self.assertIs(self.d["signed_component_intervals_retained_unchanged"], True)
        self.assertIs(self.d["V8_correlated_block_norm_retained_for_radial_upper"], True)
        self.assertIs(self.d["cartesian_component_box_norm_not_compared_to_correlated_parent"], True)
        self.assertIs(self.d["component_box_used_only_for_sign_and_radial_lower"], True)
        self.assertIs(self.d["radial_upper_formula_changed_here"], False)

    def test_source_and_shipping_contracts_unchanged(self):
        self.assertIs(self.d["source_generated_not_trajectory_fit"], True)
        self.assertIs(self.d["source_replay_used"], False)
        self.assertIs(self.d["filter_changed"], False)
        self.assertEqual(float(self.d["deployed_correction_limit_rad"]), 6.0)
        self.assertIs(self.d["deployed_correction_limit_increased"], False)

    def test_result_is_fail_closed_or_closed(self):
        st = self.d["P5_SAMPLE1_SIGNED_RADIAL_SUBCELLS_V13E"]
        self.assertIn(st, ("PASS", "NOT_ESTABLISHED"))
        f = V13E.validate(self.d)
        if self.d.get("V12D_prerequisite_passed"):
            self.assertEqual(f, [])
            if st == "PASS":
                self.assertEqual(self.d["unclosed_radial_subcells"], 0)
                self.assertLessEqual(float(self.d["max_radial_upper"]), 9.0)
            else:
                self.assertIsNotNone(self.d["first_unclosed_radial_subcell"])
        else:
            self.assertTrue(any("prerequisite" in x.lower() for x in f))

    def test_never_promotes_q8_or_word(self):
        for k in (
            "signed_cayley_q8_composed_here", "complete_sample1_branch_closed_here",
            "q8_word_promoted_here", "whole_word_promoted_here", "N_H_words_set_here",
        ):
            self.assertIs(self.d[k], False)


if __name__ == "__main__":
    unittest.main()
