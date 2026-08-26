from __future__ import annotations
import math
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import ou3_p5_sample1_structured_full_gain_v12d as V12D


class StructuredFullGainV12DTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = V12D.build(source_pieces=4, source_cell_index=0,
                           p_pieces=2, tangent_pieces=2, axial_pieces=2)

    def test_first_psd_uses_tangent_innovation_only(self):
        self.assertIs(self.d["V12C_one_plus_two_resolvent_retained"], True)
        self.assertIs(self.d["V11_first_PSD_generic_axial_noise_floor_retired"], True)
        self.assertIs(self.d["first_PSD_tangent_innovation_exact_structure_used"], True)
        self.assertIs(self.d["first_PSD_axial_innovation_perturbation_exact_zero"], True)
        if self.d.get("rows"):
            row = self.d["rows"][0]
            self.assertIs(row["first_PSD_innovation_perturbation_tangent_only"], True)
            self.assertIs(row["first_PSD_innovation_axial_row_column_exact_zero"], True)
            self.assertGreater(row["first_perturbed_tangent_innovation_lower"], 0.0)
            self.assertTrue(math.isfinite(row["first_perturbed_tangent_inverse_operator_upper"]))

    def test_coarse_fixture_remains_fail_closed(self):
        self.assertIn(self.d["P5_SAMPLE1_FIRST_PSD_TANGENT_REFINEMENT_V12D"],
                      ("PASS", "NOT_ESTABLISHED"))
        vf = V12D.validate(self.d)
        v10_missing = any("V10 prerequisite did not pass" in x for x in vf)
        if v10_missing:
            self.assertEqual(self.d["P5_SAMPLE1_FIRST_PSD_TANGENT_REFINEMENT_V12D"],
                             "NOT_ESTABLISHED")
        else:
            self.assertEqual(vf, [])
            if self.d["P5_SAMPLE1_FIRST_PSD_TANGENT_REFINEMENT_V12D"] == "PASS":
                self.assertEqual(self.d["unclosed_joint_cells"], 0)
            else:
                self.assertIsNotNone(self.d["first_unclosed_joint_cell"])

    def test_no_filter_or_theorem_widening(self):
        self.assertEqual(float(self.d["deployed_correction_limit_rad"]), 6.0)
        for k in (
            "source_replay_used", "filter_changed", "deployed_correction_limit_increased",
            "complete_sample1_branch_closed_here", "signed_cayley_q8_composed_here",
            "q8_word_promoted_here", "whole_word_promoted_here", "N_H_words_set_here",
        ):
            self.assertIs(self.d[k], False)


if __name__ == "__main__":
    unittest.main()
