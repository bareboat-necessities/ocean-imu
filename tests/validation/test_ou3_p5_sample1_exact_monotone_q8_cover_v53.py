import copy
import pathlib
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import ou3_p5_sample1_exact_monotone_cover_lift_v52 as V52
import ou3_p5_sample1_exact_monotone_q8_cover_v53 as V53


def _cover(*, unclosed, first_cell, first_q, worst_q):
    return {
        "P5_SAMPLE1_V40_FULL_SOURCE_CELL0_Q8_LIFT_V41": "NOT_ESTABLISHED",
        "evaluated_signed_cayley_cells": 461376,
        "unclosed_q8_cells": unclosed,
        "geodesic_bound_newly_closed_cells": 19444,
        "current_yz_support_newly_closed_cells": 148,
        "max_sample1_q_after_product_tightening": 2.6625654378541017,
        "first_unclosed_q8_cell": {
            "p_cell": first_cell[0],
            "tangent_residual_cell": first_cell[1],
            "axial_residual_cell": first_cell[2],
            "post_sample1_cayley_norm_upper": first_q,
            "sample1_current_cayley_norm_upper": 0.74,
            "correction_radial_lower_rad": 1.19,
            "correction_radial_upper_rad": 1.96,
        },
        "worst_q8_cell": {
            "p_cell": 23, "tangent_residual_cell": 18,
            "axial_residual_cell": 4,
            "post_sample1_cayley_norm_upper": worst_q,
        },
    }


PARENT = _cover(unclosed=235738, first_cell=(0, 0, 23),
                first_q=V53.V41_Q_POST, worst_q=525593.677323)
REFINED = _cover(unclosed=219574, first_cell=(0, 2, 23),
                 first_q=8.475205389989586, worst_q=499303.8238549043)


def _ok_artifact():
    parent = V53._summary(PARENT)
    refined = V53._summary(REFINED)
    newly = parent["unclosed_q8_cells"] - refined["unclosed_q8_cells"]
    return {
        "schema": V53.SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_EXACT_MONOTONE_Q8_COVER_V53",
        "source_generated_not_trajectory_fit": True,
        "exact_monotone_corner_enclosure_used": True,
        "both_covers_regenerated_here": True,
        "temporary_V8_block_hook_restored": True,
        "source_replay_used": False,
        "filter_changed": False,
        "archived_parent_cover_quoted_instead_of_regenerated": False,
        "deployed_correction_limit_increased": False,
        "q8_composed_here": False,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_established_here": False,
        "V41_first_survivor_row": list(V53.WITNESS),
        "archived_V41_post_sample1_q_reference": V53.V41_Q_POST,
        "deployed_correction_limit_rad": 6.0,
        "q_target": V53.Q_TARGET,
        "cover_comparison": {
            "parent": parent,
            "refined": refined,
            "additional_cells_closed": newly,
            "open_cell_reduction_fraction": newly / parent["unclosed_q8_cells"],
            "archived_first_survivor_closed_by_refinement": True,
            "parent_worst_q_upper": V53._worst_q(parent),
            "refined_worst_q_upper": V53._worst_q(refined),
            "q_target": V53.Q_TARGET,
            "cover_fully_closed": False,
        },
        "source_cell0_q8_cover_fully_closed": False,
        "P5_SAMPLE1_EXACT_MONOTONE_Q8_COVER_V53": "PASS",
        "next_obligation": "REFINE_NOMINAL_GEOMETRY_AT_THE_REMAINING_FIRST_UNCLOSED_Q8_CELL",
        "failures": [],
    }


class Sample1ExactMonotoneQ8CoverV53Tests(unittest.TestCase):
    def test_targets_and_archived_reference_are_frozen(self):
        self.assertEqual(V53.SCHEMA, 5300)
        self.assertEqual(V53.WITNESS, (0, 0, 23))
        self.assertEqual(V53.Q_TARGET, 8.0)
        self.assertEqual(V53.V41_Q_POST, 8.344528951460543)

    def test_summary_extracts_counts_and_cells(self):
        s = V53._summary(PARENT)
        self.assertEqual(s["evaluated_signed_cayley_cells"], 461376)
        self.assertEqual(s["unclosed_q8_cells"], 235738)
        self.assertEqual(s["closed_q8_cells"], 461376 - 235738)
        self.assertEqual(V53._first_key(s), (0, 0, 23))
        self.assertEqual(V53._first_q(s), V53.V41_Q_POST)
        self.assertEqual(V53._worst_q(s), 525593.677323)

    def test_missing_cells_do_not_crash_the_summary(self):
        cover = dict(PARENT)
        cover["first_unclosed_q8_cell"] = None
        cover["worst_q8_cell"] = None
        s = V53._summary(cover)
        self.assertIsNone(s["first_unclosed_q8_cell"])
        self.assertIsNone(V53._first_key(s))
        self.assertEqual(V53._worst_q(s), float("inf"))

    def test_v53_uses_the_same_hook_as_the_cover_lift(self):
        self.assertIs(V53.V52._exact_block, V52._exact_block)

    def test_validate_accepts_a_well_formed_artifact(self):
        self.assertEqual(V53.validate(_ok_artifact()), [])

    def test_validate_rejects_a_regressed_or_unanchored_artifact(self):
        more_open = copy.deepcopy(_ok_artifact())
        more_open["cover_comparison"]["refined"]["unclosed_q8_cells"] = 240000
        self.assertTrue(V53.validate(more_open))

        miscounted = copy.deepcopy(_ok_artifact())
        miscounted["cover_comparison"]["additional_cells_closed"] = 99
        self.assertTrue(V53.validate(miscounted))

        worse_q = copy.deepcopy(_ok_artifact())
        worse_q["cover_comparison"]["refined_worst_q_upper"] = 1.0e9
        self.assertTrue(V53.validate(worse_q))

        resized = copy.deepcopy(_ok_artifact())
        resized["cover_comparison"]["refined"]["evaluated_signed_cayley_cells"] = 1000
        self.assertTrue(V53.validate(resized))

        survivor_open = copy.deepcopy(_ok_artifact())
        survivor_open["cover_comparison"][
            "archived_first_survivor_closed_by_refinement"] = False
        self.assertTrue(V53.validate(survivor_open))

        quoted = copy.deepcopy(_ok_artifact())
        quoted["archived_parent_cover_quoted_instead_of_regenerated"] = True
        self.assertTrue(V53.validate(quoted))

        not_regenerated = copy.deepcopy(_ok_artifact())
        not_regenerated["both_covers_regenerated_here"] = False
        self.assertTrue(V53.validate(not_regenerated))

        leaked = copy.deepcopy(_ok_artifact())
        leaked["temporary_V8_block_hook_restored"] = False
        self.assertTrue(V53.validate(leaked))

        promoted = copy.deepcopy(_ok_artifact())
        promoted["N_H_words_set_here"] = True
        self.assertTrue(V53.validate(promoted))

        drifted = copy.deepcopy(_ok_artifact())
        drifted["archived_V41_post_sample1_q_reference"] = 8.0
        self.assertTrue(V53.validate(drifted))

        mismatched = copy.deepcopy(_ok_artifact())
        mismatched["source_cell0_q8_cover_fully_closed"] = True
        self.assertTrue(V53.validate(mismatched))


if __name__ == "__main__":
    unittest.main()
