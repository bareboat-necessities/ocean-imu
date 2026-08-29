import math
import pathlib
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from ou3_interval import Interval
import ou3_p5_sample1_structured_full_gain_v8 as V8
import ou3_p5_sample1_exact_monotone_source_gain_v51 as V51
import ou3_p5_sample1_exact_monotone_cover_lift_v52 as V52

G = 9.80665


def _ok_artifact():
    cover = {
        "evaluated_joint_cells": 12816,
        "unclosed_joint_cells": 0,
        "status": "PASS",
        "max_full_residual_norm_upper_mps2": 46.01460061569009,
        "max_Ktheta_perpendicular_block_upper": 1.3986770467171177,
        "max_Ktheta_parallel_block_upper": 0.19289137244367335,
        "max_combined_directional_correction_norm_upper_rad": 7.016940736774492,
    }
    refined = dict(cover)
    refined.update({
        "max_full_residual_norm_upper_mps2": 42.67399431981034,
        "max_Ktheta_perpendicular_block_upper": 1.2187023949259914,
        "max_Ktheta_parallel_block_upper": 0.19082905441803827,
        "max_combined_directional_correction_norm_upper_rad": 6.7576910288356276,
    })
    return {
        "schema": V52.SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_EXACT_MONOTONE_COVER_LIFT_V52",
        "source_generated_not_trajectory_fit": True,
        "exact_monotone_corner_enclosure_used": True,
        "parent_enclosure_retained_as_intersection": True,
        "temporary_V8_block_hook_restored": True,
        "complete_source_cell0_cover_evaluated": True,
        "source_replay_used": False,
        "filter_changed": False,
        "deployed_correction_limit_increased": False,
        "V41_signed_chart_q8_composition_rerun_here": False,
        "q8_composed_here": False,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_established_here": False,
        "V41_first_survivor_row": list(V52.WITNESS),
        "deployed_correction_limit_rad": 6.0,
        "q_target": V52.Q_TARGET,
        "cover_comparison": {
            "cells_evaluated": 12816,
            "cells_narrowed": 12816,
            "cells_widened": 0,
            "every_cell_inside_parent": True,
            "narrowing_ratio_min": 1.0106991820188718,
            "narrowing_ratio_max": 1.495825177084859,
            "parent_cover": cover,
            "refined_cover": refined,
        },
        "authoritative_witness": {
            "V51_status": "PASS",
            "authoritative_witness_closed": True,
        },
        "P5_SAMPLE1_EXACT_MONOTONE_COVER_LIFT_V52": "PASS",
        "next_obligation": "RERUN_V41_SIGNED_CHART_Q8_COMPOSITION_WITH_EXACT_MONOTONE_BLOCK",
        "failures": [],
    }


class Sample1ExactMonotoneCoverLiftV52Tests(unittest.TestCase):
    def test_targets_and_archived_cover_are_frozen(self):
        self.assertEqual(V52.SCHEMA, 5200)
        self.assertEqual(V52.WITNESS, (0, 0, 23))
        self.assertEqual(V52.Q_TARGET, 8.0)
        self.assertEqual(V52.PARENT_COVER["evaluated_joint_cells"], 12816)
        self.assertEqual(
            V52.PARENT_COVER["max_combined_directional_correction_norm_upper_rad"],
            7.016940736774492)

    def test_hook_signature_matches_the_v8_block_helper(self):
        t = Interval.outward_bounds(1.0e-3, 2.0e-3)
        p = Interval.outward_bounds(2.0e-3, 4.0e-3)
        r = Interval.outward_bounds(8.0e-4, 9.0e-4)
        parent = V8._first_block_quantities(t=t, p=p, r=r, g=G)
        exact = V52._exact_block(t=t, p=p, r=r, g=G)
        self.assertEqual(set(parent), set(exact))
        for key in parent:
            self.assertGreaterEqual(exact[key].lo, parent[key].lo, key)
            self.assertLessEqual(exact[key].hi, parent[key].hi, key)

    def test_installing_and_removing_the_hook_restores_the_parent(self):
        root = V8._first_block_quantities
        V8._first_block_quantities = V52._exact_block
        try:
            self.assertIs(V8._first_block_quantities, V52._exact_block)
        finally:
            V8._first_block_quantities = root
        self.assertIs(V8._first_block_quantities, root)

    def test_v51_compares_against_the_unhooked_parent_block(self):
        # V51 binds the parent helper at import time, so a stage that installs a
        # refined block into V8 cannot make V51's comparison vacuous.
        self.assertIsNot(V51._PARENT_FIRST_BLOCK, V52._exact_block)
        root = V8._first_block_quantities
        V8._first_block_quantities = V52._exact_block
        try:
            t = Interval.outward_bounds(1.0e-3, 2.0e-3)
            p = Interval.outward_bounds(2.0e-3, 4.0e-3)
            r = Interval.outward_bounds(8.0e-4, 9.0e-4)
            parent = V51._first_block(t=t, p=p, r=r, g=G, exact=False)
            exact = V51._first_block(t=t, p=p, r=r, g=G, exact=True)
            self.assertLess(exact["kz"].hi, parent["kz"].hi)
        finally:
            V8._first_block_quantities = root

    def test_cell_key_and_correction_extraction(self):
        rows = [{"p_cell": 1, "tangent_residual_cell": 2,
                 "axial_residual_cell": 3, V52.CORRECTION_KEY: 1.5}]
        self.assertEqual(V52._cell_key(rows[0]), (1, 2, 3))
        self.assertEqual(V52._corrections({"rows": rows}), {(1, 2, 3): 1.5})

    def test_validate_accepts_a_well_formed_artifact(self):
        self.assertEqual(V52.validate(_ok_artifact()), [])

    def test_validate_rejects_a_widened_or_promoted_artifact(self):
        widened = _ok_artifact()
        widened["cover_comparison"] = dict(
            widened["cover_comparison"], cells_widened=1,
            every_cell_inside_parent=False)
        self.assertTrue(V52.validate(widened))

        grew = _ok_artifact()
        refined = dict(grew["cover_comparison"]["refined_cover"])
        refined["max_combined_directional_correction_norm_upper_rad"] = 9.0
        grew["cover_comparison"] = dict(grew["cover_comparison"],
                                        refined_cover=refined)
        self.assertTrue(V52.validate(grew))

        shrunk_cover = _ok_artifact()
        shrunk_cover["cover_comparison"] = dict(
            shrunk_cover["cover_comparison"], cells_evaluated=12000)
        self.assertTrue(V52.validate(shrunk_cover))

        unclosed = _ok_artifact()
        refined = dict(unclosed["cover_comparison"]["refined_cover"])
        refined["unclosed_joint_cells"] = 3
        unclosed["cover_comparison"] = dict(unclosed["cover_comparison"],
                                            refined_cover=refined)
        self.assertTrue(V52.validate(unclosed))

        leaked = _ok_artifact()
        leaked["temporary_V8_block_hook_restored"] = False
        self.assertTrue(V52.validate(leaked))

        promoted = _ok_artifact()
        promoted["q8_composed_here"] = True
        self.assertTrue(V52.validate(promoted))

        claimed_rerun = _ok_artifact()
        claimed_rerun["V41_signed_chart_q8_composition_rerun_here"] = True
        self.assertTrue(V52.validate(claimed_rerun))

        open_witness = _ok_artifact()
        open_witness["authoritative_witness"] = {
            "V51_status": "PASS", "authoritative_witness_closed": False}
        self.assertTrue(V52.validate(open_witness))

        limit = _ok_artifact()
        limit["deployed_correction_limit_rad"] = 9.0
        self.assertTrue(V52.validate(limit))


if __name__ == "__main__":
    unittest.main()
