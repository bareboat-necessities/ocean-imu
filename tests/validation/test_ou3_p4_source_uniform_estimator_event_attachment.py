#!/usr/bin/env python3
from __future__ import annotations

import sys
import unittest
from pathlib import Path

TOOLS = Path(__file__).resolve().parents[2] / "tools" / "stability"
sys.path.insert(0, str(TOOLS))

import ou3_p4_source_uniform_estimator_event_attachment as ATTACH


class SourceUniformEstimatorEventAttachmentTest(unittest.TestCase):
    def test_relation_attachment_closes_without_promoting_p4(self) -> None:
        d = ATTACH.build()
        self.assertEqual(ATTACH.validate(d), [])
        self.assertTrue(d["source_uniform_estimator_owned_event_attachment_relation_closed"])
        self.assertTrue(d["kernel_joint_children_matched_by_exact_frontend_state_not_ordinal"])
        self.assertTrue(d["trusted_event_local_P_H_R_imported_without_reconstruction"])
        self.assertTrue(d["exact_nonlinear_state_succession_materialized_between_literal_events"])
        self.assertTrue(d["same_finite_map_generates_state_and_Jacobian"])
        self.assertTrue(d["actual_applied_RS_owned_by_same_joint_image"])
        self.assertFalse(d["favorable_successor_selected"])
        self.assertFalse(d["branch_ordinal_correspondence_assumed"])
        self.assertFalse(d["independent_P_H_R_K_reconstruction_used"])
        self.assertFalse(d["numeric_COMPLETE_BRMM_enumeration_used"])
        self.assertFalse(d["production_augmented_PrefixInput_assembled_here"])
        self.assertFalse(d["P4_MOTION_PASS"])
        self.assertFalse(d["P4_PASS"])
        self.assertFalse(d["P5_MAY_START"])

    def test_smoke_retains_both_modes_and_literal_order(self) -> None:
        s = ATTACH._smoke()
        self.assertGreater(s["compatible_pairs"], 0)
        self.assertTrue(s["all_H_cells_bound"])
        self.assertTrue(s["all_A_cells_bound"])
        self.assertTrue(s["all_event_orders_equal"])
        self.assertTrue(s["all_frontend_matches_exact"])


if __name__ == "__main__":
    unittest.main()
