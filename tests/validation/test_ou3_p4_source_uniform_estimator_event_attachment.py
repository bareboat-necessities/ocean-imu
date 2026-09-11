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
        self.assertTrue(d["kernel_authoritative_for_current_Riccati_slice_only"])
        self.assertTrue(d["JOINT_image_authoritative_for_next_frontend_state"])
        self.assertTrue(d["current_schedule_committed_before_post_measurement_adaptation"])
        self.assertTrue(d["every_joint_successor_retained"])
        self.assertTrue(d["all_kernel_frontend_children_share_current_post_Riccati_state"])
        self.assertTrue(d["trusted_event_local_P_H_R_imported_without_reconstruction"])
        self.assertTrue(d["exact_nonlinear_state_succession_materialized_between_literal_events"])
        self.assertTrue(d["same_finite_map_generates_state_and_Jacobian"])
        self.assertTrue(d["actual_applied_RS_owned_by_same_joint_image"])
        self.assertFalse(d["kernel_frontend_successor_used_to_select_next_theorem_state"])
        self.assertFalse(d["kernel_joint_child_frontend_equality_required"])
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
        self.assertGreater(s["joint_successors_retained"], 0)
        self.assertTrue(s["all_H_cells_bound"])
        self.assertTrue(s["all_A_cells_bound"])
        self.assertTrue(s["all_event_orders_equal"])
        self.assertTrue(s["authoritative_next_frontend_is_joint"])


if __name__ == "__main__":
    unittest.main()
