#!/usr/bin/env python3
from __future__ import annotations

import sys
import unittest
from dataclasses import replace
from pathlib import Path

TOOLS=Path(__file__).resolve().parents[2]/"tools"/"stability"
sys.path.insert(0,str(TOOLS))

import ou3_p4_source_uniform_event_lineage_induction as IND


class SourceUniformEventLineageInductionTest(unittest.TestCase):
    def test_relation_closes_without_promoting_p4(self):
        d=IND.build()
        self.assertEqual(IND.validate(d),[])
        self.assertTrue(d["literal_P_H_R_state_attachment_closed_by_induction"])
        self.assertTrue(d["source_uniform_local_event_coefficient_cover_closed"])
        self.assertTrue(d["primitive_prefix_structural_premises_consumed"])
        self.assertFalse(d["numeric_D_S_qualification_closed_here"])
        self.assertFalse(d["endpoint_augmented_LDLT_closed_here"])
        self.assertFalse(d["every_prefix_augmented_LDLT_closed_here"])
        self.assertFalse(d["P4_PASS"])
        self.assertFalse(d["P5_MAY_START"])

    def test_two_step_successor_is_valid_and_radial_origin_is_retained(self):
        s=IND._smoke()
        self.assertGreater(s["first_successors"],0)
        self.assertGreater(s["second_successors"],0)
        self.assertTrue(s["all_first_prefix_lengths_one"])
        self.assertTrue(s["all_second_prefix_lengths_two"])
        self.assertTrue(s["all_successors_valid"])
        self.assertTrue(s["all_radial_cells_preserved"])

    def test_detached_frontend_is_rejected(self):
        root,_=IND._root_and_sample()
        bad=replace(root,joint_state=replace(root.joint_state,frontend=IND.JOINT._smoke_state().frontend))
        # Replace with a structurally different frontend by advancing its token-bearing state once.
        # If the smoke states happen to compare equal, detach the branch frontend instead.
        if not IND.validate_node(bad):
            bad=replace(root,branch=replace(root.branch,frontend=IND.ATTACH.next_execution_branch(*IND.ATTACH.synchronize_sample(
                branch=root.branch,joint_state=root.joint_state,sample=IND._root_and_sample()[1],
                state_in_H=root.H_error_state,state_in_A=root.A_error_state,radial_scale=root.radial_scale,
                true_bias=root.true_bias,bias_projection_limit=.4,tau_ba=IND.I(1800),sample_index=0,next_cell_prefix="mut")[0]).frontend))
        self.assertTrue(IND.validate_node(bad))


if __name__=="__main__":
    unittest.main()
