#!/usr/bin/env python3
from __future__ import annotations

import copy
import sys
import unittest
from pathlib import Path

TOOLS = Path(__file__).resolve().parents[2] / "tools" / "stability"
sys.path.insert(0, str(TOOLS))

from ou3_interval import Interval
import ou3_p4_joint_brmm_frontend_transition as JOINT
import ou3_p4_source_uniform_production_lineage as LINEAGE


class SourceUniformProductionLineageTest(unittest.TestCase):
    def test_transition_operator_closes_without_promoting_p4(self) -> None:
        d = LINEAGE.build()
        self.assertEqual(LINEAGE.validate(d), [])
        self.assertTrue(d["branch_complete_recursive_transition_operator_materialized"])
        self.assertTrue(d["source_cover_all_BRMM_continuations_covered_by_relation_plus_operator"])
        self.assertTrue(d["all_BIAS0_BIAS1_BIAS2_absolute_prefix_cells_consumed"])
        self.assertTrue(d["one_radial_scale_retained_across_entire_lineage"])
        self.assertFalse(d["favorable_successor_selected"])
        self.assertFalse(d["branch_ordinal_correspondence_assumed"])
        self.assertFalse(d["production_augmented_PrefixInput_assembled_here"])
        self.assertFalse(d["P4_MOTION_PASS"])
        self.assertFalse(d["P4_PASS"])
        self.assertFalse(d["P5_MAY_START"])

    def test_two_sample_smoke_retains_every_family(self) -> None:
        s = LINEAGE._smoke()
        self.assertEqual(s["families"], ["BIAS0", "BIAS1", "BIAS2"])
        self.assertTrue(s["all_nonempty"])
        self.assertTrue(s["all_two_sample_lineages"])
        self.assertTrue(s["all_literal_cells_nonempty"])
        self.assertTrue(s["all_endpoint_frontends_joint_owned"])

    def test_radial_origin_cannot_be_restarted_outside_unit_interval(self) -> None:
        js = JOINT._smoke_state()
        samples = (LINEAGE._sample(True, False),)
        with self.assertRaises(ValueError):
            LINEAGE.execute_lineages(
                joint_entry=js, P0_H=LINEAGE._identity(18), P0_A=LINEAGE._identity(21),
                H_state0=[LINEAGE.I(0) for _ in range(18)],
                A_state0=[LINEAGE.I(0) for _ in range(21)], samples=samples,
                family="BIAS0", radial_scale=Interval(-0.1, 1.0),
            )

    def test_branch_limit_fails_closed_instead_of_selecting_successor(self) -> None:
        js = JOINT._smoke_state()
        samples = (LINEAGE._sample(True, False),)
        with self.assertRaises(RuntimeError):
            LINEAGE.execute_lineages(
                joint_entry=js, P0_H=LINEAGE._identity(18), P0_A=LINEAGE._identity(21),
                H_state0=[LINEAGE.I(0) for _ in range(18)],
                A_state0=[LINEAGE.I(0) for _ in range(21)], samples=samples,
                family="BIAS1", branch_limit=0,
            )


if __name__ == "__main__":
    unittest.main()
