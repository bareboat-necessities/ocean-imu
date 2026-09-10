from __future__ import annotations

import copy
import pathlib
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools" / "stability"))

from ou3_interval import Interval
import ou3_p4_complete_brmm_same_history_prefix_selectors as SELECTORS
import ou3_p4_source_reachable_selector_family as FAMILY


class SourceReachableSelectorFamilyTests(unittest.TestCase):
    def test_relation_cover_closes_without_promoting_p4(self):
        d = FAMILY.build()
        self.assertEqual(FAMILY.validate(d), [])
        self.assertTrue(d["source_reachable_COMPLETE_BRMM_selector_family_relation_closed"])
        self.assertTrue(d["complete_BRMM_left_inclusion_closed"])
        self.assertTrue(d["frontend_predecessor_family_closed"])
        self.assertTrue(d["same_history_estimator_coefficient_relation_closed"])
        self.assertTrue(d["branch_complete_prefix_selector_relation_available"])
        self.assertTrue(d["full_zero_to_one_hard_entry_radial_relation_attached"])
        self.assertTrue(d["all_bias_absolute_prefix_relations_attached"])
        self.assertFalse(d["numeric_601_sample_interval_hull_materialized_here"])
        self.assertFalse(d["P4_PASS"])
        self.assertEqual(d["P3_delta"], 1e-18)

    def test_false_shortcuts_are_rejected(self):
        d = FAMILY.build()
        d["favorable_branch_selected"] = True
        d["independent_sample_boxes_generate_estimator_history"] = True
        d["finite_frequency_or_direction_grid_used"] = True
        d["P4_PASS"] = True
        f = FAMILY.validate(d)
        self.assertIn("favorable_branch_selected not false", f)
        self.assertIn("independent_sample_boxes_generate_estimator_history not false", f)
        self.assertIn("finite_frequency_or_direction_grid_used not false", f)
        self.assertIn("P4_PASS not false", f)

    def test_endpoint_attachment_keeps_full_radial_and_all_bias_families(self):
        frontend = SELECTORS.FRONTEND._point_state()
        P0_H, P0_A, _ = SELECTORS._live_structured_point_covariance_fixture(
            frontend, FAMILY.OUTER.DEFAULT_DOMAIN
        )
        sample = SELECTORS._point_sample()
        endpoints, selectors, _ = SELECTORS.execute_with_prefix_selectors(
            frontend_entry=frontend,
            P0_H=P0_H,
            P0_A=P0_A,
            samples=[sample, sample],
            branch_limit=128,
        )
        endpoint = endpoints[0].source_cell_id
        a = FAMILY.attach_endpoint_family(selectors, endpoint)
        self.assertEqual(a.endpoint_source_cell_id, endpoint)
        self.assertEqual(len(a.selector_lineage), 2)
        self.assertEqual((a.radial_scale.lo, a.radial_scale.hi), (0.0, 1.0))
        self.assertEqual(set(a.bias_lineages), {"BIAS0", "BIAS1", "BIAS2"})
        for cert in a.bias_lineages.values():
            self.assertTrue(cert.source_uniform)
            self.assertTrue(cert.recurrence_outer_relation)
            self.assertEqual(set(cert.prefix_boxes), {s.source_cell_id for s in a.selector_lineage})

    def test_broken_selector_endpoint_cannot_be_attached(self):
        with self.assertRaises(ValueError):
            FAMILY.attach_endpoint_family([], "missing")


if __name__ == "__main__":
    unittest.main()
