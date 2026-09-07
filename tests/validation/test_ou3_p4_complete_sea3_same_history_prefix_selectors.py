#!/usr/bin/env python3
import copy
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "tools" / "stability"))

import ou3_p4_complete_sea3_same_history_prefix_selectors as SELECTORS
import ou3_sea3_complete_window_execution_kernel as KERNEL
import ou3_sea3_frontend_state_step as FRONTEND


class P4CompleteSea3SameHistoryPrefixSelectorsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.payload = SELECTORS.build()

    def test_status_is_non_promoting_and_preserves_frozen_contract(self):
        self.assertEqual(SELECTORS.validate(self.payload), [])
        self.assertEqual(
            self.payload["canonical_source"], "COMPLETE_SEA3_NORMAL_LIVE_WORD"
        )
        self.assertEqual(self.payload["P3_delta_preserved"], 1.0e-18)
        self.assertTrue(
            self.payload["branch_correlated_every_prefix_selector_available"]
        )
        self.assertTrue(
            self.payload["actual_applied_anisotropic_RS_retained_per_transition"]
        )
        self.assertFalse(self.payload["source_generator"])
        self.assertFalse(self.payload["trajectory_replay_used"])
        self.assertFalse(self.payload["source_uniform_provider_family_closed_here"])
        self.assertFalse(self.payload["P4_promoted_here"])

    def test_smoke_retains_every_prefix_and_both_full_state_modes(self):
        smoke = self.payload["smoke"]
        self.assertEqual(smoke["samples_executed"], 2)
        self.assertGreaterEqual(smoke["prefix_selectors"], 2)
        self.assertTrue(smoke["selector_graph_valid"])
        self.assertEqual(smoke["selector_graph_failures"], [])
        self.assertTrue(smoke["all_endpoint_lineages_cover_every_prefix"])
        self.assertTrue(smoke["all_selectors_retain_actual_RS"])
        self.assertTrue(smoke["all_selectors_retain_H18_A21_before_after"])
        self.assertTrue(smoke["all_selectors_retain_exact_event_slices"])
        self.assertFalse(smoke["favorable_frontend_successor_selected"])
        self.assertFalse(smoke["shipping_transition_reimplemented"])

    def test_selector_event_slice_matches_exact_shipping_sample(self):
        sample = SELECTORS._point_sample()
        endpoints, selectors, meta = SELECTORS.execute_with_prefix_selectors(
            frontend_entry=FRONTEND._point_state(),
            P0_H=KERNEL._diag_P(18, 2.0),
            P0_A=KERNEL._diag_P(21, 2.0),
            samples=[sample],
            branch_limit=64,
        )
        self.assertGreaterEqual(len(endpoints), 1)
        self.assertEqual(meta["samples_executed"], 1)
        self.assertEqual(len(selectors), len(endpoints))
        self.assertEqual(
            SELECTORS.validate_selector_graph(
                selectors,
                samples_executed=1,
                endpoint_source_cell_ids=[b.source_cell_id for b in endpoints],
            ),
            [],
        )
        for selector in selectors:
            self.assertEqual(selector.parent_source_cell_id, "root")
            self.assertEqual(selector.prefix_length, 1)
            self.assertEqual(selector.sample_index, 0)
            self.assertEqual(
                selector.H_events_this_sample,
                ("prediction", "aw_floor", "S_zero", "accelerometer", "magnetometer"),
            )
            self.assertEqual(
                selector.A_events_this_sample,
                ("prediction", "aw_floor", "S_zero", "accelerometer", "magnetometer"),
            )
            self.assertEqual(len(selector.actual_rs_std_xyz), 3)
            self.assertEqual(selector.H_before.mode, "H")
            self.assertEqual(selector.H_after.mode, "H")
            self.assertEqual(selector.A_before.mode, "A")
            self.assertEqual(selector.A_after.mode, "A")
            self.assertEqual(selector.H_after.S_updates, selector.H_before.S_updates + 1)
            self.assertEqual(selector.A_after.S_updates, selector.A_before.S_updates + 1)

    def test_lineage_is_explicit_across_two_prefixes(self):
        sample = SELECTORS._point_sample()
        endpoints, selectors, _ = SELECTORS.execute_with_prefix_selectors(
            frontend_entry=FRONTEND._point_state(),
            P0_H=KERNEL._diag_P(18, 2.0),
            P0_A=KERNEL._diag_P(21, 2.0),
            samples=[copy.deepcopy(sample), copy.deepcopy(sample)],
            branch_limit=128,
        )
        self.assertGreaterEqual(len(endpoints), 1)
        for endpoint in endpoints:
            lineage = SELECTORS.lineage_for_endpoint(
                selectors, endpoint.source_cell_id
            )
            self.assertEqual([s.prefix_length for s in lineage], [1, 2])
            self.assertEqual(lineage[0].parent_source_cell_id, "root")
            self.assertEqual(
                lineage[1].parent_source_cell_id, lineage[0].source_cell_id
            )
            self.assertEqual(lineage[-1].source_cell_id, endpoint.source_cell_id)

    def test_broken_or_duplicate_ancestry_is_rejected(self):
        sample = SELECTORS._point_sample()
        endpoints, selectors, _ = SELECTORS.execute_with_prefix_selectors(
            frontend_entry=FRONTEND._point_state(),
            P0_H=KERNEL._diag_P(18, 2.0),
            P0_A=KERNEL._diag_P(21, 2.0),
            samples=[sample],
            branch_limit=64,
        )
        self.assertGreaterEqual(len(selectors), 1)
        duplicate = list(selectors) + [selectors[0]]
        failures = SELECTORS.validate_selector_graph(
            duplicate,
            samples_executed=1,
            endpoint_source_cell_ids=[b.source_cell_id for b in endpoints],
        )
        self.assertTrue(any("duplicate source cell id" in x for x in failures))


if __name__ == "__main__":
    unittest.main()
