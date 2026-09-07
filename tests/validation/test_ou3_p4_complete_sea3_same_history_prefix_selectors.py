#!/usr/bin/env python3
import copy
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "tools" / "stability"))

import ou3_p4_complete_sea3_same_history_prefix_selectors as SELECTORS
import ou3_sea3_frontend_state_step as FRONTEND


class P4CompleteSea3SameHistoryPrefixSelectorsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.payload = SELECTORS.build()

    @staticmethod
    def _source_seeded_frontend_and_covariance():
        frontend = FRONTEND._point_state()
        P0_H, P0_A, meta = SELECTORS._source_generated_point_covariance_seed(
            frontend, SELECTORS.DEFAULT_DOMAIN
        )
        if meta.get("arbitrary_P0_used") is not False:
            raise AssertionError("selector test reverted to arbitrary covariance")
        return frontend, P0_H, P0_A

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
        self.assertTrue(
            self.payload["event_local_same_P_H_R_cells_materialized_on_typed_execution"]
        )
        self.assertTrue(
            self.payload["event_local_due_S_cells_use_exact_actual_committed_RS"]
        )
        self.assertFalse(self.payload["source_generator"])
        self.assertFalse(self.payload["trajectory_replay_used"])
        self.assertFalse(self.payload["source_uniform_provider_family_closed_here"])
        self.assertFalse(
            self.payload["source_uniform_event_local_same_P_H_R_cells_closed_here"]
        )
        self.assertFalse(self.payload["nonlinear_residual_history_graph_materialized_here"])
        self.assertFalse(self.payload["A21_projection_graph_attached_here"])
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
        self.assertTrue(smoke["all_measurement_event_cells_retain_same_P_H_R"])
        self.assertTrue(smoke["all_due_S_cells_retain_actual_committed_RS"])
        self.assertTrue(smoke["event_local_cells_captured_inside_same_transition"])
        seed = smoke["source_generated_covariance_seed"]
        self.assertTrue(seed["live_seed_contract_consumed"])
        self.assertFalse(seed["arbitrary_P0_used"])
        self.assertTrue(seed["aw_seed_uses_same_committed_sigma"])
        self.assertTrue(seed["A21_ba_release_floor_attached"])
        self.assertFalse(smoke["favorable_frontend_successor_selected"])
        self.assertFalse(smoke["shipping_transition_reimplemented"])

    def test_selector_event_cells_match_exact_shipping_sample(self):
        sample = SELECTORS._point_sample()
        frontend, P0_H, P0_A = self._source_seeded_frontend_and_covariance()
        endpoints, selectors, meta = SELECTORS.execute_with_prefix_selectors(
            frontend_entry=frontend,
            P0_H=P0_H,
            P0_A=P0_A,
            samples=[sample],
            branch_limit=64,
        )
        self.assertGreaterEqual(len(endpoints), 1)
        self.assertEqual(meta["samples_executed"], 1)
        self.assertTrue(meta["event_local_cells_captured_inside_same_transition"])
        self.assertEqual(len(selectors), len(endpoints))
        self.assertEqual(
            SELECTORS.validate_selector_graph(
                selectors,
                samples_executed=1,
                endpoint_source_cell_ids=[b.source_cell_id for b in endpoints],
            ),
            [],
        )
        expected_events = (
            "prediction",
            "aw_floor",
            "S_zero",
            "accelerometer",
            "magnetometer",
        )
        for selector in selectors:
            self.assertEqual(selector.parent_source_cell_id, "root")
            self.assertEqual(selector.prefix_length, 1)
            self.assertEqual(selector.sample_index, 0)
            self.assertEqual(selector.H_events_this_sample, expected_events)
            self.assertEqual(selector.A_events_this_sample, expected_events)
            self.assertEqual(
                tuple(cell.kind for cell in selector.H_event_cells), expected_events
            )
            self.assertEqual(
                tuple(cell.kind for cell in selector.A_event_cells), expected_events
            )
            self.assertEqual(len(selector.actual_rs_std_xyz), 3)
            self.assertEqual(selector.H_before.mode, "H")
            self.assertEqual(selector.H_after.mode, "H")
            self.assertEqual(selector.A_before.mode, "A")
            self.assertEqual(selector.A_after.mode, "A")
            self.assertEqual(selector.H_after.S_updates, selector.H_before.S_updates + 1)
            self.assertEqual(selector.A_after.S_updates, selector.A_before.S_updates + 1)

            for mode, cells in (
                ("H", selector.H_event_cells),
                ("A", selector.A_event_cells),
            ):
                self.assertEqual([c.event_index_in_sample for c in cells], list(range(5)))
                self.assertIsNotNone(cells[0].F)
                self.assertIsNotNone(cells[0].Q)
                self.assertIsNotNone(cells[1].floor_increment)
                s_cell = cells[2]
                self.assertTrue(s_cell.actual_rs_from_committed_schedule)
                self.assertEqual(
                    s_cell.R, SELECTORS.WORD.R_S_zero(selector.actual_rs_std_xyz)
                )
                self.assertEqual(s_cell.H, SELECTORS.WORD.H_S_zero(mode))
                for cell in cells[2:]:
                    self.assertIsNotNone(cell.H)
                    self.assertIsNotNone(cell.R)
                    self.assertTrue(cell.P_before)
                    self.assertTrue(cell.P_after)
                final_P = (
                    selector.H_after.riccati.P
                    if mode == "H"
                    else selector.A_after.riccati.P
                )
                self.assertEqual(cells[-1].P_after, final_P)

    def test_lineage_is_explicit_across_two_prefixes(self):
        sample = SELECTORS._point_sample()
        frontend, P0_H, P0_A = self._source_seeded_frontend_and_covariance()
        endpoints, selectors, _ = SELECTORS.execute_with_prefix_selectors(
            frontend_entry=frontend,
            P0_H=P0_H,
            P0_A=P0_A,
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
            self.assertTrue(all(s.H_event_cells and s.A_event_cells for s in lineage))

    def test_broken_or_duplicate_ancestry_is_rejected(self):
        sample = SELECTORS._point_sample()
        frontend, P0_H, P0_A = self._source_seeded_frontend_and_covariance()
        endpoints, selectors, _ = SELECTORS.execute_with_prefix_selectors(
            frontend_entry=frontend,
            P0_H=P0_H,
            P0_A=P0_A,
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
