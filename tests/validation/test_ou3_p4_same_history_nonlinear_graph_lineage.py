#!/usr/bin/env python3
import sys
from dataclasses import replace
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "tools" / "stability"))

from ou3_interval import Interval
import ou3_p4_complete_brmm_same_history_prefix_selectors as SELECTORS
import ou3_p4_same_history_nonlinear_graph_lineage as GRAPH
import ou3_brmm_frontend_state_step as FRONTEND


class P4SameHistoryNonlinearGraphLineageTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.payload = GRAPH.build()

    @staticmethod
    def _fixture():
        frontend = FRONTEND._point_state()
        P0_H, P0_A, _ = SELECTORS._live_structured_point_covariance_fixture(
            frontend, SELECTORS.DEFAULT_DOMAIN
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
        lineage = SELECTORS.lineage_for_endpoint(selectors, endpoint)
        bias = GRAPH._zero_bias_lineage(lineage)
        return endpoint, selectors, lineage, bias

    def test_status_is_connected_and_nonpromoting(self):
        self.assertEqual(GRAPH.validate(self.payload), [])
        self.assertEqual(
            self.payload["canonical_source"], "COMPLETE_BRMM_NORMAL_LIVE_WORD"
        )
        self.assertEqual(self.payload["P3_delta_preserved"], 1.0e-18)
        self.assertTrue(self.payload["event_local_same_P_H_R_cells_consumed"])
        self.assertTrue(self.payload["actual_applied_RS_provenance_retained"])
        self.assertTrue(self.payload["A21_absolute_bias_history_required"])
        self.assertFalse(self.payload["independent_true_bias_event_boxes_allowed"])
        self.assertFalse(self.payload["source_uniform_BRMM_window_family_materialized_here"])
        self.assertFalse(self.payload["joint_graph_sectors_assembled_here"])
        self.assertFalse(self.payload["endpoint_augmented_LDLT_closed_here"])
        self.assertFalse(self.payload["P4_promoted_here"])

    def test_two_prefix_zero_point_composes_H18_A21_and_projection(self):
        endpoint, selectors, _, bias = self._fixture()
        H, A = GRAPH.consume_endpoint_lineage(
            selectors,
            endpoint,
            initial_H_state=[Interval.point(0.0) for _ in range(18)],
            initial_A_state=[Interval.point(0.0) for _ in range(21)],
            A21_bias_lineage=bias,
        )
        self.assertEqual(H.prefixes, 2)
        self.assertEqual(A.prefixes, 2)
        self.assertEqual(len(H.J_word), 18)
        self.assertEqual(len(A.J_word), 21)
        self.assertEqual(H.source_token, A.source_token)
        self.assertEqual(
            [r.kind for r in H.event_records], [r.kind for r in A.event_records]
        )
        self.assertTrue(all(
            r.same_P_H_R_cell
            for r in H.event_records + A.event_records
            if r.kind in ("S_zero", "accelerometer", "magnetometer")
        ))
        self.assertTrue(all(
            r.actual_rs_provenance
            for r in H.event_records + A.event_records if r.kind == "S_zero"
        ))
        self.assertTrue(all(
            r.projection_branch == "inactive"
            for r in A.event_records
            if r.kind in ("S_zero", "accelerometer", "magnetometer")
        ))

    def test_A21_fails_closed_without_same_history_absolute_bias(self):
        endpoint, selectors, _, _ = self._fixture()
        with self.assertRaisesRegex(RuntimeError, "same-history absolute bias"):
            GRAPH._consume_mode_lineage(
                mode="A",
                lineage=SELECTORS.lineage_for_endpoint(selectors, endpoint),
                initial_state=[Interval.point(0.0) for _ in range(21)],
                source_token="COMPLETE_BRMM_NORMAL_LIVE_WORD:test",
                constants=GRAPH.KERNEL._process_constants(GRAPH.DEFAULT_DOMAIN),
                bias_lineage=None,
                projection_limit=0.4,
            )

    def test_bias_history_must_match_exact_endpoint_lineage(self):
        endpoint, selectors, lineage, bias = self._fixture()
        z = (Interval.point(0.0), Interval.point(0.0), Interval.point(0.0))
        missing = GRAPH.SameHistoryBiasLineage(
            endpoint_source_cell_id=endpoint,
            bias_true_by_source_cell_id={lineage[0].source_cell_id: z},
        )
        with self.assertRaisesRegex(RuntimeError, "every and only selector prefix"):
            GRAPH.consume_endpoint_lineage(
                selectors,
                endpoint,
                initial_H_state=[Interval.point(0.0) for _ in range(18)],
                initial_A_state=[Interval.point(0.0) for _ in range(21)],
                A21_bias_lineage=missing,
            )

        wrong_endpoint = GRAPH.SameHistoryBiasLineage(
            endpoint_source_cell_id="wrong-endpoint",
            bias_true_by_source_cell_id=bias.bias_true_by_source_cell_id,
        )
        with self.assertRaisesRegex(RuntimeError, "different endpoint lineage"):
            GRAPH.consume_endpoint_lineage(
                selectors,
                endpoint,
                initial_H_state=[Interval.point(0.0) for _ in range(18)],
                initial_A_state=[Interval.point(0.0) for _ in range(21)],
                A21_bias_lineage=wrong_endpoint,
            )

    def test_matching_prefix_ids_do_not_replace_the_GM_history(self):
        _, _, lineage, _ = self._fixture()
        root = (Interval.point(0.1), Interval.point(-0.05), Interval.point(0.02))
        bias = GRAPH.homogeneous_bias_lineage(lineage, root)
        constants = GRAPH.KERNEL._process_constants(GRAPH.DEFAULT_DOMAIN)
        bias.validate_homogeneous(lineage, constants, 0.4)
        detached = dict(bias.bias_true_by_source_cell_id)
        detached[lineage[-1].source_cell_id] = root
        with self.assertRaisesRegex(RuntimeError, "common GM root"):
            replace(bias, bias_true_by_source_cell_id=detached).validate_homogeneous(lineage, constants, 0.4)
        with self.assertRaisesRegex(RuntimeError, "tau mismatch"):
            replace(bias, tau_s=Interval.point(4000.0)).validate_homogeneous(lineage, constants, 0.4)
        outside = GRAPH.homogeneous_bias_lineage(lineage, (Interval.point(0.5),) * 3)
        with self.assertRaisesRegex(RuntimeError, "zero-error invariance"):
            outside.validate_homogeneous(lineage, constants, 0.4)


if __name__ == "__main__":
    unittest.main()
