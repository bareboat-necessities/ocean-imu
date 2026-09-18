"""Regression checks for real CI producer wiring; no production promotion."""
from __future__ import annotations

import copy
import dataclasses
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'tools' / 'stability'))
import ou3_p4_brmm_a21_projection_prefix as A21
import ou3_p4_brmm_bias_projection_prediction_chain as CHAIN
import ou3_p4_brmm_event_lineage_cover as LINEAGE
import ou3_p4_brmm_nonlinear_correction_domain as GRAPH
import ou3_p4_brmm_same_Dtheta_reset_binding as BIND


class CIEventPrefixTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cell = A21._smoke_cell()
        cls.event = A21.EVENT.build_event_master(cls.cell)
        cls.chart = GRAPH.certify_event(
            cls.event, delta_candidates=(1.5,), multiplier_grid=(1.0,))

    def test_actual_a21_graph_has_a_reproducible_strict_certificate(self):
        self.assertTrue(self.chart['closed'])
        cert = self.chart['certificate']
        ok, pivots = GRAPH.verify_event_certificate(self.event, cert)
        self.assertTrue(ok)
        self.assertGreater(min(pivots), 0)
        self.assertEqual(cert['delta'], 1.5)
        self.assertEqual(self.chart['hard_premise_count'], 7)
        self.assertEqual(self.chart['nonlinear_premise_count'], 11)
        self.assertEqual(self.chart['original_dimension'] - self.chart['reduced_dimension'], 3)
        weights = cert['premise_multipliers']
        self.assertNotEqual(weights['chord_orthogonality_plus'], weights['chord_orthogonality_minus'])
        self.assertNotEqual(weights['integral_displacement_norm_m_s'], weights['attitude_cayley_norm'])

    def test_stored_pivots_do_not_override_a_failed_matrix_certificate(self):
        cert = copy.deepcopy(self.chart['certificate'])
        cert['delta'] = 0.25
        cert['pivot_lowers'] = [1e9] * len(cert['pivot_lowers'])
        with self.assertRaisesRegex(ValueError, 'failed outward LDLT'):
            BIND.bind_event_certified_graph(self.event, cert)
        amplified = dataclasses.replace(
            self.event, D_theta=[[A21.I(100)*x for x in row] for row in self.event.D_theta])
        with self.assertRaisesRegex(ValueError, 'failed outward LDLT'):
            BIND.bind_event_certified_graph(amplified, self.chart['certificate'])

    def test_detached_tokens_and_missing_premises_are_rejected(self):
        for key in ('event_source_token', 'estimator_source_token'):
            with self.subTest(key=key):
                cert = copy.deepcopy(self.chart['certificate'])
                cert[key] = 'detached'
                with self.assertRaisesRegex(ValueError, 'token detached'):
                    BIND.bind_event_certified_graph(self.event, cert)
        cert = copy.deepcopy(self.chart['certificate'])
        del cert['premise_multipliers']['attitude_cayley_norm']
        with self.assertRaisesRegex(ValueError, 'premise set changed'):
            BIND.bind_event_certified_graph(self.event, cert)

    def test_local_certificate_does_not_widen_the_canonical_first_exit_chart(self):
        bound = BIND.bind_event_certified_graph(self.event, self.chart['certificate'])
        self.assertTrue(bound['closed'])
        self.assertEqual(bound['declared_correction_chart_delta'], 0.25)
        self.assertFalse(bound['within_declared_correction_chart'])
        self.assertFalse(bound['production_correction_domain_target_closed_here'])
        with self.assertRaisesRegex(ValueError, 'outside declared'):
            BIND.bind_event_first_exit(self.event, 1.5)
        first_exit = BIND.bind_event_first_exit(self.event)
        self.assertEqual(first_exit['delta'], 0.25)
        self.assertFalse(first_exit['correction_domain_target_proved_here'])

    def test_empty_or_invalid_candidate_search_cannot_promote(self):
        result = GRAPH.certify_event(self.event, delta_candidates=(0.25,), multiplier_grid=(0.0,))
        self.assertFalse(result['closed'])
        self.assertIsNone(result['certificate'])
        for scale in (-1.0, float('nan'), float('inf')):
            with self.subTest(scale=scale), self.assertRaises(ValueError):
                GRAPH.certify_event(self.event, multiplier_grid=(scale,))

    def test_joseph_fixture_tokens_continue_the_literal_event_namespace(self):
        radial = LINEAGE.root_radial()
        image, initial = LINEAGE._smoke_lineage(radial)
        szero, accel = LINEAGE._smoke_joseph_cells(image, radial)
        self.assertEqual(szero.source_token, image.source_token + ':e2')
        self.assertEqual(accel.source_token, image.source_token + ':e3')
        self.assertEqual(szero.predecessor_token, initial.cells[-1].source_token)
        self.assertEqual(accel.predecessor_token, szero.source_token)
        self.assertEqual((szero.event_ordinal, accel.event_ordinal), (2, 3))
        for cell in (szero, accel):
            self.assertEqual(LINEAGE.COVER.validate_cell(cell, require_estimator_provenance=True), [])
            prefix = LINEAGE.certify_joseph_prefix(cell, radial)
            self.assertEqual(LINEAGE.validate_certified_prefix(prefix, cell, radial), [])
            self.assertEqual(prefix.within_declared_correction_chart, cell.kind == 'S_zero')

    def test_successor_carries_the_full_outward_covariance_image(self):
        _, lineage = LINEAGE._smoke_lineage(LINEAGE.root_radial())
        before, after = lineage.cells
        F, Q = LINEAGE._identity(18), LINEAGE._zeros(18)
        self.assertTrue(LINEAGE.COV.prediction_transition_closed(before, after, F, Q))
        stale = dataclasses.replace(after, P=before.P)
        self.assertFalse(LINEAGE.COV.prediction_transition_closed(before, stale, F, Q))

    def test_bias_chain_uses_the_current_same_graph_prefix_api(self):
        report = CHAIN.build()
        self.assertEqual(CHAIN.validate(report), [])
        self.assertTrue(report['projected_ba_is_next_prediction_error_input'])
        self.assertTrue(report['same_persistent_beta_is_next_prediction_true_bias_input'])
        self.assertFalse(report['production_full_24state_prefix_transport_closed_here'])
        self.assertFalse(report['P4_PASS'])
        self.assertFalse(report['production_correction_chart_qualified_here'])
        self.assertFalse(report['local_projection_within_declared_correction_chart'])


if __name__ == '__main__':
    unittest.main()
