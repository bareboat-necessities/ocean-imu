from __future__ import annotations

import pathlib
import sys
import unittest

ROOT=pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools'/'stability'))

from ou3_interval import Interval
import ou3_p4_all_bias_24state_prefix_cocycle as PREFIX
import ou3_p4_complete_brmm_same_history_prefix_selectors as SELECTORS

I=Interval.point


class AllBias24StatePrefixCocycleTests(unittest.TestCase):
    def _fixture(self):
        frontend=SELECTORS.FRONTEND._point_state()
        P0_H,P0_A,_=SELECTORS._live_structured_point_covariance_fixture(frontend,PREFIX.GRAPH.DEFAULT_DOMAIN)
        sample=SELECTORS._point_sample()
        endpoints,selectors,_=SELECTORS.execute_with_prefix_selectors(
            frontend_entry=frontend,P0_H=P0_H,P0_A=P0_A,
            samples=[sample,sample],domain_path=PREFIX.GRAPH.DEFAULT_DOMAIN,branch_limit=128)
        return selectors,endpoints[0].source_cell_id

    def test_status_is_prefix_complete_and_nonpromoting(self):
        d=PREFIX.build()
        self.assertEqual(PREFIX.validate(d),[])
        self.assertTrue(d['every_literal_A21_event_prefix_snapshot_available'])
        self.assertTrue(d['all_previous_supply_blocks_suffix_propagated_through_later_events'])
        self.assertFalse(d['every_prefix_augmented_LDLT_closed_here'])
        self.assertFalse(d['P4_PASS'])

    def test_every_event_is_retained_and_endpoint_matches(self):
        selectors,endpoint=self._fixture()
        state=[I(0.0)]*21
        for family in ('BIAS0','BIAS1','BIAS2'):
            p=PREFIX.materialize_A21_joint_prefixes(selectors,endpoint,family=family,initial_error_state=state)
            lineage=SELECTORS.lineage_for_endpoint(selectors,endpoint)
            expected=sum(len(s.A_event_cells) for s in lineage)
            self.assertEqual(len(p),expected)
            self.assertEqual([x.literal_prefix_ordinal for x in p],list(range(1,expected+1)))
            self.assertEqual(p[-1].prediction_count,2)
            self.assertEqual(len(p[-1].A_prefix),24)
            self.assertEqual(len(p[-1].B_prefix[0]),12)
            for x in p:
                self.assertEqual(len(x.A_prefix),24)
                self.assertEqual(len(x.A_prefix[0]),24)
                self.assertEqual(len(x.B_prefix[0]),6*x.prediction_count)
                if x.kind=='S_zero': self.assertTrue(x.actual_rs_provenance)
                if x.kind in ('S_zero','accelerometer','magnetometer'):
                    self.assertNotEqual(x.projection_branch,'not_applicable')

    def test_unknown_family_and_bad_state_rejected(self):
        selectors,endpoint=self._fixture()
        with self.assertRaises(ValueError):
            PREFIX.materialize_A21_joint_prefixes(selectors,endpoint,family='BAD',initial_error_state=[I(0)]*21)
        with self.assertRaises(ValueError):
            PREFIX.materialize_A21_joint_prefixes(selectors,endpoint,family='BIAS0',initial_error_state=[I(0)]*18)


if __name__=='__main__': unittest.main()
