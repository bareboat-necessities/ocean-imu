#!/usr/bin/env python3
from __future__ import annotations
import sys,unittest
from pathlib import Path
TOOLS=Path(__file__).resolve().parents[2]/'tools'/'stability';sys.path.insert(0,str(TOOLS))
import ou3_p4_a21_post_prediction_word as W
import ou3_p4_a21_source_indexed_prefix_transport as A21
import ou3_p4_source_uniform_bias_prefix_lineage as BIAS

class A21PostPredictionWordTest(unittest.TestCase):
    def test_all_families_share_post_prediction_boundary_topology(self):
        d=W.build();self.assertEqual(W.validate(d),[])
        self.assertTrue(d['all_BIAS0_BIAS1_BIAS2_post_prediction_words_materialized'])
        self.assertTrue(d['A21_post_prediction_word_transport_closed'])
        self.assertTrue(d['first_prediction_driver_not_double_counted'])
        self.assertEqual(d['metric_comparison_mu'],1.0)
        self.assertFalse(d['endpoint_augmented_LDLT_closed_here'])
        self.assertFalse(d['P4_PASS']);self.assertFalse(d['P5_MAY_START'])
    def test_word_supply_restarts_after_entrance_prediction(self):
        for family in BIAS.FAMILIES:
            w=W.build_word(A21.smoke_objects(family),family)
            self.assertEqual(w['events'][-1]['kind'],'prediction')
            self.assertTrue(w['first_prediction_supply_absorbed_into_boundary_state'])
            self.assertTrue(w['one_new_6D_supply_block_per_in_word_prediction'])
            self.assertEqual(len(w['source_prefix_maps']),len(w['events']))
            # Before the following prediction no new bias driver belongs to this word.
            first_pred=next(i for i,e in enumerate(w['events']) if e['kind']=='prediction')
            self.assertEqual(first_pred,len(w['events'])-1)
            for B in w['source_prefix_maps'][:first_pred]:self.assertEqual(A21.shape(B)[1],0)
            self.assertEqual(A21.shape(w['source_prefix_maps'][first_pred])[1],6)
    def test_projection_is_literal_prefix_and_covariance_identity_event(self):
        w=W.build_word(A21.smoke_objects('BIAS2'),'BIAS2')
        kinds=[e['kind'] for e in w['events']]
        self.assertIn('bias_projection',kinds)
        for i,k in enumerate(kinds):
            if k=='bias_projection':self.assertEqual(w['covariance_nodes'][i],w['covariance_nodes'][i+1])

if __name__=='__main__':unittest.main()
