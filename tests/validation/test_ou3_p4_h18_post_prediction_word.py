#!/usr/bin/env python3
from __future__ import annotations
import sys,unittest
from pathlib import Path
TOOLS=Path(__file__).resolve().parents[2]/'tools'/'stability';sys.path.insert(0,str(TOOLS))
import ou3_p4_h18_post_prediction_word as W
import ou3_p4_h18_source_indexed_prefix_transport as HTR

class H18PostPredictionWordTest(unittest.TestCase):
    def test_storage_and_transport_share_shipping_prediction_boundaries(self):
        d=W.build();self.assertEqual(W.validate(d),[])
        self.assertTrue(d['H18_post_prediction_word_transport_closed'])
        self.assertEqual(d['metric_comparison_mu'],1.0)
        self.assertTrue(d['word_entrance_is_actual_shipping_post_prediction_boundary'])
        self.assertTrue(d['word_endpoint_is_actual_following_prediction_boundary'])
        self.assertFalse(d['P4_PASS']);self.assertFalse(d['P5_MAY_START'])
    def test_smoke_slice_stops_at_following_prediction(self):
        w=W.build_word(HTR.smoke_objects())
        self.assertGreater(len(w['events']),0)
        self.assertNotEqual(w['events'][0]['kind'],'prediction')
        self.assertEqual(w['events'][-1]['kind'],'prediction')
        self.assertEqual(len(w['covariance_nodes']),len(w['events'])+1)
        self.assertEqual(w['entrance_P'],w['covariance_nodes'][0])
        self.assertEqual(w['terminal_P'],w['covariance_nodes'][-1])
        self.assertTrue(w['all_prefix_identity_residuals_contain_zero'])
    def test_one_sample_cannot_fake_complete_boundary_word(self):
        with self.assertRaises(ValueError):W.build_word(HTR.smoke_objects()[:1])

if __name__=='__main__':unittest.main()
