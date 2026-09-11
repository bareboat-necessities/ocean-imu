#!/usr/bin/env python3
from __future__ import annotations
import sys
import unittest
from pathlib import Path

TOOLS=Path(__file__).resolve().parents[2]/'tools'/'stability'
sys.path.insert(0,str(TOOLS))
import ou3_p4_h18_sequential_word_search as SEARCH
import ou3_p4_h18_source_indexed_prefix_transport as HTR

class H18SequentialWordSearchTest(unittest.TestCase):
    def test_search_runs_fail_closed(self):
        d=SEARCH.build()
        self.assertEqual(SEARCH.validate(d),[])
        self.assertTrue(d['actual_post_prediction_shipping_word_consumed'])
        self.assertTrue(d['coupled_a0_a1_J012_sectors_consumed_at_prediction'])
        self.assertTrue(d['conditional_binary32_event_channel_consumed'])
        self.assertGreater(d['search']['attempts'],0)
        self.assertFalse(d['production_endpoint_augmented_LDLT_closed_here'])
        self.assertFalse(d['production_every_prefix_augmented_LDLT_closed_here'])
        self.assertFalse(d['P4_PASS'])

    def test_invalid_supply_or_multiplier_rejected(self):
        samples=HTR.smoke_objects()
        with self.assertRaises(ValueError):
            SEARCH.backward_required_form(samples,gamma_s=-1,gamma_n=1,multipliers=(1,1,1))
        with self.assertRaises(ValueError):
            SEARCH.backward_required_form(samples,gamma_s=1,gamma_n=0,multipliers=(1,1,1))
        with self.assertRaises(ValueError):
            SEARCH.backward_required_form(samples,gamma_s=1,gamma_n=1,multipliers=(-1,1,1))

    def test_rho_must_be_strict(self):
        q=SEARCH.zero(SEARCH.NP)
        j=SEARCH.zero(SEARCH.NX)
        with self.assertRaises(ValueError): SEARCH.entrance_margin(q,j,1.0)
        with self.assertRaises(ValueError): SEARCH.entrance_margin(q,j,0.0)

if __name__=='__main__':
    unittest.main()
