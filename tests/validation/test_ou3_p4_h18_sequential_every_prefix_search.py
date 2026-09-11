#!/usr/bin/env python3
from __future__ import annotations
import sys,unittest
from pathlib import Path
TOOLS=Path(__file__).resolve().parents[2]/'tools'/'stability';sys.path.insert(0,str(TOOLS))
import ou3_p4_h18_sequential_every_prefix_search as P

class TestH18SequentialEveryPrefixSearch(unittest.TestCase):
    def test_every_prefix_search_executes_fail_closed(self):
        d=P.build();self.assertEqual(P.validate(d),[])
        self.assertTrue(d['finite_Gamma_search_executed_for_every_literal_prefix'])
        self.assertGreater(d['diagnostic']['prefix_count'],0)
        self.assertEqual(len(d['diagnostic']['prefixes']),d['diagnostic']['prefix_count'])
        self.assertFalse(d['production_every_prefix_augmented_LDLT_closed_here'])
        self.assertFalse(d['P4_PASS'])
    def test_gain_requires_positive_finite_Gamma(self):
        q=P.zero(P.NP);j=P.zero(P.NX)
        with self.assertRaises(ValueError):P.gain_margin(q,j,0.0)
        with self.assertRaises(ValueError):P.gain_margin(q,j,float('inf'))
if __name__=='__main__':unittest.main()
