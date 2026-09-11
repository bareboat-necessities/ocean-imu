#!/usr/bin/env python3
from __future__ import annotations
import sys,unittest
from pathlib import Path
TOOLS=Path(__file__).resolve().parents[2]/'tools'/'stability';sys.path.insert(0,str(TOOLS))
import ou3_p4_a21_sequential_every_prefix_search as P
import ou3_p4_bias_family_joint_iss_supply as BIAS

class TestA21SequentialEveryPrefixSearch(unittest.TestCase):
    def test_all_families_all_prefixes_execute_fail_closed(self):
        d=P.build();self.assertEqual(P.validate(d),[])
        self.assertEqual(set(d['family_prefix_search']),set(BIAS.REQUIRED_BIAS_FAMILIES))
        self.assertTrue(d['preterminal_prefixes_do_not_double_charge_future_prediction_source'])
        for row in d['family_prefix_search'].values():
            self.assertGreater(row['prefix_count'],0)
            self.assertEqual(len(row['prefixes']),row['prefix_count'])
        self.assertFalse(d['production_every_prefix_augmented_LDLT_closed_here'])
        self.assertFalse(d['P4_PASS'])
    def test_Gamma_must_be_positive_finite(self):
        q=P.zero(P.NP);j=P.zero(P.NE)
        with self.assertRaises(ValueError):P.gain_margin(q,j,0.0)
        with self.assertRaises(ValueError):P.gain_margin(q,j,float('inf'))
if __name__=='__main__':unittest.main()
