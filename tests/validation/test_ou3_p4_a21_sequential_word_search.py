#!/usr/bin/env python3
from __future__ import annotations
import sys,unittest
from pathlib import Path
TOOLS=Path(__file__).resolve().parents[2]/'tools'/'stability';sys.path.insert(0,str(TOOLS))
import ou3_p4_a21_sequential_word_search as S
import ou3_p4_bias_family_joint_iss_supply as BIAS

class TestA21SequentialWordSearch(unittest.TestCase):
    def test_all_families_search_fail_closed(self):
        d=S.build();self.assertEqual(S.validate(d),[])
        self.assertEqual(set(d['family_search']),set(BIAS.REQUIRED_BIAS_FAMILIES))
        self.assertTrue(d['local_Kuu_and_persistent_supply_search_staged_without_changing_inequality'])
        for row in d['family_search'].values():
            self.assertGreater(row['attempts'],0)
            self.assertGreaterEqual(row['local_feasible_count'],0)
        self.assertFalse(d['production_endpoint_augmented_LDLT_closed_here'])
        self.assertFalse(d['P4_PASS'])

    def test_invalid_values_rejected(self):
        samples=__import__('ou3_p4_a21_source_indexed_prefix_transport').smoke_objects('BIAS1')
        with self.assertRaises(ValueError):
            S.backward_required_form(samples,'BIAS1',gamma_phys=-1,gamma_bias=1,gamma_n=1,multipliers=(1,1,1,1,1))
        with self.assertRaises(ValueError):
            S.backward_required_form(samples,'BIAS1',gamma_phys=1,gamma_bias=1,gamma_n=0,multipliers=(1,1,1,1,1))
        with self.assertRaises(ValueError):
            S.entrance_margin(S.zero(S.NP),S.zero(S.NE),1.0)

if __name__=='__main__':unittest.main()
