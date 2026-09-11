#!/usr/bin/env python3
from __future__ import annotations
import sys,unittest
from pathlib import Path
TOOLS=Path(__file__).resolve().parents[2]/'tools'/'stability';sys.path.insert(0,str(TOOLS))
import ou3_p4_a21_sequential_source_bias_fp_schur as S

class TestA21SequentialSchur(unittest.TestCase):
    def test_operator_contract(self):
        d=S.build();self.assertEqual(S.validate(d),[])
        self.assertEqual(d['persistent_dimension'],29)
        self.assertTrue(d['same_true_bias_state_carried_through_projection'])
        self.assertTrue(d['same_w_mtau_6D_supply_map_consumed'])
        self.assertTrue(d['physical_and_bias_source_scales_distinct'])
        self.assertFalse(d['P4_PASS'])
    def test_negative_multiplier_rejected(self):
        Q=S.zero(S.NP,S.NP);T=S.ordinary_transition(S.eye(S.NJOINT),1e-6)
        sec=S.ordinary_sectors(len(T[0]),.4)
        with self.assertRaises(ValueError):
            S.local_master(Q,T,gamma_phys=0,gamma_bias=0,gamma_n=1,sectors=sec,multipliers=(-1,),prediction=False)
if __name__=='__main__':unittest.main()
