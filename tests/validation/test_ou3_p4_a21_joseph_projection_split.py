#!/usr/bin/env python3
from __future__ import annotations
import sys,unittest
from pathlib import Path
TOOLS=Path(__file__).resolve().parents[2]/'tools'/'stability';sys.path.insert(0,str(TOOLS))
import ou3_p4_a21_joseph_projection_split as S
class A21JosephProjectionSplitTest(unittest.TestCase):
    def test_split_is_structural_and_nonpromoting(self):
        d=S.build();self.assertEqual(S.validate(d),[])
        self.assertTrue(d['A21_Joseph_projection_split_closed'])
        self.assertTrue(d['smooth_Joseph_C_minus_L_is_minus_GKH'])
        self.assertTrue(d['projection_C_equals_L_generalized_Jacobian'])
        self.assertTrue(d['same_physical_true_bias_source_coordinate_retained'])
        self.assertTrue(d['projection_defect_is_not_independent_source_port'])
        self.assertFalse(d['all_BIAS0_BIAS1_BIAS2_recurrences_composed_here'])
        self.assertFalse(d['A21_augmented_LDLT_closed_here'])
        self.assertFalse(d['P4_PASS'])
if __name__=='__main__':unittest.main()
