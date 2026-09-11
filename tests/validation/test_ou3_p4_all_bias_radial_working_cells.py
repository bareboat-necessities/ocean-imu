#!/usr/bin/env python3
from __future__ import annotations
import sys,unittest
from pathlib import Path
TOOLS=Path(__file__).resolve().parents[2]/'tools'/'stability';sys.path.insert(0,str(TOOLS))
import ou3_p4_all_bias_radial_working_cells as W
import ou3_p4_source_uniform_bias_prefix_lineage as BIAS

class TestAllBiasRadialWorkingCells(unittest.TestCase):
    def test_all_bias_cells_materialized_fail_closed(self):
        d=W.build();self.assertEqual(W.validate(d),[])
        self.assertTrue(d['A21_BIAS0_BIAS1_BIAS2_cells_materialized_separately'])
        self.assertEqual(set(d['A21_family_counts']),set(BIAS.FAMILIES))
        self.assertEqual(len(set(d['A21_family_counts'].values())),1)
        self.assertTrue(d['all_Szero_cells_use_actual_applied_RS'])
        self.assertFalse(d['captured_source_image_promoted_to_COMPLETE_BRMM_family'])
        self.assertFalse(d['P4_PASS'])
    def test_bias_boxes_are_family_qualified(self):
        for family in BIAS.FAMILIES:
            box=W.true_bias_box(family);self.assertEqual(len(box),3)
            cap=BIAS.family_parameters(family)['cap']
            self.assertTrue(all(x.lo<=-cap and x.hi>=cap for x in box))
if __name__=='__main__':unittest.main()
