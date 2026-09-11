#!/usr/bin/env python3
from __future__ import annotations
import sys,unittest
from pathlib import Path
TOOLS=Path(__file__).resolve().parents[2]/'tools'/'stability';sys.path.insert(0,str(TOOLS))
import ou3_p4_all_bias_correction_chart_cover as C
import ou3_p4_bias_family_joint_iss_supply as BIAS

class TestAllBiasCorrectionChartCover(unittest.TestCase):
    def test_all_families_materialized_fail_closed(self):
        d=C.build();self.assertEqual(C.validate(d),[])
        self.assertTrue(d['BIAS0_BIAS1_BIAS2_certified_separately'])
        self.assertTrue(d['shared_component_envelope_lemma_consumed'])
        self.assertFalse(d['family_driver_recurrences_collapsed_into_one'])
        self.assertFalse(d['BIAS0_or_BIAS2_inferred_from_BIAS1'])
        self.assertGreater(d['record_count'],0)
        self.assertEqual(d['event_count'],3*d['record_count'])
        self.assertFalse(d['P4_PASS'])
    def test_family_boxes_use_authoritative_component_envelopes(self):
        b=BIAS.build()
        for family in BIAS.REQUIRED_BIAS_FAMILIES:
            box=C.family_box(family,b);self.assertEqual(len(box),3)
            cap=b['true_bias_component_envelope_mps2'][family]
            self.assertTrue(all(x.lo<=-cap and x.hi>=cap for x in box))
if __name__=='__main__':unittest.main()
