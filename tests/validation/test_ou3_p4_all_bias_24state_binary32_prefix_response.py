#!/usr/bin/env python3
from __future__ import annotations

import sys
import unittest
from pathlib import Path

TOOLS=Path(__file__).resolve().parents[2]/"tools"/"stability"
sys.path.insert(0,str(TOOLS))

import ou3_p4_all_bias_24state_binary32_prefix_response as FP


class Binary32PrefixResponseTest(unittest.TestCase):
    def test_contract_stays_conditional_and_fail_closed(self):
        d=FP.build()
        self.assertEqual(FP.validate(d),[])
        self.assertTrue(d["A21_every_prefix_binary32_response_materializable"])
        self.assertTrue(d["conditional_platform_execution_premise_consumed"])
        self.assertTrue(d["previous_roundoff_inputs_suffix_propagated_through_later_event_Jacobians"])
        self.assertFalse(d["deployment_toolchain_qualification_required_for_mathematical_P4"])
        self.assertFalse(d["deployment_toolchain_qualification_closed_here"])
        self.assertFalse(d["proof_true_bias_coordinate_receives_roundoff"])
        self.assertFalse(d["packet_count_times_worst_roundoff_used"])
        self.assertFalse(d["P4_PASS"])

    def test_injection_has_no_true_bias_rows(self):
        G=FP._fp_injection(0.25)
        self.assertEqual((len(G),len(G[0])),(24,21))
        for i in range(21):
            self.assertEqual(G[i][i].lo,0.25)
            self.assertEqual(G[i][i].hi,0.25)
        for row in G[21:24]:
            self.assertTrue(all(x.lo==0.0 and x.hi==0.0 for x in row))

    def test_suffix_propagation_appends_new_event_block(self):
        I=FP.I
        A=[[I(1.0 if i==j else 0.0) for j in range(24)] for i in range(24)]
        C=FP._zero(24,0)
        C1=FP._propagate(A,C,0.1)
        C2=FP._propagate(A,C1,0.2)
        self.assertEqual(FP._shape(C1),(24,21))
        self.assertEqual(FP._shape(C2),(24,42))
        self.assertEqual(C2[0][0].lo,0.1)
        self.assertEqual(C2[0][21].lo,0.2)


if __name__=="__main__":
    unittest.main()
