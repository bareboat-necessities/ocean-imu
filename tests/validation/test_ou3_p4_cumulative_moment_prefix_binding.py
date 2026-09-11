#!/usr/bin/env python3
from __future__ import annotations
import sys,unittest
from pathlib import Path
TOOLS=Path(__file__).resolve().parents[2]/'tools'/'stability';sys.path.insert(0,str(TOOLS))
import ou3_p4_joint_brmm_prefix_ldlt as P

I=P.Interval.point

def mat(r,c):return [[I(0.0) for _ in range(c)] for _ in range(r)]

class CumulativeMomentPrefixBindingTest(unittest.TestCase):
    def binding(self,w,n=7):
        return P.MomentSectorBinding(w,mat(9,n),mat(1,n),1.0)
    def prefix(self,bindings,n=7):
        return P.PrefixInput(None,None,mat(n,n),[mat(n,n)],[0.0],source_map=mat(1,n),fp_map=mat(1,n),moment_sector_bindings=bindings)
    def test_multiple_distinct_transition_sectors_are_retained(self):
        q=self.prefix([self.binding('tr0'),self.binding('tr1')])
        self.assertEqual(P.moment_binding_count(q),2)
        self.assertEqual(P._moment_binding_failures(q,7,0),[])
    def test_duplicate_transition_witness_is_rejected(self):
        q=self.prefix([self.binding('tr0'),self.binding('tr0')])
        self.assertTrue(any('duplicate moment witness' in x for x in P._moment_binding_failures(q,7,0)))
    def test_legacy_and_cumulative_forms_cannot_be_mixed(self):
        q=P.PrefixInput(None,None,mat(7,7),[mat(7,7)],[0.0],source_map=mat(1,7),fp_map=mat(1,7),moment_sector_bindings=[self.binding('tr0')],normalized_acceleration_moment_map=mat(9,7),source_radial_scale_map=mat(1,7),moment_multiplier=1.0)
        self.assertTrue(any('both supplied' in x for x in P._moment_binding_failures(q,7,0)))
    def test_legacy_single_transition_remains_accepted_for_diagnostics(self):
        q=P.PrefixInput(None,None,mat(7,7),[mat(7,7)],[0.0],source_map=mat(1,7),fp_map=mat(1,7),normalized_acceleration_moment_map=mat(9,7),source_radial_scale_map=mat(1,7),moment_multiplier=1.0)
        self.assertEqual(P.moment_binding_count(q),1)
        self.assertEqual(P._moment_binding_failures(q,7,0),[])

if __name__=='__main__':unittest.main()
