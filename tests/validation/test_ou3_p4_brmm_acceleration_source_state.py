#!/usr/bin/env python3
import sys,unittest
from pathlib import Path
TOOLS=Path(__file__).resolve().parents[2]/'tools'/'stability';sys.path.insert(0,str(TOOLS))
import ou3_p4_brmm_acceleration_source_state as S
class AccelerationSourceStateTest(unittest.TestCase):
    def test_contract(self):
        d=S.build();self.assertEqual(S.validate(d),[])
        self.assertTrue(d['a1_k_equals_a0_kplus1_enforced'])
        self.assertTrue(d['sequential_source_state_elimination_architecture_available'])
        self.assertFalse(d['physical_acceleration_endpoint_reselected_independently_each_prediction'])
        self.assertFalse(d['P4_PASS'])
    def test_detached_endpoint_rejected(self):
        I=S.I;z=(I(0),)*3;a=(I(.1),I(.2),I(.3));b=(I(.2),I(.1),I(.25));j=z
        good=(S.SourceStep('r','k0',z,a,j,j,j),S.SourceStep('k0','k1',a,b,j,j,j))
        self.assertEqual(S.validate_chain(good),[])
        bad=(good[0],S.SourceStep('k0','k1',z,b,j,j,j))
        self.assertIn('step 1: a1/a0 same-history continuity detached',S.validate_chain(bad))
    def test_maps_have_expected_shapes(self):
        Q=S.q_assembly_map();N=S.next_source_state_map();self.assertEqual((len(Q),len(Q[0])),(15,15));self.assertEqual((len(N),len(N[0])),(3,15))
if __name__=='__main__':unittest.main()
