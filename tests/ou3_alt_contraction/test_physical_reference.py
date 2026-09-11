"""Theorem algebra regressions for ALT physical-reference forcing."""
import sys,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT),str(ROOT/'tools/stability')]
from ou3_interval import Interval
import ou3_brmm_complete_window_execution_kernel as K
from tools.stability.ou3_alt_contraction import physical_reference as P
I=Interval.point

def primitive(S=(7,-2,1)):
    z=(I(0),I(0),I(0));return K.WavePrimitivePayload('generator','qin','qout','one-live-origin',z,z,tuple(map(I,S)),z,z,z)

class PhysicalReferenceTests(unittest.TestCase):
    def test_s_zero_is_error_minus_same_physical_S(self):
        r=P.s_zero_residual((I(8),I(-1),I(4)),primitive())
        self.assertTrue(all(x.contains(v) for x,v in zip(r,[1,1,3])))
        M=P.s_zero_residual_operator();q=[I(8),I(-1),I(4),I(7),I(-2),I(1)]
        mr=P._matvec(M,q);self.assertTrue(all(mr[i].contains([1,1,3][i]) for i in range(3)))
    def test_bias_error_uses_one_truth_history(self):
        e=(I(.1),I(-.2),I(.3));b0=(I(.04),I(-.03),I(.02));b1=(I(.041),I(-.029),I(.018));ph=I(.99)
        direct=P.bias_error_prediction(e,b0,b1,ph);M=P.bias_prediction_operator(ph);lift=P._matvec(M,[*e,*b0,*b1])
        self.assertTrue(all(direct[i].contains_interval(lift[i]) or lift[i].contains_interval(direct[i]) for i in range(3)))
    def test_H18_held_bias_is_not_dropped(self):
        r=P.h18_accelerometer_residual((I(1),I(2),I(3)),(I(.1),I(-.2),I(.3)))
        self.assertTrue(all(x.contains(v) for x,v in zip(r,[1.1,1.8,3.3])))
    def test_shared_physical_forcing_and_sector_are_consumed(self):
        d=P.build();self.assertEqual(P.validate(d),[])
        self.assertTrue(d['shared_exact_BRMM_vs_OU_prediction_forcing_consumed'])
        self.assertTrue(d['shared_joint_15D_acceleration_witness_sector_consumed'])
        self.assertTrue(d['shared_one_time_S_origin_prefix_invariant_consumed'])
        self.assertFalse(d['complete_word_master_closed']);self.assertFalse(d['ALT_LIVE_PASS'])
    def test_detached_or_missing_primitive_is_rejected(self):
        with self.assertRaises(TypeError):P.validate_primitive_reference(None)
        w=primitive();object.__setattr__(w,'live_origin_id','')
        with self.assertRaises(ValueError):P.validate_primitive_reference(w)
if __name__=='__main__':unittest.main()
