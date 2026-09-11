import sys
from pathlib import Path
import unittest
ROOT=Path(__file__).resolve().parents[2];sys.path[:0]=[str(ROOT),str(ROOT/'tools/stability')]
from ou3_interval import Interval
import ou3_brmm_complete_window_execution_kernel as K
import ou3_p4_complete_brmm_source_cover_contract as C
import ou3_p4_complete_brmm_differential_events as E
from tools.stability.ou3_alt_contraction import physical_events as P

def I(x):return Interval.point(float(x))
def eye(n):return [[I(1 if i==j else 0) for j in range(n)] for i in range(n)]
def primitive(s=(0,0,0)):
    v=tuple(I(x) for x in (0,0,0));p=v;S=tuple(I(x) for x in s)
    return K.WavePrimitivePayload('generator','p0','p1','live0',v,p,S,p,p,p)
def cell(mode='H',s=(0,0,0)):
    n=18 if mode=='H' else 21
    return C.SourceCoverCell('est:e0','est','H' if mode=='H' else 'A',0,0,'S_zero',[I(0) for _ in range(n)],eye(n),I(.005),I(1.5),I(.4),I(.1),R=eye(3),R_provenance=E.ACTUAL_RS_PROVENANCE,true_bias=[I(0),I(0),I(0)] if mode=='A' else None,bias_projection_limit=.4 if mode=='A' else None,radial_scale=I(1),estimator_source_token='est',estimator_predecessor_token='root',estimator_generated_coefficients=True,wave_primitive=primitive(s))

class PhysicalEventTests(unittest.TestCase):
    def test_nonzero_physical_S_changes_actual_finite_event(self):
        zero=cell('H',(0,0,0));nonzero=cell('H',(2,-1,.5))
        old=E.source_joseph_event(**C.joseph_event_kwargs(nonzero));new=P.s_zero_event(nonzero)
        self.assertTrue(any(not (x.lo==0 and x.hi==0) for x in new['residual']))
        self.assertNotEqual([(x.lo,x.hi) for x in old['state_out']],[(x.lo,x.hi) for x in new['state_out']])
        parity=P.s_zero_event(zero);baseline=E.source_joseph_event(**C.joseph_event_kwargs(zero))['state_out']
        self.assertLess(max(max(abs(x.lo-y.lo),abs(x.hi-y.hi)) for x,y in zip(baseline,parity['state_out'])),1e-320)
    def test_A21_projection_is_composed(self):
        d=P.s_zero_event(cell('A',(1,0,0)));self.assertEqual(len(d['state_out']),21);self.assertEqual(len(d['J_state']),21);self.assertTrue(d['same_history_physical_reference_attached'])
    def test_missing_or_detached_physical_ancestry_fails(self):
        c=cell();object.__setattr__(c,'wave_primitive',None)
        with self.assertRaises(ValueError):P.s_zero_event(c)
        c=cell();object.__setattr__(c,'R_provenance','fake')
        with self.assertRaises(ValueError):P.s_zero_event(c)
    def test_contract_remains_fail_closed(self):
        d=P.build();self.assertEqual(P.validate(d),[]);self.assertFalse(d['complete_physical_word_closed']);self.assertFalse(d['ALT_LIVE_PASS'])
if __name__=='__main__':unittest.main()
