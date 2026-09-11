import sys,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path[:0]=[str(ROOT),str(ROOT/'tools/stability')]
from ou3_interval import Interval
import ou3_brmm_complete_window_execution_kernel as K
import ou3_p4_complete_brmm_source_cover_contract as C
import ou3_p4_complete_brmm_differential_events as E
from tools.stability.ou3_alt_contraction import joint24_events as J
I=Interval.point
def eye(n):return [[I(1 if i==j else 0) for j in range(n)] for i in range(n)]
def prim(S=(1,0,0)):
 z=(I(0),)*3;return K.WavePrimitivePayload('g','p0','p1','live',z,z,tuple(map(I,S)),z,z,z)
def cell(mode,kind):
 n=18 if mode=='H' else 21
 kw=dict(R=eye(3),radial_scale=I(1),estimator_source_token='est',estimator_predecessor_token='root',estimator_generated_coefficients=True,wave_primitive=prim())
 if kind=='S_zero':kw['R_provenance']=E.ACTUAL_RS_PROVENANCE
 if kind=='accelerometer':kw.update(f_hat=[I(0),I(0),I(-9.8)],R_hat=eye(3))
 if mode=='A':kw.update(true_bias=[I(.1),I(0),I(0)],bias_projection_limit=.4)
 return C.SourceCoverCell('est:e0','est',mode,0,0,kind,[I(0) for _ in range(n)],eye(n),I(.005),I(1.5),I(.4),I(.1),**kw)
class Joint24Tests(unittest.TestCase):
 def test_H18_accel_has_held_bias_columns(self):
  d=J.h18_accelerometer_event(cell('H','accelerometer'));M=d['J_joint24'];self.assertEqual((len(M),len(M[0])),(24,24));self.assertTrue(d['same_history_held_bias_retained']);self.assertTrue(any(not (M[i][18+j].lo==0 and M[i][18+j].hi==0) for i in range(18) for j in range(3)))
 def test_H18_prediction_shares_one_physical_driver(self):
  A,B=J.h18_prediction_lift(eye(18),I(.99))
  for i in range(3):self.assertTrue(B[18+i][i].contains(1) and B[21+i][i].contains(1))
  self.assertTrue(any(not (A[18+i][21+i].lo==0 and A[18+i][21+i].hi==0) for i in range(3)))
 def test_physical_S_columns_exist_in_H18_and_A21(self):
  for mode in ('H','A'):
   d=J.s_zero_event_joint24(cell(mode,'S_zero'));self.assertEqual((len(d['J_joint24']),len(d['J_joint24'][0])),(24,24));self.assertEqual((len(d['J_S_phys']),len(d['J_S_phys'][0])),(24,3));self.assertTrue(any(not (x.lo==0 and x.hi==0) for r in d['J_S_phys'] for x in r))
 def test_contract_still_blocks_storage(self):
  d=J.build();self.assertEqual(J.validate(d),[]);self.assertFalse(d['storage_search_allowed']);self.assertFalse(d['ALT_LIVE_PASS'])
if __name__=='__main__':unittest.main()
