"""Startup goLive/enterLive exact real-arithmetic regressions."""
from fractions import Fraction as F
from pathlib import Path
import sys,unittest
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))

from tools.stability.ou3_alt_contraction import finite_startup_live_entry as X
from tools.stability.ou3_alt_contraction import finite_startup_handoff_init as I
from tools.stability.ou3_alt_contraction import finite_startup_handoff_seed as S
from tools.stability.ou3_alt_contraction import finite_runtime_parameters as A
from tools.stability.ou3_alt_contraction import deployment_scope as D


def active():
    sig=((F(1),0,0),(0,F(2),0),(0,0,F(3)))
    rs=((F(4),0,0),(0,F(5),0),(0,0,F(6)))
    return A.ActiveParameters(F(11,10),sig,F(3,200),rs)

def handoff(q=(1,0,0,0)):
    seed=S.Result(q,S.TILT_SIGMA,S.YAW_SIGMA_GAUGED,False,True)
    x=tuple(F(i+1,100) for i in range(21))
    P=[[F((i+1)*100+(j+1),100000) for j in range(21)] for i in range(21)]
    return I.initialize_from_gauged_seed_zero_heel(seed,x,P,scope=D.certified_scope())

class Tests(unittest.TestCase):
    def test_fresh_entry_is_H18_and_seats_aw_covariance(self):
        h=handoff();a=active();o=X.enter_live(h,a,scope=D.certified_scope())
        self.assertEqual(o.startup_stage,'Live');self.assertEqual(o.startup_stage_t,0)
        self.assertFalse(o.acc_bias_updates_enabled);self.assertTrue(o.gauged);self.assertTrue(o.zero_wind_heel)
        for i in range(21):
            if i<15 or i>=18:
                for j in range(15,18):
                    self.assertEqual(o.P[i][j],0);self.assertEqual(o.P[j][i],0)
        self.assertEqual(tuple(tuple(r[15:18]) for r in o.P[15:18]),a.Sigma_aw)
        self.assertEqual(o.x,h.x)

    def test_handoff_body_to_world_seed_is_stored_as_world_to_body_qhat(self):
        h=handoff((F(3,5),F(4,5),0,0))
        o=X.enter_live(h,active(),scope=D.certified_scope())
        self.assertEqual(o.q_hat,(F(3,5),F(-4,5),0,0))

    def test_other_covariance_entries_survive_aw_seating(self):
        h=handoff();o=X.enter_live(h,active(),scope=D.certified_scope())
        for i in range(21):
            for j in range(21):
                aw=(15<=i<18) or (15<=j<18)
                if not aw:self.assertEqual(o.P[i][j],h.P[i][j],(i,j))

    def test_requires_locked_gauged_zeroheel_and_live_RS(self):
        h=handoff()
        with self.assertRaisesRegex(ValueError,'lock retained'):
            X.enter_live(h,active(),scope=D.certified_scope(),accel_bias_locked=False)
        a=active();bad=A.ActiveParameters(a.tau,a.Sigma_aw,a.pseudo_period,None)
        with self.assertRaisesRegex(ValueError,'Live R_S'): X.enter_live(h,bad,scope=D.certified_scope())
        ung=I.Result(h.x,h.P,h.q_seed,h.tilt_sigma,S.YAW_SIGMA_FREE,False)
        with self.assertRaisesRegex(ValueError,'gauged'): X.enter_live(ung,a,scope=D.certified_scope())

    def test_readiness_remains_fail_closed(self):
        r=X.readiness();self.assertTrue(r['fresh_real_arithmetic_H18_entry_composed'])
        self.assertTrue(r['internal_world_to_body_handoff_quaternion_retained'])
        self.assertFalse(r['complete_same_history_startup_to_Live_word'])
        self.assertFalse(r['ALT_STARTUP_PASS']);self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__':unittest.main()
