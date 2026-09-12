"""Startup proxy-to-MEKF handoff seed/control regressions."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))

from tools.stability.ou3_alt_contraction import finite_startup_handoff_seed as X
from tools.stability.ou3_alt_contraction import finite_mag_startup_ready as R
from tools.stability.ou3_alt_contraction import finite_mag_tilt_frame as T
from tools.stability.ou3_alt_contraction import finite_vertical_complementary_runtime as V


def root(v): return T.SqrtWitness(v*v,v)


class Tests(unittest.TestCase):
    def test_ungauged_handoff_uses_proxy_directly_and_free_yaw_sigma(self):
        p=V.State((F(3,5),0,0,F(4,5)),(0,0,0),True,F(2),0)
        out=X.seed(p,None)
        self.assertEqual(out.q_seed,p.q); self.assertFalse(out.gauged)
        self.assertEqual(out.tilt_sigma,F(35,1000)); self.assertEqual(out.yaw_sigma,F(15708,10000))
        self.assertFalse(out.allow_acc_bias)
        with self.assertRaisesRegex(ValueError,'consumes no yaw-strip'):
            X.seed(p,None,proxy_q_norm=root(1))

    def test_gauged_handoff_replaces_proxy_yaw_but_preserves_tilt(self):
        # proxy = yaw(3/5,4/5) * roll(5/13,12/13)
        q_yaw=(F(3,5),0,0,F(4,5)); q_roll=(F(5,13),F(12,13),0,0)
        from tools.stability.ou3_alt_contraction import finite_prediction_graph as P
        q=tuple(P.quat_mul(q_yaw,q_roll))
        proxy=V.State(q,(0,0,0),True,F(2),0)
        proxy_yaw=T.YawHalfWitness(F(-7,25),F(24,25),root(1),F(3,5),F(4,5))
        pending=R.PendingYaw((1,0,0,0),(1,0))
        out=X.seed(proxy,pending,proxy_q_norm=root(1),proxy_yaw_half=proxy_yaw)
        self.assertEqual(out.q_seed,q_roll); self.assertTrue(out.gauged)
        self.assertEqual(out.yaw_sigma,F(87,1000)); self.assertFalse(out.allow_acc_bias)

    def test_nonzero_pending_absolute_yaw_composes_after_yaw_strip(self):
        p=V.State((1,0,0,0),(0,0,0),True,F(2),0)
        pending=R.PendingYaw((F(3,5),0,0,F(-4,5)),(F(-7),24))
        zero=T.YawHalfWitness(1,0,root(1),1,0)
        out=X.seed(p,pending,proxy_q_norm=root(1),proxy_yaw_half=zero)
        self.assertEqual(out.q_seed,pending.q_abs)

    def test_readiness_keeps_initialize_and_capture_open(self):
        r=X.readiness()
        self.assertTrue(r['gauged_yaw_sigma_0087_materialized'])
        self.assertTrue(r['handoff_allow_acc_bias_false'])
        self.assertFalse(r['initialize_from_attitude_state_covariance_attached'])
        self.assertFalse(r['startup_capture_closed']); self.assertFalse(r['ALT_STARTUP_PASS'])

if __name__=='__main__': unittest.main()
