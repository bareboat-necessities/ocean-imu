"""Literal startup MEKF initialize_from_attitude reset regressions."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))

from tools.stability.ou3_alt_contraction import finite_startup_handoff_init as X
from tools.stability.ou3_alt_contraction import finite_startup_handoff_seed as S
from tools.stability.ou3_alt_contraction import deployment_scope as DS


def seed(gauged=True,q=(1,0,0,0)):
    return S.Result(q,S.TILT_SIGMA,
                    S.YAW_SIGMA_GAUGED if gauged else S.YAW_SIGMA_FREE,
                    False,gauged)
def state(): return tuple(F(i+1,100) for i in range(21))
def cov(): return [[F((i+1)*100+(j+1),100000) for j in range(21)] for i in range(21)]

class Tests(unittest.TestCase):
    def test_down_z_gives_diagonal_tilt_tilt_yaw_split(self):
        s=seed(True); out=X.initialize(s,state(),cov(),down_body=X.UnitDownWitness((0,0,1)))
        tv=S.TILT_SIGMA*S.TILT_SIGMA; yv=S.YAW_SIGMA_GAUGED*S.YAW_SIGMA_GAUGED
        self.assertEqual([list(r[:3]) for r in out.P[:3]],[[tv,0,0],[0,tv,0],[0,0,yv]])
        self.assertEqual(out.x[:3],(0,0,0))

    def test_arbitrary_unit_axis_keeps_projector_formula_exact(self):
        u=X.UnitDownWitness((F(3,5),F(4,5),0)); s=seed(True)
        A=X.attitude_covariance(u,s.tilt_sigma,s.yaw_sigma)
        tv=s.tilt_sigma*s.tilt_sigma; yv=s.yaw_sigma*s.yaw_sigma
        self.assertEqual(A[0][0],tv*F(16,25)+yv*F(9,25))
        self.assertEqual(A[1][1],tv*F(9,25)+yv*F(16,25))
        self.assertEqual(A[0][1],(yv-tv)*F(12,25)); self.assertEqual(A[2][2],tv)

    def test_zero_heel_axis_is_directly_derived_from_unit_boat_quaternion(self):
        scope=DS.certified_scope()
        self.assertEqual(X.down_from_unit_boat_quaternion_zero_heel((1,0,0,0),scope=scope).vector,(0,0,1))
        # 90 deg roll: q=(sqrt(2)/2,sqrt(2)/2,0,0); use rational unit quaternion 3-4-5 half-angle.
        q=(F(4,5),F(3,5),0,0)
        u=X.down_from_unit_boat_quaternion_zero_heel(q,scope=scope).vector
        self.assertEqual(u,(0,F(24,25),F(7,25)))
        out=X.initialize_from_gauged_seed_zero_heel(seed(True,q),state(),cov(),scope=scope)
        self.assertTrue(out.gauged)

    def test_nonzero_or_dynamic_wind_heel_is_outside_ALT_scope(self):
        with self.assertRaisesRegex(ValueError,'nonzero wind heel'): DS.Scope(F(1,100),0)
        with self.assertRaisesRegex(ValueError,'update_wind_heel'): DS.Scope(0,1)

    def test_only_attitude_block_and_attitude_gyro_cross_are_covariance_mutated(self):
        P0=cov(); x0=state(); out=X.initialize(seed(True),x0,P0,down_body=X.UnitDownWitness((0,0,1)))
        for i in range(21):
            for j in range(21):
                targeted=(i<3 and j<3) or (i<3 and 3<=j<6) or (j<3 and 3<=i<6)
                if not targeted:self.assertEqual(out.P[i][j],P0[i][j],(i,j))
        self.assertEqual(out.x[3:],x0[3:])
        for i in range(3):
            for j in range(3,6): self.assertEqual(out.P[i][j],0); self.assertEqual(out.P[j][i],0)

    def test_gauged_and_free_handoffs_use_shipping_yaw_sigmas(self):
        down=X.UnitDownWitness((0,0,1))
        a=X.initialize(seed(True),state(),cov(),down_body=down)
        b=X.initialize(seed(False),state(),cov(),down_body=down)
        self.assertEqual(a.P[2][2],S.YAW_SIGMA_GAUGED**2); self.assertEqual(b.P[2][2],S.YAW_SIGMA_FREE**2)

    def test_nonunit_down_axis_and_acc_bias_unlock_fail_closed(self):
        with self.assertRaisesRegex(ValueError,'unit'): X.UnitDownWitness((1,1,0))
        bad=S.Result((1,0,0,0),S.TILT_SIGMA,S.YAW_SIGMA_GAUGED,True,True)
        with self.assertRaisesRegex(ValueError,'accelerometer-bias'):
            X.initialize(bad,state(),cov(),down_body=X.UnitDownWitness((0,0,1)))

    def test_readiness_does_not_overpromote_startup(self):
        r=X.readiness()
        self.assertTrue(r['wind_heel_branch_excluded_by_declared_ALT_scope'])
        self.assertTrue(r['zero_heel_body_equals_body_prime'])
        self.assertTrue(r['gauged_unit_seed_to_down_axis_attached_exactly'])
        self.assertFalse(r['startup_capture_closed']); self.assertFalse(r['ALT_STARTUP_PASS'])

if __name__=='__main__': unittest.main()
