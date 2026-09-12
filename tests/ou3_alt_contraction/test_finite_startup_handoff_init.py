"""Literal startup MEKF initialize_from_attitude reset regressions."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))

from tools.stability.ou3_alt_contraction import finite_startup_handoff_init as X
from tools.stability.ou3_alt_contraction import finite_startup_handoff_seed as S


def seed(gauged=True):
    return S.Result((1,0,0,0),S.TILT_SIGMA,
                    S.YAW_SIGMA_GAUGED if gauged else S.YAW_SIGMA_FREE,
                    False,gauged)

def state(): return tuple(F(i+1,100) for i in range(21))
def cov():
    # Dense symmetric matrix with unique rational entries so unintended reset
    # of any non-target entry is immediately visible.
    return [[F((i+1)*100+(j+1),100000) for j in range(21)] for i in range(21)]


class Tests(unittest.TestCase):
    def test_down_z_gives_diagonal_tilt_tilt_yaw_split(self):
        s=seed(True); out=X.initialize(s,state(),cov(),down_body_prime=X.UnitDownWitness((0,0,1)))
        tv=S.TILT_SIGMA*S.TILT_SIGMA; yv=S.YAW_SIGMA_GAUGED*S.YAW_SIGMA_GAUGED
        self.assertEqual([list(r[:3]) for r in out.P[:3]],
                         [[tv,0,0],[0,tv,0],[0,0,yv]])
        self.assertEqual(out.x[:3],(0,0,0))

    def test_arbitrary_unit_axis_keeps_projector_formula_exact(self):
        # (3/5,4/5,0) is exactly unit.
        u=X.UnitDownWitness((F(3,5),F(4,5),0)); s=seed(True)
        A=X.attitude_covariance(u,s.tilt_sigma,s.yaw_sigma)
        tv=s.tilt_sigma*s.tilt_sigma; yv=s.yaw_sigma*s.yaw_sigma
        self.assertEqual(A[0][0],tv*F(16,25)+yv*F(9,25))
        self.assertEqual(A[1][1],tv*F(9,25)+yv*F(16,25))
        self.assertEqual(A[0][1],(yv-tv)*F(12,25))
        self.assertEqual(A[2][2],tv)
        self.assertEqual(A[1][0],A[0][1])

    def test_only_attitude_block_and_attitude_gyro_cross_are_covariance_mutated(self):
        P0=cov(); x0=state(); out=X.initialize(seed(True),x0,P0,down_body_prime=X.UnitDownWitness((0,0,1)))
        for i in range(21):
            for j in range(21):
                targeted=(i<3 and j<3) or (i<3 and 3<=j<6) or (j<3 and 3<=i<6)
                if not targeted:
                    self.assertEqual(out.P[i][j],P0[i][j],(i,j))
        self.assertEqual(out.x[3:],x0[3:])
        for i in range(3):
            for j in range(3,6):
                self.assertEqual(out.P[i][j],0); self.assertEqual(out.P[j][i],0)

    def test_gauged_and_free_handoffs_use_shipping_yaw_sigmas(self):
        down=X.UnitDownWitness((0,0,1))
        a=X.initialize(seed(True),state(),cov(),down_body_prime=down)
        b=X.initialize(seed(False),state(),cov(),down_body_prime=down)
        self.assertEqual(a.P[2][2],S.YAW_SIGMA_GAUGED**2)
        self.assertEqual(b.P[2][2],S.YAW_SIGMA_FREE**2)
        self.assertLess(a.P[2][2],b.P[2][2])

    def test_nonunit_down_axis_and_acc_bias_unlock_fail_closed(self):
        with self.assertRaisesRegex(ValueError,'unit'):
            X.UnitDownWitness((1,1,0))
        bad=S.Result((1,0,0,0),S.TILT_SIGMA,S.YAW_SIGMA_GAUGED,True,True)
        with self.assertRaisesRegex(ValueError,'accelerometer-bias'):
            X.initialize(bad,state(),cov(),down_body_prime=X.UnitDownWitness((0,0,1)))

    def test_readiness_does_not_overpromote_startup(self):
        r=X.readiness()
        self.assertTrue(r['attitude_error_zero_reset_materialized'])
        self.assertTrue(r['anisotropic_tilt_yaw_covariance_formula_materialized'])
        self.assertTrue(r['all_other_mean_and_covariance_entries_retained'])
        self.assertFalse(r['body_prime_down_axis_from_boat_quaternion_and_wind_heel_attached'])
        self.assertFalse(r['startup_capture_closed']); self.assertFalse(r['ALT_STARTUP_PASS'])

if __name__=='__main__': unittest.main()
