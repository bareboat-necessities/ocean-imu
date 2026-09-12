"""Same-operand Live tilt-reset runtime regressions."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))

from tools.stability.ou3_alt_contraction import finite_tilt_reset_runtime as X
from tools.stability.ou3_alt_contraction import finite_tilt_watchdog as W
from tools.stability.ou3_alt_contraction import finite_mag_tilt_frame as T
from tools.stability.ou3_alt_contraction import finite_core as C
from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_prediction_graph as P
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as S
from tools.stability.ou3_alt_contraction import finite_accel_guard_runtime as G

GRAV=F(196133,20000)


def core(qhat=(1,0,0,0)):
    ref=C.Reference(time=F(1),live_origin=0,q_world_to_body=(1,0,0,0),
        acceleration=(0,0,0),velocity=(0,0,0),position=(0,0,0),centered_S=(0,0,0),
        beta=(0,0,0),gyro_bias=(0,0,0),history_id='h',bias_root='b',bias_family='BIAS0')
    z=[F(0)]*24
    z[:3]=C.cayley(P.quat_mul(ref.q_world_to_body,P.quat_conj(qhat)))
    cov=M.scaled(M.eye(21),F(1,100))
    return C.State('H',tuple(z),cov,qhat,ref)


def guarded(st, acc=(None,None,None)):
    if acc[0] is None:
        acc=(F(0),F(0),-GRAV)
    acc=tuple(F(x) for x in acc)
    # Physical acceleration is zero, so the no-bias/no-noise specific-force
    # baseline is (0,0,-g). Choose the source residual required to make this
    # exact raw packet, then let the zero-cutoff guard pass the same packet.
    residual=(acc[0],acc[1],acc[2]+GRAV)
    raw=S.RawImuSample(st.reference,(0,0,0),(0,0,0),residual,(0,0,0),acc,(0,0,GRAV))
    return S.guarded_sample(raw,G.State(),G.Config(cutoff_hz=0),dt=F(1,200))


def yaw_half(c=F(1),s=F(0),ch=F(1),sh=F(0)):
    return T.YawHalfWitness(c,s,T.SqrtWitness(c*c+s*s,F(1)),ch,sh)


def atan_zero():
    return yaw_half()


class Tests(unittest.TestCase):
    def test_watchdog_predicate_reads_current_nominal_attitude(self):
        self.assertFalse(X.watchdog_over_limit(core()))
        self.assertTrue(X.watchdog_over_limit(core((F(3,5),F(4,5),0,0))))
        self.assertLess(X.COS70_HI-X.COS70_LO,F(1,10**12))

    def test_preserve_yaw_output_is_derived_from_same_accel_and_predecessor(self):
        st=core((F(3,5),F(4,5),0,0)); sample=guarded(st)
        aw=X.AccTiltWitness(T.SqrtWitness(GRAV*GRAV,GRAV))
        out=X.preserve_yaw_witness(st,sample,
            old_q_norm=T.SqrtWitness(1,1),old_yaw_half=atan_zero(),acc_tilt=aw,
            pitch_cos=T.SqrtWitness(1,1),pitch_half=X.SignedHalfWitness(0,1,1,0),
            roll_half=atan_zero())
        self.assertEqual(out.q_new_hat,(1,0,0,0))
        self.assertEqual(out.down_body_unit,(0,0,1))
        reset=W.preserve_yaw_reset(st,sample,out)
        self.assertEqual(reset.q_hat,(1,0,0,0))
        self.assertEqual(reset.z[:3],(0,0,0))

    def test_covariance_axis_is_accel_intermediate_before_yaw_restore(self):
        # Rational 3-4-5 half angles avoid numerical/trigonometric oracles.
        # Old BODY->WORLD is a yaw with half-angle (4/5,3/5). The guarded
        # accelerometer yields accel-only qref=(4/5,-3/5,0,0), hence shipping
        # reseeds P about down_body=(0,24/25,7/25) *before* restoring yaw.
        st=core((F(4,5),0,0,F(-3,5)))
        acc=(F(0),-F(24,25)*GRAV,-F(7,25)*GRAV)
        sample=guarded(st,acc)
        old_yaw=yaw_half(F(7,25),F(24,25),F(4,5),F(3,5))
        aw=X.AccTiltWitness(
            T.SqrtWitness(GRAV*GRAV,GRAV),
            T.SqrtWitness(F(576,625),F(24,25)),F(4,5),F(3,5))
        pitch=X.SignedHalfWitness(F(0),F(1),F(1),F(0))
        roll=yaw_half(F(7,25),F(24,25),F(4,5),F(3,5))
        out=X.preserve_yaw_witness(st,sample,
            old_q_norm=T.SqrtWitness(1,1),old_yaw_half=old_yaw,acc_tilt=aw,
            pitch_cos=T.SqrtWitness(1,1),pitch_half=pitch,roll_half=roll)
        self.assertEqual(out.down_body_unit,(0,F(24,25),F(7,25)))
        self.assertEqual(out.q_new_hat,(F(16,25),F(-12,25),F(-9,25),F(-12,25)))
        reset=W.preserve_yaw_reset(st,sample,out)
        u=out.down_body_unit
        ts,ys=F(35,1000),F(15708,10000)
        expected=[[ts*ts*((1 if i==j else 0)-u[i]*u[j])+ys*ys*u[i]*u[j]
                   for j in range(3)] for i in range(3)]
        self.assertEqual([list(r[:3]) for r in reset.covariance[:3]],expected)

    def test_near_parallel_axis_uses_shipping_strict_real_cutoff(self):
        # Rational unit vector parameterization. Its transverse component is
        # about 2e-9, hence norm_axis < 1e-8 although it is not exactly zero.
        t=F(1,10**9); d=1+t*t
        y=2*t/d; z=(1-t*t)/d
        acc=(F(0),-y,-z)  # norm exactly one, so target=(0,y,z)
        self.assertEqual(M.dot(acc,acc),F(1))
        self.assertGreater(M.dot((y,0,0),(y,0,0)),0)
        self.assertLess(y*y,X.AXIS_NORM_CUTOFF2)
        near=X.AccTiltWitness(T.SqrtWitness(F(1),F(1)))
        self.assertEqual(X._qref_from_acc(acc,near),(F(1),F(0),F(0),F(0)))
        with self.assertRaises(ValueError):
            X._qref_from_acc(acc,X.AccTiltWitness(
                T.SqrtWitness(F(1),F(1)),T.SqrtWitness(y*y,y),F(1),F(0)))

    def test_wrong_old_yaw_or_acc_norm_is_rejected(self):
        st=core((F(3,5),F(4,5),0,0)); sample=guarded(st)
        aw=X.AccTiltWitness(T.SqrtWitness(GRAV*GRAV,GRAV))
        with self.assertRaises(ValueError):
            X.preserve_yaw_witness(st,sample,old_q_norm=T.SqrtWitness(1,1),old_yaw_half=None,
                acc_tilt=aw,pitch_cos=T.SqrtWitness(1,1),pitch_half=X.SignedHalfWitness(0,1,1,0),
                roll_half=atan_zero())
        with self.assertRaises(ValueError):
            X.preserve_yaw_witness(st,sample,old_q_norm=T.SqrtWitness(1,1),old_yaw_half=atan_zero(),
                acc_tilt=X.AccTiltWitness(T.SqrtWitness(F(1),F(1))),pitch_cos=T.SqrtWitness(1,1),
                pitch_half=X.SignedHalfWitness(0,1,1,0),roll_half=atan_zero())

    def test_watchdog_predicate_recurrence_is_identical(self):
        s=W.State()
        self.assertEqual(W.step(s,dt=F(1,10),tilt_deg=71),
                         W.step_over_limit(s,dt=F(1,10),over_limit=True))
        self.assertEqual(W.step(s,dt=F(1,10),tilt_deg=70),
                         W.step_over_limit(s,dt=F(1,10),over_limit=False))

    def test_readiness_stays_fail_closed(self):
        r=X.readiness()
        self.assertTrue(r['free_final_preserve_yaw_quaternion_removed_from_theorem_entry'])
        self.assertTrue(r['watchdog_threshold_rigorous_cosine_enclosure'])
        self.assertTrue(r['preserve_yaw_axis_norm_strict_cutoff_materialized_real'])
        self.assertTrue(r['reset_covariance_axis_comes_from_accel_only_intermediate_before_yaw_restore'])
        self.assertTrue(r['final_yaw_restore_does_not_reseed_covariance'])
        self.assertFalse(r['near_parallel_cutoff_binary32_closed'])
        self.assertFalse(r['watchdog_boundary_sliver_and_binary32_libm_closed'])
        self.assertFalse(r['complete_word_finite_identity']); self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
