"""Private measurement-only vertical observer finite runtime regressions."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))
from tools.stability.ou3_alt_contraction import finite_vertical_complementary_runtime as V


def cfg(ki=0): return V.Config(F(1,5),ki,1,20)


class Tests(unittest.TestCase):
    def test_initialized_level_rest_sample_stays_level_and_outputs_zero_vertical(self):
        s=V.State(initialized=True)
        out=V.step(s,cfg(),dt=F(1,100),gyro=(0,0,0),acc=(0,0,-1),
                   accel_invnorm=V.InvSqrtWitness(1,1),quat_invnorm=V.InvSqrtWitness(1,1))
        self.assertEqual(out.state.q,(1,0,0,0)); self.assertEqual(out.state.integral,(0,0,0))
        self.assertEqual(out.vertical_accel,0); self.assertEqual(out.state.elapsed,F(1,100))

    def test_first_valid_sample_seeds_then_runs_same_mahony_step(self):
        out=V.step(V.State(),cfg(),dt=F(1,100),gyro=(0,0,0),acc=(0,0,-1),
                   seed=V.SeedWitness((1,0,0,0)),accel_invnorm=V.InvSqrtWitness(1,1),
                   quat_invnorm=V.InvSqrtWitness(1,1))
        self.assertTrue(out.state.initialized); self.assertEqual(out.state.q,(1,0,0,0)); self.assertEqual(out.vertical_accel,0)

    def test_first_degenerate_accel_is_literal_identity(self):
        s=V.State()
        out=V.step(s,cfg(),dt=F(1,100),gyro=(1,2,3),acc=(0,0,0))
        self.assertEqual(out.state,s); self.assertEqual(out.vertical_accel,0)

    def test_accel_normalization_must_belong_to_same_raw_sample(self):
        with self.assertRaises(ValueError):
            V.step(V.State(initialized=True),cfg(),dt=F(1,100),gyro=(0,0,0),acc=(0,0,-1),
                   accel_invnorm=V.InvSqrtWitness(4,F(1,2)),quat_invnorm=V.InvSqrtWitness(1,1))

    def test_quaternion_normalization_must_belong_to_same_euler_successor(self):
        with self.assertRaises(ValueError):
            V.step(V.State(initialized=True),cfg(),dt=F(1,100),gyro=(0,0,0),acc=(0,0,-1),
                   accel_invnorm=V.InvSqrtWitness(1,1),quat_invnorm=V.InvSqrtWitness(4,F(1,2)))

    def test_positive_integral_gain_uses_same_error_and_persists_integral(self):
        # identity attitude with measured gravity tilted +x after Mahony's -acc
        s=V.State(initialized=True)
        # raw acc=(-1,0,0) -> Mahony normalized measured vector=(1,0,0),
        # identity estimated gravity=(0,0,1/2), giving halfey=-1/2.
        # Choose gyro y to cancel proportional+integral feedback so q stays identity.
        ki=F(1,10); kp=F(1,5); dt=F(1,10)
        integ_y=ki*F(-1,2)*dt
        gyro_y=-(integ_y+kp*F(-1,2))
        out=V.step(s,V.Config(kp,ki,1,20),dt=dt,gyro=(0,gyro_y,0),acc=(-1,0,0),
                   accel_invnorm=V.InvSqrtWitness(1,1),quat_invnorm=V.InvSqrtWitness(1,1))
        self.assertEqual(out.state.integral,(0,integ_y,0)); self.assertEqual(out.state.q,(1,0,0,0))

    def test_readiness_is_fail_closed_at_seed_fast_inverse_sqrt_and_source(self):
        r=V.readiness(); self.assertTrue(r['raw_gyro_accel_same_sample_inputs_retained']); self.assertTrue(r['mahony_feedback_and_quaternion_recurrence_materialized']); self.assertTrue(r['vertical_levelled_output_materialized']); self.assertFalse(r['first_sample_FromTwoVectors_seed_attached']); self.assertFalse(r['fast_inv_sqrt_binary32_attached']); self.assertFalse(r['raw_sensor_BRMM_disturbance_relation_attached']); self.assertFalse(r['complete_word_finite_identity']); self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
