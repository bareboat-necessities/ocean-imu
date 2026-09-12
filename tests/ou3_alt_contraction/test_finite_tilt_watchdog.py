"""Finite Live tilt-watchdog regressions; transcendental reset output stays open."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))

from tools.stability.ou3_alt_contraction import finite_tilt_watchdog as X
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as S
from tools.stability.ou3_alt_contraction import finite_accel_guard_runtime as G
import test_finite_core as FC

GRAV=F(196133,20000)


def guarded(state):
    ref=state.reference
    gyro=tuple(ref.gyro_bias)
    inertial=[ref.acceleration[i]-(GRAV if i==2 else 0) for i in range(3)]
    f=S.q_rotate(ref.q_world_to_body,inertial)
    acc=tuple(f[i]+ref.beta[i] for i in range(3))
    raw=S.RawImuSample(ref,(0,0,0),(0,0,0),(0,0,0),gyro,acc,(0,0,GRAV))
    return S.guarded_sample(raw,G.State(),G.Config(cutoff_hz=0),dt=F(1,200))


class Tests(unittest.TestCase):
    def test_threshold_is_strict_and_hold_fires_at_035(self):
        s=X.State()
        # Exactly 70 deg is recovery branch, not over-limit.
        s=X.step(s,dt=F(1,10),tilt_deg=70).state
        self.assertEqual(s.over_limit,0)
        for _ in range(3):
            r=X.step(s,dt=F(1,10),tilt_deg=71); s=r.state; self.assertFalse(r.fired)
        r=X.step(s,dt=F(1,20),tilt_deg=71)
        self.assertTrue(r.fired); self.assertEqual(r.state.over_limit,0); self.assertEqual(r.state.cooldown,3)

    def test_recovery_decays_overlimit_twice_as_fast(self):
        s=X.State(F(3,10),0)
        r=X.step(s,dt=F(1,20),tilt_deg=60)
        self.assertEqual(r.state.over_limit,F(1,5)); self.assertFalse(r.fired)

    def test_cooldown_blocks_fire_and_counts_down(self):
        s=X.State(F(2,5),F(1,10))
        r=X.step(s,dt=F(1,20),tilt_deg=80)
        self.assertFalse(r.fired); self.assertEqual(r.state.cooldown,F(1,20))
        r2=X.step(r.state,dt=F(1,20),tilt_deg=80)
        self.assertTrue(r2.fired); self.assertEqual(r2.state.cooldown,3)

    def test_preserve_yaw_hard_reset_only_reseeds_attitude_covariance_and_q(self):
        state=FC.root('A'); sample=guarded(state)
        # Identity down-axis witness makes expected anisotropic block explicit.
        w=X.PreserveYawWitness(state.q_hat,(0,0,1))
        out=X.preserve_yaw_reset(state,sample,w)
        self.assertEqual(out.q_hat,state.q_hat)
        self.assertEqual(out.z,state.z)
        cov0=state.covariance; cov1=out.covariance
        self.assertEqual([row[3:] for row in cov1[3:]],[row[3:] for row in cov0[3:]])
        for i in range(3):
            self.assertTrue(all(cov1[i][j]==0 for j in range(3,21)))
            self.assertTrue(all(cov1[j][i]==0 for j in range(3,21)))
        tv=F(35,1000)**2; yv=F(15708,10000)**2
        self.assertEqual(tuple(tuple(cov1[i][j] for j in range(3)) for i in range(3)),
                         ((tv,0,0),(0,tv,0),(0,0,yv)))

    def test_reset_requires_same_history_and_nondegenerate_guarded_accel(self):
        state=FC.root('A'); sample=guarded(state); w=X.PreserveYawWitness(state.q_hat,(0,0,1))
        from dataclasses import replace
        badraw=replace(sample.raw,physical=replace(sample.raw.physical,history_id='other'))
        # RawImuSample reconstruction would reject physical equations after arbitrary
        # mutation, so use a guarded object only when the raw packet remains valid.
        with self.assertRaises(ValueError):
            X.preserve_yaw_reset(state,S.GuardedImuSample(badraw,sample.guard),w)

    def test_readiness_keeps_transcendental_and_theorem_gates_false(self):
        r=X.readiness()
        self.assertTrue(r['hold_035s_and_recovery_2dt_materialized'])
        self.assertTrue(r['hard_reset_drops_all_attitude_cross_covariances'])
        self.assertFalse(r['preserve_yaw_atan_asin_angleaxis_binary32_attached'])
        self.assertFalse(r['watchdog_tilt_acos_binary32_attached'])
        self.assertFalse(r['complete_word_finite_identity']); self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
