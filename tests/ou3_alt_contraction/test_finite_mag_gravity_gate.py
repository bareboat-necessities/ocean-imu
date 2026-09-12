"""Startup magnetic gravity-alignment gate regressions."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))

from tools.stability.ou3_alt_contraction import finite_mag_gravity_gate as X


def root(v): return X.SqrtWitness(v*v,v)


class Tests(unittest.TestCase):
    def test_first_world_lpf_sample_is_exact_rotated_accel_and_persists_branch(self):
        s=X.State(); cfg=X.Config(world_warmup=F(5))
        out=X.imu_step(s,cfg,q_proxy_bw=(1,0,0,0),acc_body=(0,0,-10),gyro_body=(0,0,0),dt=F(1,200),
                       lpf_exp=None,gyro_norm=root(0))
        self.assertEqual(out.world_accel,(0,0,-10)); self.assertEqual(out.state.lpf,(0,0,-10))
        self.assertTrue(out.state.aligned_branch); self.assertEqual(out.align_sin,1)
        self.assertFalse(out.gravity_good_now)

    def test_warm_world_average_closes_gravity_good_and_caps_at_ten(self):
        # Start close enough that one 200 Hz good sample genuinely crosses 10 s.
        s=X.State((0,0,-10),True,F(5),F(9999,1000),True)
        cfg=X.Config(world_warmup=5,max_align_sin=F(3,40),extreme_gyro_dps=30)
        exp=X.LpfExpWitness(F(1,200),12,1)  # alpha=0 keeps same exact average
        out=X.imu_step(s,cfg,q_proxy_bw=(1,0,0,0),acc_body=(0,0,-10),gyro_body=(0,0,0),dt=F(1,200),
                       lpf_exp=exp,lpf_norm=root(10),horizontal_norm=root(0),gyro_norm=root(0))
        self.assertEqual(out.align_sin,0); self.assertTrue(out.gravity_good_now)
        self.assertEqual(out.state.gravity_good,10)

    def test_bad_alignment_and_extreme_motion_decay_good_timer_twice_dt(self):
        cfg=X.Config(world_warmup=0,extreme_gyro_dps=30)
        # Horizontal averaged force fails branch/alignment.
        s=X.State((10,0,0),True,1,1,False)
        exp=X.LpfExpWitness(F(1,10),12,1)
        out=X.imu_step(s,cfg,q_proxy_bw=(1,0,0,0),acc_body=(10,0,0),gyro_body=(0,0,0),dt=F(1,10),
                       lpf_exp=exp,lpf_norm=root(10),horizontal_norm=root(10),gyro_norm=root(0))
        self.assertFalse(out.gravity_good_now); self.assertEqual(out.state.gravity_good,F(4,5))
        # Good gravity direction but violent gyro independently vetoes it.
        s=X.State((0,0,-10),True,1,1,True)
        one=F(1)
        out=X.imu_step(s,cfg,q_proxy_bw=(1,0,0,0),acc_body=(0,0,-10),gyro_body=(one,0,0),dt=F(1,10),
                       lpf_exp=exp,lpf_norm=root(10),horizontal_norm=root(0),gyro_norm=root(1))
        self.assertGreater(out.gyro_dps,30); self.assertFalse(out.gravity_good_now)
        self.assertEqual(out.state.gravity_good,F(4,5))

    def test_eligible_origin_latches_only_after_delay_and_settle(self):
        cfg=X.Config(mag_delay=7,proxy_settle=9,fallback_sec=30)
        s=X.State()
        self.assertIsNone(X.startup_mag_admission(s,cfg,wrapper_time=6,begun=True,have_last_imu=True).state.eligible_t0)
        self.assertIsNone(X.startup_mag_admission(s,cfg,wrapper_time=8,begun=True,have_last_imu=True).state.eligible_t0)
        out=X.startup_mag_admission(s,cfg,wrapper_time=9,begun=True,have_last_imu=True)
        self.assertEqual(out.state.eligible_t0,9); self.assertFalse(out.admitted)
        self.assertEqual(out.reason,'gravity_or_fallback_wait')

    def test_gravity_hold_admits_and_fallback_is_measured_from_latched_origin(self):
        cfg=X.Config(mag_delay=0,hold_sec=2,fallback_sec=30)
        trusted=X.State(gravity_good=2,aligned_branch=True,eligible_t0=5)
        out=X.startup_mag_admission(trusted,cfg,wrapper_time=6,begun=True,have_last_imu=True)
        self.assertTrue(out.admitted); self.assertTrue(out.gravity_trusted); self.assertFalse(out.fallback_ok)
        waiting=X.State(gravity_good=0,eligible_t0=5)
        self.assertFalse(X.startup_mag_admission(waiting,cfg,wrapper_time=F(349,10),begun=True,have_last_imu=True).admitted)
        out=X.startup_mag_admission(waiting,cfg,wrapper_time=35,begun=True,have_last_imu=True)
        self.assertTrue(out.admitted); self.assertFalse(out.gravity_trusted); self.assertTrue(out.fallback_ok)

    def test_have_last_imu_is_required_even_after_gate_opens(self):
        cfg=X.Config(mag_delay=0,hold_sec=2)
        s=X.State(gravity_good=3,eligible_t0=0)
        out=X.startup_mag_admission(s,cfg,wrapper_time=10,begun=True,have_last_imu=False)
        self.assertFalse(out.admitted); self.assertEqual(out.reason,'no_last_imu')

    def test_readiness_keeps_roundoff_source_and_schedule_open(self):
        r=X.readiness()
        self.assertTrue(r['gravity_hold_or_fallback_admission_materialized'])
        self.assertTrue(r['world_acceleration_uses_same_startup_proxy_attitude'])
        self.assertFalse(r['LPF_exp_norm_binary32_attached'])
        self.assertFalse(r['raw_IMU_source_bounds_attached'])
        self.assertFalse(r['complete_word_finite_identity']); self.assertFalse(r['ALT_STARTUP_PASS'])

if __name__=='__main__': unittest.main()
