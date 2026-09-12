"""Persistent Live tilt-watchdog composition regressions."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))

from tools.stability.ou3_alt_contraction import finite_live_tilt_prefix as X
from tools.stability.ou3_alt_contraction import finite_tilt_watchdog as WATCH
from tools.stability.ou3_alt_contraction import finite_accel_guard_runtime as GUARD
from tools.stability.ou3_alt_contraction import finite_vertical_complementary_runtime as V
from tools.stability.ou3_alt_contraction import finite_band_variance_runtime as B
from tools.stability.ou3_alt_contraction import finite_frontend_runtime as FRONT
from tools.stability.ou3_alt_contraction import finite_racc_runtime as RACC
from tools.stability.ou3_alt_contraction import finite_measurement_runtime as MR
from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_wpe_runtime as W
from tools.stability.ou3_alt_contraction import finite_stillness_runtime as SF
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as C
import test_finite_live_imu_prefix as BASE

G=BASE.G


def call(watchdog,*,tilt_deg,reset_witness=None):
    live,raw,segment,angular,ou,bias,qaxis=BASE.fixture()
    state=X.State(live,watchdog)
    kw=dict(
      dt=segment.h,commit_cfg=BASE.commit_cfg(),boundary_bench_noise_sigma=0,
      guard_cfg=GUARD.Config(),guard_decay=GUARD.DecayWitness(1,1,0,0),guard_rms=GUARD.RmsWitness(0),
      vertical_cfg=V.Config(0,0,G,20),accel_invnorm=V.InvSqrtWitness(G*G,1/G),quat_invnorm=V.InvSqrtWitness(1,1),seed=None,
      band_cfg=B.BandConfig(F(1,2),4,F(1,100),6,F(3,100),F(6,5),F(1,5)),
      racc_cfg=RACC.Config(),nominal_racc_std=(1,1,1),
      angular=angular,Qbase=M.zeros(6,6),ou=ou,bias=bias,qaxis=qaxis,
      accel_conditioning=BASE.SENSOR.AccelConditioning(0,(0,0,0)),
      accel_ldlt=MR.SafeLDLT(False,False,1,F(1,10**7)),
      wpe_cfg=W.WPEConfig(1,4,F(1,2),1,180),wpe_decay=W.ExpWitness(1),
      stats_cfg=B.StatsConfig(4,F(3,10),60,F(1,20),5),
      band_decay=B.BandDecayWitness(1,1),variance_decay=B.VarianceDecayWitness(1),
      bench_noise_sigma=0,noise_sqrt=B.NoiseSqrtWitness(0),
      tracker_lpf_decay=FRONT.LPFDecayWitness(1),still_cfg=SF.Config(),
      still_attenuation=SF.AttenuationWitness(1),candidate_cfg=BASE.candidate_cfg(),
      sigma_wave_sqrt=1,spectral=C.SpectralWitness(1,1),ema=C.EmaWitness(1,1))
    return X.step(state,raw,segment,tilt_deg=tilt_deg,reset_witness=reset_witness,**kw)


class Tests(unittest.TestCase):
    def test_nonfiring_watchdog_persists_timer_and_cannot_consume_reset_witness(self):
        out=call(WATCH.State(F(1,10),0),tilt_deg=60)
        self.assertFalse(out.reset_applied)
        self.assertEqual(out.state.watchdog.over_limit,F(1,10)-2*out.live.prediction.state.reference.time)
        self.assertEqual(out.state.live.mekf,out.live.state.mekf)
        with self.assertRaisesRegex(ValueError,'nonfiring'):
            call(WATCH.State(),tilt_deg=0,
                 reset_witness=WATCH.PreserveYawWitness(BASE.fixture()[0].mekf.q_hat,(0,0,1)))

    def test_firing_watchdog_requires_reset_witness(self):
        # The fixture step is 1/25 s.  0.34 + 0.04 crosses the 0.35 s hold.
        with self.assertRaisesRegex(ValueError,'requires preserve-yaw'):
            call(WATCH.State(F(17,50),0),tilt_deg=80)

    def test_firing_watchdog_applies_covariance_reset_and_sets_cooldown(self):
        live,_,_,_,_,_,_=BASE.fixture()
        witness=WATCH.PreserveYawWitness(live.mekf.q_hat,(0,0,1))
        out=call(WATCH.State(F(17,50),0),tilt_deg=80,reset_witness=witness)
        self.assertTrue(out.reset_applied); self.assertTrue(out.watchdog.fired)
        self.assertEqual(out.state.watchdog.cooldown,3)
        for i in range(3):
            self.assertTrue(all(out.state.live.mekf.covariance[i][j]==0 for j in range(3,21)))
        self.assertEqual(out.state.live.tuner,out.live.state.tuner)
        self.assertEqual(out.state.live.aw_sync,out.live.state.aw_sync)

    def test_readiness_keeps_native_tilt_transcendentals_and_theorem_open(self):
        r=X.readiness()
        self.assertTrue(r['persistent_tilt_overlimit_and_cooldown_carried_in_Live_word'])
        self.assertTrue(r['firing_edge_applies_preserve_yaw_hard_reset_structure'])
        self.assertFalse(r['watchdog_tilt_acos_binary32_attached'])
        self.assertFalse(r['preserve_yaw_atan_asin_angleaxis_binary32_attached'])
        self.assertFalse(r['complete_word_finite_identity']); self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
