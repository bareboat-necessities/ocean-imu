"""Adaptive band + variance finite runtime regressions."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))
from tools.stability.ou3_alt_contraction import finite_band_variance_runtime as B
from tools.stability.ou3_alt_contraction import finite_wpe_runtime as W


def bcfg(): return B.BandConfig(F(1,2),4,F(1,100),6,F(3,100),F(6,5))
def scfg(): return B.StatsConfig(4,F(3,10),60,F(1,20),5)
def wpe(freq=F(1,2)): return W.UpdateResult(W.WPEState(log_period=1),True,freq,1/freq,None,None)


class Tests(unittest.TestCase):
    def test_band_state_and_white_noise_covariance_use_same_coefficients(self):
        s=B.band_step(B.BandState(),bcfg(),x=2,dt=F(1,10),f_ref=F(1,2),decay=B.BandDecayWitness(F(1,2),F(1,2)))
        self.assertEqual(s.lowpass_low,1); self.assertEqual(s.band,F(1,2))
        self.assertEqual(s.p00,F(1,4)); self.assertEqual(s.p01,F(1,8)); self.assertEqual(s.p11,F(1,16)); self.assertTrue(s.ready)

    def test_degenerate_corner_branch_is_literal_identity(self):
        cfg=B.BandConfig(F(1,2),4,1,2,F(1,10),1); s=B.BandState(band=3,p00=1,p11=1,ready=True)
        out=B.band_step(s,cfg,x=9,dt=1,f_ref=F(1,2),decay=B.BandDecayWitness(F(1,2),F(1,2)))
        self.assertEqual(out,s)

    def test_previous_tuner_frequency_drives_band_current_wpe_drives_stats(self):
        prev_stats=B.StatsState(frequency=F(1,4))
        out=B.frontend_step(B.BandState(),prev_stats,wpe(F(1,2)),vertical_accel=2,dt=F(1,10),
            band_cfg=bcfg(),stats_cfg=scfg(),band_decay=B.BandDecayWitness(F(1,2),F(1,2)),variance_decay=B.VarianceDecayWitness(F(1,2)),bench_noise_sigma=2,noise_sqrt=B.NoiseSqrtWitness(F(1,4)))
        self.assertEqual(out.band_reference_frequency,F(1,4)); self.assertEqual(out.current_tuner_frequency,F(1,2)); self.assertEqual(out.filtered_accel,F(1,2)); self.assertTrue(out.variance_ready); self.assertEqual(out.accel_variance,0); self.assertEqual(out.band_noise_sigma,F(1,2))

    def test_unready_tuner_uses_current_wpe_for_band_reference(self):
        out=B.frontend_step(B.BandState(),B.StatsState(),wpe(F(1,2)),vertical_accel=2,dt=F(1,10),band_cfg=bcfg(),stats_cfg=scfg(),band_decay=B.BandDecayWitness(F(1,2),F(1,2)),variance_decay=B.VarianceDecayWitness(F(1,2)),bench_noise_sigma=2,noise_sqrt=B.NoiseSqrtWitness(F(1,4)))
        self.assertEqual(out.band_reference_frequency,F(1,2))

    def test_unready_band_uses_raw_bench_noise_and_no_sqrt_witness(self):
        cfg=B.BandConfig(F(1,2),4,1,2,F(1,10),1)
        out=B.frontend_step(B.BandState(),B.StatsState(),wpe(F(1,2)),vertical_accel=2,dt=1,band_cfg=cfg,stats_cfg=scfg(),band_decay=B.BandDecayWitness(F(1,2),F(1,2)),variance_decay=B.VarianceDecayWitness(F(1,2)),bench_noise_sigma=F(3,25))
        self.assertFalse(out.band_state.ready); self.assertEqual(out.band_noise_sigma,F(3,25)); self.assertEqual(out.filtered_accel,0)
        with self.assertRaises(ValueError):
            B.frontend_step(B.BandState(),B.StatsState(),wpe(F(1,2)),vertical_accel=2,dt=1,band_cfg=cfg,stats_cfg=scfg(),band_decay=B.BandDecayWitness(F(1,2),F(1,2)),variance_decay=B.VarianceDecayWitness(F(1,2)),bench_noise_sigma=F(3,25),noise_sqrt=B.NoiseSqrtWitness(0))

    def test_ready_band_requires_same_covariance_sqrt(self):
        with self.assertRaises(ValueError):
            B.frontend_step(B.BandState(),B.StatsState(),wpe(),vertical_accel=2,dt=F(1,10),band_cfg=bcfg(),stats_cfg=scfg(),band_decay=B.BandDecayWitness(F(1,2),F(1,2)),variance_decay=B.VarianceDecayWitness(F(1,2)),bench_noise_sigma=2,noise_sqrt=B.NoiseSqrtWitness(F(1,3)))
        with self.assertRaises(ValueError):
            B.frontend_step(B.BandState(),B.StatsState(),wpe(),vertical_accel=2,dt=F(1,10),band_cfg=bcfg(),stats_cfg=scfg(),band_decay=B.BandDecayWitness(F(1,2),F(1,2)),variance_decay=B.VarianceDecayWitness(F(1,2)),bench_noise_sigma=2)

    def test_debiased_variance_is_exact_second_central_moment(self):
        s=B.StatsState(frequency=F(1,2),mean_value=1,mean_weight=1,sq_value=5,sq_weight=1); self.assertEqual(B.variance(s),4)

    def test_readiness_fail_closed_at_frontend_and_transcendentals(self):
        r=B.readiness(); self.assertTrue(r['adaptive_band_state_and_covariance_materialized']); self.assertTrue(r['adaptive_band_identity_branches_materialized']); self.assertTrue(r['unready_band_noise_floor_branch_materialized']); self.assertTrue(r['previous_tuner_frequency_drives_band_corner']); self.assertFalse(r['band_exp_and_noise_sqrt_binary32_ancestry_attached']); self.assertFalse(r['vertical_accel_frontend_same_history_attached']); self.assertFalse(r['complete_word_finite_identity']); self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
