"""Raw-IMU -> tracker-free tuner candidate temporal prefix regressions."""
from dataclasses import replace
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))

from tools.stability.ou3_alt_contraction import finite_tuner_frontend_prefix as X
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as RAW
from tools.stability.ou3_alt_contraction import finite_vertical_complementary_runtime as V
from tools.stability.ou3_alt_contraction import finite_wpe_runtime as W
from tools.stability.ou3_alt_contraction import finite_band_variance_runtime as B
from tools.stability.ou3_alt_contraction import finite_frontend_runtime as FRONT
from tools.stability.ou3_alt_contraction import finite_stillness_runtime as SF
from tools.stability.ou3_alt_contraction import finite_tuner_stillness_projection as SP
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as C
from tools.stability.ou3_alt_contraction.finite_tuner_commit import TuneState
import test_finite_core as FC

G=F(196133,20000); DT=F(1,200)


def candidate_cfg():
    # Shipping-shaped branch with exact unit spectral witnesses at f=0.2:
    # tau_target=(2/5)*(1/2)/(1/5)=1, sigma_target=1, T_S=1, u=1.
    return C.CandidateConfig(F(1,10),2,F(2,5),1,F(1,10),2,2,1,F(1,100),2,
                             F(1,10),2,1,F(1,2),1,1,F(2,5),F(3,2),0,F(1,10),True)


def state():
    # Keep an already-ready zero-gain adaptive band and exact variance=1 so the
    # candidate algebra can be checked without irrational sqrt/power witnesses.
    band=B.BandState(band=0,p00=0,p01=0,p11=0,ready=True)
    stats=B.StatsState(frequency=F(1,2),mean_value=0,mean_weight=1,sq_value=1,sq_weight=1)
    return X.State(V.State(initialized=True),W.WPEState(),band,stats,FRONT.LPFState(),SP.State(),TuneState(1,1,1))


def packet():
    root=FC.root('A'); phys=replace(root.reference,acceleration=(0,0,0),beta=(0,0,0),gyro_bias=(0,0,0))
    return RAW.RawImuSample(phys,(0,0,0),(0,0,0),(0,0,0),(0,0,0),(0,0,-G),(0,0,G))


def do_step(s):
    return X.step(s,packet(),dt=DT,
        vertical_cfg=V.Config(0,0,G,20),accel_invnorm=V.InvSqrtWitness(G*G,1/G),
        quat_invnorm=V.InvSqrtWitness(1,1),seed=None,
        wpe_cfg=W.WPEConfig(1,4,F(1,2),1,180),wpe_decay=W.ExpWitness(1),
        band_cfg=B.BandConfig(F(1,2),4,F(1,100),6,F(3,100),F(6,5),F(1,5)),
        stats_cfg=B.StatsConfig(4,F(3,10),60,F(1,20),5),
        band_decay=B.BandDecayWitness(1,1),variance_decay=B.VarianceDecayWitness(1),
        bench_noise_sigma=0,noise_sqrt=B.NoiseSqrtWitness(0),
        tracker_lpf_decay=FRONT.LPFDecayWitness(1),still_cfg=SF.Config(),
        still_attenuation=SF.AttenuationWitness(1),candidate_cfg=candidate_cfg(),
        sigma_wave_sqrt=1,spectral=C.SpectralWitness(1,1),ema=C.EmaWitness(1,1))


class Tests(unittest.TestCase):
    def test_same_raw_sample_reaches_exact_tuner_candidate_without_tracker(self):
        out=do_step(state())
        self.assertEqual(out.vertical.vertical_accel,0)
        self.assertEqual(out.band.external_tuner_frequency,F(1,5))
        self.assertTrue(out.stillness.state.last_is_still)
        self.assertEqual(out.candidate.frequency,F(1,5))
        self.assertEqual(out.candidate.variance_wave,1)
        self.assertEqual((out.candidate.tau_target,out.candidate.sigma_target,out.candidate.RS_target),(1,1,1))
        self.assertEqual(out.state.tune,TuneState(1,1,1))

    def test_two_samples_preserve_all_adaptation_memory(self):
        first=do_step(state()); second=do_step(first.state)
        self.assertEqual(second.state.sample_index,2); self.assertEqual(second.state.time,2*DT)
        self.assertEqual(second.state.vertical.elapsed,2*DT); self.assertEqual(second.state.wpe.elapsed,2*DT)
        self.assertEqual(second.state.stats.mean_weight,1); self.assertEqual(second.state.stats.sq_weight,1)
        self.assertGreater(second.state.stillness.still_time,first.state.stillness.still_time)
        self.assertEqual(second.candidate.frequency,F(1,5))

    def test_candidate_commit_clock_is_carried_not_restarted(self):
        # At dt=5 ms and 100 ms cadence, first samples do not fire a commit.
        first=do_step(state()); second=do_step(first.state)
        self.assertFalse(first.candidate.pending_after); self.assertFalse(second.candidate.pending_after)
        self.assertEqual(first.state.last_adapt_time,0); self.assertEqual(second.state.last_adapt_time,0)

    def test_readiness_fail_closed_at_next_boundary_and_precision(self):
        r=X.readiness()
        self.assertTrue(r['dominant_frequency_tracker_absent_from_OU_tuner_prefix'])
        self.assertTrue(r['TuneState_and_adapt_clock_persist_across_samples'])
        self.assertFalse(r['next_boundary_staged_commit_composed'])
        self.assertFalse(r['sensor_residual_source_bounds_attached'])
        self.assertFalse(r['transcendental_binary32_attached'])
        self.assertFalse(r['complete_word_finite_identity']); self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
