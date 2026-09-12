"""Guard-persistent raw-IMU -> tuner prefix regressions."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))

from tools.stability.ou3_alt_contraction import finite_guarded_tuner_prefix as X
from tools.stability.ou3_alt_contraction import finite_accel_guard_runtime as G
from tools.stability.ou3_alt_contraction import finite_tuner_frontend_prefix as T
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

GRAV=F(196133,20000); DT=F(1,200)


def candidate_cfg():
    return C.CandidateConfig(F(1,10),2,F(2,5),1,F(1,10),2,2,1,F(1,100),2,
                             F(1,10),2,1,F(1,2),1,1,F(2,5),F(3,2),0,F(1,10),True)


def tuner_state():
    band=B.BandState(band=0,p00=0,p01=0,p11=0,ready=True)
    stats=B.StatsState(frequency=F(1,2),mean_value=0,mean_weight=1,sq_value=1,sq_weight=1)
    return T.State(V.State(initialized=True),W.WPEState(),band,stats,FRONT.LPFState(),SP.State(),TuneState(1,1,1),stage='Live')


def packet(acc=(0,0,GRAV)):
    root=FC.root('A'); from dataclasses import replace
    phys=replace(root.reference,acceleration=(0,0,0),beta=(0,0,0),gyro_bias=(0,0,0))
    # For the level reference, physical specific force is -g internally; the
    # body convention in this proof fixture maps the supplied +g packet through
    # the existing source constructor exactly as the lower-level tests do.
    return RAW.RawImuSample(phys,(0,0,0),(0,0,0),(0,0,0),(0,0,0),(0,0,-GRAV),(0,0,GRAV))


def kwargs():
    return dict(
        vertical_cfg=V.Config(0,0,GRAV,20),accel_invnorm=V.InvSqrtWitness(GRAV*GRAV,1/GRAV),
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
    def test_first_guard_sample_is_identity_but_state_is_persisted(self):
        s=X.State(G.State(),tuner_state())
        out=X.step(s,packet(),dt=DT,guard_cfg=G.Config(),**kwargs())
        self.assertTrue(out.state.guard.initialized)
        self.assertEqual(out.guarded_sample.conditioned_accel_body,out.guarded_sample.raw.raw_accel_body)
        self.assertIs(out.tuner.raw_sample,out.guarded_sample)
        self.assertEqual(out.state.tuner.sample_index,1)

    def test_second_sample_cannot_restart_guard_and_same_descendant_reaches_mahony(self):
        first=X.step(X.State(G.State(),tuner_state()),packet(),dt=DT,guard_cfg=G.Config(),**kwargs())
        second=X.step(first.state,packet(),dt=DT,guard_cfg=G.Config(),
                      guard_decay=G.DecayWitness(1,1,0,0),guard_rms=G.RmsWitness(0),**kwargs())
        self.assertEqual(second.state.tuner.sample_index,2)
        self.assertTrue(second.state.guard.initialized)
        self.assertIs(second.tuner.raw_sample,second.guarded_sample)
        self.assertEqual(second.guarded_sample.guard.state,second.state.guard)
        self.assertEqual(second.tuner.vertical.vertical_accel,0)

    def test_guarded_packet_reserved_for_same_MEKF_accel_event(self):
        out=X.step(X.State(G.State(),tuner_state()),packet(),dt=DT,guard_cfg=G.Config(),**kwargs())
        self.assertEqual(RAW.assert_guarded_acc_measurement_input(out.guarded_sample,
                         out.guarded_sample.conditioned_accel_body),out.guarded_sample.internal_accel)

    def test_readiness_stays_fail_closed_at_Racc_MEKF_and_precision(self):
        r=X.readiness()
        self.assertTrue(r['raw_accel_guard_state_persists_across_samples'])
        self.assertTrue(r['one_guard_successor_feeds_private_vertical_and_future_MEKF_event'])
        self.assertFalse(r['Racc_inflation_from_same_guard_excess_attached'])
        self.assertFalse(r['finite_MEKF_event_interleaved_in_same_sample_word'])
        self.assertFalse(r['guard_exp_sqrt_binary32_ancestry_attached'])
        self.assertFalse(r['complete_word_finite_identity']); self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
