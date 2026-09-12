"""Wave-band/stillness runtime -> tuner candidate regressions; no source admission."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as C
from tools.stability.ou3_alt_contraction import finite_wpe_runtime as W
from tools.stability.ou3_alt_contraction import finite_band_variance_runtime as B
from tools.stability.ou3_alt_contraction import finite_stillness_runtime as S
from tools.stability.ou3_alt_contraction.finite_tuner_commit import TuneState


def cfg():
    return C.CandidateConfig(F(1,10),2,1,1,F(1,10),2,2,1,F(1,100),2,F(1,10),2,1,F(1,2),1,1,F(2,5),F(3,2),0,F(1,10),True)


def sample(**kw):
    d=dict(frequency_hz=F(1,2),variance_ready=True,accel_variance=1,band_noise_sigma=0,still=False,still_time=0,still_attenuation=1,sigma_wave_sqrt=1)
    d.update(kw); return C.WaveBandSample(**d)


def wpe(freq=F(1,2)):
    return W.UpdateResult(W.WPEState(log_period=1,usable_period=True),True,freq,1/freq,None,None)


def frontend(freq=F(1,2),var=1,noise=0,ready=True):
    return B.FrontendResult(B.BandState(),B.StatsState(frequency=freq),freq,freq,freq,0,ready,F(var),F(noise))


def moving_stillness():
    return S.Result(S.State(energy_ema=1,still_time=0,freq_init=True,freq_state=F(1,2),last_is_still=False),F(1,2),1)


def still_result():
    return S.Result(S.State(energy_ema=0,still_time=1,freq_init=True,freq_state=F(1,2),last_is_still=True),F(1,2),F(1,4))


class Tests(unittest.TestCase):
    def test_same_sample_generates_tau_sigma_spectral_RS_targets(self):
        prev=TuneState(F(1,2),F(1,2),F(1,2)); out=C.step(prev,sample(),cfg(),dt=F(1,100),time=F(1,5),last_adapt_time=0,spectral=C.SpectralWitness(1,1),ema=C.EmaWitness(F(9,10),F(4,5)))
        self.assertEqual(out.frequency,F(1,2)); self.assertEqual(out.variance_wave,1)
        self.assertEqual((out.tau_target,out.sigma_target,out.RS_target),(1,1,1))
        self.assertEqual(out.tune_next.tau_applied,F(11,20)); self.assertEqual(out.tune_next.sigma_applied,F(11,20)); self.assertEqual(out.tune_next.RS_applied,F(3,5))
        self.assertTrue(out.pending_after); self.assertEqual(out.last_adapt_time_after,F(1,5))

    def test_shipping_sample_fields_come_from_frontend_and_stillness_successors(self):
        fs=C.sample_from_runtime(frontend(),moving_stillness(),sigma_wave_sqrt=1)
        self.assertEqual(fs.frequency_hz,F(1,2)); self.assertTrue(fs.variance_ready)
        self.assertEqual(fs.accel_variance,1); self.assertEqual(fs.band_noise_sigma,0)
        self.assertFalse(fs.still); self.assertEqual(fs.still_time,0); self.assertEqual(fs.still_attenuation,1)
        t=C.targets(fs,cfg()); self.assertEqual((t.frequency,t.variance_wave,t.tau_target,t.sigma_target),(F(1,2),1,1,1))

    def test_runtime_stillness_attenuation_is_used_before_wave_floor(self):
        fs=C.sample_from_runtime(frontend(var=4),still_result(),sigma_wave_sqrt=1)
        self.assertTrue(fs.still); self.assertEqual(fs.still_time,1); self.assertEqual(fs.still_attenuation,F(1,4))
        self.assertEqual(C.targets(fs,cfg()).variance_wave,1)

    def test_step_from_runtime_has_no_free_frequency_variance_or_stillness_fields(self):
        prev=TuneState(F(1,2),F(1,2),F(1,2))
        out=C.step_from_runtime(prev,frontend(),moving_stillness(),cfg(),sigma_wave_sqrt=1,dt=F(1,100),time=F(1,5),last_adapt_time=0,spectral=C.SpectralWitness(1,1),ema=C.EmaWitness(1,1))
        self.assertEqual((out.frequency,out.variance_wave),(F(1,2),1))

    def test_direct_wpe_shortcut_rejects_active_or_detached_tuner_clamp(self):
        prev=TuneState(F(1,2),F(1,2),F(1,2))
        out=C.step_from_wpe(prev,wpe(),sample(),cfg(),dt=F(1,100),time=F(1,5),last_adapt_time=0,spectral=C.SpectralWitness(1,1),ema=C.EmaWitness(1,1))
        self.assertEqual(out.frequency,F(1,2))
        with self.assertRaises(ValueError):
            C.step_from_wpe(prev,wpe(F(1,3)),sample(),cfg(),dt=F(1,100),time=0,last_adapt_time=0,spectral=C.SpectralWitness(1,1),ema=C.EmaWitness(1,1))

    def test_not_ready_variance_uses_noise_only_then_wave_floor(self):
        s=sample(variance_ready=False,accel_variance=100,band_noise_sigma=F(1,5),sigma_wave_sqrt=F(1,1000))
        out=C.targets(s,cfg()); self.assertEqual(out.variance_wave,F(1,10**6))

    def test_frequency_and_target_clamps_are_literal(self):
        s=sample(frequency_hz=F(1,100),accel_variance=4,sigma_wave_sqrt=2)
        out=C.targets(s,cfg()); self.assertEqual(out.frequency,F(1,10)); self.assertEqual(out.tau_target,2); self.assertEqual(out.sigma_target,2)

    def test_spectral_witnesses_must_share_tau_sigma_and_cadence(self):
        target=C.targets(sample(),cfg())
        with self.assertRaises(ValueError): C.spectral_RS(cfg(),target,C.SpectralWitness(2,1))
        with self.assertRaises(ValueError): C.spectral_RS(cfg(),target,C.SpectralWitness(1,2))

    def test_commit_cadence_does_not_throttle_ema_inputs(self):
        prev=TuneState(F(1,2),F(1,2),F(1,2)); out=C.step(prev,sample(),cfg(),dt=F(1,100),time=F(1,20),last_adapt_time=0,spectral=C.SpectralWitness(1,1),ema=C.EmaWitness(F(9,10),F(9,10)))
        self.assertFalse(out.pending_after); self.assertNotEqual(out.tune_next,prev); self.assertEqual(out.last_adapt_time_after,0)

    def test_bad_qeff_cache_and_sqrt_variance_are_rejected(self):
        c=cfg()
        with self.assertRaises(ValueError): C.CandidateConfig(c.min_freq,c.max_freq,c.tau_coeff,c.sigma_coeff,c.min_tau,c.max_tau,c.max_sigma,c.pseudo_tau_ratio,c.pseudo_min,c.pseudo_max,c.min_RS,c.max_RS,c.rs_mse_coeff,c.accel_noise_density,2,c.adapt_tau_sec,c.adapt_tau_sea_periods,c.adapt_RS_mult,c.adapt_RS_slew_log,c.adapt_every_sec,True)
        with self.assertRaises(ValueError): C.targets(sample(sigma_wave_sqrt=2),c)

    def test_readiness_stays_fail_closed_upstream(self):
        r=C.readiness(); self.assertTrue(r['frequency_variance_to_tau_sigma_targets']); self.assertTrue(r['default_SpectralMSE_target_same_tau_sigma_cadence']); self.assertTrue(r['band_variance_tuner_frequency_same_history_attached']); self.assertTrue(r['stillness_state_same_history_attached']); self.assertFalse(r['frontend_tracker_and_vertical_measurement_attached']); self.assertFalse(r['exp_sqrt_pow_binary32_enclosed']); self.assertFalse(r['complete_word_finite_identity']); self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
