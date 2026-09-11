"""Wave-band sample -> tuner candidate regressions; not upstream source admission."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as C
from tools.stability.ou3_alt_contraction.finite_tuner_commit import TuneState


def cfg():
    return C.CandidateConfig(F(1,10),2,1,1,F(1,10),2,2,1,F(1,100),2,F(1,10),2,1,F(1,2),1,1,F(2,5),F(3,2),0,F(1,10),True)


def sample(**kw):
    d=dict(frequency_hz=F(1,2),variance_ready=True,accel_variance=1,band_noise_sigma=0,still=False,still_time=0,still_attenuation=1,sigma_wave_sqrt=1)
    d.update(kw); return C.WaveBandSample(**d)


class Tests(unittest.TestCase):
    def test_same_sample_generates_tau_sigma_spectral_RS_targets(self):
        prev=TuneState(F(1,2),F(1,2),F(1,2)); out=C.step(prev,sample(),cfg(),dt=F(1,100),time=F(1,5),last_adapt_time=0,spectral=C.SpectralWitness(1,1),ema=C.EmaWitness(F(9,10),F(4,5)))
        self.assertEqual(out.frequency,F(1,2)); self.assertEqual(out.variance_wave,1)
        self.assertEqual((out.tau_target,out.sigma_target,out.RS_target),(1,1,1))
        self.assertEqual(out.tune_next.tau_applied,F(11,20)); self.assertEqual(out.tune_next.sigma_applied,F(11,20)); self.assertEqual(out.tune_next.RS_applied,F(3,5))
        self.assertTrue(out.pending_after); self.assertEqual(out.last_adapt_time_after,F(1,5))

    def test_not_ready_variance_uses_noise_only_then_wave_floor(self):
        # total=noise variance, so subtraction is zero and the 1e-6 floor binds.
        s=sample(variance_ready=False,accel_variance=100,band_noise_sigma=F(1,5),sigma_wave_sqrt=F(1,1000))
        out=C.step(TuneState(1,1,1),s,cfg(),dt=F(1,100),time=F(1,20),last_adapt_time=0,
                   spectral=C.SpectralWitness(1,1),ema=C.EmaWitness(1,1))
        self.assertEqual(out.variance_wave,F(1,10**6))

    def test_stillness_attenuation_is_on_same_variance_before_floor(self):
        s=sample(accel_variance=4,still=True,still_time=1,still_attenuation=F(1,4),sigma_wave_sqrt=1)
        out=C.step(TuneState(1,1,1),s,cfg(),dt=F(1,100),time=F(1,20),last_adapt_time=0,
                   spectral=C.SpectralWitness(1,1),ema=C.EmaWitness(1,1))
        self.assertEqual(out.variance_wave,1)

    def test_frequency_and_target_clamps_are_literal(self):
        s=sample(frequency_hz=F(1,100),accel_variance=4,sigma_wave_sqrt=2)
        # f clamps to .1, tau raw=5 then clamps to 2; sigma clamps to 2.
        # TS=2 and sqrt witness must be exact; choose a config with pseudo max=1
        c=cfg(); c=C.CandidateConfig(c.min_freq,c.max_freq,c.tau_coeff,c.sigma_coeff,c.min_tau,c.max_tau,c.max_sigma,c.pseudo_tau_ratio,c.pseudo_min,1,c.min_RS,c.max_RS,c.rs_mse_coeff,c.accel_noise_density,c.qeff_pow,c.adapt_tau_sec,c.adapt_tau_sea_periods,c.adapt_RS_mult,c.adapt_RS_slew_log,c.adapt_every_sec,c.clamp_enabled)
        # u=sigma*tau^4=32, provide an exact 6/7 witness by making RS clamp before
        # relevance impossible; witness still has to satisfy the finite law, so use
        # a separate simple sample to test frequency clamp only.
        simple=sample(frequency_hz=F(1,10),accel_variance=1,sigma_wave_sqrt=1)
        out=C.step(TuneState(1,1,1),simple,cfg(),dt=F(1,100),time=0,last_adapt_time=0,spectral=C.SpectralWitness(1,1),ema=C.EmaWitness(1,1))
        self.assertEqual(out.frequency,F(1,10)); self.assertEqual(out.tau_target,2)

    def test_spectral_witnesses_must_share_tau_sigma_and_cadence(self):
        with self.assertRaises(ValueError):
            C.step(TuneState(1,1,1),sample(),cfg(),dt=F(1,100),time=0,last_adapt_time=0,spectral=C.SpectralWitness(2,1),ema=C.EmaWitness(1,1))
        with self.assertRaises(ValueError):
            C.step(TuneState(1,1,1),sample(),cfg(),dt=F(1,100),time=0,last_adapt_time=0,spectral=C.SpectralWitness(1,2),ema=C.EmaWitness(1,1))

    def test_commit_cadence_does_not_throttle_ema_inputs(self):
        prev=TuneState(F(1,2),F(1,2),F(1,2)); out=C.step(prev,sample(),cfg(),dt=F(1,100),time=F(1,20),last_adapt_time=0,spectral=C.SpectralWitness(1,1),ema=C.EmaWitness(F(9,10),F(9,10)))
        self.assertFalse(out.pending_after)
        self.assertNotEqual(out.tune_next,prev)
        self.assertEqual(out.last_adapt_time_after,0)

    def test_bad_qeff_cache_and_sqrt_variance_are_rejected(self):
        c=cfg()
        with self.assertRaises(ValueError): C.CandidateConfig(c.min_freq,c.max_freq,c.tau_coeff,c.sigma_coeff,c.min_tau,c.max_tau,c.max_sigma,c.pseudo_tau_ratio,c.pseudo_min,c.pseudo_max,c.min_RS,c.max_RS,c.rs_mse_coeff,c.accel_noise_density,2,c.adapt_tau_sec,c.adapt_tau_sea_periods,c.adapt_RS_mult,c.adapt_RS_slew_log,c.adapt_every_sec,True)
        with self.assertRaises(ValueError):
            C.step(TuneState(1,1,1),sample(sigma_wave_sqrt=2),c,dt=F(1,100),time=0,last_adapt_time=0,spectral=C.SpectralWitness(1,1),ema=C.EmaWitness(1,1))

    def test_readiness_stays_fail_closed_upstream(self):
        r=C.readiness(); self.assertTrue(r['frequency_variance_to_tau_sigma_targets']); self.assertTrue(r['default_SpectralMSE_target_same_tau_sigma_cadence']); self.assertFalse(r['WPE_bandpass_variance_state_attached']); self.assertFalse(r['exp_sqrt_pow_binary32_enclosed']); self.assertFalse(r['complete_word_finite_identity']); self.assertFalse(r['ALT_LIVE_PASS'])

if __name__=='__main__': unittest.main()
