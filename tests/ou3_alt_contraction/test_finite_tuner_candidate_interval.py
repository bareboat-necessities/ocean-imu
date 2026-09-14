"""Interval-valued exact-real SpectralMSE tuner-candidate regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_shipping_tau_target_binary32 as TARGET
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as C
from tools.stability.ou3_alt_contraction import finite_tuner_candidate_interval as X
from tools.stability.ou3_alt_contraction import finite_tuner_deployment_config as D
from tools.stability.ou3_alt_contraction.finite_tuner_commit import TuneState


def cfg(*,shipping_cadence=True):
    ratio=B.div(B.rn32(F(3,200)),B.rn32(F(11,10))) if shipping_cadence else F(1)
    pmin=B.rn32(F(1,200)) if shipping_cadence else F(1,100)
    pmax=B.rn32(F(3,20)) if shipping_cadence else F(12)
    return C.CandidateConfig(
        TARGET.FLOOR,TARGET.CEIL,F(1),F(1),TARGET.TAU_MIN,TARGET.TAU_MAX,F(4),
        ratio,pmin,pmax,F(1,100),F(100),F(1),F(1,2),F(1),
        B.rn32(F(9,5)),B.rn32(F(2,5)),F(3,2),F(0),F(1,10),True)


def sample(freq=F(1,5),sigma=F(1)):
    return C.WaveBandSample(freq,True,sigma*sigma,0,False,0,1,sigma)


class Tests(unittest.TestCase):
    def test_actual_prior_frequency_cell_executes_without_fake_spectral_roots(self):
        c=cfg(); prev=X.IntervalTuneState.point(TuneState(F(11,10),F(1,2),F(1,2)))
        ema=C.EmaWitness(F(99,100),F(49,50))
        out=X.step(prev,sample(),c,dt=F(1,200),time=F(1,5),last_adapt_time=0,ema=ema)
        self.assertEqual(out.target.frequency,F(1,5))
        self.assertEqual(out.target.tau_target,F(5,2))
        self.assertLess(out.spectral.sqrt_TS_lo,out.spectral.sqrt_TS_hi)
        self.assertLess(out.spectral.u_pow_6_7_lo,out.spectral.u_pow_6_7_hi)
        self.assertLessEqual(out.tune_next.RS_lo,out.tune_next.RS_hi)

    def test_deployment_path_uses_exact_qeff_root_interval_not_machine_cache_value(self):
        prev=X.IntervalTuneState.point(TuneState(F(11,10),F(1,2),F(1,2)))
        ema=C.EmaWitness(F(99,100),F(49,50))
        a=D.shipping_defaults(qeff_pow_result=B.rn32(1))
        b=D.shipping_defaults(qeff_pow_result=B.rn32(2))
        oa=X.step_deployment(prev,sample(),a,dt=F(1,200),time=F(1,5),last_adapt_time=0,ema=ema)
        ob=X.step_deployment(prev,sample(),b,dt=F(1,200),time=F(1,5),last_adapt_time=0,ema=ema)
        self.assertEqual(oa.spectral.qeff_root_lo,a.qeff_cache.true_lo)
        self.assertEqual(oa.spectral.qeff_root_hi,a.qeff_cache.true_hi)
        self.assertEqual(oa.spectral.raw_RS_lo,ob.spectral.raw_RS_lo)
        self.assertEqual(oa.spectral.raw_RS_hi,ob.spectral.raw_RS_hi)
        self.assertNotEqual(a.qeff_pow,b.qeff_pow)

    def test_RS_interval_ema_is_exact_monotone_image(self):
        c=cfg(); prev=X.IntervalTuneState(F(1),F(1),F(1,4),F(3,4))
        ema=C.EmaWitness(F(9,10),F(3,4)); out=X.step(prev,sample(),c,
            dt=F(1,200),time=F(1,5),last_adapt_time=0,ema=ema,bits=64)
        a=F(1,4)
        self.assertEqual(out.tune_next.RS_lo,(1-a)*prev.RS_lo+a*out.spectral.target_RS_lo)
        self.assertEqual(out.tune_next.RS_hi,(1-a)*prev.RS_hi+a*out.spectral.target_RS_hi)

    def test_legacy_exact_cell_is_contained(self):
        c=cfg(shipping_cadence=False); prev=TuneState(F(1),F(1),F(1))
        ema=C.EmaWitness(F(9,10),F(4,5))
        self.assertTrue(X.contains_legacy_exact(prev,sample(F(1,2),F(1)),c,
            dt=F(1,200),time=F(1,5),last_adapt_time=0,
            spectral=C.SpectralWitness(1,1),ema=ema))

    def test_readiness_closes_real_candidate_relation_not_machine_correspondence(self):
        r=X.readiness()
        for k in ('same_frequency_variance_tau_sigma_targets_as_exact_tuner',
                  'SpectralMSE_irrational_target_enclosure_consumed',
                  'deployment_exact_qeff_root_interval_consumed',
                  'deployment_cached_qeff_machine_value_not_used_as_exact_real_coefficient',
                  'RS_interval_propagated_through_literal_nonexpansive_EMA',
                  'tau_sigma_EMA_and_commit_pending_logic_remain_exact',
                  'legacy_rational_exact_candidate_embeds_in_interval_relation',
                  'actual_0p2Hz_prior_cell_representable_without_fake_rational_roots'):
            self.assertTrue(r[k])
        self.assertFalse(r['binary32_qeff_sqrt_pow_correspondence_closed'])
        self.assertFalse(r['binary32_tuner_exp_correspondence_closed'])
        self.assertFalse(r['startup_frontend_composed_with_deployment_interval_RS_state'])
        self.assertFalse(r['source_uniform_complete_600_step_word_qualified'])
        self.assertFalse(r['storage_search_allowed'])


if __name__=='__main__': unittest.main()
