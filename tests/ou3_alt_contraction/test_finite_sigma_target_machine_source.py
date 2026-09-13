"""Machine frontend + stillness -> binary32 sigma-target regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_sigma_target_machine_source as X
from tools.stability.ou3_alt_contraction import finite_stillness_runtime as SR
from tools.stability.ou3_alt_contraction import finite_stillness_sigma_binary32 as S
from tools.stability.ou3_alt_contraction import finite_tuner_deployment_config as D
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_binary32 as T
from tools.stability.ou3_alt_contraction import finite_tuner_spectral_real_enclosure as ROOT
import test_finite_machine_frontend_sigma_source as FB


def sqrtw(x):
    lo,hi=ROOT.sqrt_enclosure(x); return B.rn32((lo+hi)/2)

def expw(x):
    lo,hi=T.exp_minus_enclosure(x); return B.rn32((lo+hi)/2)
def calm_step(state:S.State):
    cfg=SR.Config(); a=B.rn32(0); alpha=B.rn32(cfg.energy_alpha)
    vals=S._sum_products(B.sub(B.rn32(1),alpha),state.energy,alpha,B.rn32(0))
    h=B.rn32(F(1,200)); st=min(B.add(state.still_time,h),S.SIXTY)
    return S.step(state,cfg,vertical_lp=a,dt=h,energy_successor=vals[0],attenuation_exp=expw(st))
def frontend_one():
    state=FB.X.initial(); x=B.rn32(F(1,2)); bc,bn,sc,sn,sg=FB.build_successors(state,x)
    return FB.X.step(state,band_coefficients=bc,band_input=x,band_successor=bn,
                     stats_coefficients=sc,stats_successor=sn,
                     bench_noise_sigma=B.rn32(F(3,100)),noise_sqrt_gain=sg)
def still_one(): return calm_step(S.State())


class Tests(unittest.TestCase):
    def test_first_machine_sample_owns_all_nonlibm_sigma_inputs(self):
        front=frontend_one(); still=still_one(); cfg=D.shipping_defaults(qeff_pow_result=B.rn32(1))
        # DebiasedEMA is already ready after the first default 5 ms update
        # (weight > 1e-6f), but a single-sample central moment is exactly zero.
        # Shipping therefore still reaches the same 1e-6 sigma variance floor,
        # now through the ready branch rather than the unready-noise branch.
        out=X.derive(front,still,cfg,sqrt_result=sqrtw(T.VAR_FLOOR))
        self.assertIs(out.frontend,front); self.assertIs(out.stillness,still)
        self.assertEqual(out.target.accel_variance,front.accel_variance)
        self.assertEqual(out.target.band_noise_sigma,front.band_noise_sigma)
        self.assertEqual(out.target.still,still.state.is_still)
        self.assertEqual(out.target.still_time,still.state.still_time)
        self.assertEqual(out.target.attenuation,still.attenuation)
        self.assertTrue(out.target.var_ready)
        self.assertEqual(out.target.var_wave,T.VAR_FLOOR)

    def test_different_sample_ordinals_cannot_splice(self):
        front=frontend_one()
        second=calm_step(still_one().state)  # valid certificate at ordinal 2
        self.assertEqual(second.state.samples,2)
        with self.assertRaisesRegex(ValueError,'same sample ordinal'):
            X.derive(front,second,D.shipping_defaults(qeff_pow_result=B.rn32(1)),sqrt_result=sqrtw(T.VAR_FLOOR))

    def test_readiness_closes_nonlibm_sigma_source_not_platform_or_master(self):
        r=X.readiness()
        for k in ('same_sample_machine_frontend_and_stillness_required','sigma_var_ready_derived_from_machine_stats_weights',
                  'sigma_accel_variance_derived_from_machine_stats_successor','sigma_band_noise_derived_from_machine_band_successor',
                  'sigma_still_flag_time_and_attenuation_derived_from_machine_stillness_successor',
                  'no_exact_real_sigma_source_inputs_substituted','separate_and_FMA_tuner_tracks_may_consume_distinct_machine_frontend_sources'):
            self.assertTrue(r[k])
        for k in ('final_sigma_sqrt_target_libm_correspondence_closed','upstream_vertical_LP_and_WPE_machine_ancestry_closed',
                  'target_compiler_contraction_membership_closed','startup_frontend_machine_history_attached',
                  'Live_600_step_machine_history_attached','source_uniform_complete_600_step_word_qualified','storage_search_allowed','ALT_LIVE_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
