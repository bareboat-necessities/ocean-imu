"""Conditional source-uniform R_S alpha supply-bound regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_source_bound_exp_enclosure as EXP
from tools.stability.ou3_alt_contraction import finite_tuner_deployment_config as D
from tools.stability.ou3_alt_contraction import finite_tuner_rs_alpha_binary32 as M
from tools.stability.ou3_alt_contraction import finite_tuner_rs_alpha_machine_real_join as J
from tools.stability.ou3_alt_contraction import finite_tuner_rs_alpha_conditional_bound as X


def joined(tau):
    cfg=D.shipping_defaults(qeff_pow_result=B.rn32(1)); tau=B.rn32(tau); dt=X.DT
    safe=min(max(tau,M.TIME_MIN),M.TIME_MAX); requested=B.mul(cfg.adapt_RS_mult,safe)
    lo=min(max(dt,M.HORIZON_MIN),M.HORIZON_MAX); rssec=min(max(requested,lo),M.HORIZON_MAX)
    x=B.div(dt,rssec); elo,ehi,_,_=EXP.enclosure(x); e=B.rn32((elo+ehi)/2)
    return J.join(cfg,tau_target=tau,dt=dt,exp_decay=e)


class Tests(unittest.TestCase):
    def test_derived_uniform_bound_has_margin_inside_simple_certificate(self):
        self.assertLess(X.DERIVED_ABS_BOUND,X.CERTIFIED_ABS_BOUND)
        self.assertEqual(X.CERTIFIED_ABS_BOUND,F(1,40000))

    def test_low_tau_clamp_endpoint_is_contained(self):
        out=joined(F(1,50)); X.validate(out)
        self.assertGreaterEqual(out.machine_minus_exact_alpha_lo,-X.CERTIFIED_ABS_BOUND)
        self.assertLessEqual(out.machine_minus_exact_alpha_hi,X.CERTIFIED_ABS_BOUND)

    def test_high_tau_clamp_endpoint_is_contained(self):
        out=joined(12); X.validate(out)
        self.assertGreaterEqual(out.machine_minus_exact_alpha_lo,-X.CERTIFIED_ABS_BOUND)
        self.assertLessEqual(out.machine_minus_exact_alpha_hi,X.CERTIFIED_ABS_BOUND)

    def test_noncanonical_dt_is_not_promoted_by_uniform_lemma(self):
        cfg=D.shipping_defaults(qeff_pow_result=B.rn32(1)); tau=B.rn32(1); dt=B.rn32(F(1,100))
        safe=min(max(tau,M.TIME_MIN),M.TIME_MAX); rssec=B.mul(cfg.adapt_RS_mult,safe)
        x=B.div(dt,rssec); elo,ehi,_,_=EXP.enclosure(x); e=B.rn32((elo+ehi)/2)
        j=J.join(cfg,tau_target=tau,dt=dt,exp_decay=e)
        with self.assertRaisesRegex(ValueError,'canonical compiled 5 ms'):
            X.validate(j)

    def test_readiness_is_explicitly_conditional_on_libm_qualification(self):
        r=X.readiness()
        for k in ('canonical_dt_binary32_fixed','deployed_slew_zero_horizon_domain_reduced_to_safe_tau_interval',
                  'binary32_horizon_multiply_and_divide_roundoff_uniformly_bounded',
                  'second_order_exp_enclosure_width_uniformly_bounded',
                  'conditional_abs_alpha_supply_le_2p5e_minus5','derived_bound_strictly_inside_certified_simple_bound'):
            self.assertTrue(r[k])
        for k in ('target_libm_exp_correspondence_closed','unconditional_machine_execution_alpha_bound_closed',
                  'source_uniform_complete_600_step_word_qualified','storage_search_allowed','ALT_LIVE_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
