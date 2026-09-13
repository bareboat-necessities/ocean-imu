"""Same-source exact-real/binary32 R_S alpha regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_source_bound_exp_enclosure as EXP
from tools.stability.ou3_alt_contraction import finite_tuner_deployment_config as D
from tools.stability.ou3_alt_contraction import finite_tuner_rs_alpha_binary32 as M
from tools.stability.ou3_alt_contraction import finite_tuner_rs_alpha_machine_real_join as X


def witness(cfg,tau,dt):
    safe=min(max(tau,M.TIME_MIN),M.TIME_MAX)
    requested=B.mul(cfg.adapt_RS_mult,safe)
    lo=min(max(dt,M.HORIZON_MIN),M.HORIZON_MAX)
    rssec=min(max(requested,lo),M.HORIZON_MAX)
    x=B.div(dt,rssec); elo,ehi,_,_=EXP.enclosure(x)
    return B.rn32((elo+ehi)/2)


class Tests(unittest.TestCase):
    def test_same_config_tau_dt_root_both_alpha_relations(self):
        cfg=D.shipping_defaults(qeff_pow_result=B.rn32(1))
        tau=B.rn32(F(5,2)); dt=B.rn32(F(1,200)); e=witness(cfg,tau,dt)
        out=X.join(cfg,tau_target=tau,dt=dt,exp_decay=e)
        self.assertEqual(out.machine.mult,cfg.adapt_RS_mult)
        self.assertEqual(out.machine.tau_target,tau)
        self.assertEqual(out.machine.dt,dt)
        self.assertLessEqual(out.exact_alpha_lo,out.exact_alpha_hi)
        self.assertEqual(out.machine_minus_exact_alpha_lo,out.machine.alpha-out.exact_alpha_hi)
        self.assertEqual(out.machine_minus_exact_alpha_hi,out.machine.alpha-out.exact_alpha_lo)

    def test_exact_horizon_does_not_reuse_rounded_machine_product(self):
        cfg=D.shipping_defaults(qeff_pow_result=B.rn32(1))
        tau=B.rn32(F(7,3)); dt=B.rn32(F(1,200)); out=X.join(cfg,tau_target=tau,dt=dt,exp_decay=witness(cfg,tau,dt))
        self.assertEqual(out.exact_requested_horizon,F(cfg.adapt_RS_mult)*out.exact_safe_tau)
        self.assertEqual(out.machine.requested_horizon,B.mul(cfg.adapt_RS_mult,tau))
        # Equality is allowed on special cells but is never an assumption.
        self.assertEqual(out.exact_exp_argument,dt/out.exact_RS_sec)
        self.assertEqual(out.machine.exp_argument,B.div(dt,out.machine.RS_sec))

    def test_nonbinary_join_operands_fail_closed(self):
        cfg=D.shipping_defaults(qeff_pow_result=B.rn32(1))
        with self.assertRaisesRegex(ValueError,'actual binary32'):
            X.join(cfg,tau_target=F(1,3),dt=B.rn32(F(1,200)),exp_decay=B.rn32(F(99,100)))

    def test_readiness_removes_free_exact_alpha_but_not_libm_gap(self):
        r=X.readiness()
        for k in ('exact_real_RS_horizon_rooted_at_same_config_tau_dt',
                  'exact_real_exp_decay_rigorously_enclosed_at_same_horizon_argument',
                  'exact_real_alpha_interval_not_free_witness',
                  'machine_alpha_path_rooted_at_same_config_tau_dt',
                  'machine_minus_exact_alpha_interval_exposed'):
            self.assertTrue(r[k])
        for k in ('target_libm_exp_correspondence_closed','source_uniform_tight_alpha_supply_bound_closed',
                  'RS_EMA_interval_consumes_alpha_interval','source_uniform_complete_600_step_word_qualified',
                  'storage_search_allowed','ALT_LIVE_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
