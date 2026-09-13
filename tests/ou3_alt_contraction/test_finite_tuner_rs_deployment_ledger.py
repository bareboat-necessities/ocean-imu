"""Persistent dual-compiler R_S deployment-ledger regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_source_bound_exp_enclosure as EXP
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as C
from tools.stability.ou3_alt_contraction import finite_tuner_deployment_config as D
from tools.stability.ou3_alt_contraction import finite_tuner_rs_alpha_binary32 as AB
from tools.stability.ou3_alt_contraction import finite_tuner_rs_alpha_machine_real_join as A
from tools.stability.ou3_alt_contraction import finite_tuner_rs_deployment_ledger as X
from tools.stability.ou3_alt_contraction import finite_tuner_spectral_machine_real_join as T


def cfg(): return D.shipping_defaults(qeff_pow_result=B.rn32(1))


def target(c,tau,sigma,freq):
    t=C.TargetState(F(freq),F(1),B.rn32(F(tau)),B.rn32(F(sigma)))
    return T.join(c,t,pow_result=B.rn32(16),sqrt_result=B.rn32(F(1,5)))


def alpha(t):
    c=t.cfg; tau=t.target.tau_target; dt=B.rn32(F(1,200))
    safe=min(max(tau,AB.TIME_MIN),AB.TIME_MAX); requested=B.mul(c.adapt_RS_mult,safe)
    lo=min(max(dt,AB.HORIZON_MIN),AB.HORIZON_MAX); rssec=min(max(requested,lo),AB.HORIZON_MAX)
    x=B.div(dt,rssec); elo,ehi,_,_=EXP.enclosure(x); e=B.rn32((elo+ehi)/2)
    return A.join(c,tau_target=tau,dt=dt,exp_decay=e)


class Tests(unittest.TestCase):
    def test_initial_state_is_literal_shipping_half(self):
        s=X.initial()
        self.assertEqual(s.separate,B.rn32(F(1,2)))
        self.assertEqual(s.fma,B.rn32(F(1,2)))
        self.assertEqual(s.updates,0)
        self.assertIs(X.hold(s),s)

    def test_two_global_histories_may_consume_distinct_targets_and_alphas(self):
        c=cfg()
        st=target(c,F(5,2),F(9,10),F(1,5))
        ft=target(c,F(2),F(4,5),F(1,4))
        sa=alpha(st); fa=alpha(ft)
        out=X.step(X.initial(),separate_target=st,separate_alpha=sa,
                   fma_target=ft,fma_alpha=fa)
        self.assertEqual(out.state.separate,out.separate.ema.next_separate)
        self.assertEqual(out.state.fma,out.fma.ema.next_fma)
        self.assertEqual(out.state.updates,1)
        self.assertNotEqual(st.target.tau_target,ft.target.tau_target)
        self.assertNotEqual(sa.machine.RS_sec,fa.machine.RS_sec)

    def test_successive_updates_use_same_track_predecessor(self):
        c=cfg(); st=target(c,F(5,2),F(9,10),F(1,5)); ft=target(c,F(2),F(4,5),F(1,4))
        first=X.step(X.initial(),separate_target=st,separate_alpha=alpha(st),
                     fma_target=ft,fma_alpha=alpha(ft))
        second=X.step(first.state,separate_target=st,separate_alpha=alpha(st),
                      fma_target=ft,fma_alpha=alpha(ft))
        self.assertEqual(second.separate.ema.previous,first.state.separate)
        self.assertEqual(second.fma.ema.previous,first.state.fma)
        self.assertEqual(second.state.updates,2)

    def test_crossed_or_detached_alpha_is_rejected(self):
        c=cfg(); st=target(c,F(5,2),F(9,10),F(1,5)); ft=target(c,F(2),F(4,5),F(1,4))
        with self.assertRaisesRegex(ValueError,'same target/config'):
            X.step(X.initial(),separate_target=st,separate_alpha=alpha(ft),
                   fma_target=ft,fma_alpha=alpha(ft))

    def test_different_deployment_configs_cannot_share_one_ledger_step(self):
        c1=cfg(); c2=D.shipping_defaults(qeff_pow_result=B.rn32(F(9,10)))
        st=target(c1,F(5,2),F(9,10),F(1,5)); ft=target(c2,F(2),F(4,5),F(1,4))
        with self.assertRaisesRegex(ValueError,'common deployment config'):
            X.step(X.initial(),separate_target=st,separate_alpha=alpha(st),
                   fma_target=ft,fma_alpha=alpha(ft))

    def test_readiness_closes_persistence_not_commit_or_master_word(self):
        r=X.readiness()
        for k in ('shipping_RS_seed_and_update_source_shape_matches','source_locked_binary32_RS_seed_0p5',
                  'persistent_separate_and_FMA_RS_histories_carried','compiler_histories_may_consume_distinct_SpectralMSE_targets',
                  'compiler_histories_may_consume_distinct_alpha_RS_values','each_track_alpha_bound_to_its_same_target_tau_and_config',
                  'Cold_or_no_candidate_is_literal_RS_identity'):
            self.assertTrue(r[k])
        for k in ('pending_next_sample_commit_snapshot_attached','startup_frontend_RS_machine_history_attached',
                  'Live_600_step_RS_machine_history_attached','target_libm_correspondence_closed',
                  'source_uniform_complete_600_step_word_qualified','storage_search_allowed','ALT_LIVE_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
