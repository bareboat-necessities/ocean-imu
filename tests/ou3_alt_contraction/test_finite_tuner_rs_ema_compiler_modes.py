"""Dual compiler-history R_S EMA regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_source_bound_exp_enclosure as EXP
from tools.stability.ou3_alt_contraction import finite_tuner_rs_alpha_binary32 as A
from tools.stability.ou3_alt_contraction import finite_tuner_rs_ema_compiler_modes as X


def alpha(tau=B.rn32(F(5,2))):
    mult=B.rn32(F(3,2)); dt=B.rn32(F(1,200))
    safe=min(max(tau,A.TIME_MIN),A.TIME_MAX); requested=B.mul(mult,safe)
    lo=min(max(dt,A.HORIZON_MIN),A.HORIZON_MAX); rssec=min(max(requested,lo),A.HORIZON_MAX)
    x=B.div(dt,rssec); elo,ehi,_,_=EXP.enclosure(x); e=B.rn32((elo+ehi)/2)
    return A.step(mult=mult,tau_target=tau,dt=dt,exp_decay=e)


class Tests(unittest.TestCase):
    def test_same_operands_feed_separate_and_fma_histories(self):
        p=B.rn32(F(49,100)); t=B.rn32(F(13,10)); a=alpha()
        out=X.step(p,t,a)
        d=B.sub(t,p)
        self.assertEqual(out.delta,d)
        self.assertEqual(out.next_separate,B.add(p,B.mul(a.alpha,d)))
        self.assertEqual(out.next_fma,B.fma(a.alpha,d,p))
        self.assertEqual(X.committed(out,contracted=False),out.next_separate)
        self.assertEqual(X.committed(out,contracted=True),out.next_fma)

    def test_source_owned_alpha_is_mandatory(self):
        with self.assertRaisesRegex(TypeError,'source-owned alpha'):
            X.step(B.rn32(F(1,2)),B.rn32(1),B.rn32(F(1,100)))

    def test_mode_must_be_literal_when_selecting_result(self):
        out=X.step(B.rn32(F(1,2)),B.rn32(1),alpha())
        with self.assertRaisesRegex(TypeError,'literal compiler contraction'):
            X.committed(out,contracted=1)

    def test_readiness_refuses_to_guess_toolchain_mode(self):
        r=X.readiness()
        for k in ('shipping_RS_EMA_source_shape_matches','required_target_minus_previous_subtraction_materialized',
                  'separate_multiply_add_history_materialized','contracted_FMA_history_materialized',
                  'both_histories_share_same_predecessor_target_and_source_owned_alpha'):
            self.assertTrue(r[k])
        for k in ('shipping_compiler_FP_contraction_mode_qualified','source_uniform_RS_roundoff_supply_bound_closed',
                  'complete_word_RS_compiler_history_attached','source_uniform_complete_600_step_word_qualified',
                  'storage_search_allowed','ALT_LIVE_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
