"""Adaptive-band binary32 contraction-set regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_band_binary32_contraction as X


def q(x): return B.rn32(F(x))


def separate_successor(s,x,ql,qh):
    one=B.rn32(1); two=B.rn32(2); al=B.sub(one,ql); ah=B.sub(one,qh)
    low=B.add(B.mul(ql,s.lowpass_low),B.mul(al,x))
    hp=B.mul(ql,B.sub(x,s.lowpass_low))
    band=B.add(B.mul(qh,s.band),B.mul(ah,hp))
    a00=ql; a10=-B.mul(ah,ql); a11=qh; b0=al; b1=B.mul(ah,ql)
    p00=B.add(B.mul(B.mul(a00,a00),s.p00),B.mul(b0,b0)); p00=max(B.rn32(0),p00)
    inner=B.add(B.mul(a10,s.p00),B.mul(a11,s.p01))
    p01=B.add(B.mul(a00,inner),B.mul(b0,b1))
    t1=B.mul(B.mul(a10,a10),s.p00)
    t2=B.mul(B.mul(B.mul(two,a10),a11),s.p01)
    t3=B.mul(B.mul(a11,a11),s.p11); t4=B.mul(b1,b1)
    p11=B.add(B.add(B.add(t1,t2),t3),t4); p11=max(B.rn32(0),p11)
    return X.State(low,band,p00,p01,p11,True)


class Tests(unittest.TestCase):
    def test_zero_seed_step_contains_literal_separate_shipping_evaluation(self):
        s=X.State(); x=q(F(1,2)); ql=q(F(9,10)); qh=q(F(4,5))
        env=X.step(s,x=x,q_low=ql,q_high=qh); sep=separate_successor(s,x,ql,qh)
        self.assertTrue(env.accepts(sep)); self.assertEqual(X.choose(env,lowpass_low=sep.lowpass_low,band=sep.band,p00=sep.p00,p01=sep.p01,p11=sep.p11),sep)

    def test_nonzero_predecessor_enumerates_multiple_local_contraction_outcomes(self):
        s=X.State(q(F(1,7)),q(F(-1,9)),q(F(1,5)),q(F(-1,20)),q(F(3,10)),True)
        env=X.step(s,x=q(F(2,5)),q_low=q(F(91,100)),q_high=q(F(83,100)))
        sep=separate_successor(s,q(F(2,5)),q(F(91,100)),q(F(83,100)))
        self.assertTrue(env.accepts(sep))
        self.assertGreaterEqual(len(env.p01_values),1); self.assertGreaterEqual(len(env.p11_values),1)
        self.assertTrue(len(env.lowpass_values)>1 or len(env.band_values)>1 or len(env.p00_values)>1 or len(env.p01_values)>1 or len(env.p11_values)>1)

    def test_outside_contraction_set_fails_closed(self):
        env=X.step(X.State(),x=q(F(1,2)),q_low=q(F(9,10)),q_high=q(F(4,5)))
        with self.assertRaisesRegex(ValueError,'outside all legal'):
            X.choose(env,lowpass_low=q(7),band=env.band_values[0],p00=env.p00_values[0],p01=env.p01_values[0],p11=env.p11_values[0])

    def test_readiness_does_not_promote_toolchain_or_persistence(self):
        r=X.readiness()
        self.assertTrue(r['shipping_adaptive_band_signal_and_covariance_source_shape_matches'])
        self.assertTrue(r['one_step_all_local_no_reassociation_contraction_outcomes_enumerated'])
        self.assertTrue(r['no_single_global_FMA_bit_assumed_for_band_polynomial'])
        for k in ('band_corner_and_exp_binary32_coefficient_production_closed','target_compiler_no_reassociation_contract_qualified',
                  'persistent_band_machine_history_composed','band_noise_floor_sqrt_machine_relation_composed',
                  'source_uniform_complete_600_step_word_qualified','storage_search_allowed','ALT_LIVE_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
