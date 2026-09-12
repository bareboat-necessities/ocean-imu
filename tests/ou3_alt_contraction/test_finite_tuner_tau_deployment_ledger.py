"""Persistent binary32 tau deployment-ledger regressions."""
from dataclasses import replace
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_frequency_binary32 as FREQ
from tools.stability.ou3_alt_contraction import finite_tuner_tau_deployment_ledger as X
from tools.stability.ou3_alt_contraction import finite_shipping_tau_target_binary32 as TARGET
import test_finite_tuner_tau_binary32 as BASE


def cfg():
    c=BASE.cfg()
    return replace(c,min_freq=TARGET.FLOOR,max_freq=TARGET.CEIL,tau_coeff=B.rn32(1),
                   min_tau=TARGET.TAU_MIN,max_tau=TARGET.TAU_MAX,
                   adapt_tau_sec=X.ADAPT_SEC,adapt_tau_sea_periods=X.ADAPT_PERIODS,
                   clamp_enabled=True)


def stored(freq=F(1,2)):
    return FREQ.store(B.rn32(freq),B.rn32(F(1,20)),B.rn32(5))


def decay():
    # f=.5 -> sea=1, adapt=.4, dt=.005 -> x=.0125.
    return B.rn32(F(98755,100000))


class Tests(unittest.TestCase):
    def test_initial_state_is_literal_shipping_seed_and_hold_is_identity(self):
        s=X.initial()
        self.assertEqual((s.separate,s.fma),(B.rn32(F(11,10)),)*2)
        self.assertEqual(s.updates,0)
        self.assertIs(X.hold(s),s)

    def test_one_update_advances_two_global_compiler_tracks_and_retains_certificates(self):
        s=X.initial(); out=X.step(s,stored(),cfg(),dt=F(1,200),exp_decay=decay())
        self.assertEqual(out.state.updates,1)
        self.assertEqual(out.state.separate,out.separate_step.next_separate)
        self.assertEqual(out.state.fma,out.fma_step.next_fma)
        self.assertEqual(out.target, TARGET.evaluate(stored().stored_hz))
        self.assertLessEqual(abs(out.separate_certificate.supply.residual_separate),out.separate_certificate.bound)
        self.assertLessEqual(abs(out.fma_certificate.supply.residual_fma),out.fma_certificate.bound)

    def test_successive_updates_use_previous_value_from_same_compiler_track(self):
        first=X.step(X.initial(),stored(),cfg(),dt=F(1,200),exp_decay=decay())
        second=X.step(first.state,stored(),cfg(),dt=F(1,200),exp_decay=decay())
        self.assertEqual(second.separate_step.previous,first.state.separate)
        self.assertEqual(second.fma_step.previous,first.state.fma)
        self.assertEqual(second.state.updates,2)

    def test_wrong_config_dt_or_frequency_object_is_fail_closed(self):
        with self.assertRaisesRegex(ValueError,'tau_coeff'):
            X.step(X.initial(),stored(),replace(cfg(),tau_coeff=F(2)),dt=F(1,200),exp_decay=decay())
        with self.assertRaisesRegex(ValueError,'5 ms'):
            X.step(X.initial(),stored(),cfg(),dt=F(1,100),exp_decay=decay())
        with self.assertRaises(TypeError):
            X.step(X.initial(),F(1,2),cfg(),dt=F(1,200),exp_decay=decay())

    def test_readiness_keeps_upstream_and_master_attachment_open(self):
        r=X.readiness()
        self.assertTrue(r['global_separate_and_FMA_compiler_tracks_persist_without_branch_explosion'])
        self.assertTrue(r['each_update_retains_source_locked_target_and_both_roundoff_certificates'])
        self.assertFalse(r['upstream_WPE_to_StoredFrequency_binary32_correspondence_closed'])
        self.assertFalse(r['tuner_exp_libm_binary32_correspondence_closed'])
        self.assertFalse(r['startup_master_product_carries_tau_ledger'])
        self.assertFalse(r['Live_master_product_carries_tau_ledger'])
        self.assertFalse(r['storage_search_allowed'])


if __name__=='__main__': unittest.main()
