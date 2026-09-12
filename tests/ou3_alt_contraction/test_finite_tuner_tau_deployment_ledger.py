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
    return FREQ.store(B.rn32(freq),TARGET.FLOOR,TARGET.CEIL)


def decay(freq=F(1,2)):
    # Build a conservative binary32 exp witness inside the rigorous enclosure
    # for the shipping 5 ms period-scaled tau EMA horizon.
    f=B.rn32(freq); sea=B.div(B.rn32(F(1,2)),f)
    safe=max(B.rn32(F(1,2)),min(B.rn32(6),sea))
    adapt=max(B.rn32(F(1,20)),min(B.rn32(35),B.mul(X.ADAPT_PERIODS,safe)))
    x=B.div(X.DT,adapt)
    return B.rn32(F(1)-x+x*x/F(4))


class Tests(unittest.TestCase):
    def test_initial_state_is_literal_shipping_seed_and_hold_is_identity(self):
        s=X.initial()
        self.assertEqual((s.separate,s.fma),(B.rn32(F(11,10)),)*2)
        self.assertEqual(s.updates,0)
        self.assertIs(X.hold(s),s)

    def test_common_input_helper_advances_two_global_tracks_and_retains_certificates(self):
        s=X.initial(); q=stored(); out=X.step(s,q,cfg(),dt=F(1,200),exp_decay=decay())
        self.assertEqual(out.state.updates,1)
        self.assertEqual(out.state.separate,out.separate_step.next_separate)
        self.assertEqual(out.state.fma,out.fma_step.next_fma)
        self.assertEqual(out.target,TARGET.evaluate(q.stored_hz))
        self.assertLessEqual(abs(out.separate_certificate.supply.residual_separate),out.separate_certificate.bound)
        self.assertLessEqual(abs(out.fma_certificate.supply.residual_fma),out.fma_certificate.bound)

    def test_mode_coherent_tracks_may_consume_different_WPE_frequencies_and_exp_results(self):
        sf=stored(F(1,2)); ff=stored(F(51,100))
        out=X.step_tracks(X.initial(),separate_frequency=sf,fma_frequency=ff,cfg=cfg(),dt=X.DT,
                          separate_exp_decay=decay(F(1,2)),fma_exp_decay=decay(F(51,100)))
        self.assertNotEqual(out.separate_target,out.fma_target)
        self.assertEqual(out.separate_step.frequency,sf.stored_hz)
        self.assertEqual(out.fma_step.frequency,ff.stored_hz)
        self.assertEqual(out.separate_target,TARGET.evaluate(sf.stored_hz))
        self.assertEqual(out.fma_target,TARGET.evaluate(ff.stored_hz))
        with self.assertRaisesRegex(ValueError,'no common target'):
            _=out.target

    def test_successive_updates_use_previous_value_from_same_compiler_track(self):
        first=X.step_tracks(X.initial(),separate_frequency=stored(F(1,2)),fma_frequency=stored(F(51,100)),
                            cfg=cfg(),dt=X.DT,separate_exp_decay=decay(F(1,2)),fma_exp_decay=decay(F(51,100)))
        second=X.step_tracks(first.state,separate_frequency=stored(F(49,100)),fma_frequency=stored(F(52,100)),
                             cfg=cfg(),dt=X.DT,separate_exp_decay=decay(F(49,100)),fma_exp_decay=decay(F(52,100)))
        self.assertEqual(second.separate_step.previous,first.state.separate)
        self.assertEqual(second.fma_step.previous,first.state.fma)
        self.assertEqual(second.state.updates,2)

    def test_wrong_config_dt_or_frequency_object_is_fail_closed(self):
        with self.assertRaisesRegex(ValueError,'tau_coeff'):
            X.step(X.initial(),stored(),replace(cfg(),tau_coeff=F(2)),dt=X.DT,exp_decay=decay())
        with self.assertRaisesRegex(ValueError,'5 ms'):
            X.step(X.initial(),stored(),cfg(),dt=F(1,100),exp_decay=decay())
        with self.assertRaises(TypeError):
            X.step_tracks(X.initial(),separate_frequency=F(1,2),fma_frequency=stored(),cfg=cfg(),dt=X.DT,
                          separate_exp_decay=decay(),fma_exp_decay=decay())

    def test_readiness_keeps_upstream_and_master_attachment_open(self):
        r=X.readiness()
        self.assertTrue(r['global_separate_and_FMA_compiler_tracks_persist_without_branch_explosion'])
        self.assertTrue(r['global_compiler_tracks_accept_mode_coherent_distinct_WPE_frequencies'])
        self.assertTrue(r['global_compiler_tracks_accept_mode_coherent_distinct_exp_results'])
        self.assertFalse(r['common_frequency_across_compiler_tracks_assumed'])
        self.assertTrue(r['each_update_retains_source_locked_target_and_roundoff_certificate_per_track'])
        self.assertFalse(r['global_compiler_track_coherence_includes_WPE_log_state'])
        self.assertFalse(r['upstream_WPE_to_StoredFrequency_binary32_correspondence_closed'])
        self.assertFalse(r['tuner_exp_libm_binary32_correspondence_closed'])
        self.assertFalse(r['storage_search_allowed'])


if __name__=='__main__': unittest.main()
