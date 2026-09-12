"""Admitted Live tau-ledger product regressions; upstream libm remains open."""
from dataclasses import replace
import unittest

from tools.stability.ou3_alt_contraction import finite_admitted_tau_interleaved_prefix as X
from tools.stability.ou3_alt_contraction import finite_tuner_tau_deployment_ledger as LEDGER
import test_finite_admitted_iss_interleaved_prefix as BASE
import test_finite_complete_word_tau_qualification as QBASE
import test_finite_source_bound_live_word as LBASE


def qualified_prefix():
    p=BASE.begin()
    source=p.prefix.live.live_word
    source=replace(source,runtime=QBASE.shipping_runtime())
    admitted=replace(p.prefix.live,live_word=source)
    interleaved=replace(p.prefix,live=admitted)
    return replace(p,prefix=interleaved)


class Tests(unittest.TestCase):
    def test_component_toy_runtime_cannot_enter_tau_theorem_product(self):
        with self.assertRaisesRegex(ValueError,'detached from shipping binary32 default'):
            X.begin(BASE.begin(),LEDGER.initial())

    def test_source_qualified_live_product_accepts_startup_carried_ledger(self):
        p=qualified_prefix(); tau=LEDGER.State(updates=17)
        s=X.begin(p,tau)
        self.assertIs(s.prefix,p)
        self.assertEqual(s.tau,tau)
        self.assertEqual(s.live_entry_tau_updates,17)

    def test_mag_and_hold_are_literal_tau_ledger_identities(self):
        s=X.begin(qualified_prefix(),LEDGER.State(updates=17)); tau=s.tau
        s,_=X.mag_step(s,**LBASE.mag_kwargs(s.prefix.prefix.live.live_word))
        self.assertIs(s.tau,tau); self.assertEqual(s.live_entry_tau_updates,17)
        s,_=X.set_hold(s,hold=False)
        self.assertIs(s.tau,tau); self.assertEqual(s.live_entry_tau_updates,17)

    def test_live_entry_must_leave_room_for_full_600_step_word(self):
        too_late=LEDGER.State(updates=LEDGER.MAX_UPDATES-599)
        with self.assertRaisesRegex(ValueError,'no certified room'):
            X.begin(qualified_prefix(),too_late)

    def test_short_prefix_cannot_claim_tau_complete_word(self):
        s=X.begin(qualified_prefix(),LEDGER.State(updates=17))
        with self.assertRaisesRegex(ValueError,'exactly 600 physical source transitions'):
            X.complete(s)

    def test_readiness_closes_live_attachment_not_startup_or_libm(self):
        r=X.readiness()
        self.assertTrue(r['Live_product_carries_persistent_dual_compiler_tau_ledger'])
        self.assertTrue(r['each_Live_IMU_requires_same_candidate_frequency_target_and_decay_as_tau_ledger'])
        self.assertTrue(r['MAG_and_HOLD_preserve_tau_ledger_exactly'])
        self.assertFalse(r['startup_master_product_derives_Live_entry_tau_ledger'])
        self.assertFalse(r['upstream_WPE_to_StoredFrequency_binary32_correspondence_closed'])
        self.assertFalse(r['tuner_exp_libm_binary32_correspondence_closed'])
        self.assertFalse(r['source_uniform_complete_600_step_word_qualified'])
        self.assertFalse(r['storage_search_allowed'])


if __name__=='__main__': unittest.main()
