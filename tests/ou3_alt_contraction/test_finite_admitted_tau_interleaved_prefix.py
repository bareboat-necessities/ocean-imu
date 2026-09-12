"""Admitted Live tau-ledger product regressions; upstream libm remains open."""
from dataclasses import replace
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_admitted_tau_interleaved_prefix as X
from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B32
from tools.stability.ou3_alt_contraction import finite_shipping_tau_target_binary32 as TARGET
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as CAND
from tools.stability.ou3_alt_contraction import finite_tuner_tau_deployment_ledger as LEDGER
from tools.stability.ou3_alt_contraction import finite_wpe_frequency_binary32 as WPEF
import test_finite_admitted_iss_interleaved_prefix as BASE
import test_finite_admitted_source_imu_word as IBASE
import test_finite_complete_word_tau_qualification as QBASE
import test_finite_source_bound_live_word as LBASE
import test_finite_source_bound_prediction_word as PBASE


def qualified_prefix():
    p=BASE.begin()
    source=p.prefix.live.live_word
    source=replace(source,runtime=QBASE.shipping_runtime())
    admitted=replace(p.prefix.live,live_word=source)
    interleaved=replace(p.prefix,live=admitted)
    return replace(p,prefix=interleaved)


def wpe_source_for(state,dynamic):
    entry=X._entry_wpe(state)
    if not entry.usable_period:
        return WPEF.tuner_frequency(entry,min_hz=TARGET.FLOOR,max_hz=TARGET.CEIL)
    exact_f=F(dynamic['wpe_current_frequency'])
    frequency=B32.rn32(exact_f)
    log=WPEF.bind_log_state(entry,B32.rn32(entry.log_period))
    getter=WPEF.getters(log,period_exp=B32.rn32(F(1,frequency)),frequency_exp=frequency)
    return WPEF.tuner_frequency(entry,min_hz=TARGET.FLOOR,max_hz=TARGET.CEIL,getter=getter)


class Tests(unittest.TestCase):
    def test_component_toy_runtime_cannot_enter_tau_theorem_product(self):
        with self.assertRaisesRegex(ValueError,'detached from shipping binary32 default'):
            X.begin(BASE.begin(),LEDGER.initial())

    def test_source_qualified_live_product_accepts_startup_carried_ledger_for_component_algebra(self):
        p=qualified_prefix(); tau=LEDGER.State(updates=17)
        s=X.begin(p,tau)
        self.assertIs(s.prefix,p)
        self.assertEqual(s.tau,tau)
        self.assertEqual(s.live_entry_tau_updates,17)

    def test_strong_IMU_frequency_source_is_same_preupdate_WPE_state(self):
        s=X.begin(qualified_prefix(),LEDGER.State(updates=17))
        witness,segment,raw,r,b,dynamic=IBASE.operands(s.prefix.prefix.live)
        dynamic=dict(dynamic)
        # The tau deployment witness and exact tuner shadow consume the same
        # carried decay.  0.995 lies inside the 5ms/prior-horizon exp enclosure.
        e=B32.rn32(F(199,200))
        dynamic['ema']=CAND.EmaWitness(e,dynamic['ema'].decay_RS)
        source=wpe_source_for(s,dynamic)
        out=X.imu_step_from_wpe(s,tuner_frequency=source,tau_exp_decay=e,
            restricted=r,bias_restricted=b,witness=witness,raw=raw,
            packet_id='imu-wpe-source',**PBASE.root_args(),**dynamic)
        self.assertEqual(out.state.tau.updates,18)
        self.assertEqual(out.tau_step.target.binary32_target,out.tau_step.separate_step.tau_target)
        self.assertEqual(source.shadow,X._entry_wpe(s))

    def test_detached_preupdate_WPE_frequency_source_rejected_before_execution(self):
        s=X.begin(qualified_prefix(),LEDGER.State(updates=17))
        entry=X._entry_wpe(s)
        other=replace(entry,elapsed=entry.elapsed+1)
        if other.usable_period:
            log=WPEF.bind_log_state(other,B32.rn32(other.log_period))
            getter=WPEF.getters(log,period_exp=B32.rn32(2),frequency_exp=B32.rn32(F(1,2)))
            source=WPEF.tuner_frequency(other,min_hz=TARGET.FLOOR,max_hz=TARGET.CEIL,getter=getter)
        else:
            source=WPEF.tuner_frequency(other,min_hz=TARGET.FLOOR,max_hz=TARGET.CEIL)
        with self.assertRaisesRegex(ValueError,'sample-entry WPE'):
            X.imu_step_from_wpe(s,tuner_frequency=source,tau_exp_decay=B32.rn32(F(199,200)))

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

    def test_readiness_closes_goLive_handoff_but_not_admitted_startup_or_libm(self):
        r=X.readiness()
        self.assertTrue(r['Live_product_carries_persistent_dual_compiler_tau_ledger'])
        self.assertTrue(r['strong_Live_constructor_requires_exact_goLive_filter_frontend_state_and_tau_ledger'])
        self.assertTrue(r['goLive_tau_ledger_identity_bridge_available'])
        self.assertTrue(r['each_Live_IMU_requires_same_candidate_frequency_target_and_decay_as_tau_ledger'])
        self.assertTrue(r['strong_Live_IMU_frequency_source_bound_to_sample_entry_WPE_or_prior'])
        self.assertTrue(r['WPE_getter_to_tuner_store_topology_available'])
        self.assertTrue(r['MAG_and_HOLD_preserve_tau_ledger_exactly'])
        self.assertFalse(r['admitted_startup_reachability_with_tau_ledger_closed'])
        self.assertFalse(r['WPE_binary32_log_period_production_closed'])
        self.assertFalse(r['upstream_WPE_to_StoredFrequency_binary32_correspondence_closed'])
        self.assertFalse(r['tuner_exp_libm_binary32_correspondence_closed'])
        self.assertFalse(r['source_uniform_complete_600_step_word_qualified'])
        self.assertFalse(r['storage_search_allowed'])


if __name__=='__main__': unittest.main()
