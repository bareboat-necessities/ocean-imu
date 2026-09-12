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
    p=BASE.begin(); source=p.prefix.live.live_word
    source=replace(source,runtime=QBASE.shipping_runtime())
    admitted=replace(p.prefix.live,live_word=source)
    return replace(p,prefix=replace(p.prefix,live=admitted))


def wpe_source_for(state,dynamic):
    entry=X._entry_wpe(state)
    if not entry.usable_period:
        return WPEF.tuner_frequency(entry,min_hz=TARGET.FLOOR,max_hz=TARGET.CEIL)
    exact_f=F(dynamic['wpe_current_frequency']); frequency=B32.rn32(exact_f)
    log=WPEF.bind_log_state(entry,B32.rn32(entry.log_period))
    getter=WPEF.getters(log,period_exp=B32.rn32(F(1)/frequency),frequency_exp=frequency)
    return WPEF.tuner_frequency(entry,min_hz=TARGET.FLOOR,max_hz=TARGET.CEIL,
                                getter=getter,shadow_frequency=exact_f)


class Tests(unittest.TestCase):
    def test_component_toy_runtime_cannot_enter_tau_theorem_product(self):
        with self.assertRaisesRegex(ValueError,'detached from shipping binary32 default'):
            X.begin(BASE.begin(),LEDGER.initial())

    def test_source_qualified_live_product_accepts_startup_carried_ledger_for_component_algebra(self):
        p=qualified_prefix(); tau=LEDGER.State(updates=17); s=X.begin(p,tau)
        self.assertIs(s.prefix,p); self.assertEqual(s.tau,tau); self.assertEqual(s.live_entry_tau_updates,17)

    def test_WPE_source_edge_retains_exact_vs_machine_input_supplies(self):
        s=X.begin(qualified_prefix(),LEDGER.State(updates=17))
        witness,segment,raw,r,b,dynamic=IBASE.operands(s.prefix.prefix.live); dynamic=dict(dynamic)
        source=wpe_source_for(s,dynamic)
        # Keep the machine decay inside the shipping enclosure; it need not be
        # exactly identical to the exact-real candidate decay.
        e=B32.rn32(F(dynamic['ema'].decay_tau_sigma))
        out=X.imu_step_from_wpe(s,tuner_frequency=source,tau_exp_decay=e,
            restricted=r,bias_restricted=b,witness=witness,raw=raw,
            packet_id='imu-wpe-source',**PBASE.root_args(),**dynamic)
        self.assertEqual(out.state.tau.updates,18)
        self.assertEqual(out.frequency_supply,source.machine_minus_shadow)
        self.assertEqual(out.decay_input_supply,e-F(dynamic['ema'].decay_tau_sigma))
        machine_exact=out.tau_step.separate_target.exact_target
        cand=X._candidate(out.event)
        self.assertEqual(out.target_input_supply,machine_exact-F(cand.tau_target))
        self.assertEqual(source.shadow,X._entry_wpe(s))

    def test_detached_preupdate_WPE_frequency_source_rejected_before_execution(self):
        s=X.begin(qualified_prefix(),LEDGER.State(updates=17)); entry=X._entry_wpe(s)
        other=replace(entry,elapsed=entry.elapsed+1)
        if other.usable_period:
            log=WPEF.bind_log_state(other,B32.rn32(other.log_period))
            getter=WPEF.getters(log,period_exp=B32.rn32(2),frequency_exp=B32.rn32(F(1,2)))
            source=WPEF.tuner_frequency(other,min_hz=TARGET.FLOOR,max_hz=TARGET.CEIL,
                                        getter=getter,shadow_frequency=F(1,2))
        else:
            source=WPEF.tuner_frequency(other,min_hz=TARGET.FLOOR,max_hz=TARGET.CEIL)
        with self.assertRaisesRegex(ValueError,'sample-entry WPE'):
            X.imu_step_from_wpe(s,tuner_frequency=source,tau_exp_decay=B32.rn32(F(199,200)))

    def test_mag_and_hold_are_literal_tau_ledger_identities(self):
        s=X.begin(qualified_prefix(),LEDGER.State(updates=17)); tau=s.tau
        s,_=X.mag_step(s,**LBASE.mag_kwargs(s.prefix.prefix.live.live_word))
        self.assertIs(s.tau,tau); s,_=X.set_hold(s,hold=False); self.assertIs(s.tau,tau)

    def test_live_entry_must_leave_room_for_full_600_step_word(self):
        with self.assertRaisesRegex(ValueError,'no certified room'):
            X.begin(qualified_prefix(),LEDGER.State(updates=LEDGER.MAX_UPDATES-599))

    def test_short_prefix_cannot_claim_tau_complete_word(self):
        with self.assertRaisesRegex(ValueError,'exactly 600 physical source transitions'):
            X.complete(X.begin(qualified_prefix(),LEDGER.State(updates=17)))

    def test_readiness_is_precise_about_remaining_global_WPE_join(self):
        r=X.readiness()
        for k in ('Live_product_carries_persistent_dual_compiler_tau_ledger',
                  'strong_Live_constructor_requires_exact_goLive_filter_frontend_state_and_tau_ledger',
                  'WPE_source_edge_keeps_exact_candidate_and_machine_frequency_distinct',
                  'WPE_machine_minus_shadow_frequency_supply_attached',
                  'machine_exact_target_minus_exact_candidate_target_supply_attached',
                  'machine_tau_decay_minus_exact_candidate_decay_supply_attached',
                  'global_compiler_tau_tracks_accept_distinct_WPE_inputs'):
            self.assertTrue(r[k])
        self.assertFalse(r['global_compiler_WPE_log_tracks_composed_into_Live_product'])
        self.assertFalse(r['WPE_binary32_log_period_production_closed'])
        self.assertFalse(r['source_uniform_WPE_frequency_supply_bound_closed'])
        self.assertFalse(r['source_uniform_complete_600_step_word_qualified'])
        self.assertFalse(r['storage_search_allowed'])


if __name__=='__main__': unittest.main()
