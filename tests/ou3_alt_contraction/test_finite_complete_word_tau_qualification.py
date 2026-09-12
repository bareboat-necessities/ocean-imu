"""Tau-specific complete-word source-qualification regressions."""
from dataclasses import replace
import unittest

from tools.stability.ou3_alt_contraction import finite_complete_word_tau_qualification as X
from tools.stability.ou3_alt_contraction import finite_shipping_tau_target_binary32 as TARGET
import test_finite_source_bound_live_word as BASE


def shipping_runtime():
    r=BASE.root_state().runtime; c=r.candidate_cfg
    c=replace(c,min_freq=TARGET.FLOOR,max_freq=TARGET.CEIL,tau_coeff=X.B.rn32(1),
              min_tau=TARGET.TAU_MIN,max_tau=TARGET.TAU_MAX,
              adapt_tau_sec=X.ADAPT_SEC,adapt_tau_sea_periods=X.ADAPT_PERIODS,
              clamp_enabled=True)
    return replace(r,candidate_cfg=c)


class Tests(unittest.TestCase):
    def test_actual_shipping_tau_fields_qualify_persistent_runtime(self):
        r=shipping_runtime()
        self.assertIs(X.qualify_runtime(r),r)

    def test_component_toy_candidate_config_cannot_enter_shipping_tau_theorem(self):
        r=BASE.root_state().runtime
        with self.assertRaisesRegex(ValueError,'detached from shipping binary32 default'):
            X.qualify_runtime(r)

    def test_each_tau_relevant_field_is_fail_closed(self):
        r=shipping_runtime(); c=r.candidate_cfg
        for name in ('min_freq','max_freq','tau_coeff','min_tau','max_tau',
                     'adapt_tau_sec','adapt_tau_sea_periods'):
            bad=replace(c,**{name:getattr(c,name)+1})
            with self.assertRaisesRegex(ValueError,name):
                X.qualify_runtime(replace(r,candidate_cfg=bad))
        with self.assertRaisesRegex(ValueError,'clamp-enabled'):
            X.qualify_runtime(replace(r,candidate_cfg=replace(c,clamp_enabled=False)))

    def test_readiness_exposes_available_bound_but_not_per_event_attachment(self):
        q=X.readiness()
        self.assertTrue(q['tau_relevant_persistent_runtime_config_source_qualified'])
        self.assertTrue(q['shipping_tau_target_binary32_cell_bound_available'])
        self.assertTrue(q['shipping_tau_predecessor_domain_inductively_closed_for_600_step_word'])
        self.assertTrue(q['source_uniform_tau_roundoff_supply_bound_available_for_600_step_word'])
        self.assertFalse(q['every_IMU_event_retains_binary32_tau_step_and_supply_certificate'])
        self.assertFalse(q['source_uniform_complete_600_step_word_qualified'])
        self.assertFalse(q['storage_search_allowed'])


if __name__=='__main__': unittest.main()
