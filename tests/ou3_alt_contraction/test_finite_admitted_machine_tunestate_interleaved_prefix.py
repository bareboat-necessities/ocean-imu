"""Admitted Live whole-machine TuneState interleaver regressions."""
from dataclasses import replace
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_admitted_machine_tunestate_interleaved_prefix as X
from tools.stability.ou3_alt_contraction import finite_admitted_wpe_tau_interleaved_prefix as LOWER
from tools.stability.ou3_alt_contraction import finite_tuner_machine_tunestate_product as PRODUCT
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_deployment_ledger as SIG
from tools.stability.ou3_alt_contraction import finite_tuner_rs_deployment_ledger as RS
from tools.stability.ou3_alt_contraction import finite_tuner_deployment_config as D
from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
import test_finite_admitted_wpe_tau_interleaved_prefix as BASE
import test_finite_source_bound_live_word as LBASE


def dcfg(): return D.shipping_defaults(qeff_pow_result=B.rn32(1))

def state():
    b0=BASE.base_state(); b=LOWER.begin(b0,BASE.machine_wpe_for(b0))
    n=b.base.tau.updates
    exact=X._entry_live(b)
    m=PRODUCT.State(b.base.tau,SIG.State(updates=n),RS.State(updates=n),exact.tuner.pending)
    return X.begin(b,m,dcfg(),separate_active=exact.active,fma_active=exact.active)


class Tests(unittest.TestCase):
    def test_product_rejects_tau_pending_and_config_splices(self):
        s=state()
        with self.assertRaisesRegex(ValueError,'tau ledger detached'):
            X.State(s.base,replace(s.machine,tau=replace(s.machine.tau,separate=B.rn32(1))),
                    s.deployment_cfg,s.separate_active,s.fma_active,s.live_entry_machine_updates)
        with self.assertRaisesRegex(ValueError,'pending bit detached'):
            X.State(s.base,replace(s.machine,pending=not s.machine.pending),
                    s.deployment_cfg,s.separate_active,s.fma_active,s.live_entry_machine_updates)
        bad=replace(s.deployment_cfg,sigma_coeff=B.rn32(F(4,5)))
        with self.assertRaisesRegex(ValueError,'config detached'):
            X.State(s.base,s.machine,bad,s.separate_active,s.fma_active,s.live_entry_machine_updates)

    def test_applied_parameters_are_distinct_state_from_candidate_memory(self):
        s=state()
        self.assertEqual(s.separate_active_join.supply.tau,0)
        self.assertEqual(s.fma_active_join.supply.tau,0)
        changed=replace(s.machine,sigma=replace(s.machine.sigma,separate=B.rn32(F(1,5))))
        # A candidate-memory change alone does not rewrite already-applied active
        # parameters; it would only take effect after a later pending boundary.
        altered=X.State(s.base,changed,s.deployment_cfg,s.separate_active,s.fma_active,
                        s.live_entry_machine_updates)
        self.assertEqual(altered.separate_active,s.separate_active)

    def test_MAG_and_HOLD_preserve_candidate_and_applied_machine_state(self):
        s=state(); m=s.machine; a=s.separate_active; w=s.base.wpe
        s,_=X.mag_step(s,**LBASE.mag_kwargs(s.base.base.prefix.prefix.live.live_word))
        self.assertIs(s.machine,m); self.assertIs(s.separate_active,a); self.assertIs(s.base.wpe,w)
        s,_=X.set_hold(s,hold=False)
        self.assertIs(s.machine,m); self.assertIs(s.separate_active,a); self.assertIs(s.base.wpe,w)

    def test_complete_word_cannot_be_claimed_without_600_common_machine_updates(self):
        s=state()
        with self.assertRaises((ValueError,TypeError)):
            X.complete(s)

    def test_readiness_attaches_full_tuner_state_but_keeps_live_coefficients_open(self):
        r=X.readiness()
        self.assertTrue(r['whole_tau_sigma_RS_machine_TuneState_carried_in_same_Live_product'])
        self.assertTrue(r['candidate_memory_and_applied_machine_parameters_carried_separately'])
        self.assertTrue(r['exact_vs_machine_applied_parameter_supplies_exposed_every_Live_state'])
        self.assertTrue(r['Live_600_step_machine_TuneState_product_attached'])
        self.assertTrue(r['MAG_and_HOLD_preserve_WPE_TuneState_and_applied_machine_parameters_by_identity'])
        self.assertFalse(r['machine_active_parameter_displacement_injected_into_Live_coefficients'])
        self.assertFalse(r['source_uniform_machine_supply_bounds_closed'])
        self.assertFalse(r['source_uniform_complete_600_step_word_qualified'])
        self.assertFalse(r['storage_search_allowed'])
        self.assertFalse(r['ALT_LIVE_PASS']); self.assertFalse(r['ALT_END_TO_END_PASS'])


if __name__=='__main__': unittest.main()
