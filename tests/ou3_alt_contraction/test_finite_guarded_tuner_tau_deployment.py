"""Startup guarded-frontend tau deployment regressions."""
from dataclasses import replace
import unittest

from tools.stability.ou3_alt_contraction import finite_accel_guard_runtime as G
from tools.stability.ou3_alt_contraction import finite_guarded_tuner_prefix as FRONT
from tools.stability.ou3_alt_contraction import finite_guarded_tuner_tau_deployment as X
from tools.stability.ou3_alt_contraction import finite_tuner_tau_deployment_ledger as LEDGER
import test_finite_guarded_tuner_prefix as BASE


def cold_state():
    t=replace(BASE.tuner_state(),stage='Cold',stage_time=0,warmup_sec=5)
    return X.State(FRONT.State(G.State(),t),LEDGER.initial())


def cold_kwargs():
    k=BASE.kwargs()
    for name in ('candidate_cfg','sigma_wave_sqrt','spectral','ema'):
        k.pop(name,None)
    return k


class Tests(unittest.TestCase):
    def test_cold_frontend_advances_but_tau_ledger_is_identity(self):
        s=cold_state(); tau=s.tau
        out=X.step(s,BASE.packet(),dt=BASE.DT,guard_cfg=G.Config(),**cold_kwargs())
        self.assertEqual(out.state.frontend.tuner.sample_index,1)
        self.assertIs(out.state.tau,tau)
        self.assertIsNone(out.tau_step)

    def test_cold_branch_rejects_free_tau_deployment_witnesses(self):
        s=cold_state()
        with self.assertRaisesRegex(ValueError,'consumes no tau deployment witnesses'):
            X.step(s,BASE.packet(),dt=BASE.DT,guard_cfg=G.Config(),
                   stored_frequency='detached',tau_exp_decay=1,**cold_kwargs())

    def test_initial_wrapper_uses_shipping_tau_seed(self):
        raw=FRONT.State(G.State(),replace(BASE.tuner_state(),stage='Cold'))
        s=X.initial(raw)
        self.assertEqual(s.tau,LEDGER.initial())

    def test_readiness_keeps_upstream_machine_facts_open(self):
        r=X.readiness()
        self.assertTrue(r['Cold_noncandidate_branch_preserves_tau_ledger'])
        self.assertTrue(r['postCold_candidate_advances_both_global_compiler_tau_tracks'])
        self.assertFalse(r['upstream_WPE_to_StoredFrequency_binary32_correspondence_closed'])
        self.assertFalse(r['tuner_exp_libm_binary32_correspondence_closed'])
        self.assertFalse(r['admitted_startup_source_history_carries_this_frontend_product'])
        self.assertFalse(r['storage_search_allowed'])


if __name__=='__main__': unittest.main()
