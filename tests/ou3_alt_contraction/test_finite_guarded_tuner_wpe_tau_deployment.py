"""Coherent guarded-startup WPE/tau deployment regressions."""
from dataclasses import replace
import unittest

from tools.stability.ou3_alt_contraction import finite_accel_guard_runtime as G
from tools.stability.ou3_alt_contraction import finite_guarded_tuner_prefix as FRONT
from tools.stability.ou3_alt_contraction import finite_guarded_tuner_wpe_tau_deployment as X
from tools.stability.ou3_alt_contraction import finite_tuner_tau_deployment_ledger as TAU
import test_finite_guarded_tuner_prefix as BASE


def construction_frontend():
    t=BASE.tuner_state()
    t=replace(t,tune=replace(t.tune,tau_applied=TAU.INITIAL),stage='Cold',stage_time=0,
              warmup_sec=5,sample_index=0)
    return FRONT.State(G.State(),t)


def cold_kwargs():
    k=BASE.kwargs()
    for name in ('candidate_cfg','sigma_wave_sqrt','spectral','ema'):
        k.pop(name,None)
    return k


class Tests(unittest.TestCase):
    def test_construction_root_owns_both_machine_ledgers_before_sample_one(self):
        s=X.initial(construction_frontend())
        self.assertEqual(s.tau,TAU.initial())
        self.assertEqual(s.wpe.samples,0)
        self.assertIsNone(s.wpe.separate.log_period); self.assertIsNone(s.wpe.fma.log_period)

    def test_cold_sample_advances_exact_frontend_and_WPE_sample_but_not_tau(self):
        s=X.initial(construction_frontend()); tau=s.tau
        out=X.step(s,BASE.packet(),dt=BASE.DT,guard_cfg=G.Config(),**cold_kwargs())
        self.assertEqual(out.state.frontend.tuner.sample_index,1)
        self.assertEqual(out.state.wpe.samples,1)
        self.assertIs(out.state.tau,tau)
        self.assertIsNone(out.tau_step)
        self.assertFalse(out.wpe_step.produced_period)
        self.assertIsNone(out.separate_supply); self.assertIsNone(out.fma_supply)

    def test_cold_sample_rejects_tau_machine_operands(self):
        s=X.initial(construction_frontend())
        with self.assertRaisesRegex(ValueError,'consumes no tau deployment witnesses'):
            X.step(s,BASE.packet(),dt=BASE.DT,guard_cfg=G.Config(),
                   separate_tau_exp_decay=1,fma_tau_exp_decay=1,**cold_kwargs())

    def test_arbitrary_initialized_machine_WPE_cannot_be_spliced_at_construction(self):
        from tools.stability.ou3_alt_contraction import finite_wpe_log_binary32 as W
        from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
        from fractions import Fraction as F
        front=construction_frontend(); t=W.Track(B.rn32(F(1,2)),1)
        with self.assertRaisesRegex(ValueError,'initialization detached'):
            X.State(front,TAU.initial(),W.State(t,t,0))

    def test_readiness_keeps_reachability_and_libm_open(self):
        r=X.readiness()
        for k in ('construction_seed_roots_both_tau_and_WPE_machine_ledgers',
                  'Cold_samples_preserve_tau_but_advance_WPE_machine_sample_history',
                  'postCold_sample_entry_WPE_frequency_precedes_tau_candidate_and_current_WPE_update',
                  'startup_global_compiler_track_coherence_includes_WPE_log_and_tau_states'):
            self.assertTrue(r[k])
        self.assertFalse(r['WPE_log_std_log_target_libm_correspondence_closed'])
        self.assertFalse(r['tuner_exp_libm_binary32_correspondence_closed'])
        self.assertFalse(r['every_admitted_startup_history_reaches_TunerReady_with_this_product'])
        self.assertFalse(r['goLive_carries_this_exact_WPE_tau_product'])
        self.assertFalse(r['storage_search_allowed'])


if __name__=='__main__': unittest.main()
