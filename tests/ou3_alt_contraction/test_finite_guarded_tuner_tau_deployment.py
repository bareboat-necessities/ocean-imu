"""Startup guarded-frontend tau deployment regressions."""
from dataclasses import replace
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_accel_guard_runtime as G
from tools.stability.ou3_alt_contraction import finite_guarded_tuner_prefix as FRONT
from tools.stability.ou3_alt_contraction import finite_guarded_tuner_tau_deployment as X
from tools.stability.ou3_alt_contraction import finite_tuner_tau_deployment_ledger as LEDGER
from tools.stability.ou3_alt_contraction import finite_tuner_frequency_binary32 as FREQ
from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B32
from tools.stability.ou3_alt_contraction import finite_shipping_tau_target_binary32 as TARGET
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as C
import test_finite_guarded_tuner_prefix as BASE


def cold_state():
    t=replace(BASE.tuner_state(),stage='Cold',stage_time=0,warmup_sec=5)
    return X.State(FRONT.State(G.State(),t),LEDGER.initial())


def cold_kwargs():
    k=BASE.kwargs()
    for name in ('candidate_cfg','sigma_wave_sqrt','spectral','ema'):
        k.pop(name,None)
    return k


def live_kwargs(decay):
    k=BASE.kwargs()
    # Make the pre-WPE external frequency .5 Hz.  The carried stats state is
    # already ready at .5 Hz, so the same .5 reaches SeaStateAutoTuner and the
    # candidate.  Source tau target is then exactly .5/.5 = 1 s.
    k['band_cfg']=replace(k['band_cfg'],tune_freq_prior=F(1,2))
    c=k['candidate_cfg']
    k['candidate_cfg']=replace(c,min_freq=TARGET.FLOOR,max_freq=TARGET.CEIL,
        tau_coeff=B32.rn32(1),min_tau=TARGET.TAU_MIN,max_tau=TARGET.TAU_MAX,
        adapt_tau_sec=LEDGER.ADAPT_SEC,adapt_tau_sea_periods=LEDGER.ADAPT_PERIODS,
        clamp_enabled=True)
    k['ema']=C.EmaWitness(decay,1)
    k['spectral']=C.SpectralWitness(1,1)
    return k


def stored_half():
    return FREQ.store(B32.rn32(F(1,2)),B32.rn32(F(1,20)),B32.rn32(5))


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

    def test_source_qualified_live_candidate_advances_both_tau_tracks(self):
        # x=.005/.4=.0125; this binary32 value lies in the rigorous exp enclosure.
        decay=B32.rn32(F(98755,100000))
        s=X.State(FRONT.State(G.State(),BASE.tuner_state()),LEDGER.initial())
        out=X.step(s,BASE.packet(),dt=BASE.DT,guard_cfg=G.Config(),
                   stored_frequency=stored_half(),tau_exp_decay=decay,
                   **live_kwargs(decay))
        self.assertIsNotNone(out.frontend.tuner.candidate)
        self.assertEqual(out.frontend.tuner.candidate.frequency,F(1,2))
        self.assertEqual(out.frontend.tuner.candidate.tau_target,F(1))
        self.assertEqual(out.state.tau.updates,1)
        self.assertEqual(out.state.tau.separate,out.tau_step.separate_step.next_separate)
        self.assertEqual(out.state.tau.fma,out.tau_step.fma_step.next_fma)

    def test_postcold_exp_witness_must_equal_exact_candidate_decay(self):
        decay=B32.rn32(F(98755,100000))
        s=X.State(FRONT.State(G.State(),BASE.tuner_state()),LEDGER.initial())
        with self.assertRaisesRegex(ValueError,'detached from exact candidate'):
            X.step(s,BASE.packet(),dt=BASE.DT,guard_cfg=G.Config(),
                   stored_frequency=stored_half(),tau_exp_decay=B32.rn32(F(99,100)),
                   **live_kwargs(decay))

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
