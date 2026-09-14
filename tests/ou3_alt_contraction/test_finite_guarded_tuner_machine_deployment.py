"""Guarded startup whole-machine TuneState regressions."""
from dataclasses import replace
import unittest

from tools.stability.ou3_alt_contraction import finite_accel_guard_runtime as G
from tools.stability.ou3_alt_contraction import finite_guarded_tuner_prefix as FRONT
from tools.stability.ou3_alt_contraction import finite_guarded_tuner_machine_deployment as X
from tools.stability.ou3_alt_contraction import finite_tuner_tau_deployment_ledger as TAU
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_deployment_ledger as SIG
from tools.stability.ou3_alt_contraction import finite_tuner_rs_deployment_ledger as RS
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
    def test_construction_roots_whole_machine_state_before_sample_one(self):
        s=X.initial(construction_frontend())
        self.assertEqual(s.machine.tau,TAU.initial())
        self.assertEqual(s.machine.sigma,SIG.initial())
        self.assertEqual(s.machine.rs,RS.initial())
        self.assertFalse(s.machine.pending)
        self.assertEqual(s.base.tau,s.machine.tau)

    def test_cold_sample_advances_WPE_but_preserves_sigma_RS_pending(self):
        s=X.initial(construction_frontend())
        out=X.step(s,BASE.packet(),dt=BASE.DT,guard_cfg=G.Config(),**cold_kwargs())
        self.assertEqual(out.state.base.wpe.samples,1)
        self.assertEqual(out.state.machine.tau,s.machine.tau)
        self.assertEqual(out.state.machine.sigma,s.machine.sigma)
        self.assertEqual(out.state.machine.rs,s.machine.rs)
        self.assertEqual(out.state.machine.pending,s.machine.pending)
        self.assertIsNone(out.machine_result)

    def test_cold_sample_rejects_detached_sigma_or_RS_result(self):
        s=X.initial(construction_frontend())
        with self.assertRaisesRegex(ValueError,'consumes no sigma/R_S'):
            X.step(s,BASE.packet(),dt=BASE.DT,guard_cfg=G.Config(),
                   sigma_result=object(),**cold_kwargs())

    def test_machine_tau_cannot_detach_from_lower_startup_history(self):
        s=X.initial(construction_frontend())
        bad=replace(s.machine,tau=TAU.State(TAU.INITIAL,TAU.INITIAL,1))
        with self.assertRaisesRegex(ValueError,'tau detached'):
            X.State(s.base,bad)

    def test_readiness_closes_product_topology_not_sources_or_reachability(self):
        r=X.readiness()
        for k in ('construction_roots_whole_tau_sigma_RS_machine_TuneState',
                  'machine_tau_is_exact_same_history_as_startup_WPE_tau_product',
                  'Cold_noncandidate_samples_preserve_sigma_RS_and_pending',
                  'candidate_samples_require_tau_sigma_RS_to_advance_from_one_predecessor',
                  'exact_candidate_pending_bit_controls_whole_machine_successor',
                  'machine_minus_exact_sigma_and_RS_target_supplies_retained_per_mode',
                  'startup_global_compiler_track_coherence_includes_WPE_tau_sigma_RS'):
            self.assertTrue(r[k])
        for k in ('upstream_sigma_binary32_source_production_closed','upstream_RS_binary32_source_production_closed',
                  'all_target_libm_correspondence_closed','every_admitted_startup_history_reaches_TunerReady_with_this_product',
                  'goLive_carries_this_exact_WPE_machine_TuneState_product','storage_search_allowed','ALT_STARTUP_PASS','ALT_LIVE_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
