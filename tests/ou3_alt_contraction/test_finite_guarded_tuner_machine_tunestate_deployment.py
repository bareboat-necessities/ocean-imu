"""Guarded startup whole-machine TuneState product regressions."""
from dataclasses import replace
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_accel_guard_runtime as G
from tools.stability.ou3_alt_contraction import finite_tuner_commit as COMMIT
from tools.stability.ou3_alt_contraction import finite_tuner_deployment_config as D
from tools.stability.ou3_alt_contraction import finite_tuner_machine_tunestate_product as PRODUCT
from tools.stability.ou3_alt_contraction import finite_guarded_tuner_machine_tunestate_deployment as X
import test_finite_guarded_tuner_wpe_tau_deployment as BASE
import test_finite_guarded_tuner_prefix as RAW
from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B


def dcfg(): return D.shipping_defaults(qeff_pow_result=B.rn32(1))

def commit_cfg():
    return COMMIT.CommitConfig(F(3,22),F(1,200),F(3,20),F(3,200),True,False,
                               F(3,20),100,F(4,5),F(6,5),F(7,5))


def cold_machine(s):
    from test_finite_machine_frontend_sigma_source import source_step
    k=BASE.cold_kwargs()
    a=dict(band_cfg=k['band_cfg'],stats_cfg=k['stats_cfg'],frequency=B.rn32(F(1,5)),dt=RAW.DT,bench_noise_sigma=k['bench_noise_sigma'])
    return dict(separate_frontend=source_step(s.frontends.separate,**a),fma_frontend=source_step(s.frontends.fma,**a,last=True))


class Tests(unittest.TestCase):
    def test_consecutive_Cold_samples_do_not_restart_machine_frontend(self):
        s=X.initial(BASE.construction_frontend())
        first=X.step(s,RAW.packet(),dt=RAW.DT,deployment_cfg=dcfg(),guard_cfg=G.Config(),**cold_machine(s),**BASE.cold_kwargs())
        second=X.step(first.state,RAW.packet(),dt=RAW.DT,deployment_cfg=dcfg(),guard_cfg=G.Config(),
            guard_decay=G.DecayWitness(1,1,0,0),guard_rms=G.RmsWitness(0),**cold_machine(first.state),**BASE.cold_kwargs())
        self.assertEqual(second.state.frontends.samples,2)
        self.assertIs(second.separate_frontend.before,first.state.frontends.separate)
        self.assertEqual(second.state.machine,s.machine)
        with self.assertRaisesRegex(ValueError,'sample count detached'):
            replace(second.state,frontends=s.frontends)

    def test_construction_roots_whole_machine_product_without_synthetic_state(self):
        s=X.initial(BASE.construction_frontend())
        self.assertEqual(s.machine,PRODUCT.initial())
        self.assertEqual(s.machine.tau,s.lower.tau)
        self.assertFalse(s.machine.pending)
        self.assertEqual(s.machine.pending,s.lower.frontend.tuner.pending)

    def test_cold_sample_advances_lower_WPE_but_holds_whole_TuneState(self):
        s=X.initial(BASE.construction_frontend()); before=s.machine
        out=X.step(s,RAW.packet(),dt=RAW.DT,deployment_cfg=dcfg(),guard_cfg=G.Config(),**cold_machine(s),**BASE.cold_kwargs())
        self.assertEqual(out.state.machine,before)
        self.assertEqual(out.state.lower.wpe.samples,1)
        self.assertEqual(out.state.frontends.samples,1)
        self.assertTrue(out.state.frontends.separate.band.machine.ready)
        self.assertIsNone(out.machine)
        self.assertIsNone(out.separate_sigma_join); self.assertIsNone(out.fma_sigma_join)

    def test_cold_sample_rejects_sigma_or_RS_machine_witnesses(self):
        s=X.initial(BASE.construction_frontend())
        with self.assertRaisesRegex(ValueError,'consumes no sigma/R_S machine witnesses'):
            X.step(s,RAW.packet(),dt=RAW.DT,deployment_cfg=dcfg(),guard_cfg=G.Config(),
                   separate_spectral_pow=B.rn32(1),**cold_machine(s),**BASE.cold_kwargs())

    def test_product_state_rejects_pending_or_tau_detachment(self):
        s=X.initial(BASE.construction_frontend())
        with self.assertRaisesRegex(ValueError,'pending bits differ'):
            X.State(s.lower,replace(s.machine,pending=True))
        with self.assertRaisesRegex(ValueError,'tau history detached'):
            X.State(s.lower,replace(s.machine,tau=replace(s.machine.tau,separate=B.rn32(1))))

    def test_no_pending_combined_boundary_is_identity_on_both_products(self):
        s=X.initial(BASE.construction_frontend())
        out=X.boundary(s,commit_cfg(),bench_noise_sigma=F(1,100))
        self.assertEqual(out.state,s)
        self.assertIsNone(out.exact.commit)
        self.assertFalse(out.machine.consumed)

    def test_preLive_pending_boundary_consumes_binary32_outputs_and_preserves_ledgers(self):
        from test_finite_startup_live_machine_tunestate_bridge import startup
        from tools.stability.ou3_alt_contraction import finite_band_variance_runtime as BAND
        s,*_=startup(True)
        from test_finite_machine_frontend_sigma_source import ready_pair
        s=replace(s,frontends=ready_pair())
        bench=B.rn32(F(1,10)); sb,fb=B.mul(bench,11),B.mul(bench,12)
        out=X.boundary(s,commit_cfg(),bench_noise_sigma=bench,
            exact_noise_sqrt=BAND.NoiseSqrtWitness(0),
            separate_noise_sqrt_gain=11,fma_noise_sqrt_gain=12)
        self.assertFalse(out.state.machine.pending)
        self.assertFalse(out.state.lower.frontend.tuner.pending)
        self.assertIs(out.state.machine.tau,s.machine.tau)
        self.assertIs(out.state.machine.sigma,s.machine.sigma)
        self.assertIs(out.state.machine.rs,s.machine.rs)
        self.assertIs(out.state.lower.wpe,s.lower.wpe)
        self.assertEqual(out.machine.separate_commit.Sigma_aw[2][2],B.mul(sb,sb))
        self.assertEqual(out.machine.fma_commit.Sigma_aw[2][2],B.mul(fb,fb))
        self.assertIsNone(out.machine.separate_commit.R_S)
        self.assertIsNone(out.machine.fma_commit.R_S)
        self.assertNotEqual(out.machine.separate_commit.Sigma_aw,out.exact.commit.Sigma_aw)

    def test_readiness_attaches_startup_topology_but_not_reachability_or_storage(self):
        r=X.readiness()
        self.assertTrue(r['startup_frontend_machine_TuneState_product_attached'])
        self.assertTrue(r['combined_next_IMU_boundary_clears_exact_and_machine_pending_together'])
        for k in ('upstream_frontend_binary32_correspondence_closed','all_target_libm_correspondence_closed',
                  'source_uniform_machine_supply_bounds_closed','every_admitted_startup_history_reaches_TunerReady_with_this_product',
                  'goLive_carries_whole_machine_TuneState_product','source_uniform_complete_600_step_word_qualified',
                  'storage_search_allowed','ALT_STARTUP_PASS','ALT_LIVE_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
