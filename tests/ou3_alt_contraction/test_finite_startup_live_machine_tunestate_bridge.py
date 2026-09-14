"""Whole-machine TuneState goLive bridge regressions."""
from dataclasses import replace
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_guarded_tuner_machine_tunestate_deployment as START
from tools.stability.ou3_alt_contraction import finite_guarded_tuner_wpe_tau_deployment as LOWER
from tools.stability.ou3_alt_contraction import finite_tuner_tau_deployment_ledger as TAU
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_deployment_ledger as SIG
from tools.stability.ou3_alt_contraction import finite_tuner_rs_deployment_ledger as RS
from tools.stability.ou3_alt_contraction import finite_tuner_machine_tunestate_product as PRODUCT
from tools.stability.ou3_alt_contraction import finite_wpe_log_binary32 as WPE
from tools.stability.ou3_alt_contraction import finite_racc_runtime as RACC
from tools.stability.ou3_alt_contraction import finite_periodic_aw_sync as AWSYNC
from tools.stability.ou3_alt_contraction import finite_band_variance_runtime as BAND
from tools.stability.ou3_alt_contraction import finite_startup_live_machine_tunestate_bridge as X
from tools.stability.ou3_alt_contraction import deployment_scope as SCOPE
import test_finite_startup_live_runtime_bridge as BASE


def startup(pending=True):
    front,entry,fresh,active,scheduler=BASE.objects(pending)
    # This fixture isolates the goLive identity edge: all three machine values
    # are chosen equal to the exact fixture's TuneState.  Construction-to-ready
    # provenance is tested in the guarded startup module, not fabricated here as
    # a claim about universal startup reachability.
    tau=TAU.State(F(1),F(1),0); sig=SIG.State(F(1),F(1),0); rs=RS.State(F(1),F(1),0)
    lower=LOWER.State(front,tau,WPE.initial())
    strong=START.State(lower,PRODUCT.State(tau,sig,rs,pending))
    return strong,entry,fresh,active,scheduler


class Tests(unittest.TestCase):
    def test_goLive_preserves_whole_machine_history_and_pending_identity(self):
        s,entry,fresh,active,scheduler=startup(True)
        out=X.bridge(s,entry,fresh,scope=SCOPE.certified_scope(),commit_cfg=BASE.cfg(),
            bench_noise_sigma=0,noise_sqrt=BAND.NoiseSqrtWitness(0),
            scheduler=scheduler,racc=RACC.State(),aw_sync=AWSYNC.State())
        self.assertIs(out.machine,s.machine)
        self.assertIs(out.frontends,s.frontends)
        self.assertEqual(out.wpe,s.lower.wpe)
        self.assertTrue(out.live.frontend_live.tuner.pending)
        self.assertTrue(out.machine.pending)
        self.assertEqual(out.separate_active.tau,active.tau)
        self.assertEqual(out.separate_active.pseudo_period,B.rn32(active.pseudo_period))
        self.assertEqual(out.separate_active,out.fma_active)
        self.assertEqual(out.separate_commit.aw_floor_target,out.separate_active.Sigma_aw)
        self.assertIs(out.arithmetic.arithmetic.before.tau,s.machine.tau)

    def test_goLive_applies_even_when_online_pending_is_false_and_preserves_false(self):
        s,entry,fresh,active,scheduler=startup(False)
        out=X.bridge(s,entry,fresh,scope=SCOPE.certified_scope(),commit_cfg=BASE.cfg(),
            bench_noise_sigma=0,noise_sqrt=BAND.NoiseSqrtWitness(0),
            scheduler=scheduler,racc=RACC.State(),aw_sync=AWSYNC.State())
        self.assertFalse(out.machine.pending)
        self.assertFalse(out.live.frontend_live.tuner.pending)
        self.assertIsNotNone(out.separate_commit); self.assertIsNotNone(out.fma_commit)

    def test_goLive_ready_band_requires_sqrt_of_its_own_covariance(self):
        s,entry,fresh,active,scheduler=startup(False)
        from test_finite_machine_frontend_sigma_source import ready_pair
        s=replace(s,frontends=ready_pair())
        with self.assertRaisesRegex(ValueError,'requires sqrt'):
            X.bridge(s,entry,fresh,scope=SCOPE.certified_scope(),commit_cfg=BASE.cfg(),
                bench_noise_sigma=0,noise_sqrt=BAND.NoiseSqrtWitness(0),
                scheduler=scheduler,racc=RACC.State(),aw_sync=AWSYNC.State())

    def test_goLive_retains_roundoff_and_distinct_mode_noise_floors(self):
        s,entry,fresh,active,scheduler=startup(False)
        from test_finite_machine_frontend_sigma_source import ready_pair
        s=replace(s,frontends=ready_pair())
        bench=B.rn32(F(1,10)); sb=B.mul(bench,11); fb=B.mul(bench,12)
        out=X.bridge(s,entry,fresh,separate_noise_sqrt_gain=11,fma_noise_sqrt_gain=12,
            scope=SCOPE.certified_scope(),commit_cfg=BASE.cfg(),bench_noise_sigma=bench,
            noise_sqrt=BAND.NoiseSqrtWitness(0),scheduler=scheduler,racc=RACC.State(),aw_sync=AWSYNC.State())
        self.assertIs(out.machine,s.machine); self.assertFalse(out.machine.pending)
        self.assertEqual(out.separate_active.Sigma_aw[2][2],B.mul(sb,sb))
        self.assertEqual(out.fma_active.Sigma_aw[2][2],B.mul(fb,fb))
        self.assertNotEqual(out.separate_active.Sigma_aw[2][2],sb*sb)
        self.assertNotEqual(out.separate_active,out.fma_active)
        with self.assertRaisesRegex(ValueError,'binary32 transaction'):
            replace(out,separate_commit=replace(out.separate_commit,Sigma_aw=active.Sigma_aw))

    def test_readiness_does_not_identify_machine_and_exact_active_parameters_generally(self):
        r=X.readiness()
        self.assertTrue(r['goLive_carries_whole_machine_TuneState_product'])
        self.assertTrue(r['goLive_preserves_tau_sigma_RS_and_common_pending_bit_by_identity'])
        self.assertFalse(r['machine_active_parameters_identified_with_exact_shadow_active_parameters'])
        self.assertFalse(r['machine_active_parameter_displacement_bound_closed'])
        self.assertFalse(r['every_admitted_startup_history_reaches_this_goLive_product'])
        self.assertFalse(r['Live_600_step_machine_TuneState_product_attached'])
        self.assertFalse(r['storage_search_allowed'])
        self.assertFalse(r['ALT_STARTUP_PASS']); self.assertFalse(r['ALT_LIVE_PASS'])


if __name__=='__main__': unittest.main()
