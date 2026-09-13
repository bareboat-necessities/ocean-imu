"""Whole-machine TuneState goLive bridge regressions."""
from fractions import Fraction as F
import unittest

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
        self.assertEqual(out.wpe,s.lower.wpe)
        self.assertTrue(out.live.frontend_live.tuner.pending)
        self.assertTrue(out.machine.pending)
        self.assertEqual(out.separate_active,active)
        self.assertEqual(out.fma_active,active)

    def test_goLive_applies_even_when_online_pending_is_false_and_preserves_false(self):
        s,entry,fresh,active,scheduler=startup(False)
        out=X.bridge(s,entry,fresh,scope=SCOPE.certified_scope(),commit_cfg=BASE.cfg(),
            bench_noise_sigma=0,noise_sqrt=BAND.NoiseSqrtWitness(0),
            scheduler=scheduler,racc=RACC.State(),aw_sync=AWSYNC.State())
        self.assertFalse(out.machine.pending)
        self.assertFalse(out.live.frontend_live.tuner.pending)
        self.assertIsNotNone(out.separate_commit); self.assertIsNotNone(out.fma_commit)

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
