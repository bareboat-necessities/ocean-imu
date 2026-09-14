"""Local binary32 versus exact-real SpectralMSE residual regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as C
from tools.stability.ou3_alt_contraction import finite_tuner_deployment_config as D
from tools.stability.ou3_alt_contraction import finite_tuner_spectral_machine_real_join as X


def target():
    return C.TargetState(F(1,5),F(1),B.rn32(F(5,2)),B.rn32(F(9,10)))


class Tests(unittest.TestCase):
    def test_join_exposes_exact_machine_minus_real_interval(self):
        cfg=D.shipping_defaults(qeff_pow_result=B.rn32(1))
        out=X.join(cfg,target(),pow_result=B.rn32(16),sqrt_result=B.rn32(F(1,5)))
        self.assertEqual(out.machine.rs_qeff_pow,cfg.qeff_cache.result)
        self.assertEqual(out.machine.tau,out.target.tau_target)
        self.assertEqual(out.machine.sigma,out.target.sigma_target)
        self.assertEqual(out.residual_lo,out.machine_target_RS-out.exact.target_RS_hi)
        self.assertEqual(out.residual_hi,out.machine_target_RS-out.exact.target_RS_lo)
        self.assertLessEqual(out.residual_lo,out.residual_hi)

    def test_changing_only_machine_cache_moves_residual_not_exact_relation(self):
        a=X.join(D.shipping_defaults(qeff_pow_result=B.rn32(1)),target(),
                 pow_result=B.rn32(16),sqrt_result=B.rn32(F(1,5)))
        b=X.join(D.shipping_defaults(qeff_pow_result=B.rn32(2)),target(),
                 pow_result=B.rn32(16),sqrt_result=B.rn32(F(1,5)))
        self.assertEqual(a.exact.target_RS_lo,b.exact.target_RS_lo)
        self.assertEqual(a.exact.target_RS_hi,b.exact.target_RS_hi)
        self.assertNotEqual(a.machine.rs_qeff_pow,b.machine.rs_qeff_pow)
        self.assertNotEqual((a.residual_lo,a.residual_hi),(b.residual_lo,b.residual_hi))

    def test_nonbinary_candidate_inputs_fail_closed(self):
        cfg=D.shipping_defaults(qeff_pow_result=B.rn32(1))
        bad=C.TargetState(F(1,5),F(1),F(1,3),B.rn32(1))
        with self.assertRaisesRegex(ValueError,'tau target must be actual binary32'):
            X.join(cfg,bad,pow_result=B.rn32(1),sqrt_result=B.rn32(1))

    def test_readiness_is_local_supply_not_platform_or_commit_proof(self):
        r=X.readiness()
        for k in ('same_deployment_config_feeds_exact_real_and_binary32_targets',
                  'same_binary32_tau_sigma_candidate_feeds_both_relations',
                  'shipping_RS_clamp_composed_on_machine_target',
                  'machine_minus_exact_real_RS_residual_interval_exposed',
                  'residual_includes_qeff_cache_pow_and_per_candidate_pow_sqrt_supplies'):
            self.assertTrue(r[k])
        for k in ('target_libm_correspondence_closed','source_uniform_machine_real_RS_residual_bound_closed',
                  'RS_commit_correspondence_closed','source_uniform_complete_startup_reachability_closed',
                  'storage_search_allowed','ALT_STARTUP_PASS','ALT_LIVE_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
