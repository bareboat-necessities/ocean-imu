"""Deployment-faithful SpectralMSE configuration regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_deployment_config as D
from tools.stability.ou3_alt_contraction import finite_tuner_qeff_cache_binary32 as Q


class Tests(unittest.TestCase):
    def test_shipping_scalar_defaults_are_binary32_and_cache_is_source_bound(self):
        c=D.shipping_defaults(qeff_pow_result=B.rn32(1))
        for name in ('min_freq','max_freq','tau_coeff','sigma_coeff','min_tau','max_tau','max_sigma',
                     'pseudo_tau_ratio','pseudo_min','pseudo_max','min_RS','max_RS','rs_mse_coeff',
                     'accel_noise_density','adapt_tau_sec','adapt_tau_sea_periods','adapt_RS_mult',
                     'adapt_RS_slew_log','adapt_every_sec'):
            self.assertTrue(B.is_binary32(getattr(c,name)),name)
        self.assertEqual(c.accel_noise_density,Q.default_r_a_binary32())
        self.assertEqual(c.qeff_pow,c.qeff_cache.result)
        self.assertEqual(c.qeff_exact_root_interval,(c.qeff_cache.true_lo,c.qeff_cache.true_hi))

    def test_machine_cache_value_is_distinct_from_exact_root_relation(self):
        a=D.shipping_defaults(qeff_pow_result=B.rn32(1))
        b=D.shipping_defaults(qeff_pow_result=B.rn32(2))
        self.assertNotEqual(a.qeff_pow,b.qeff_pow)
        # Both caches are rooted at the same exact shipping 2*r_a argument, so
        # changing only the unqualified libm witness must not change the exact root.
        self.assertEqual(a.qeff_exact_root_interval,b.qeff_exact_root_interval)
        self.assertNotEqual(a.qeff_cache.error_lo,b.qeff_cache.error_lo)

    def test_detached_noise_density_cache_is_rejected(self):
        cache=Q.produce(r_a=B.rn32(F(1,1000)),pow_result=B.rn32(1))
        with self.assertRaisesRegex(ValueError,'r_a detached'):
            D.DeploymentConfig(D.MIN_FREQ,D.MAX_FREQ,D.TAU_COEFF,D.SIGMA_COEFF,
                D.MIN_TAU,D.MAX_TAU,D.MAX_SIGMA,D.PSEUDO_RATIO,D.PSEUDO_MIN,D.PSEUDO_MAX,
                D.MIN_RS,D.MAX_RS,D.RS_MSE_COEFF,Q.default_r_a_binary32(),cache,
                D.ADAPT_TAU_SEC,D.ADAPT_TAU_SEA_PERIODS,D.ADAPT_RS_MULT,
                D.ADAPT_RS_SLEW_LOG,D.ADAPT_EVERY_SEC,True)

    def test_readiness_keeps_all_libm_and_theorem_gates_false(self):
        r=D.readiness()
        self.assertTrue(r['shipping_default_source_shape_matches'])
        self.assertTrue(r['legacy_exact_rational_qeff_identity_not_required'])
        self.assertTrue(r['configured_r_a_and_cached_qeff_machine_value_carried_distinctly'])
        self.assertTrue(r['exact_qeff_fourteenth_root_interval_carried_separately_from_cache'])
        for k in ('qeff_cache_target_libm_correspondence_closed',
                  'per_sample_sqrt_pow_target_libm_correspondence_closed',
                  'source_uniform_complete_startup_reachability_closed',
                  'storage_search_allowed','ALT_STARTUP_PASS','ALT_LIVE_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
