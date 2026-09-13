"""Cached q_eff^(1/14) binary32 production regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_qeff_cache_binary32 as X


class Tests(unittest.TestCase):
    def test_default_r_a_follows_literal_source_order_binary32_products(self):
        noise=B.rn32(F(148,10000)); dt=B.div(B.rn32(1),B.rn32(200))
        expected=B.mul(B.mul(noise,noise),dt)
        self.assertEqual(X.default_r_a_binary32(),expected)
        self.assertTrue(B.is_binary32(expected)); self.assertGreater(expected,0)

    def test_cache_witness_retains_same_argument_exponent_and_exact_root_supply(self):
        ra=X.default_r_a_binary32()
        # Explicit witness only; no target-libm correctness is asserted here.
        out=X.produce(r_a=ra,pow_result=B.rn32(F(2,5)))
        self.assertEqual(out.qeff_argument,B.mul(B.rn32(2),ra))
        self.assertEqual(out.exponent,B.div(B.rn32(1),B.rn32(14)))
        self.assertLessEqual(out.true_lo**14,out.qeff_argument)
        self.assertGreaterEqual(out.true_hi**14,out.qeff_argument)
        self.assertEqual(out.error_lo,out.result-out.true_hi)
        self.assertEqual(out.error_hi,out.result-out.true_lo)

    def test_candidate_cache_splice_is_rejected(self):
        out=X.produce(r_a=X.default_r_a_binary32(),pow_result=B.rn32(F(2,5)))
        self.assertIs(X.qualify_candidate_cache(out.result,out),out)
        with self.assertRaisesRegex(ValueError,'detached'):
            X.qualify_candidate_cache(B.rn32(F(1,2)),out)

    def test_readiness_keeps_libm_and_storage_fail_closed(self):
        r=X.readiness()
        for k in ('shipping_qeff_cache_source_shape_matches',
                  'default_r_a_source_order_binary32_graph_materialized',
                  'qeff_argument_binary32_2_times_r_a_materialized',
                  'compiled_binary32_one_over_fourteen_exponent_materialized',
                  'cached_pow_result_bound_to_same_machine_argument',
                  'cached_pow_error_interval_against_exact_fourteenth_root_exposed',
                  'candidate_qeff_cache_can_be_qualified_against_produced_cache'):
            self.assertTrue(r[k])
        self.assertFalse(r['target_libm_qeff_pow_correspondence_closed'])
        self.assertFalse(r['source_uniform_qeff_cache_supply_bound_closed'])
        self.assertFalse(r['binary32_SpectralMSE_target_correspondence_closed'])
        self.assertFalse(r['storage_search_allowed'])
        self.assertFalse(r['ALT_STARTUP_PASS'])
        self.assertFalse(r['ALT_LIVE_PASS'])


if __name__=='__main__': unittest.main()
