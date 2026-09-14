"""Binary32 R_S EMA update-order regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_rs_ema_binary32 as X


class Tests(unittest.TestCase):
    def test_source_order_subtract_multiply_add_is_retained(self):
        p=B.rn32(F(1,2)); t=B.rn32(F(3,4)); a=B.rn32(F(1,100))
        out=X.step(p,t,a)
        self.assertEqual(out.delta,B.sub(t,p))
        self.assertEqual(out.increment,B.mul(a,out.delta))
        self.assertEqual(out.next,B.add(p,out.increment))

    def test_nonbinary_or_invalid_alpha_fails_closed(self):
        with self.assertRaisesRegex(ValueError,'previous RS must be actual binary32'):
            X.step(F(1,3),B.rn32(1),B.rn32(F(1,2)))
        with self.assertRaisesRegex(ValueError,'alpha_RS'):
            X.step(B.rn32(1),B.rn32(1),B.rn32(2))

    def test_readiness_keeps_alpha_libm_and_commit_open(self):
        r=X.readiness()
        self.assertTrue(r['shipping_RS_EMA_source_shape_matches'])
        self.assertTrue(r['RS_target_minus_previous_binary32_subtraction_materialized'])
        self.assertTrue(r['alpha_times_delta_binary32_multiply_materialized'])
        self.assertTrue(r['previous_plus_increment_binary32_add_materialized'])
        for k in ('alpha_RS_target_libm_production_closed','source_uniform_alpha_RS_supply_bound_closed',
                  'machine_RS_EMA_to_exact_interval_join_closed','RS_commit_correspondence_closed',
                  'source_uniform_complete_600_step_word_qualified','storage_search_allowed','ALT_LIVE_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
