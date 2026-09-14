"""Source-locked shipping tau-target binary32 regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_shipping_tau_target_binary32 as X


class Tests(unittest.TestCase):
    def test_prior_and_frequency_extremes_share_one_source_locked_target_cell(self):
        for raw in (F(1,1000),F(1,5),F(1,2),F(6,5),F(5)):
            q=X.evaluate(B.rn32(raw))
            self.assertGreaterEqual(q.clamped_frequency,X.FLOOR)
            self.assertLessEqual(q.clamped_frequency,X.CEIL)
            self.assertGreaterEqual(q.binary32_target,X.TAU_MIN)
            self.assertLessEqual(q.binary32_target,X.TAU_MAX)
            self.assertLessEqual(abs(q.error),X.TARGET_ERROR_MAX)

    def test_upper_tau_clamp_makes_large_raw_division_roundoff_irrelevant(self):
        q=X.evaluate(B.rn32(F(1,1000)))
        self.assertEqual(q.clamped_frequency,X.FLOOR)
        self.assertEqual(q.exact_target,X.TAU_MAX)
        self.assertEqual(q.binary32_target,X.TAU_MAX)
        self.assertEqual(q.error,0)

    def test_nonbinary_or_nonpositive_source_frequency_is_rejected(self):
        with self.assertRaisesRegex(ValueError,'stored binary32'):
            X.evaluate(F(1,10))
        with self.assertRaisesRegex(ValueError,'positive stored'):
            X.evaluate(F(0))

    def test_readiness_closes_target_cell_not_upstream_or_predecessor_domain(self):
        r=X.readiness()
        self.assertTrue(r['shipping_tau_frequency_and_target_constant_source_shape_matches'])
        self.assertTrue(r['compiled_binary32_constants_used_as_exact_shadow_operands'])
        self.assertTrue(r['shipping_tau_target_exact_vs_binary32_cell_bound_source_locked'])
        self.assertFalse(r['upstream_WPE_binary32_frequency_production_closed'])
        self.assertFalse(r['shipping_tau_predecessor_domain_inductively_closed'])
        self.assertFalse(r['source_uniform_tau_roundoff_supply_bound_closed'])
        self.assertFalse(r['storage_search_allowed'])


if __name__=='__main__': unittest.main()
