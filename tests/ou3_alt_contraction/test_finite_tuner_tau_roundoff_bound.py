"""Conditional source-uniform tau EMA roundoff-bound regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as C
from tools.stability.ou3_alt_contraction import finite_tuner_tau_binary32 as T
from tools.stability.ou3_alt_contraction import finite_tuner_tau_roundoff_bound as X
import test_finite_tuner_tau_binary32 as BASE


class Tests(unittest.TestCase):
    def edge(self,previous=F(1,2)):
        e=B.rn32(F(9753,10000)); p=B.rn32(previous)
        binary=T.step(p,BASE.sample(),BASE.cfg(),dt=F(1,100),exp_decay=e)
        exact=BASE.exact_candidate(p,e)
        return binary,exact

    def test_both_compiler_shapes_fit_one_history_independent_bound(self):
        binary,exact=self.edge()
        c=X.certify(binary,exact,exact_target=exact.tau_target)
        self.assertLessEqual(abs(c.supply.residual_fma),X.FMA_DERIVED_MAX)
        self.assertLessEqual(abs(c.supply.residual_separate),X.SEPARATE_DERIVED_MAX)
        self.assertLessEqual(abs(c.supply.residual_fma),X.UNIFORM_RESIDUAL_MAX)
        self.assertLessEqual(abs(c.supply.residual_separate),X.UNIFORM_RESIDUAL_MAX)
        self.assertEqual(X.UNIFORM_RESIDUAL_MAX,F(1,1<<14))

    def test_target_cell_and_predecessor_domain_are_fail_closed(self):
        binary,exact=self.edge()
        with self.assertRaisesRegex(ValueError,'target argument detached'):
            X.certify(binary,exact,exact_target=exact.tau_target+F(1,100))
        # A fabricated TauStep outside the predecessor theorem domain is rejected
        # even when all fields remain representable binary32 values.
        bad=T.TauStep(B.rn32(16),binary.frequency,binary.tau_target,binary.sea_time,
                      binary.adapt_sec,binary.exp_decay,binary.alpha,
                      binary.next_separate,binary.next_fma)
        with self.assertRaisesRegex(ValueError,'previous tau outside'):
            X.certify(bad,exact,exact_target=exact.tau_target)

    def test_readiness_is_conditional_not_complete_source_uniform_promotion(self):
        r=X.readiness()
        self.assertTrue(r['tau_FMA_roundoff_supply_uniform_bound_conditional_on_scalar_domain'])
        self.assertTrue(r['tau_separate_mul_add_roundoff_supply_uniform_bound_conditional_on_scalar_domain'])
        self.assertFalse(r['shipping_tau_target_exact_vs_binary32_cell_bound_source_locked'])
        self.assertFalse(r['shipping_tau_predecessor_domain_inductively_closed'])
        self.assertFalse(r['source_uniform_tau_roundoff_supply_bound_closed'])
        self.assertFalse(r['storage_search_allowed'])
        self.assertFalse(r['ALT_LIVE_PASS'])


if __name__=='__main__': unittest.main()
