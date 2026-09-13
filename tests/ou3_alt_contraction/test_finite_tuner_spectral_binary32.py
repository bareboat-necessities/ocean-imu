"""Binary32 SpectralMSE graph and libm-supply regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_spectral_binary32 as X
import test_finite_complete_word_tau_qualification as QBASE


class Tests(unittest.TestCase):
    def cfg(self): return QBASE.shipping_runtime().candidate_cfg

    def test_shipping_graph_materializes_ordinary_float_ops_around_libm(self):
        c=self.cfg()
        # The values below are explicit binary32 libm WITNESSES only.  This test
        # validates graph ancestry and supply accounting; it does not claim they
        # are target-platform powf/sqrtf results.
        out=X.step_from_config(c,tau=B.rn32(F(5,2)),sigma=B.rn32(F(9,10)),
            pow_result=B.rn32(F(16)),sqrt_result=B.rn32(F(1,5)))
        self.assertEqual(out.pow_witness.argument,out.u)
        self.assertEqual(out.sqrt_witness.argument,out.TS)
        self.assertEqual(out.pow_witness.exponent,X.POW_EXPONENT)
        self.assertEqual(out.pow_supply.error_lo,out.pow_witness.result-out.pow_supply.true_hi)
        self.assertEqual(out.pow_supply.error_hi,out.pow_witness.result-out.pow_supply.true_lo)
        self.assertEqual(out.sqrt_supply.error_lo,out.sqrt_witness.result-out.sqrt_supply.true_hi)
        self.assertEqual(out.sqrt_supply.error_hi,out.sqrt_witness.result-out.sqrt_supply.true_lo)
        self.assertTrue(B.is_binary32(out.raw_RS))

    def test_source_order_for_u_is_left_associative_binary32(self):
        c=self.cfg(); tau=B.rn32(F(7,3)); sigma=B.rn32(F(13,10))
        out=X.step_from_config(c,tau=tau,sigma=sigma,
            pow_result=B.rn32(F(8)),sqrt_result=B.rn32(F(1,4)))
        expected_tau2=B.mul(tau,tau)
        expected_div=B.div(sigma,c.sigma_coeff)
        expected_sab=max(expected_div,X.SIGMA_AB_MIN)
        expected_u=B.mul(B.mul(expected_sab,expected_tau2),expected_tau2)
        self.assertEqual(out.tau2,expected_tau2)
        self.assertEqual(out.u,expected_u)

    def test_pow_exponent_is_compiled_binary32_six_over_seven(self):
        self.assertEqual(X.POW_EXPONENT,B.div(B.rn32(6),B.rn32(7)))
        with self.assertRaisesRegex(ValueError,'pow exponent detached'):
            X.LibmWitness('pow',B.rn32(1),B.rn32(1),B.rn32(F(3,4)))
        with self.assertRaisesRegex(ValueError,'sqrt witness consumes no exponent'):
            X.LibmWitness('sqrt',B.rn32(1),B.rn32(1),X.POW_EXPONENT)

    def test_nonbinary32_or_nonpositive_machine_operands_fail_closed(self):
        c=self.cfg()
        with self.assertRaisesRegex(ValueError,'tau must be an actual binary32'):
            X.step_from_config(c,tau=F(1,3),sigma=B.rn32(1),pow_result=B.rn32(1),sqrt_result=B.rn32(1))
        with self.assertRaisesRegex(ValueError,'libm result must be positive'):
            X.step_from_config(c,tau=B.rn32(1),sigma=B.rn32(1),pow_result=0,sqrt_result=B.rn32(1))

    def test_readiness_closes_graph_not_target_libm_or_master_word(self):
        r=X.readiness()
        for k in ('shipping_SpectralMSE_source_shape_matches',
                  'pseudo_cadence_binary32_mul_and_clamp_materialized',
                  'sigma_div_floor_tau2_u_binary32_graph_materialized',
                  'pow_and_sqrt_kept_as_distinct_same_argument_binary32_witnesses',
                  'final_coefficient_multiply_multiply_divide_graph_materialized',
                  'pow_witness_error_interval_against_exact_root_of_same_machine_u_exposed',
                  'sqrt_witness_error_interval_against_exact_root_of_same_machine_TS_exposed'):
            self.assertTrue(r[k])
        for k in ('target_libm_pow_correspondence_closed','target_libm_sqrt_correspondence_closed',
                  'cached_qeff_pow_target_libm_correspondence_closed','source_uniform_root_supply_bounds_closed',
                  'binary32_SpectralMSE_target_correspondence_closed',
                  'source_uniform_complete_startup_reachability_closed','storage_search_allowed',
                  'ALT_STARTUP_PASS','ALT_LIVE_PASS'):
            self.assertFalse(r[k])


if __name__=='__main__': unittest.main()
