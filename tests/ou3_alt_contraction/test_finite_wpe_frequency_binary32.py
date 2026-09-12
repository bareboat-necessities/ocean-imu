"""WPE binary32 getter/tuner-store deployment regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_frequency_binary32 as STORE
from tools.stability.ou3_alt_contraction import finite_wpe_frequency_binary32 as X
from tools.stability.ou3_alt_contraction import finite_wpe_runtime as WPE


class Tests(unittest.TestCase):
    def shadow(self):
        return WPE.WPEState(log_period=F(7,10),usable_period=True)

    def test_two_shipping_exp_calls_share_one_stored_log_but_not_bit_reciprocity(self):
        log=X.bind_log_state(self.shadow(),B.rn32(F(7,10)))
        period=B.rn32(F(201,100))
        frequency=B.rn32(F(497,1000))
        out=X.getters(log,period_exp=period,frequency_exp=frequency)
        self.assertEqual(out.period_argument,log.stored_log_period)
        self.assertEqual(out.frequency_argument,-log.stored_log_period)
        self.assertEqual(out.period_result,period)
        self.assertEqual(out.frequency_result,frequency)
        # Deployment witnesses are independent libm results; forcing an exact
        # reciprocal product here would be a false binary32 theorem.
        self.assertNotEqual(out.period_result*out.frequency_result,1)

    def test_frequency_getter_result_is_exact_input_to_tuner_store(self):
        log=X.bind_log_state(self.shadow(),B.rn32(F(7,10)))
        out=X.getters(log,period_exp=B.rn32(2),frequency_exp=B.rn32(F(1,2)))
        stored=X.store_frequency(out,B.rn32(F(3,100)),B.rn32(F(6,5)))
        self.assertIsInstance(stored,STORE.StoredFrequency)
        self.assertEqual(stored.input_hz,out.frequency_result)
        self.assertEqual(stored.stored_hz,out.frequency_result)

    def test_log_state_and_getter_arguments_cannot_be_spliced(self):
        shadow=self.shadow()
        with self.assertRaisesRegex(ValueError,'deployment residual detached'):
            X.StoredLogPeriod(shadow.log_period,B.rn32(F(7,10)),0)
        log=X.bind_log_state(shadow,B.rn32(F(7,10)))
        with self.assertRaisesRegex(ValueError,'arguments detached'):
            X.GetterResult(log,log.stored_log_period,
                           log.stored_log_period,
                           B.rn32(2),B.rn32(F(1,2)))

    def test_getter_requires_binary32_results_without_claiming_libm_correctness(self):
        log=X.bind_log_state(self.shadow(),B.rn32(F(7,10)))
        with self.assertRaisesRegex(ValueError,'binary32 witnesses'):
            X.getters(log,period_exp=F(2),frequency_exp=F(1,3))
        r=X.readiness()
        self.assertTrue(r['shipping_WPE_dual_exp_getter_source_shape_matches'])
        self.assertTrue(r['period_and_frequency_getters_share_same_stored_binary32_log_state'])
        self.assertTrue(r['WPE_frequency_getter_to_tuner_binary32_store_topology_closed'])
        self.assertFalse(r['binary32_period_frequency_bit_reciprocity_assumed'])
        self.assertFalse(r['WPE_binary32_log_period_production_closed'])
        self.assertFalse(r['WPE_frequency_exp_target_libm_correspondence_closed'])
        self.assertFalse(r['storage_search_allowed'])
        self.assertFalse(r['ALT_LIVE_PASS'])


if __name__=='__main__': unittest.main()
