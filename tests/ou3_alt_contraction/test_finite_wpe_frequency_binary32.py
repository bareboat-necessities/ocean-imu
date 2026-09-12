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

    def getter(self):
        log=X.bind_log_state(self.shadow(),B.rn32(F(7,10)))
        return X.getters(log,period_exp=B.rn32(F(201,100)),frequency_exp=B.rn32(F(497,1000)))

    def test_two_shipping_exp_calls_share_one_stored_log_but_not_bit_reciprocity(self):
        out=self.getter()
        self.assertEqual(out.period_argument,out.log.stored_log_period)
        self.assertEqual(out.frequency_argument,-out.log.stored_log_period)
        self.assertNotEqual(out.period_result*out.frequency_result,1)

    def test_usable_preupdate_WPE_getter_is_exact_input_to_tuner_store(self):
        shadow=self.shadow(); out=self.getter()
        q=X.tuner_frequency(shadow,min_hz=B.rn32(F(3,100)),max_hz=B.rn32(F(6,5)),getter=out)
        self.assertEqual(q.branch,'wpe')
        self.assertIsInstance(q.stored,STORE.StoredFrequency)
        self.assertEqual(q.stored.input_hz,out.frequency_result)
        self.assertEqual(q.stored.stored_hz,out.frequency_result)

    def test_preusable_WPE_uses_literal_prior_and_consumes_no_getter(self):
        shadow=WPE.WPEState()
        q=X.tuner_frequency(shadow,min_hz=B.rn32(F(3,100)),max_hz=B.rn32(F(6,5)))
        self.assertEqual(q.branch,'prior')
        self.assertIsNone(q.getter)
        self.assertEqual(q.stored.input_hz,X.PRIOR)
        self.assertEqual(q.stored.stored_hz,X.PRIOR)
        with self.assertRaisesRegex(ValueError,'consumes no getter'):
            X.tuner_frequency(shadow,min_hz=B.rn32(F(3,100)),max_hz=B.rn32(F(6,5)),getter=self.getter())

    def test_log_state_and_getter_arguments_cannot_be_spliced(self):
        shadow=self.shadow()
        with self.assertRaisesRegex(ValueError,'deployment residual detached'):
            X.StoredLogPeriod(shadow.log_period,B.rn32(F(7,10)),0)
        log=X.bind_log_state(shadow,B.rn32(F(7,10)))
        with self.assertRaisesRegex(ValueError,'arguments detached'):
            X.GetterResult(log,log.stored_log_period,log.stored_log_period,
                           B.rn32(2),B.rn32(F(1,2)))
        other=WPE.WPEState(log_period=F(4,5),usable_period=True)
        with self.assertRaisesRegex(ValueError,'preupdate canonical'):
            X.tuner_frequency(other,min_hz=B.rn32(F(3,100)),max_hz=B.rn32(F(6,5)),getter=self.getter())

    def test_getter_requires_binary32_results_without_claiming_libm_correctness(self):
        log=X.bind_log_state(self.shadow(),B.rn32(F(7,10)))
        with self.assertRaisesRegex(ValueError,'binary32 witnesses'):
            X.getters(log,period_exp=F(2),frequency_exp=F(1,3))
        r=X.readiness()
        self.assertTrue(r['shipping_WPE_dual_exp_and_tuner_source_shape_matches'])
        self.assertTrue(r['period_and_frequency_getters_share_same_stored_binary32_log_state'])
        self.assertTrue(r['preusable_WPE_uses_literal_binary32_0p2_prior_without_exp'])
        self.assertTrue(r['usable_WPE_frequency_getter_bound_to_sample_entry_log_state'])
        self.assertTrue(r['WPE_frequency_getter_to_tuner_binary32_store_topology_closed'])
        self.assertFalse(r['binary32_period_frequency_bit_reciprocity_assumed'])
        self.assertFalse(r['WPE_binary32_log_period_production_closed'])
        self.assertFalse(r['WPE_frequency_exp_target_libm_correspondence_closed'])
        self.assertFalse(r['storage_search_allowed'])
        self.assertFalse(r['ALT_LIVE_PASS'])


if __name__=='__main__': unittest.main()
