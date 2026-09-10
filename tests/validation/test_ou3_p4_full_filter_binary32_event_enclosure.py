from __future__ import annotations

import pathlib
import sys
import unittest

ROOT=pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools'/'stability'))

import ou3_p4_full_filter_binary32_event_enclosure as FP


class FullFilterBinary32EventEnclosureTests(unittest.TestCase):
    def test_explicit_scalar_arithmetic_closes_but_full_filter_stays_open(self):
        d=FP.build()
        self.assertEqual(FP.validate(d),[])
        self.assertTrue(d['source_audit_closed'])
        self.assertTrue(d['projection_finite_precision_enclosure_closed'])
        self.assertTrue(d['joseph_scalar_loop_local_rounding_enclosure_available'])
        self.assertTrue(d['left_reset_scalar_loop_local_rounding_enclosure_available'])
        self.assertFalse(d['eigen_3x3_LDLT_backward_error_enclosure_closed'])
        self.assertFalse(d['covariance_roundoff_recursively_attached_to_next_event_P'])
        self.assertFalse(d['full_Kalman_reset_finite_precision_enclosure_closed'])
        self.assertFalse(d['terminal_gamma_n_may_claim_full_shipping_binary32'])
        self.assertFalse(d['P4_PASS'])

    def test_rounding_helpers_are_monotone(self):
        self.assertLess(FP.rn_abs_error(1.0),FP.rn_abs_error(2.0))
        self.assertLess(FP.dot3_abs_error(1.0,1.0),FP.dot3_abs_error(2.0,2.0))
        self.assertLess(
            FP.state_correction_component_error(1.0,1.0,1.0),
            FP.state_correction_component_error(2.0,2.0,2.0))
        self.assertLess(
            FP.joseph_entry_local_error(1.0,1.0,1.0,1.0),
            FP.joseph_entry_local_error(2.0,2.0,2.0,2.0))
        self.assertLess(
            FP.reset_covariance_entry_local_error(1.0,0.1),
            FP.reset_covariance_entry_local_error(2.0,0.2))

    def test_cancellation_does_not_remove_absolute_floor(self):
        self.assertGreater(FP.rn_abs_error(0.0),0.0)
        self.assertGreater(FP.dot3_abs_error(0.0,0.0),0.0)

    def test_false_full_arithmetic_promotion_is_rejected(self):
        d=FP.build()
        d['full_Kalman_reset_finite_precision_enclosure_closed']=True
        d['terminal_gamma_n_may_claim_full_shipping_binary32']=True
        d['P4_PASS']=True
        f=FP.validate(d)
        self.assertIn('full_Kalman_reset_finite_precision_enclosure_closed not false',f)
        self.assertIn('terminal_gamma_n_may_claim_full_shipping_binary32 not false',f)
        self.assertIn('P4_PASS not false',f)


if __name__=='__main__': unittest.main()
