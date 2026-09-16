"""Exact all-input exp/log certificates and stale-coefficient rejection."""
from fractions import Fraction as F
import unittest
from unittest.mock import patch

from tools.stability.ou3_alt_contraction import target_qaxis_exp as Q
from tools.stability.ou3_alt_contraction import target_wpe_libm as W
from tools.stability.ou3_alt_contraction import finite_qaxis_exp_binary32 as E


class Tests(unittest.TestCase):
    def test_target_Qaxis_error_fits_original_real_envelope(self):
        proof=Q.error_certificate()
        self.assertLess(proof['total_absolute_exp_error'],F(1,2**24))
        combined=E.pinned_target_graph_envelope()
        self.assertTrue(combined['target_graph_results_satisfy_original_Qaxis_envelope_under_scalar_link_premises'])
        self.assertFalse(combined['shared_scalar_and_final_link_premises_discharged_here'])

    def test_altered_Qaxis_polynomial_cannot_reuse_certificate(self):
        with patch.object(Q,'P',(Q.P[0]+F(1,1000),)+Q.P[1:]):
            with self.assertRaises(AssertionError):
                Q.error_certificate()

    def test_target_WPE_bounds_cover_tolerant_runtime_domains(self):
        exp=W.exp_certificate(); log=W.log_certificate()
        self.assertEqual(exp['argument_interval'],(F(-60),F(60)))
        self.assertEqual(log['argument_interval'],(F(1,2**62),F(2**16)))
        self.assertLess(exp['total_relative_exp_error'],F(1,2**20))
        self.assertLess(log['total_absolute_log_error'],F(1,2**14))
        self.assertFalse(exp['transcendental_correct_rounding_assumed'])
        self.assertFalse(log['transcendental_correct_rounding_assumed'])

    def test_log_polynomial_error_is_derived_from_actual_coefficients(self):
        with patch.object(W,'LG',(W.LG[0]+F(1,10),)+W.LG[1:]):
            with self.assertRaises(AssertionError):
                W.log_certificate()

    def test_named_profile_closes_library_step_and_retains_firmware_obligation(self):
        profile=W.profile_correspondence()
        self.assertTrue(profile['pinned_WPE_exp_log_approximation_correspondence_closed'])
        self.assertTrue(profile['pinned_WPE_sqrt_approximation_correspondence_closed'])
        self.assertTrue(E.readiness()['Qaxis_pinned_libm_approximation_correspondence_closed'])
        self.assertFalse(profile['whole_firmware_compiler_and_link_correspondence_closed'])
        self.assertEqual(profile['rounding_premise'],
                         'FCR.RM=0 throughout the admitted execution')

    def test_shared_target_premises_are_not_claimed_by_algebra(self):
        proof=W.certificate()
        self.assertFalse(proof['scalar_instruction_and_division_premises_discharged_here'])
        self.assertFalse(proof['whole_firmware_link_resolution_qualified_here'])
        self.assertFalse(proof['WPE_sqrt_error_qualification_supplied_here'])


if __name__=='__main__':
    unittest.main()
