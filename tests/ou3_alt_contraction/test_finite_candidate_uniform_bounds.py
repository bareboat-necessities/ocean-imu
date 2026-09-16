import unittest

from tools.stability.ou3_alt_contraction import finite_candidate_uniform_bounds as C


class Tests(unittest.TestCase):
    def test_configured_candidate_has_finite_source_uniform_ranges(self):
        report = C.build()
        self.assertTrue(report['candidate_arithmetic_totality_under_named_pow_range_closed'])
        self.assertTrue(report['pinned_powf_positive_finite_range_under_scalar_profile_closed'])
        self.assertTrue(report['compiled_pow_exponents_identified_with_exact_rational_exponents'])
        self.assertFalse(report['full_target_candidate_compiler_and_pow_profile_correspondence_closed'])
        self.assertFalse(report['startup_capture_and_MEKF_covariance_totality_proved_here'])
        self.assertLess(report['candidate_inductive_images']['tau'][1], 13)

    def test_candidate_state_and_inputs_are_not_detached(self):
        with self.assertRaises(TypeError):
            C.require_config(object())
        with self.assertRaises(ValueError):
            C.require_pow(C.POW_BASE[0], C.Q.EXPONENT, C.POW_OUTPUT[0])


if __name__ == '__main__':
    unittest.main()
