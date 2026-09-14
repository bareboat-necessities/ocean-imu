"""Universal primary-physics -> finite-word restriction theorem regressions."""
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_complete_brmm_restriction as X
from tools.stability.ou3_alt_contraction import finite_brmm_moment_prefix as MOMENT


class Tests(unittest.TestCase):
    def test_primary_physical_history_restricts_without_common_generator(self):
        d=X.build(); self.assertEqual(X.validate(d),[])
        self.assertTrue(d['primary_physics_is_definition'])
        self.assertFalse(d['common_generator_representation_required'])
        self.assertTrue(d['certificate_methods_sufficient_not_definitional'])
        self.assertTrue(d['same_history_restriction_required'])

    def test_finite_caps_are_exactly_the_primary_complete_brmm_caps(self):
        d=X.build(); self.assertTrue(d['finite_endpoint_caps_equal_primary_hard_bounds'])
        self.assertEqual(MOMENT.BOUNDS.acceleration,F(44,5))
        self.assertEqual(MOMENT.BOUNDS.velocity,F(11,2))
        self.assertEqual(MOMENT.BOUNDS.position,F(81,10))
        self.assertEqual(MOMENT.BOUNDS.centered_S,F(1100))
        self.assertEqual(d['finite_rate_cap_rad_s'],'11/18')

    def test_bounded_primitive_maps_to_one_live_origin_without_word_reset(self):
        d=X.build()
        self.assertTrue(d['centered_S_is_restriction_of_integral_from_Live_origin'])
        self.assertTrue(d['bounded_primitive_implies_every_prefix_S_cap'])
        self.assertTrue(d['wordwise_S_rezero_forbidden'])
        self.assertTrue(d['position_reanchor_forbidden'])

    def test_moment_iqc_is_a_consequence_of_same_acceleration_history(self):
        d=X.build()
        self.assertTrue(d['coupled_J0_J1_J2_are_moments_of_one_acceleration_history'])
        self.assertTrue(d['moment_projection_energy_le_pointwise_energy'])
        self.assertTrue(d['pointwise_acceleration_cap_implies_segment_and_prefix_moment_IQC'])

    def test_restriction_does_not_promote_unbound_runtime_parts(self):
        d=X.build()
        for k in ('arbitrary_runtime_tokens_prove_primary_history_membership',
                  'frequency_support_bound_consumed_by_current_finite_runtime_graph',
                  'full_BIAS_generating_history_attached_here',
                  'source_to_runtime_floating_arithmetic_closed_here',
                  'complete_finite_shipping_word_closed_here','storage_search_allowed',
                  'ALT_LIVE_PASS','ALT_STARTUP_PASS','ALT_END_TO_END_PASS'):
            self.assertFalse(d[k])

if __name__=='__main__': unittest.main()
