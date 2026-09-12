"""Universal restriction of primary COMPLETE-BRMM physics to the finite ALT word.

This is the theorem-level bridge that the finite source composer was missing.
COMPLETE-BRMM is defined by a physical same-history condition, not by one common
spectral or shaping-state generator.  Therefore universal finite-word ancestry
must follow by restricting an arbitrary admitted physical history to the Live
sampling grid, rather than by inventing a uniform generator certificate.

The result below proves that every admitted primary physical history restricted
to any finite 5 ms prefix satisfies the physical constraints already consumed by
``finite_source_continuation``: endpoint p/v/a/S caps, exact p/v/S primitive
recurrences, the coupled acceleration-moment IQC and the projective rotation
chord bound.  It does NOT prove that an arbitrary runtime object belongs to the
primary history merely because it carries matching string tokens.  It also does
not close the frequency/frontend relation, BIAS generating histories or any
floating-point/runtime arithmetic branch.
"""
from __future__ import annotations

from fractions import Fraction as F

import ou3_brmm_physical_wave_condition as PHYSICAL
import ou3_brmm_finite_window_primitive_qualification as PRIMITIVE
from tools.stability.ou3_alt_contraction import finite_brmm_moment_prefix as MOMENT

QUALIFICATION = 'OU3_ALT_COMPLETE_BRMM_FINITE_RESTRICTION_V1'
DT = F(1, 200)
TRANSITIONS = 600
HORIZON = TRANSITIONS * DT


def _exact_decimal(x):
    """Exact rational interpretation of declared decimal theorem constants."""
    return F(str(x))


def build():
    physical = PHYSICAL.build(); pf = PHYSICAL.validate(physical)
    primitive = PRIMITIVE.build(); qf = PRIMITIVE.validate(primitive)
    if pf or qf:
        raise RuntimeError('primary COMPLETE-BRMM prerequisite failed: '+repr({'physical':pf,'primitive':qf}))

    condition = physical_condition = physical
    hard = condition['complete_BRMM_numeric_physical_envelope']['complete_BRMM_hard_bounds']
    uniform = primitive['uniform_physical_primitives']
    relation = primitive['exact_source_relation']

    expected = {
        'position': _exact_decimal(hard['wave_position_norm_upper_m']),
        'velocity': _exact_decimal(hard['wave_velocity_norm_upper_mps']),
        'acceleration': _exact_decimal(hard['wave_acceleration_norm_upper_mps2']),
        'centered_S': _exact_decimal(hard['centered_primitive_D_S_upper_m_s']),
    }
    actual = {
        'position': MOMENT.BOUNDS.position,
        'velocity': MOMENT.BOUNDS.velocity,
        'acceleration': MOMENT.BOUNDS.acceleration,
        'centered_S': MOMENT.BOUNDS.centered_S,
    }
    if expected != actual:
        raise RuntimeError('finite physical caps detached from primary COMPLETE-BRMM envelope: '
                           +repr({'expected':expected,'actual':actual}))

    if relation['translation'] != 'p_dot=v; v_dot=a; S_dot=p in one fixed wave/world frame':
        raise RuntimeError('primary translation primitive relation changed')
    if relation['one_common_history_required'] is not True or relation['S_is_derived_not_independent'] is not True:
        raise RuntimeError('primary same-history/S ancestry changed')
    if condition['same_history_ancestry_required'] is not True:
        raise RuntimeError('primary physical condition lost same-history requirement')
    if condition['wordwise_rezero_of_S_allowed'] is not False or condition['position_reanchor_allowed'] is not False:
        raise RuntimeError('primary physical condition unexpectedly permits word restart')
    if condition['certificate_methods_are_sufficient_not_definitional'] is not True:
        raise RuntimeError('certificate representation became part of source definition')
    if condition['physical_D_S_numeric_qualification_closed_for_complete_family'] is not True:
        raise RuntimeError('primary numeric COMPLETE-BRMM envelope is not closed')

    # pi < 22/7 gives 35*pi/180 < 11/18 exactly, matching the finite chord cap.
    omega_outward = F(35, 180) * F(22, 7)
    if omega_outward != MOMENT.BOUNDS.angular_rate_upper:
        raise RuntimeError('finite rotation-rate cap detached from primary 35 deg/s bound')

    return {
        'qualification': QUALIFICATION,
        'canonical_source': condition['canonical_source'],
        'primary_physics_is_definition': condition['physics_is_primary'],
        'common_generator_representation_required': False,
        'certificate_methods_sufficient_not_definitional': True,
        'restriction_horizon_s': str(HORIZON),
        'restriction_sample_period_s': str(DT),
        'restriction_transition_count': TRANSITIONS,
        'same_history_restriction_required': True,
        'finite_endpoint_caps_equal_primary_hard_bounds': True,
        'centered_S_is_restriction_of_integral_from_Live_origin': True,
        'bounded_primitive_implies_every_prefix_S_cap': True,
        'wordwise_S_rezero_forbidden': True,
        'position_reanchor_forbidden': True,
        'exact_translation_recurrence_follows_from_primary_jet': True,
        'coupled_J0_J1_J2_are_moments_of_one_acceleration_history': True,
        'moment_projection_energy_le_pointwise_energy': True,
        'pointwise_acceleration_cap_implies_segment_and_prefix_moment_IQC': True,
        'body_rate_cap_implies_projective_rotation_chord_cap': True,
        'finite_rate_cap_rad_s': str(omega_outward),
        'arbitrary_runtime_tokens_prove_primary_history_membership': False,
        'frequency_support_bound_consumed_by_current_finite_runtime_graph': False,
        'full_BIAS_generating_history_attached_here': False,
        'source_to_runtime_floating_arithmetic_closed_here': False,
        'complete_finite_shipping_word_closed_here': False,
        'storage_search_allowed': False,
        'ALT_LIVE_PASS': False,
        'ALT_STARTUP_PASS': False,
        'ALT_END_TO_END_PASS': False,
    }


def validate(d):
    f=[]
    if d.get('qualification') != QUALIFICATION: f.append('qualification mismatch')
    for k in (
        'primary_physics_is_definition','certificate_methods_sufficient_not_definitional',
        'same_history_restriction_required','finite_endpoint_caps_equal_primary_hard_bounds',
        'centered_S_is_restriction_of_integral_from_Live_origin',
        'bounded_primitive_implies_every_prefix_S_cap','wordwise_S_rezero_forbidden',
        'position_reanchor_forbidden','exact_translation_recurrence_follows_from_primary_jet',
        'coupled_J0_J1_J2_are_moments_of_one_acceleration_history',
        'moment_projection_energy_le_pointwise_energy',
        'pointwise_acceleration_cap_implies_segment_and_prefix_moment_IQC',
        'body_rate_cap_implies_projective_rotation_chord_cap'):
        if d.get(k) is not True: f.append(k+' not true')
    for k in (
        'common_generator_representation_required','arbitrary_runtime_tokens_prove_primary_history_membership',
        'frequency_support_bound_consumed_by_current_finite_runtime_graph',
        'full_BIAS_generating_history_attached_here','source_to_runtime_floating_arithmetic_closed_here',
        'complete_finite_shipping_word_closed_here','storage_search_allowed',
        'ALT_LIVE_PASS','ALT_STARTUP_PASS','ALT_END_TO_END_PASS'):
        if d.get(k) is not False: f.append(k+' not false')
    if d.get('restriction_sample_period_s') != '1/200': f.append('sample period changed')
    if d.get('restriction_horizon_s') != '3': f.append('restriction horizon changed')
    if d.get('restriction_transition_count') != 600: f.append('transition count changed')
    if d.get('finite_rate_cap_rad_s') != '11/18': f.append('outward angular-rate cap changed')
    return f
