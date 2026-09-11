"""Primary BRMM physical-source declaration, not a numerical P3 certificate.

The deterministic theorem source is a bounded oscillatory physical wave history
about a local equilibrium.  Its centered primitive is uniformly bounded.  The
explicit padded numerical COMPLETE-BRMM envelope is owned by the primary
physical-wave condition; spectral/shaping realizations are certificate methods,
not the source definition.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import ou3_brmm_finite_window_primitive_qualification as PRIMITIVE
import ou3_brmm_gravity_direction_forcing_qualification as DIRECTION
import ou3_brmm_physical_wave_condition as PHYSICAL


def build():
    physical=PHYSICAL.build();wf=PHYSICAL.validate(physical)
    if wf:raise RuntimeError("physical wave condition failed: "+repr(wf))
    primitive=PRIMITIVE.build();pf=PRIMITIVE.validate(primitive)
    direction=DIRECTION.build();df=DIRECTION.validate(direction)
    if pf:raise RuntimeError('BRMM primitive qualification failed: '+repr(pf))
    if df:raise RuntimeError('BRMM gravity-direction qualification failed: '+repr(df))
    u=primitive['uniform_physical_primitives'];m=primitive['finite_window_moment_consequences']
    a=float(u['acceleration_norm_upper_mps2']);g=9.80665
    return {
        "qualification": "OU3_COMPLETE_BRMM_BOUNDED_BIAS_MOTION_V6",
        "canonical_source": "COMPLETE_BRMM_NORMAL_LIVE_WORD",
        "physics_is_primary":True,
        "spectral_membership_required": False,
        "physical_bounded_oscillatory_wave_condition_required": True,
        "physical_wave_condition_contract": physical,
        "certificate_methods_are_sufficient_not_definitional":True,
        "wave_displacement_DC_prohibited_by_physical_theorem": not physical["constant_nonzero_position_zero_velocity_history_admitted"],
        "constant_nonzero_position_zero_velocity_history_admitted": physical["constant_nonzero_position_zero_velocity_history_admitted"],
        "finite_window_only_definition_is_current_source": False,
        "common_source_frontend_tuner_geometry_and_filter_history_required": True,
        "every_due_actual_anisotropic_RS_required": True,
        "fixed_frame_CoG_acceleration_required": True,
        "hard_caps": {"acceleration_norm_mps2": a,
                      "body_rate_deg_s":u['body_rate_norm_upper_deg_s'],
                      "specific_force_norm_mps2": [g-a, g+a],
                      "positive_sea_height_m": {"lower_open": 0., "upper_closed": u['H_s_upper_m']},
                      "zero_height_quiet_case_admitted": True,
                      "estimated_bias_projection_radius_mps2": .4,
                      "physical_velocity_norm_mps":u['V_m_norm_upper_mps'],
                      "physical_wave_frame_position_norm_m":u['P_m_norm_upper_m'],
                      "centered_primitive_D_S_upper_m_s":u['S_centered_primitive_norm_upper_m_s'],
                      "frequency_support_hz":u['frequency_support_hz']},
        "legacy_active_bias_interior_mps2": .35,
        "primary_motion_target_covers_closed_projection_ball": True,
        "shipping_tuner_default_frequency_hz": [.03, 1.2],
        "proposed_outer_tuner_frequency_hz": [.018, 1.2],
        "tuner_frequency_is_physical_motion_bandlimit": False,
        "candidate_audit": {"recurrence_window_s": 10., "quiet_rms_mps2": [.03, .05],
                            "impulse_cap_mps": 2., "lobe_chi_for_diagnostic_only": 1.},
        "uniform_primitive_qualification":primitive,
        "gravity_direction_forcing_qualification":direction,
        "physical_constants": {"T_R": None, "E_q": None, "V_R": None,
                               "E_min": None, "J_min": None, "chi": None,
                               "V_a": None,
                               "V_m":u['V_m_norm_upper_mps'],
                               "P_m":u['P_m_norm_upper_m'],
                               "S_m":u['S_centered_primitive_norm_upper_m_s'],
                               "A_m":a},
        "S_primitive_contract":"no independent arbitrary S ball; S_L=integral_from_Live p_wave and the physical COMPLETE-BRMM family gives ||S_L||<=D_S on every continuation",
        "physical_numeric_family_qualified": physical["physical_D_S_numeric_qualification_closed_for_complete_family"],
        "S_three_second_increment_norm_upper_m_s":m['Delta_S_norm_upper_from_uniform_position_m_s'],
        "source_applicability_refined_to_uniform_primitives":True,
        "gravity_direction_low_frequency_qualification_required":True,
        "gravity_direction_same_history_decomposition":direction['direction_forcing_decomposition'],
        "gravity_direction_mean_chord_norm_upper":direction['mean_direction_chord_norm_upper'],
        "gravity_direction_primitive_norm_upper_s":direction['direction_primitive_norm_upper_s'],
        "pointwise_padded_acceleration_cap_retained":direction['pointwise_padded_acceleration_cap_retained'],
        "mean_impulse_bound_applies_to_Q_and_O": True,
        "Q_is_not_a_global_minimum_excitation_requirement": True,
        "Q_AC_energy_alone_excludes_constant_acceleration": False,
        "fixed_window_impulse_cap_alone_excludes_all_DC": False,
        "nonzero_fixed_physical_acceleration_DC_prohibited_in_Q_and_O": True,
        "bounded_velocity_primitive_required": True,
        "zero_DC_definition": "a_m=dv_m/dt and the physical COMPLETE-BRMM family has bounded v_m; hence a persistent nonzero acceleration DC is excluded",
        "constant_measurement_bias_is_marine_motion": False,
        "bounded_same_history_primitives_or_explicit_forcing_budget_required": True,
        "equivalent_hard_finite_window_dynamic_constraint_available":True,
        "Q_and_O_mode_mixes_and_vector_PE_require_separate_coverage": True,
        "bias_error_decay_required": False,
        "zero_navigation_error_floor_required": False,
        "P3_delta": 1e-18,
        "P3_matrix_implication_requires_its_actual_premises": True,
        "source_declaration_promotes_P3": False,
        "BRMM_PRIMITIVE_QUALIFICATION_PASS":True,
        "BRMM_DIRECTION_FORCING_QUALIFICATION_PASS":True,
        "BRMM_SOURCE_UNIFORM_PASS": False,
        "source_declaration_is_P3_certificate": False,
        "BRMM_P4_MOTION_PASS": False,
        "BRMM_P5_MOTION_MAY_START": False,
        "sampled_audit_sets_physical_theorem_constants": False,
    }


def validate(d):
    expected=build()
    return [f"{k} changed or falsely certified" for k,v in expected.items() if d.get(k)!=v]


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", required=True, type=Path)
    args = p.parse_args()
    args.output.write_text(json.dumps(build(), indent=2)+"\n")
