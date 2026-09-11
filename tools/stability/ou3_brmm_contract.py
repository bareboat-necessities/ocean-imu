"""Primary BRMM physical-source declaration, not a numerical P3 certificate.

The deterministic theorem source requires a bounded oscillatory physical wave
generator. Its potential derivative is p; the centered S bound is a derived
generator consequence, not an independent S_m ball. Finite-window p/v/a
constraints remain necessary outer constraints, not sufficient admission. The widened
8 m/s^2 acceleration cap is paired with an explicit same-history low-frequency
gravity-direction qualification so large oscillatory acceleration is admitted
without silently allowing arbitrary persistent false-gravity forcing.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import ou3_brmm_finite_window_primitive_qualification as PRIMITIVE
import ou3_brmm_gravity_direction_forcing_qualification as DIRECTION
import ou3_brmm_physical_wave_source as WAVE


def build():
    wave=WAVE.build();wf=WAVE.validate(wave)
    if wf:raise RuntimeError("physical wave source failed: "+repr(wf))
    primitive=PRIMITIVE.build();pf=PRIMITIVE.validate(primitive)
    direction=DIRECTION.build();df=DIRECTION.validate(direction)
    if pf:raise RuntimeError('BRMM primitive qualification failed: '+repr(pf))
    if df:raise RuntimeError('BRMM gravity-direction qualification failed: '+repr(df))
    u=primitive['uniform_physical_primitives'];m=primitive['finite_window_moment_consequences']
    a=float(u['acceleration_norm_upper_mps2']);g=9.80665
    return {
        "qualification": "OU3_COMPLETE_BRMM_BOUNDED_BIAS_MOTION_V5",
        "canonical_source": "COMPLETE_BRMM_NORMAL_LIVE_WORD",
        "spectral_membership_required": False,  # bounded shaping realization is also allowed
        "physical_wave_generator_required": True,
        "physical_wave_generator_contract": wave,
        "wave_displacement_DC_prohibited_by_generator_theorem": wave["old_witness_excluded_by_corrected_physical_theorem"],
        "constant_nonzero_position_zero_velocity_history_admitted": wave["constant_nonzero_position_zero_velocity_history_admitted"],
        "finite_window_only_definition_is_current_source": False,
        "common_source_frontend_tuner_geometry_and_filter_history_required": True,
        "every_due_actual_anisotropic_RS_required": True,
        "fixed_frame_CoG_acceleration_required": True,
        "hard_caps": {"acceleration_norm_mps2": a, "body_rate_deg_s": 30.,
                      "specific_force_norm_mps2": [g-a, g+a],
                      "positive_sea_height_m": {"lower_open": 0., "upper_closed": 8.5},
                      "zero_height_quiet_case_admitted": True,
                      "estimated_bias_projection_radius_mps2": .4,
                      "physical_velocity_norm_mps":u['V_m_norm_upper_mps'],
                      "physical_wave_frame_position_norm_m":u['P_m_norm_upper_m']},
        "legacy_active_bias_interior_mps2": .35,
        "primary_motion_target_covers_closed_projection_ball": True,
        "shipping_tuner_default_frequency_hz": [.03, 1.2],
        "proposed_outer_tuner_frequency_hz": [.02, 1.2],
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
                               "S_m":None,
                               "A_m":a},
        "S_primitive_contract":"no independent S_m; S_L=phi(x_s)-phi(x_s_at_Live), dphi/dt=p; the same bounded physical generator owns p/v/a, all moments and every continuation",
        "physical_generator_numeric_family_qualified": wave["physical_D_S_numeric_qualification_closed"],
        "S_three_second_increment_norm_upper_m_s":m['Delta_S_norm_upper_from_uniform_position_m_s'],
        "source_applicability_refined_to_uniform_primitives":True,
        "gravity_direction_low_frequency_qualification_required":True,
        "gravity_direction_same_history_decomposition":direction['direction_forcing_decomposition'],
        "gravity_direction_mean_chord_norm_upper":direction['mean_direction_chord_norm_upper'],
        "gravity_direction_primitive_norm_upper_s":direction['direction_primitive_norm_upper_s'],
        "pointwise_8mps2_cap_retained":direction['pointwise_8mps2_cap_retained'],
        "mean_impulse_bound_applies_to_Q_and_O": True,
        "Q_is_not_a_global_minimum_excitation_requirement": True,
        "Q_AC_energy_alone_excludes_constant_acceleration": False,
        "fixed_window_impulse_cap_alone_excludes_all_DC": False,
        "nonzero_fixed_physical_acceleration_DC_prohibited_in_Q_and_O": True,
        "bounded_velocity_primitive_required": True,
        "zero_DC_definition": "a_m=dv_m/dt and sup_t ||v_m(t)||<=V_m=5 m/s; hence uniform long-time mean acceleration is zero",
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
