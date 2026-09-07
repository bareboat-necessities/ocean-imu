"""Primary BRMM physical-source target; legacy SEA3 certificates stay scoped.

This declaration does not turn sampled recurrence statistics into uniform
source coverage or transfer a SEA3 P3 certificate to a different source family.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def build():
    return {
        "qualification": "OU3_COMPLETE_BRMM_BOUNDED_BIAS_MOTION_V1",
        "canonical_source": "COMPLETE_BRMM_NORMAL_LIVE_WORD",
        "spectral_SEA3_membership_required": False,
        "common_source_frontend_tuner_geometry_and_filter_history_required": True,
        "every_due_actual_anisotropic_RS_required": True,
        "fixed_frame_CoG_acceleration_required": True,
        "hard_caps": {"acceleration_norm_mps2": 4., "body_rate_deg_s": 30.,
                      "specific_force_norm_mps2": [5.80665, 13.80665],
                      "positive_sea_height_m": {"lower_open": 0., "upper_closed": 8.5},
                      "zero_height_quiet_case_admitted": True,
                      "estimated_bias_projection_radius_mps2": .4},
        "legacy_active_bias_interior_mps2": .35,
        "primary_motion_target_covers_closed_projection_ball": True,
        "shipping_tuner_default_frequency_hz": [.03, 1.2],
        "proposed_outer_tuner_frequency_hz": [.02, 1.2],
        "tuner_frequency_is_physical_motion_bandlimit": False,
        "candidate_audit": {"recurrence_window_s": 10., "quiet_rms_mps2": [.03, .05],
                            "impulse_cap_mps": 2., "lobe_chi_for_diagnostic_only": 1.},
        "unfrozen_physical_constants": {"T_R": None, "E_q": None, "V_R": None,
                                        "E_min": None, "J_min": None, "chi": None,
                                        "V_a": None, "V_m": None, "P_m": None, "S_m": None},
        "mean_impulse_bound_applies_to_Q_and_O": True,
        "Q_is_not_a_global_minimum_excitation_requirement": True,
        "Q_AC_energy_alone_excludes_constant_acceleration": False,
        "fixed_window_impulse_cap_alone_excludes_all_DC": False,
        "nonzero_fixed_physical_acceleration_DC_prohibited_in_Q_and_O": True,
        "bounded_velocity_primitive_required": True,
        "zero_DC_definition": "a_m=dv_m/dt and sup_t ||v_m(t)||<=V_m<infinity; hence uniform long-time mean acceleration is zero",
        "constant_measurement_bias_is_marine_motion": False,
        "bounded_same_history_primitives_or_explicit_forcing_budget_required": True,
        "Q_and_O_mode_mixes_and_vector_PE_require_separate_coverage": True,
        "bias_error_decay_required": False,
        "zero_navigation_error_floor_required": False,
        "P3_delta": 1e-18,
        "legacy_P3_matrix_implication_reusable_with_its_actual_premises": True,
        "legacy_SEA3_P3_flags_renamed_to_BRMM": False,
        "BRMM_SOURCE_UNIFORM_PASS": False,
        "BRMM_P3_PASS": False,
        "BRMM_P4_MOTION_PASS": False,
        "BRMM_P5_MOTION_MAY_START": False,
        "sampled_audit_sets_physical_theorem_constants": False,
    }


def validate(d):
    return [f"{k} changed or falsely certified" for k, v in build().items() if d.get(k) != v]


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", required=True, type=Path)
    args = p.parse_args()
    args.output.write_text(json.dumps(build(), indent=2)+"\n")
