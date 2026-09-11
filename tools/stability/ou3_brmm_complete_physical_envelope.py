#!/usr/bin/env python3
"""Hard numerical physical envelope for the COMPLETE-BRMM proof family.

This is a theorem-domain design envelope, not a replay fit.  The existing
Hs=8.5 m reference sea and its already-declared physical bounds are used as the
engineering anchor, then padded outward before being admitted to the theorem.
The padding is fixed here and is not selected by P4 success/failure.

Reference anchors:
  Hs <= 8.5 m
  ||p_wave|| <= 7.361215932167728 m
  ||v_wave|| <= 5.0 m/s
  ||a_wave|| <= 8.0 m/s^2
  ||omega_body|| <= 30 deg/s
  finite-harmonic 28-ft RAO D_S <= 863.7793367602501 m*s
  reference harmonic support f in [0.02, 0.8] Hz

The complete proof family applies a 10% amplitude/kinematic pad and a 10%
outward support pad.  Acceleration remains below standard gravity so the
existing gravity-direction proof retains a strictly positive specific-force
lower bound.  D_S is derived outwardly from the padded amplitude and padded
low-frequency edge, then rounded upward to a declared hard theorem constant.
"""
from __future__ import annotations

import json
import math

QUALIFICATION = "OU3_COMPLETE_BRMM_PADDED_PHYSICAL_ENVELOPE_V1"
REFERENCE_HS_M = 8.5
REFERENCE_POSITION_NORM_M = 7.361215932167728
REFERENCE_VELOCITY_NORM_MPS = 5.0
REFERENCE_ACCELERATION_NORM_MPS2 = 8.0
REFERENCE_BODY_RATE_DEG_S = 30.0
REFERENCE_F_MIN_HZ = 0.02
REFERENCE_F_MAX_HZ = 0.8
REFERENCE_D_S_M_S = 863.7793367602501
AMPLITUDE_PAD = 1.10
FREQUENCY_LOW_FACTOR = 0.90
FREQUENCY_HIGH_FACTOR = 1.10
GRAVITY_MPS2 = 9.80665

# Deliberately simple, reviewable outward theorem constants.
HS_MAX_M = 9.35
POSITION_NORM_MAX_M = 8.10
VELOCITY_NORM_MAX_MPS = 5.50
ACCELERATION_NORM_MAX_MPS2 = 8.80
BODY_RATE_NORM_MAX_DEG_S = 35.0
F_MIN_HZ = 0.018
F_MAX_HZ = 0.88
D_S_MAX_M_S = 1100.0


def build() -> dict:
    derived_ds = REFERENCE_D_S_M_S * AMPLITUDE_PAD / FREQUENCY_LOW_FACTOR
    specific_force_lower = GRAVITY_MPS2 - ACCELERATION_NORM_MAX_MPS2
    specific_force_upper = GRAVITY_MPS2 + ACCELERATION_NORM_MAX_MPS2
    return {
        "qualification": QUALIFICATION,
        "theorem_domain_design_not_replay_fit": True,
        "reference_case": {
            "Hs_m": REFERENCE_HS_M,
            "position_norm_m": REFERENCE_POSITION_NORM_M,
            "velocity_norm_mps": REFERENCE_VELOCITY_NORM_MPS,
            "acceleration_norm_mps2": REFERENCE_ACCELERATION_NORM_MPS2,
            "body_rate_norm_deg_s": REFERENCE_BODY_RATE_DEG_S,
            "frequency_support_hz": [REFERENCE_F_MIN_HZ, REFERENCE_F_MAX_HZ],
            "centered_primitive_D_S_upper_m_s": REFERENCE_D_S_M_S,
        },
        "padding": {
            "amplitude_and_kinematic_factor": AMPLITUDE_PAD,
            "frequency_lower_edge_factor": FREQUENCY_LOW_FACTOR,
            "frequency_upper_edge_factor": FREQUENCY_HIGH_FACTOR,
            "chosen_independently_of_P4_certificate": True,
        },
        "complete_BRMM_hard_bounds": {
            "Hs_upper_m": HS_MAX_M,
            "wave_position_norm_upper_m": POSITION_NORM_MAX_M,
            "wave_velocity_norm_upper_mps": VELOCITY_NORM_MAX_MPS,
            "wave_acceleration_norm_upper_mps2": ACCELERATION_NORM_MAX_MPS2,
            "body_rate_norm_upper_deg_s": BODY_RATE_NORM_MAX_DEG_S,
            "frequency_support_hz": [F_MIN_HZ, F_MAX_HZ],
            "centered_primitive_D_S_upper_m_s": D_S_MAX_M_S,
            "specific_force_norm_lower_mps2": specific_force_lower,
            "specific_force_norm_upper_mps2": specific_force_upper,
        },
        "D_S_derivation": {
            "formula": "D_S,padded <= D_S,8.5 * amplitude_pad / frequency_low_factor",
            "unrounded_upper_m_s": derived_ds,
            "declared_outward_rounded_upper_m_s": D_S_MAX_M_S,
            "legacy_300_m_s_used": False,
        },
        "full_COMPLETE_BRMM_numeric_physical_envelope_closed": True,
        "histories_outside_this_envelope_covered_by_theorem": False,
        "P3_delta": 1e-18,
        "P4_PASS": False,
        "P5_MAY_START": False,
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("qualification") != QUALIFICATION:
        f.append("qualification mismatch")
    if d.get("theorem_domain_design_not_replay_fit") is not True:
        f.append("envelope became replay-fit")
    pad = d.get("padding", {})
    if pad.get("chosen_independently_of_P4_certificate") is not True:
        f.append("padding tied to P4 result")
    b = d.get("complete_BRMM_hard_bounds", {})
    expected = {
        "Hs_upper_m": HS_MAX_M,
        "wave_position_norm_upper_m": POSITION_NORM_MAX_M,
        "wave_velocity_norm_upper_mps": VELOCITY_NORM_MAX_MPS,
        "wave_acceleration_norm_upper_mps2": ACCELERATION_NORM_MAX_MPS2,
        "body_rate_norm_upper_deg_s": BODY_RATE_NORM_MAX_DEG_S,
        "frequency_support_hz": [F_MIN_HZ, F_MAX_HZ],
        "centered_primitive_D_S_upper_m_s": D_S_MAX_M_S,
    }
    for k, v in expected.items():
        if b.get(k) != v:
            f.append("hard bound changed: " + k)
    if not (float(b.get("specific_force_norm_lower_mps2", -1.0)) > 0.0):
        f.append("padded acceleration destroyed gravity-direction lower bound")
    deriv = d.get("D_S_derivation", {})
    unrounded = float(deriv.get("unrounded_upper_m_s", math.inf))
    declared = float(deriv.get("declared_outward_rounded_upper_m_s", -math.inf))
    if not (math.isfinite(unrounded) and unrounded <= declared == D_S_MAX_M_S):
        f.append("D_S was not rounded outward from padded 8.5 reference")
    if deriv.get("legacy_300_m_s_used") is not False:
        f.append("legacy 300 m*s contaminated physical bound")
    if d.get("full_COMPLETE_BRMM_numeric_physical_envelope_closed") is not True:
        f.append("numeric physical envelope not closed")
    if d.get("histories_outside_this_envelope_covered_by_theorem") is not False:
        f.append("theorem overclaims histories outside declared envelope")
    if d.get("P4_PASS") or d.get("P5_MAY_START"):
        f.append("physical envelope improperly promoted P4/P5")
    return f


if __name__ == "__main__":
    d = build(); failures = validate(d)
    print(json.dumps({**d, "validation_pass": not failures,
                      "validation_failures": failures}, indent=2, sort_keys=True))
    raise SystemExit(0 if not failures else 2)
