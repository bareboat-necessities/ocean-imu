"""Explicit premises of the BRMM P3 matrix implication.

This manifest is not an execution admission test and does not assume P3's
conclusion. Physical BRMM recurrence alone does not guarantee vector PE,
accepted packets, transparent frontend operation or successful factorizations.
The quantitative builders prove their conclusions from these named premises.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import ou3_brmm_contract as BRMM

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = ROOT / "tools/stability/ou3_proof_operating_domain.json"
MEKF = ROOT / "src/kalman_ou_iii/Kalman3D_Wave_OU_III.h"
BRANCHES = ("Q", "O", "MIXED")
SCOPE = "BRMM_AND_EXPLICIT_EXECUTION_PREMISES_IMPLY_P3_NOT_PHYSICAL_ADMISSION"


def projection_covariance_parity(text: str) -> bool:
    """Match the complete shipping projection body; unknown changes fail closed.

    The body modifies only the nominal bias. It has no P, Psi or Omega write,
    so the P3 matrix comparison is unchanged by this event. Its nonlinear
    error map and influence on future coefficients remain P4 obligations.
    """
    signature = "void Kalman3D_Wave_OU_III<T, with_gyro_bias, with_accel_bias>::project_acc_bias_()"
    if text.count(signature) != 1:
        return False
    tail = text.split(signature, 1)[1].lstrip()
    depth = 0
    end = None
    for i, c in enumerate(tail):
        depth += (c == "{") - (c == "}")
        if c == "}" and depth == 0:
            end = i + 1
            break
    expected = """{
        if constexpr (with_accel_bias) {
            if (!(acc_bias_limit_ > T(0))) return;
            auto b = xext.template segment<3>(OFF_BA);
            if (!b.allFinite()) { b.setZero(); return; }
            const T n = b.norm();
            if (n > acc_bias_limit_) { b *= (acc_bias_limit_ / n); }
        }
    }"""
    return end is not None and "".join(tail[:end].split()) == "".join(expected.split())


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    raw = Path(domain_path).read_bytes()
    domain = json.loads(raw)
    live = domain["normal_live"]
    brmm = BRMM.build()
    return {
        "qualification": "OU3_BRMM_P3_EXPLICIT_EXECUTION_PREMISES_V1",
        "scope": SCOPE,
        "domain_sha256": hashlib.sha256(raw).hexdigest(),
        "physical_BRMM_contract": brmm,
        "execution_premises": {
            "same_complete_physical_frontend_tuner_geometry_history": True,
            "declared_Normal_Live_branch_not_just_runtime_Live_flag": True,
            "lever_arm_disabled_and_vibration_guard_dormant_transparent": True,
            "accepted_accelerometer_at_every_valid_IMU_sample": live[
                "accelerometer_update_required_each_valid_imu_sample_after_live_entry"],
            "no_accelerometer_rejection_in_word": not live["accelerometer_rejection_in_normal_live_scope"],
            "vector_PE": {k: live[k] for k in (
                "vector_pe_recurrence_window_s", "specific_force_norm_lower_mps2",
                "specific_force_norm_upper_mps2", "magnetic_vector_norm_lower_uT",
                "magnetic_vector_norm_upper_uT", "vector_sine_separation_lower",
                "body_rate_norm_upper_deg_s")},
            "PE_geometry_attached_to_actual_measurement_jacobians": True,
            "committed_OU_parameters_obey_shipping_clamp_and_commit_invariants": True,
            "every_due_S_with_actual_anisotropic_RS_and_full_cross_covariance": True,
            "shipping_configured_R_full_Q_and_PSD_aw_floors": True,
            "successful_Joseph_updates_and_immediate_covariance_resets": True,
            "no_hard_attitude_rewrite_inside_same_mode_word": not live["hard_attitude_rewrite_inside_word"],
            "source_generated_Live_seed_and_shipping_H_to_A_release": True,
            "premises_continue_after_every_projection_and_hybrid_entry": True,
        },
        "claim": "Omega_W - 1e-18 P_W >= 0 for full H18 and A21 matrix recursion",
        "Psi_scope": "formal covariance-word transport; not the derivative of the projected nonlinear error map",
        "P3_conclusion_is_an_assumption": False,
        "runtime_Live_flag_implies_all_premises": False,
        "BRMM_alone_implies_vector_PE": False,
        "spectral_membership_required": False,
        "bias_error_decay_required": False,
        "bias_estimate_0p35_interior_required_for_covariance_algebra": False,
        "closed_bias_projection_radius_mps2": live["active_accelerometer_bias_projection_limit_mps2"],
        "projection_leaves_covariance_unchanged_source_parity": projection_covariance_parity(MEKF.read_text()),
        "projection_error_sector_or_prefix_retention_proved_here": False,
        "physical_execution_admission_proved_here": False,
        "nonlinear_P4_proved_here": False,
        "delta": 1e-18,
    }


def validate(d: dict, domain_path: Path = DEFAULT_DOMAIN) -> list[str]:
    expected = build(domain_path)
    f = [f"P3 premise binding changed: {k}" for k, v in expected.items() if d.get(k) != v]
    if set(d) != set(expected):
        f.append("P3 premise manifest has unknown or missing fields")
    if not d.get("projection_leaves_covariance_unchanged_source_parity"):
        f.append("shipping projection covariance identity no longer established")
    return f


def conditional_coverage(h18: bool, a21: bool) -> dict:
    """Branch-independent implication, never a source membership assertion."""
    return {branch: {"H18": bool(h18), "A21": bool(a21), "scope": SCOPE,
                     "physical_admission_certified": False} for branch in BRANCHES}


def coverage_closed(d: dict) -> bool:
    return d == conditional_coverage(True, True)
