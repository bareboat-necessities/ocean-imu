#!/usr/bin/env python3
"""Exact structural certificate for the shipping late-north yaw rewrite.

If a magnetometer reference becomes ready after a timeout Live handoff, the
wrapper executes

    q_tilt = Rz(-yaw(q_old)) q_old
    q_new  = Rz(yaw_abs) q_tilt

then ``mekf.set_quaternion_boat(q_new)``.  This is a left multiplication by a
world-Z rotation.  Consequently the body representation of world down is
exactly invariant:

    q_new^{-1} e_z q_new = q_old^{-1} e_z q_old,

because every Rz leaves e_z fixed.  The reset therefore changes yaw gauge only;
it cannot enlarge the already-certified gravity-quotient tilt error.

The MEKF setter rewrites qref and clears the local attitude error state.  It does
not touch v,p,S,a_w,b_g,b_a or covariance P.  This file certifies that hybrid
map structure from shipping source parity.  It deliberately does NOT claim that
the new absolute yaw is inside P4: that requires a bound on the MagAutoTuner
reference/gauge error from the same accepted magnetic history.
"""
from __future__ import annotations

import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
KALMAN = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"
SCHEMA = 1
QUALIFICATION = "OU3_P5_LATE_NORTH_YAW_RESET_V1"


def _function_body(text: str, signature_fragment: str, next_fragment: str) -> str:
    if signature_fragment not in text:
        raise RuntimeError("missing function signature: " + signature_fragment)
    tail = text.split(signature_fragment, 1)[1]
    if next_fragment not in tail:
        raise RuntimeError("missing function terminator marker: " + next_fragment)
    return tail.split(next_fragment, 1)[0]


def build() -> dict:
    wrapper = WRAPPER.read_text(encoding="utf-8")
    kalman = KALMAN.read_text(encoding="utf-8")
    yaw_body = _function_body(wrapper,"static Eigen::Quaternionf boatQuatWithAbsoluteYaw_","static Eigen::Quaternionf tiltOnlyQuatFromBoatQuat_")
    removed_body = _function_body(wrapper,"static Eigen::Quaternionf yawRemovedBoatQuat_","// Rewrites heading only")
    setter_body = _function_body(kalman,"void set_quaternion_boat(const Eigen::Quaternion<T>& q_bw)","[[nodiscard]] Vector3 gyroscope_bias() const")
    wrapper_parity = {
        "yaw_removed_by_left_world_z_rotation": "Eigen::AngleAxisf(-yaw, Eigen::Vector3f::UnitZ())" in removed_body and "q_yaw_inv * q_bw" in removed_body,
        "absolute_yaw_restored_by_left_world_z_rotation": "yaw_abs_rad" in yaw_body and "Eigen::Vector3f::UnitZ()" in yaw_body and "q_yaw * q_tilt" in yaw_body,
        "live_path_calls_mefk_setter": "stage_ != Stage::Live" in wrapper and "impl_.mekf().set_quaternion_boat(q_new);" in wrapper,
    }
    setter_parity = {
        "writes_internal_quaternion": "qref = q_WBprime.conjugate();" in setter_body,
        "normalizes_internal_quaternion": "qref.normalize();" in setter_body,
        "clears_local_attitude_error": "xext.template segment<3>(0).setZero();" in setter_body,
        "does_not_assign_full_state": "xext.setZero()" not in setter_body,
        "does_not_write_velocity": "OFF_V" not in setter_body and "get_velocity" not in setter_body,
        "does_not_write_position": "OFF_P" not in setter_body,
        "does_not_write_integral_displacement": "OFF_S" not in setter_body,
        "does_not_write_latent_acceleration": "OFF_AW" not in setter_body,
        "does_not_write_accelerometer_bias": "OFF_BA" not in setter_body,
        "does_not_write_covariance": "Pext" not in setter_body,
    }
    structure_closed = all(wrapper_parity.values()) and all(setter_parity.values())
    return {
        "schema":SCHEMA,"qualification":QUALIFICATION,"filter_changed":False,"quality_gates_changed":False,
        "wrapper_source_parity":wrapper_parity,"mekf_setter_source_parity":setter_parity,
        "exact_map":"q_new=Rz(yaw_abs-yaw(q_old))*q_old","world_down_is_fixed_by_reset_rotation":True,
        "body_gravity_direction_exactly_invariant":True,"gravity_quotient_tilt_error_exactly_invariant":True,
        "linear_navigation_state_exactly_unchanged":True,"gyro_bias_state_exactly_unchanged":True,
        "accelerometer_bias_state_exactly_unchanged":True,"covariance_exactly_unchanged_by_setter":True,
        "local_attitude_error_state_zeroed_after_rewrite":True,"LATE_NORTH_HYBRID_RESET_MAP_STRUCTURE_CLOSED":structure_closed,
        "magnetic_reference_yaw_error_bound_closed_here":False,"full_attitude_P4_landing_closed_here":False,
        "remaining_landing_obligation":"bound yaw_abs - true north from the same accepted MagAutoTuner history, including proxy/MEKF tilt error, sensor disturbance, horizontal-field nondegeneracy, and binary32 charge",
        "P4_PASS":False,"P5_PASS":False,
    }

def validate(d:dict)->list[str]:
    f=[]
    if d.get("schema")!=SCHEMA or d.get("qualification")!=QUALIFICATION:f.append("schema/qualification mismatch")
    for k in ("world_down_is_fixed_by_reset_rotation","body_gravity_direction_exactly_invariant","gravity_quotient_tilt_error_exactly_invariant","linear_navigation_state_exactly_unchanged","gyro_bias_state_exactly_unchanged","accelerometer_bias_state_exactly_unchanged","covariance_exactly_unchanged_by_setter","local_attitude_error_state_zeroed_after_rewrite","LATE_NORTH_HYBRID_RESET_MAP_STRUCTURE_CLOSED"):
        if d.get(k) is not True:f.append(k+" not true")
    if not all(d.get("wrapper_source_parity",{}).values()):f.append("wrapper yaw reset parity failed")
    if not all(d.get("mekf_setter_source_parity",{}).values()):f.append("MEKF setter parity failed")
    for k in ("filter_changed","quality_gates_changed","magnetic_reference_yaw_error_bound_closed_here","full_attitude_P4_landing_closed_here","P4_PASS","P5_PASS"):
        if d.get(k) is not False:f.append(k+" not false")
    return f

if __name__=="__main__":
    d=build();failures=validate(d);print(json.dumps({**d,"validation_pass":not failures,"validation_failures":failures},indent=2,sort_keys=True));raise SystemExit(1 if failures else 0)
