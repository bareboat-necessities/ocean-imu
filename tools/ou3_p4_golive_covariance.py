#!/usr/bin/env python3
"""Source-derived goLive covariance seed used by the active OU-III P4 route.

Only the current handoff covariance facts required by P4 are retained here.
Historical P5 entrance/prefix search bookkeeping belongs to the PR history and
is intentionally not part of this module.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import ou3_implementation_proof_manifest as MANIFEST
import ou3_source_domain_contract as SOURCE

REPO = Path(__file__).resolve().parents[1]
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
MEKF = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"
CORE = REPO / "src" / "kalman_ou_common" / "KalmanOUCoreMath.h"
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"


def _one(text: str, pattern: str, label: str) -> float:
    m = re.search(pattern, text, flags=re.MULTILINE | re.DOTALL)
    if not m:
        raise RuntimeError(f"cannot extract {label}")
    return float(m.group(1))


def _require(text: str, marker: str, label: str) -> None:
    if marker not in text:
        raise RuntimeError(f"missing source semantic {label}: {marker}")


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain_path = Path(domain_path).resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("goLive covariance seed must not be trajectory fitted")

    w = WRAPPER.read_text(encoding="utf-8")
    k = MEKF.read_text(encoding="utf-8")
    c = CORE.read_text(encoding="utf-8")
    manifest = MANIFEST.build()
    mf = MANIFEST.validate(manifest)
    if mf:
        raise RuntimeError(f"implementation manifest prerequisite failed: {mf}")

    joined = w + "\n" + k + "\n" + c
    for label, marker in {
        "go_live_initializes_attitude": "mekf_->initialize_from_attitude(q_bw, tilt_sigma_rad, yaw_sigma_rad);",
        "live_applies_ou_before_enable": "apply_ou_tune_(true);",
        "bootstrap_withholds_time_update": "impl_.updateFrontEnd(dt, gyro_body_ned, acc_body_ned);",
        "attitude_init_zeros_AL": "zero_AL_cross_cov_once_();",
        "live_resets_aw": "reset_aw_covariance_to_stationary();",
        "aw_reset_zeros_cross": "Pext.template block<3,1>(OFF_AW, i).setZero();",
    }.items():
        _require(joined, marker, label)

    sigma_v0 = _one(k, r"const\s+T\s+sigma_v0\s*=\s*T\(([0-9.eE+-]+)\)", "sigma_v0")
    sigma_p0 = _one(k, r"const\s+T\s+sigma_p0\s*=\s*T\(([0-9.eE+-]+)\)", "sigma_p0")
    sigma_S0 = _one(k, r"const\s+T\s+sigma_S0\s*=\s*T\(([0-9.eE+-]+)\)", "sigma_S0")

    startup = manifest["startup"]
    tilt_sigma = float(startup["handoff_tilt_sigma_rad"])
    yaw_sigma = float(startup["handoff_yaw_sigma_rad"])
    yaw_free_sigma = float(startup["handoff_yaw_sigma_free_rad"])
    pbox = SOURCE.build(WRAPPER)["validated_parameter_box"]["continuous_parameters"]
    sigma_aw_lo, sigma_aw_hi = map(float, pbox["sigma_aw_mps2"])

    seed = {
        "mode": "H",
        "dimension": 18,
        "mekf_was_not_propagated_during_bootstrap": True,
        "attitude_linear_cross_covariance_exact_zero": True,
        "P_vv_variance_per_axis": sigma_v0 * sigma_v0,
        "P_pp_variance_per_axis": sigma_p0 * sigma_p0,
        "P_SS_variance_per_axis": sigma_S0 * sigma_S0,
        "P_awaw_reset_to_current_stationary_covariance": True,
        "P_awaw_source_std_outward_mps2": [sigma_aw_lo, sigma_aw_hi],
        "attitude_covariance_seed": {
            "tilt_sigma_rad": tilt_sigma,
            "gauged_yaw_sigma_rad": yaw_sigma,
            "ungauged_yaw_sigma_rad": yaw_free_sigma,
            "tilt_variance": tilt_sigma * tilt_sigma,
            "gauged_yaw_variance": yaw_sigma * yaw_sigma,
            "ungauged_yaw_variance": yaw_free_sigma * yaw_free_sigma,
        },
    }
    return {
        "source_generated_not_trajectory_fit": True,
        "filter_changed": False,
        "goLive_H_covariance_seed": seed,
    }


def validate(d: dict) -> list[str]:
    failures = []
    if d.get("source_generated_not_trajectory_fit") is not True:
        failures.append("goLive seed is not source generated")
    if d.get("filter_changed") is not False:
        failures.append("goLive seed changes filter")
    seed = d.get("goLive_H_covariance_seed", {})
    for key in (
        "mekf_was_not_propagated_during_bootstrap",
        "attitude_linear_cross_covariance_exact_zero",
        "P_awaw_reset_to_current_stationary_covariance",
    ):
        if seed.get(key) is not True:
            failures.append(f"goLive seed lost {key}")
    return failures
