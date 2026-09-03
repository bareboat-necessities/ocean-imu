#!/usr/bin/env python3
"""Source-correlated signed-Joseph feasibility audit for OU-III P4.

This is deliberately a diagnostic, not a P4 theorem producer.  It answers one
specific question before we build the full 18/21-state nonlinear word: can an
accepted vector operation be scalarized locally after paying both the actual
innovation attenuation and the exact finite-angle residual defect, or must its
signed directional form be retained until complete-word accumulation?

For a pure rotational vector residual at Cayley radius q,

    y = (R-I)v,        h = [c]x v,        eta = y-h,

exact Cayley geometry gives

    y' eta = 0,        ||eta||^2 / ||y||^2 = q^2/4.

The exact Joseph identity contains

    y' S^-1 y - eta' R^-1 eta.

If S <= s_max I and R = r I, its coefficient relative to ||y||^2/r is at
least

    r/s_max - q^2/4.

The first term is evaluated from the *same-history* covariance envelopes in the
canonical P3->P4 metric attachment.  For accelerometer innovation covariance we
retain arbitrary PSD cross covariance through the safe marginal bound

    H P H' <= ( |f| sqrt(U_theta) + sqrt(U_aw) [+ sqrt(U_ba)] )^2 I.

For magnetometer only the attitude marginal enters.  We scan every one of the
800 P2 source endpoints and both retained phase envelopes (stage-boundary and
positive phase), in H=18 and A=21 modes.

A negative local coefficient is not instability and is not P4 failure: it only
proves that per-operation scalarization is too lossy and the signed directional
Joseph forms must survive through prediction/reset until the recurrent word is
scalarized.  Conversely a positive coefficient is still only a local
directional fact and cannot promote P4.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_cayley_sector_certificate as CAYLEY
import ou3_p4_p3_metric_attachment as METRIC
import ou3_p4_vector_remainder_sector as REMAINDER
import ou3_vector_uco_certificate as VECTOR

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
MODES = ("H", "A")
PHASE_ENVELOPES = (
    ("stage_boundary_0", "boundary_history_envelope"),
    ("positive_1_25", "positive_phase_history_envelope"),
)


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def sqrt_up(x: float) -> float:
    if not (math.isfinite(float(x)) and float(x) >= 0.0):
        raise ValueError("finite nonnegative square-root input required")
    return up(math.sqrt(float(x)))


def add_up(*xs: float) -> float:
    y = 0.0
    for x in xs:
        y = up(y + float(x))
    return y


def mul_up(a: float, b: float) -> float:
    return up(float(a) * float(b))


def _positive(x, label: str) -> float:
    y = float(x)
    if not (math.isfinite(y) and y > 0.0):
        raise RuntimeError(f"{label} must be finite positive")
    return y


def _measurement_variances(vector: dict) -> dict:
    c = vector["configured_measurement_bounds"]
    acc_std = _positive(c["acc_measurement_std_mps2"], "accelerometer std")
    mag_std = _positive(c["mag_measurement_std_uT"], "magnetometer std")
    return {
        "acc_lower": down(acc_std * acc_std),
        "acc_upper": up(acc_std * acc_std),
        "mag_lower": down(mag_std * mag_std),
        "mag_upper": up(mag_std * mag_std),
    }


def _innovation_bounds(env: dict, mode: str, fmax: float, mmax: float,
                       rv: dict) -> dict:
    trans = env["translation_covariance_upper_groups"]
    bias = env[f"{mode}_bias_covariance_upper"]
    utheta = _positive(bias["theta_covariance_upper"], f"{mode} theta upper")
    uaw = _positive(trans["a_w"], f"{mode} aw upper")
    uba = 0.0
    if mode == "A":
        uba = _positive(bias["accel_bias_covariance_upper"], "A accel-bias upper")

    # For a PSD joint covariance with marginal ceilings U_i I, every cross
    # covariance is bounded by sqrt(U_i U_j).  Therefore the measurement-space
    # covariance is dominated by the square of the sum of block operator radii.
    acc_radius = add_up(mul_up(fmax, sqrt_up(utheta)), sqrt_up(uaw))
    if mode == "A":
        acc_radius = add_up(acc_radius, sqrt_up(uba))
    acc_state = mul_up(acc_radius, acc_radius)
    mag_state = mul_up(mul_up(mmax, mmax), utheta)
    return {
        "theta_covariance_upper": utheta,
        "aw_covariance_upper": uaw,
        "accel_bias_covariance_upper": None if mode == "H" else uba,
        "accelerometer_state_innovation_covariance_upper": acc_state,
        "magnetometer_state_innovation_covariance_upper": mag_state,
        "accelerometer_S_lambda_max_upper": add_up(rv["acc_upper"], acc_state),
        "magnetometer_S_lambda_max_upper": add_up(rv["mag_upper"], mag_state),
    }


def evaluate(metric: dict, cayley: dict, remainder: dict, vector: dict,
             domain: dict) -> dict:
    mf = METRIC.validate(metric)
    cf = CAYLEY.validate(cayley)
    rf = REMAINDER.validate(remainder)
    vf = VECTOR.validate(vector)
    failures = [f"metric: {x}" for x in mf]
    failures += [f"Cayley: {x}" for x in cf]
    failures += [f"remainder: {x}" for x in rf]
    failures += [f"vector: {x}" for x in vf]
    if failures:
        return {"schema": SCHEMA, "failures": failures}

    if metric.get("same_history_P3_frontier_consumed") is not True:
        failures.append("metric attachment is not same-history P3")
    if metric.get("independent_cartesian_tau_sigma_R_S_extrema_used") is not False:
        failures.append("metric attachment reintroduced Cartesian tuner extrema")
    if metric.get("finite_source_nodes") != 800:
        failures.append("metric attachment does not cover 800 source endpoints")

    live = domain["normal_live"]
    fmax = _positive(live["specific_force_norm_upper_mps2"], "specific-force upper")
    mmax = _positive(live["magnetic_vector_norm_upper_uT"], "magnetic-vector upper")
    rv = _measurement_variances(vector)

    q = _positive(cayley["cayley_radius_upper"], "Cayley radius")
    eta_over_y2 = up(mul_up(q, q) / 4.0)
    rot_minus_I = _positive(remainder["rotation_minus_identity_norm_upper"], "rotation-I norm")
    aw_eta_Rinv_coeff = up(mul_up(rot_minus_I, rot_minus_I) / rv["acc_lower"])

    rows = []
    worst = {
        mode: {
            "accelerometer_force_rotation_signed_fraction_lower": math.inf,
            "magnetometer_rotation_signed_fraction_lower": math.inf,
            "accelerometer_innovation_retention_lower": math.inf,
            "magnetometer_innovation_retention_lower": math.inf,
            "limiting_accelerometer": None,
            "limiting_magnetometer": None,
        }
        for mode in MODES
    }

    for endpoint in metric.get("endpoint_rows", []):
        source = int(endpoint["source_node"])
        for phase_name, phase_key in PHASE_ENVELOPES:
            env = endpoint[phase_key]
            for mode in MODES:
                ib = _innovation_bounds(env, mode, fmax, mmax, rv)
                acc_ret = down(rv["acc_lower"] / float(ib["accelerometer_S_lambda_max_upper"]))
                mag_ret = down(rv["mag_lower"] / float(ib["magnetometer_S_lambda_max_upper"]))
                acc_signed = down(acc_ret - eta_over_y2)
                mag_signed = down(mag_ret - eta_over_y2)
                row = {
                    "source_node": source,
                    "phase_envelope": phase_name,
                    "mode": mode,
                    **ib,
                    "accelerometer_innovation_retention_R_over_S_lower": acc_ret,
                    "magnetometer_innovation_retention_R_over_S_lower": mag_ret,
                    "exact_pure_rotation_eta_squared_over_exact_residual_squared_upper": eta_over_y2,
                    "accelerometer_force_rotation_signed_fraction_lower": acc_signed,
                    "magnetometer_rotation_signed_fraction_lower": mag_signed,
                }
                rows.append(row)
                w = worst[mode]
                if acc_signed < w["accelerometer_force_rotation_signed_fraction_lower"]:
                    w["accelerometer_force_rotation_signed_fraction_lower"] = acc_signed
                    w["accelerometer_innovation_retention_lower"] = acc_ret
                    w["limiting_accelerometer"] = {
                        "source_node": source, "phase_envelope": phase_name,
                        "S_lambda_max_upper": ib["accelerometer_S_lambda_max_upper"],
                    }
                if mag_signed < w["magnetometer_rotation_signed_fraction_lower"]:
                    w["magnetometer_rotation_signed_fraction_lower"] = mag_signed
                    w["magnetometer_innovation_retention_lower"] = mag_ret
                    w["limiting_magnetometer"] = {
                        "source_node": source, "phase_envelope": phase_name,
                        "S_lambda_max_upper": ib["magnetometer_S_lambda_max_upper"],
                    }

    expected = 800 * len(PHASE_ENVELOPES) * len(MODES)
    if len(rows) != expected:
        failures.append(f"signed-Joseph scan covered {len(rows)} classes, expected {expected}")

    local_scalar = all(
        worst[m]["accelerometer_force_rotation_signed_fraction_lower"] > 0.0
        and worst[m]["magnetometer_rotation_signed_fraction_lower"] > 0.0
        for m in MODES
    )
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_SOURCE_CORRELATED_SIGNED_JOSEPH_FEASIBILITY_AUDIT",
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "same_history_P3_metric_consumed": True,
        "independent_cartesian_tau_sigma_R_S_extrema_used": False,
        "outer_angle_rad": float(cayley["outer_angle_rad"]),
        "cayley_radius_upper": q,
        "exact_pure_rotation_eta_squared_over_exact_residual_squared_upper": eta_over_y2,
        "accelerometer_aw_eta_Rinv_quadratic_coefficient_upper": aw_eta_Rinv_coeff,
        "finite_source_phase_mode_classes_scanned": len(rows),
        "expected_source_phase_mode_classes": expected,
        "worst_by_mode": worst,
        "per_operation_scalar_signed_rotation_credit_positive_everywhere": local_scalar,
        "directional_signed_forms_must_survive_to_word_scalarization": not local_scalar,
        "accelerometer_aw_remainder_requires_joint_directional_word_accounting": True,
        "complete_H18_A21_word_established_here": False,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "rows": rows,
        "next_obligation": (
            "if any local signed fraction is nonpositive, do not scalarize vector operations; "
            "transport the source-correlated signed Joseph directional forms through prediction/reset and accumulate the recurrent word before scalarization"
            if not local_scalar else
            "local vector signed credit survives the conservative innovation bound; retain it directionally and construct the complete recurrent H/A word before any P4 promotion"
        ),
        "failures": failures,
    }


def build(metric: dict, domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("signed-Joseph feasibility audit must not be trajectory fitted")
    return evaluate(metric, CAYLEY.build(path), REMAINDER.build(path), VECTOR.build(), domain)


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P4_SOURCE_CORRELATED_SIGNED_JOSEPH_FEASIBILITY_AUDIT":
        f.append("wrong qualification")
    for key in ("source_generated_not_trajectory_fit", "same_history_P3_metric_consumed",
                "accelerometer_aw_remainder_requires_joint_directional_word_accounting"):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in ("trajectory_replay_used", "filter_changed", "declared_domain_changed",
                "independent_cartesian_tau_sigma_R_S_extrema_used",
                "complete_H18_A21_word_established_here", "P4_USABLE_CERTIFICATE_PROMOTED"):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if int(d.get("finite_source_phase_mode_classes_scanned", 0)) != 800 * 2 * 2:
        f.append("signed-Joseph audit did not cover all 800 x 2 x H/A classes")
    if d.get("directional_signed_forms_must_survive_to_word_scalarization") is not (
        not bool(d.get("per_operation_scalar_signed_rotation_credit_positive_everywhere"))
    ):
        f.append("directional/scalarization verdict is inconsistent")
    for mode in MODES:
        row = d.get("worst_by_mode", {}).get(mode, {})
        for key in ("accelerometer_force_rotation_signed_fraction_lower",
                    "magnetometer_rotation_signed_fraction_lower",
                    "accelerometer_innovation_retention_lower",
                    "magnetometer_innovation_retention_lower"):
            x = row.get(key)
            if not isinstance(x, (int, float)) or isinstance(x, bool) or not math.isfinite(float(x)):
                f.append(f"{mode}: {key} is not finite")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--metric", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    metric = json.loads(a.metric.read_text(encoding="utf-8"))
    d = build(metric, a.domain)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "eta2_over_y2_upper": d.get("exact_pure_rotation_eta_squared_over_exact_residual_squared_upper"),
        "H": d.get("worst_by_mode", {}).get("H"),
        "A": d.get("worst_by_mode", {}).get("A"),
        "local_scalar_signed_credit_everywhere": d.get("per_operation_scalar_signed_rotation_credit_positive_everywhere"),
        "directional_word_required": d.get("directional_signed_forms_must_survive_to_word_scalarization"),
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
