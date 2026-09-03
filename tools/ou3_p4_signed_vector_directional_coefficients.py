#!/usr/bin/env python3
"""Signed finite-angle vector coefficients against the P3 tangent packet forms.

The first signed-Joseph audit reports innovation retention R/S on every retained
800 x 2 x H/A source/phase class.  This producer converts those numbers into
coefficients multiplying the *linear tangent measurement forms* used by P3.
That is the natural interface for the complete-word directional accumulator.

Accelerometer.
The exact co-rotated a_w coordinate makes eta pure attitude rotation.  Cayley
geometry gives

    ||y_R||^2 >= a_R ||h_R||^2,       a_R=4/(4+q^2),
    ||eta_R||^2 <= s^2 ||h_R||^2,     s=sin(theta/2).

The current Cayley certificate exposes ``a_R`` as
``chart_sigma_min_lower``: for c=q u, the exact residual-to-linear-tangent
amplitude ratio is cos(theta/2), hence its squared ratio is
cos^2(theta/2)=4/(4+q^2).  This is deliberately *not* the separate
``exact_vector_information_retention_factor_lower=cos^4(theta/2)`` used for the
Jacobian information bound.

Thus on the pure rotational tangent direction

    y' S^-1 y - eta' R^-1 eta
      >= [ (R/S_max) a_R - s^2 ] h' R^-1 h.

The linear a_w and b_a directions remain inside y and inside S; they are not
charged again as nonlinear eta.

Magnetometer.
The radial Joseph contribution cancels exactly.  With the effective tangent
coordinate d from the retained radial-cancellation primitive,

    ||d|| >= A ||c_perp||,             A=4/(4+q^2),
    ||d-c_perp|| <= beta ||c_perp||,   beta=q/sqrt(4+q^2).

Since ||H x||=||v|| ||x|| on v-perp,

    D_m >= [ (R/S_max) A^2 - beta^2 ] h_m' R^-1 h_m.

These are still local directional coefficients, not full-state scalar packet
margins.  A nonpositive coefficient means the operation must remain signed and
be accumulated with the rest of the recurrent word; it is not a stability
failure.  A positive coefficient likewise does not promote P4 by itself.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_cayley_sector_certificate as CAYLEY
import ou3_p4_magnetometer_radial_joseph as MAG
import ou3_p4_vector_remainder_sector as REMAINDER
import ou3_p4_signed_joseph_feasibility as AUDIT

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 2
MODES = ("H", "A")


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def mul_down(a: float, b: float) -> float:
    return down(float(a) * float(b))


def square_up(x: float) -> float:
    return up(float(x) * float(x))


def evaluate(audit: dict, cayley: dict, remainder: dict, mag: dict) -> dict:
    failures = [f"audit: {x}" for x in AUDIT.validate(audit)]
    failures += [f"Cayley: {x}" for x in CAYLEY.validate(cayley)]
    failures += [f"remainder: {x}" for x in REMAINDER.validate(remainder)]
    failures += [f"magnetometer: {x}" for x in MAG.validate(mag)]
    if failures:
        return {"schema": SCHEMA, "failures": failures}

    # For the exact vector residual r=(R-I)v versus h=[c]xv, the *squared*
    # norm retention is cos^2(theta/2)=4/(4+q^2).  The Cayley producer's
    # chart_sigma_min_lower is exactly that scalar.  Do not use its cos^4
    # Jacobian-information factor here.
    aR = float(cayley["chart_sigma_min_lower"])
    jacobian_info = float(cayley["exact_vector_information_retention_factor_lower"])
    s2 = float(remainder["acc_eta_force_rotation_quadratic_coefficient_upper"])
    A = float(mag["effective_tangent_coordinate_gain_lower"])
    beta = float(mag["effective_vs_linear_tangent_defect_ratio_upper"])
    if not (0.0 < aR <= 1.0 and 0.0 < jacobian_info <= aR <= 1.0):
        failures.append("Cayley residual/information retention factors invalid")
    if not (0.0 <= s2 < 1.0 and 0.0 < A <= 1.0 and 0.0 <= beta < 1.0):
        failures.append("finite-angle tangent coefficients invalid")

    A2_lo = mul_down(A, A)
    beta2_hi = square_up(beta)
    worst = {
        m: {
            "accelerometer_tangent_signed_coefficient_lower": math.inf,
            "magnetometer_tangent_signed_coefficient_lower": math.inf,
            "limiting_accelerometer": None,
            "limiting_magnetometer": None,
        }
        for m in MODES
    }
    count = 0
    for row in audit.get("rows", []):
        mode = row.get("mode")
        if mode not in MODES:
            failures.append("signed audit row has unknown mode")
            continue
        ar = float(row["accelerometer_innovation_retention_R_over_S_lower"])
        mr = float(row["magnetometer_innovation_retention_R_over_S_lower"])
        acc = down(mul_down(ar, aR) - s2)
        magc = down(mul_down(mr, A2_lo) - beta2_hi)
        count += 1
        w = worst[mode]
        if acc < w["accelerometer_tangent_signed_coefficient_lower"]:
            w["accelerometer_tangent_signed_coefficient_lower"] = acc
            w["limiting_accelerometer"] = {
                "source_node": int(row["source_node"]),
                "phase_envelope": row["phase_envelope"],
                "innovation_retention_lower": ar,
            }
        if magc < w["magnetometer_tangent_signed_coefficient_lower"]:
            w["magnetometer_tangent_signed_coefficient_lower"] = magc
            w["limiting_magnetometer"] = {
                "source_node": int(row["source_node"]),
                "phase_envelope": row["phase_envelope"],
                "innovation_retention_lower": mr,
            }

    expected = int(audit.get("expected_source_phase_mode_classes", 0))
    if count != expected or expected != 800 * 2 * 2:
        failures.append(f"directional coefficient scan covered {count}, expected {expected}")

    local_positive = all(
        worst[m]["accelerometer_tangent_signed_coefficient_lower"] > 0.0
        and worst[m]["magnetometer_tangent_signed_coefficient_lower"] > 0.0
        for m in MODES
    )
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_SOURCE_CORRELATED_SIGNED_VECTOR_TANGENT_COEFFICIENTS",
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "same_history_signed_Joseph_audit_consumed": True,
        "accelerometer_corotated_aw_eta_zero_consumed": True,
        "magnetometer_radial_Joseph_cancellation_consumed": True,
        "linear_tangent_form_interface_matches_P3": True,
        "outer_angle_rad": float(cayley["outer_angle_rad"]),
        "accelerometer_exact_residual_vs_tangent_norm_squared_lower": aR,
        "cayley_jacobian_information_retention_factor_lower_not_used_for_residual": jacobian_info,
        "accelerometer_eta_vs_tangent_norm_squared_upper": s2,
        "magnetometer_effective_tangent_gain_lower": A,
        "magnetometer_effective_tangent_gain_squared_lower": A2_lo,
        "magnetometer_effective_tangent_defect_ratio_upper": beta,
        "magnetometer_effective_tangent_defect_squared_upper": beta2_hi,
        "source_phase_mode_classes_scanned": count,
        "worst_by_mode": worst,
        "local_vector_tangent_signed_coefficients_positive_everywhere": local_positive,
        "word_level_directional_accumulation_required": not local_positive,
        "instantaneous_full_state_scalarization_attempted": False,
        "complete_H18_A21_word_established_here": False,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "next_obligation": (
            "transport these signed tangent forms with the S=0 linear form through the recurrent prediction/reset word and scalarize only after the H/A word has full rank"
        ),
        "failures": failures,
    }


def build(audit: dict, domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    return evaluate(
        audit,
        CAYLEY.build(path),
        REMAINDER.build(path),
        MAG.build(path),
    )


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P4_SOURCE_CORRELATED_SIGNED_VECTOR_TANGENT_COEFFICIENTS":
        f.append("wrong qualification")
    for key in (
        "source_generated_not_trajectory_fit", "same_history_signed_Joseph_audit_consumed",
        "accelerometer_corotated_aw_eta_zero_consumed",
        "magnetometer_radial_Joseph_cancellation_consumed",
        "linear_tangent_form_interface_matches_P3",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_changed",
        "instantaneous_full_state_scalarization_attempted",
        "complete_H18_A21_word_established_here", "P4_USABLE_CERTIFICATE_PROMOTED",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if int(d.get("source_phase_mode_classes_scanned", 0)) != 800 * 2 * 2:
        f.append("directional coefficient audit did not cover 800 x 2 x H/A")
    residual = d.get("accelerometer_exact_residual_vs_tangent_norm_squared_lower")
    jac = d.get("cayley_jacobian_information_retention_factor_lower_not_used_for_residual")
    if not isinstance(residual, (int, float)) or not (0.0 < float(residual) <= 1.0):
        f.append("invalid accelerometer exact-residual retention factor")
    if not isinstance(jac, (int, float)) or not (0.0 < float(jac) <= float(residual or 0.0)):
        f.append("invalid separated Cayley Jacobian-information factor")
    positive = True
    for mode in MODES:
        row = d.get("worst_by_mode", {}).get(mode, {})
        for key in (
            "accelerometer_tangent_signed_coefficient_lower",
            "magnetometer_tangent_signed_coefficient_lower",
        ):
            x = row.get(key)
            if not isinstance(x, (int, float)) or isinstance(x, bool) or not math.isfinite(float(x)):
                f.append(f"{mode}: {key} is not finite")
                positive = False
            else:
                positive = positive and float(x) > 0.0
    if d.get("local_vector_tangent_signed_coefficients_positive_everywhere") is not positive:
        f.append("local tangent positivity verdict inconsistent")
    if d.get("word_level_directional_accumulation_required") is not (not positive):
        f.append("word-level directional accumulation verdict inconsistent")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--signed-audit", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    audit = json.loads(a.signed_audit.read_text(encoding="utf-8"))
    d = build(audit, a.domain)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "H": d.get("worst_by_mode", {}).get("H"),
        "A": d.get("worst_by_mode", {}).get("A"),
        "local_tangent_positive": d.get("local_vector_tangent_signed_coefficients_positive_everywhere"),
        "word_directional_required": d.get("word_level_directional_accumulation_required"),
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
