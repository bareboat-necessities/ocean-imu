#!/usr/bin/env python3
"""Exact Joseph/quaternion/reset transport primitive for current OU-III P4.

After an accepted Kalman correction, write the signed nonlinear residual as

    y = H z + eta,

and choose signs so the tangent posterior is t=z-Ky.  For the implemented Joseph
posterior P+,

    z'P^-1 z - t'(P+)^-1 t = y'S^-1y - eta'R^-1eta.       (1)

Shipping then injects the attitude correction d with its normalized quaternion
and applies the left-error covariance reset

    G(d)=I+0.5[d]x,       P_r=G_ext P+ G_ext'.

Let z_exact be the physical post-injection error and

    rho = z_exact - G_ext t.

Then congruence gives exactly

    z_exact'P_r^-1 z_exact
      = (t+G_ext^-1 rho)'(P+)^-1(t+G_ext^-1 rho).          (2)

Furthermore G'G has eigenvalues 1, 1+|d|^2/4, 1+|d|^2/4, so

    ||G^-1||_2 = 1

for every finite correction.  There is no reset condition-number penalty.

The only nonlinear reset mismatch is rho.  If a is the Cayley vector of the
*deployed correction quaternion*, then a is parallel to d and exact Cayley
composition gives

    c+ = (a+c+0.5 a x c)/(1-0.25 a'c).

Against the homogeneous reset target G(d)(c+d), for |c|<=q, |d|<=delta and
D=1-|a|q/4>0,

    |rho_theta| <= [ |a-d|(1+q/2)
                    + |a|q/4 (q+delta+delta q/2) ] / D.   (3)

The functions below outward-enclose a and rho across the shipping polynomial
quaternion branch below 1e-2 and the source axis-angle branch above it.  This
producer is intentionally parameterized in the correction radius delta: the
complete-word route must obtain delta from the same source-correlated operation
cell, not invent a global correction bound here.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_cayley_sector_certificate as CAYLEY
import ou3_validated_transcendentals as VT

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
CORE = REPO / "src" / "kalman_ou_common" / "KalmanOUCoreMath.h"
MEKF = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"
SCHEMA = 1
SERIES_BRANCH_NORM = 1.0e-2
# Utility domain only.  It stays below pi so the correction Cayley vector is
# finite; this is not a claimed source correction bound.
CAYLEY_MONOTONE_NORM_MAX = 3.0


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def add_up(a: float, b: float) -> float:
    return up(float(a) + float(b))


def mul_up(a: float, b: float) -> float:
    return up(float(a) * float(b))


def div_up(a: float, b: float) -> float:
    if not b > 0.0:
        raise RuntimeError("positive denominator required")
    return up(float(a) / float(b))


def _series_cayley_scalar_bounds(delta: float) -> tuple[float, float, float]:
    """Bounds |a| and |a-d| for every 0<=|d|<=delta<=1e-2."""
    if not (0.0 <= delta <= SERIES_BRANCH_NORM):
        raise ValueError("series correction radius outside source branch")
    d2 = mul_up(delta, delta)
    d4 = mul_up(d2, d2)
    w_lo = down(1.0 - up(d2 / 8.0))
    w_hi = up(1.0 + up(d4 / 384.0))
    k_lo = down(0.5 - up(d2 / 48.0))
    k_hi = up(0.5 + up(d4 / 3840.0))
    if not (w_lo > 0.0 and k_lo > 0.0):
        raise RuntimeError("series Cayley coefficient lost positivity")
    coeff_lo = down((2.0 * k_lo) / w_hi)
    coeff_hi = up((2.0 * k_hi) / w_lo)
    a_hi = mul_up(coeff_hi, delta)
    coeff_err = max(abs(coeff_lo - 1.0), abs(coeff_hi - 1.0))
    adiff = mul_up(coeff_err, delta)
    return a_hi, adiff, coeff_hi


def _axis_cayley_scalar_bounds(delta: float) -> tuple[float, float]:
    """Bounds |a|=2 tan(|d|/2) and |a-d| on [0,delta], delta<pi."""
    if not (0.0 <= delta <= CAYLEY_MONOTONE_NORM_MAX):
        raise ValueError("axis-angle correction radius outside validated utility range")
    if delta == 0.0:
        return 0.0, 0.0
    half_lo = down(0.5 * delta)
    half_hi = up(0.5 * delta)
    if half_hi >= 0.5 * math.pi:
        raise RuntimeError("axis-angle correction reaches Cayley antipode")
    s_hi = VT.sin_point(half_hi).hi
    c_lo = VT.cos_point(half_hi).lo
    if not c_lo > 0.0:
        raise RuntimeError("validated correction quaternion scalar part is not positive")
    a_hi = div_up(mul_up(2.0, s_hi), c_lo)
    # 2 tan(x/2)-x is nonnegative and increasing on [0,pi).
    adiff = up(max(0.0, a_hi - half_lo*2.0))
    return a_hi, adiff


def correction_cayley_norm_bounds(delta: float) -> dict:
    delta = float(delta)
    if not (math.isfinite(delta) and 0.0 <= delta <= CAYLEY_MONOTONE_NORM_MAX):
        raise ValueError("correction radius must be finite in utility range")
    if delta < SERIES_BRANCH_NORM:
        a, diff, coeff = _series_cayley_scalar_bounds(delta)
        return {
            "correction_norm_upper": delta,
            "source_branch_family": "source_polynomial_series",
            "injected_cayley_norm_upper": a,
            "injected_cayley_minus_delta_norm_upper": diff,
            "series_cayley_coefficient_upper": coeff,
        }

    a_axis, diff_axis = _axis_cayley_scalar_bounds(delta)
    a_ser, diff_ser, _ = _series_cayley_scalar_bounds(SERIES_BRANCH_NORM)
    return {
        "correction_norm_upper": delta,
        "source_branch_family": (
            "threshold_hull" if delta == SERIES_BRANCH_NORM else "series_or_axis_angle"
        ),
        "injected_cayley_norm_upper": max(a_axis, a_ser),
        "injected_cayley_minus_delta_norm_upper": max(diff_axis, diff_ser),
        "series_cayley_coefficient_upper": None,
    }


def reset_defect_bound(q: float, delta: float) -> dict:
    q = float(q)
    delta = float(delta)
    if not (math.isfinite(q) and 0.0 <= q < 2.0):
        raise ValueError("state Cayley radius must be finite in [0,2)")
    corr = correction_cayley_norm_bounds(delta)
    a = float(corr["injected_cayley_norm_upper"])
    adiff = float(corr["injected_cayley_minus_delta_norm_upper"])
    denom = down(1.0 - mul_up(0.25, mul_up(a, q)))
    if not denom > 0.0:
        raise RuntimeError("correction/state Cayley composition can reach antipodal denominator")

    first = mul_up(adiff, add_up(1.0, 0.5*q))
    target = add_up(q, add_up(delta, mul_up(0.5*delta, q)))
    second = mul_up(mul_up(0.25*a, q), target)
    rho = div_up(add_up(first, second), denom)
    reset_max_sv = up(math.sqrt(up(1.0 + up(0.25*delta*delta))))
    return {
        **corr,
        "state_cayley_norm_upper": q,
        "cayley_composition_denominator_lower": denom,
        "reset_attitude_defect_norm_upper": rho,
        "reset_inverse_operator_norm_upper": 1.0,
        "reset_min_singular_value_lower": 1.0,
        "reset_max_singular_value_upper": reset_max_sv,
        "reset_determinant_lower": 1.0,
        "chart_safe": True,
    }


def _source_contract() -> list[str]:
    core = CORE.read_text(encoding="utf-8")
    mekf = MEKF.read_text(encoding="utf-8")
    checks = (
        (core, "if (theta < T(1e-2))", "quaternion series threshold"),
        (core, "q.normalize();", "correction quaternion normalization"),
        (core, "Eigen::Matrix<T,3,3>::Identity() + T(0.5)*skew(dtheta)", "left reset matrix"),
        (mekf, "xext.noalias() += K * r;", "state correction before reset"),
        (mekf, "joseph_update3_(K, S_mat, PCt);", "Joseph covariance update"),
        (mekf, "apply_error_state_reset_jacobian_(dtheta);", "left covariance reset"),
        (mekf, "qref = corr * qref;", "left quaternion injection"),
    )
    return [f"shipping reset semantic changed: {label}" for text, marker, label in checks if marker not in text]


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("exact reset transport must not be trajectory fitted")
    cayley = CAYLEY.build(path)
    failures = [f"Cayley: {x}" for x in CAYLEY.validate(cayley)] + _source_contract()
    q = float(cayley["cayley_radius_upper"])
    if not (0.0 < q < 1.0):
        failures.append("retained 0.8-rad state Cayley radius not inside q<1")

    # Sanity points exercise both source quaternion branches, but they are not
    # source correction assumptions and are not used for theorem promotion.
    sanity = []
    for delta in (0.0, 0.005, 0.01, 0.1):
        try:
            sanity.append(reset_defect_bound(q, delta))
        except Exception as exc:
            failures.append(f"reset arithmetic sanity delta={delta}: {exc}")

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_EXACT_JOSEPH_QUATERNION_RESET_TRANSPORT",
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "full_state_cross_terms_retained": True,
        "exact_joseph_information_identity": (
            "z^T P^-1 z-(z-Ky)^T(Pplus)^-1(z-Ky)=y^T S^-1 y-eta^T R^-1 eta"
        ),
        "exact_reset_congruence_identity": (
            "z_exact^T(Ge Pplus Ge^T)^-1 z_exact=(t+Ge^-1 rho)^T(Pplus)^-1(t+Ge^-1 rho)"
        ),
        "reset_defect_definition": "rho=z_exact-Ge*t; attitude component uses exact deployed Cayley composition",
        "reset_inverse_norm_proof": "G^T G eigenvalues are 1,1+||d||^2/4,1+||d||^2/4",
        "reset_inverse_operator_norm_upper": 1.0,
        "condition_number_multiplier_used_for_reset_transport": False,
        "parametric_correction_cayley_bound_available": True,
        "parametric_reset_defect_bound_available": True,
        "source_correlated_correction_norm_bound_supplied_here": False,
        "outer_angle_rad": float(cayley["outer_angle_rad"]),
        "outer_sector_covered": float(cayley["outer_angle_rad"]) >= 0.80,
        "arithmetic_sanity_points_not_theorem_bounds": sanity,
        "complete_H18_A21_word_established_here": False,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "next_obligation": (
            "obtain each accepted correction radius from the same source-correlated P,H,R,residual cell, apply the parametric rho bound, and accumulate its transported cross term with the signed Joseph word before endpoint scalarization"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P4_EXACT_JOSEPH_QUATERNION_RESET_TRANSPORT":
        f.append("wrong qualification")
    for key in (
        "source_generated_not_trajectory_fit", "full_state_cross_terms_retained",
        "parametric_correction_cayley_bound_available", "parametric_reset_defect_bound_available",
        "outer_sector_covered",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_changed",
        "condition_number_multiplier_used_for_reset_transport",
        "source_correlated_correction_norm_bound_supplied_here",
        "complete_H18_A21_word_established_here", "P4_USABLE_CERTIFICATE_PROMOTED",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("reset_inverse_operator_norm_upper") != 1.0:
        f.append("reset inverse norm is not exact one")
    sanity = d.get("arithmetic_sanity_points_not_theorem_bounds", [])
    if len(sanity) != 4:
        f.append("reset arithmetic sanity grid changed")
    for row in sanity:
        if row.get("chart_safe") is not True:
            f.append("reset arithmetic sanity point left Cayley chart")
        if not float(row.get("cayley_composition_denominator_lower", 0.0)) > 0.0:
            f.append("reset arithmetic sanity denominator not positive")
        if row.get("reset_inverse_operator_norm_upper") != 1.0:
            f.append("reset arithmetic sanity inverse norm changed")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(a.domain)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": "PASS" if not vf else "FAIL",
        "outer_angle_rad": d["outer_angle_rad"],
        "reset_inverse_norm": d["reset_inverse_operator_norm_upper"],
        "parametric_reset_bound": d["parametric_reset_defect_bound_available"],
        "correction_bound_supplied_here": d["source_correlated_correction_norm_bound_supplied_here"],
        "P4_promoted": d["P4_USABLE_CERTIFICATE_PROMOTED"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
