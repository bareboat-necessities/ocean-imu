#!/usr/bin/env python3
"""Exact Joseph/quaternion/reset transport primitive for OU-III P5.

P5 cannot transport a finite-angle correction by multiplying an independently
chosen gain norm, reset norm and metric condition number.  The shipping update
has much more structure.  Write the signed nonlinear measurement mismatch as

    y = H z + eta,

and choose the sign so that the implemented attitude injection is
``d=-E_theta K y``.  Before reset the tangent posterior is ``t=z-K y``.  For the
implemented Joseph posterior ``P+`` one has the exact information identity

    z' P^-1 z - t' (P+)^-1 t
      = y' S^-1 y - eta' R^-1 eta.                    (1)

The source then injects its normalized quaternion and applies the covariance
reset

    G(d) = I + 0.5 [d]x,
    P_r = G_ext P+ G_ext'.

Let ``z_exact`` use the *exact* post-injection Cayley coordinate and the actual
non-attitude posterior coordinates, and define

    rho = z_exact - G_ext t.

Because reset is a congruence, not a disturbance, its metric transport is
exact:

    z_exact' P_r^-1 z_exact
      = (t + G_ext^-1 rho)' (P+)^-1
        (t + G_ext^-1 rho).                            (2)

Furthermore ``G'G`` has eigenvalues ``1, 1+|d|^2/4, 1+|d|^2/4``.  Hence
``||G^-1||_2=1`` for every finite correction.  There is no reset condition
number penalty in the correction-word proof.

This producer also gives an outward finite bound for the *only* reset mismatch,
``rho``.  If ``a`` is the exact Cayley vector of the deployed correction
quaternion, then ``a`` is parallel to ``d`` and

    c+ = (a+c+0.5 a x c)/(1-0.25 a'c).

Against the homogeneous reset target ``G(d)(c+d)`` this gives, for
``|c|<=q, |d|<=delta`` and ``D=1-|a|q/4>0``,

    |rho_theta| <= [ |a-d|(1+q/2)
                    + |a|q/4 (q+delta+delta q/2) ] / D.       (3)

The Cayley magnitude of the source quaternion is enclosed from the source's
actual polynomial branch below 1e-2 and from validated sin/cos above it.  The
current first-due S correction is intentionally evaluated here on both gauged
P5 nodes.  This closes the algebra and chart denominator, but it does *not*
claim that the resulting coarse first-S defect budget proves a complete 1 s
word.  Later source subdivision/tighter staged covariance may reduce that
budget without changing this theorem route or the filter.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_group_algebra as GROUP
import ou3_p5_first_s_state_prefix_certificate as SPREFIX
import ou3_p5_heading_handoff_contract as HEADING
import ou3_validated_transcendentals as VT

REPO = Path(__file__).resolve().parents[1]
CORE = REPO / "src" / "kalman_ou_common" / "KalmanOUCoreMath.h"
MEKF = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


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
    """Bounds |a| and |a-d| for every 0<=|d|<=delta<1e-2."""
    if not (0.0 <= delta <= GROUP.SERIES_BRANCH_NORM):
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
    """Bounds |a|=2 tan(|d|/2) and |a-d| on [0,delta], delta<=3."""
    if not (0.0 <= delta <= GROUP.CAYLEY_MONOTONE_NORM_MAX):
        raise ValueError("axis-angle correction radius outside validated group range")
    if delta == 0.0:
        return 0.0, 0.0
    half = up(0.5 * delta)
    if half >= 0.5 * math.pi:
        raise RuntimeError("axis-angle correction reaches correction Cayley antipode")
    s_hi = VT.sin_point(half).hi
    c_lo = VT.cos_point(half).lo
    if not c_lo > 0.0:
        raise RuntimeError("validated correction quaternion scalar part is not positive")
    a_hi = div_up(mul_up(2.0, s_hi), c_lo)
    # 2 tan(x/2)-x is nonnegative and increasing on [0,pi).
    adiff = up(max(0.0, a_hi - delta))
    return a_hi, adiff


def correction_cayley_norm_bounds(delta: float) -> dict:
    """Source-faithful scalar enclosure for every correction norm <= delta."""
    delta = float(delta)
    if not (math.isfinite(delta) and 0.0 <= delta <= GROUP.CAYLEY_MONOTONE_NORM_MAX):
        raise ValueError("correction radius must be finite in [0,3]")
    if delta < GROUP.SERIES_BRANCH_NORM:
        a, diff, coeff = _series_cayley_scalar_bounds(delta)
        branch = "source_polynomial_series"
        return {
            "correction_norm_upper": delta,
            "source_branch_family": branch,
            "injected_cayley_norm_upper": a,
            "injected_cayley_minus_delta_norm_upper": diff,
            "series_cayley_coefficient_upper": coeff,
        }

    a_axis, diff_axis = _axis_cayley_scalar_bounds(delta)
    if delta == GROUP.SERIES_BRANCH_NORM:
        a_ser, diff_ser, _ = _series_cayley_scalar_bounds(delta)
    else:
        a_ser, diff_ser, _ = _series_cayley_scalar_bounds(GROUP.SERIES_BRANCH_NORM)
    return {
        "correction_norm_upper": delta,
        "source_branch_family": "series_or_axis_angle" if delta > GROUP.SERIES_BRANCH_NORM else "threshold_hull",
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
    dot_term = mul_up(0.25, mul_up(a, q))
    denom = down(1.0 - dot_term)
    if not denom > 0.0:
        raise RuntimeError("source correction/state Cayley composition can reach antipodal denominator")

    first = mul_up(adiff, add_up(1.0, 0.5 * q))
    target = add_up(q, add_up(delta, mul_up(0.5 * delta, q)))
    second = mul_up(mul_up(0.25 * a, q), target)
    rho = div_up(add_up(first, second), denom)
    return {
        **corr,
        "state_cayley_norm_upper": q,
        "cayley_composition_denominator_lower": denom,
        "reset_attitude_defect_norm_upper": rho,
        "reset_inverse_operator_norm_upper": 1.0,
        "reset_min_singular_value_lower": 1.0,
        "reset_max_singular_value_upper": up(math.sqrt(up(1.0 + up(0.25 * delta * delta)))),
        "reset_determinant_lower": 1.0,
        "chart_safe": denom > 0.0,
    }


def _source_markers() -> list[str]:
    core = CORE.read_text(encoding="utf-8")
    mekf = MEKF.read_text(encoding="utf-8")
    required = {
        "source_quaternion_series_threshold": "if (theta < T(1e-2))",
        "source_quaternion_normalized": "q.normalize();",
        "source_reset_helper": "apply_left_error_reset",
        "source_reset_matrix": "Eigen::Matrix<T,3,3>::Identity() + T(0.5)*skew(dtheta)",
        "source_state_then_joseph": "xext.noalias() += K * r;",
        "source_joseph_before_injection": "joseph_update3_(K, S_mat, PCt);",
        "source_immediate_quaternion_correction": "applyQuaternionCorrectionFromErrorState();",
    }
    joined = core + "\n" + mekf
    missing = [name for name, marker in required.items() if marker not in joined]
    return missing


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain_path = Path(domain_path).resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("exact correction transport domain must not be trajectory fitted")

    sp = SPREFIX.build(domain_path)
    sf = SPREFIX.validate(sp)
    heading = HEADING.build(domain_path)
    hf = HEADING.validate(heading)
    failures = [f"first-S-prefix: {x}" for x in sf] + [f"heading: {x}" for x in hf]
    failures += [f"missing source semantic: {x}" for x in _source_markers()]

    delta_s = float(sp["first_due_S_induced_attitude_correction_norm_upper_rad"])
    helper = float(sp["deployed_group_helper_correction_limit_rad"])
    if not (0.0 < delta_s < helper <= GROUP.CAYLEY_MONOTONE_NORM_MAX):
        failures.append("first-S correction is not inside validated deployed group helper")

    nodes = {}
    for name, row in (
        ("normal_gauged", heading["gauged_quality_handoff"]),
        ("timeout_gauged", heading["gauged_timeout_subbranch"]),
    ):
        q = float(row["full_attitude_cayley_norm_upper"])
        try:
            nodes[name] = reset_defect_bound(q, delta_s)
        except Exception as exc:
            failures.append(f"{name}: {exc}")

    for name, row in nodes.items():
        if row.get("chart_safe") is not True:
            failures.append(f"{name}: first-S correction does not preserve Cayley chart")
        if not float(row.get("cayley_composition_denominator_lower", -1.0)) > 0.0:
            failures.append(f"{name}: Cayley composition denominator not strict")
        if not math.isclose(float(row.get("reset_inverse_operator_norm_upper", math.inf)), 1.0):
            failures.append(f"{name}: exact reset inverse norm not one")

    passed = not failures
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_EXACT_JOSEPH_QUATERNION_RESET_TRANSPORT",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "full_state_cross_terms_retained": True,
        "exact_joseph_information_identity": (
            "z^T P^-1 z-(z-Ky)^T(Pplus)^-1(z-Ky)=y^T S^-1 y-eta^T R^-1 eta"
        ),
        "exact_reset_congruence_identity": (
            "z_exact^T(Ge Pplus Ge^T)^-1 z_exact=(t+Ge^-1 rho)^T(Pplus)^-1(t+Ge^-1 rho)"
        ),
        "reset_defect_definition": "rho=z_exact-Ge*t; attitude component uses exact deployed Cayley composition",
        "reset_inverse_norm_proof": "G^T G eigenvalues are 1, 1+||d||^2/4, 1+||d||^2/4",
        "condition_number_multiplier_used_for_reset_transport": False,
        "first_due_S_correction_norm_upper_rad": delta_s,
        "nodes": nodes,
        "complete_word_numerical_budget_closed_here": False,
        "P5_EXACT_CORRECTION_TRANSPORT_ALGEBRA_CERTIFICATE": "PASS" if passed else "FAIL",
        "next_obligation": (
            "use identities (1)-(3) at every accepted S/accelerometer/magnetometer prefix and prediction in the complete source-correlated 1 s word; tighten source-staged correction boxes/subdivide when the coarse reset defect exceeds the word information budget"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    for flag in ("source_generated_not_trajectory_fit", "full_state_cross_terms_retained"):
        if d.get(flag) is not True:
            failures.append(f"{flag} is not true")
    if d.get("source_replay_used") is not False:
        failures.append("exact correction transport uses replay")
    if d.get("filter_changed") is not False:
        failures.append("exact correction transport changes filter")
    if d.get("condition_number_multiplier_used_for_reset_transport") is not False:
        failures.append("reset transport reintroduced a condition-number multiplier")
    if d.get("complete_word_numerical_budget_closed_here") is not False:
        failures.append("correction algebra prematurely claims complete word")
    for name in ("normal_gauged", "timeout_gauged"):
        row = d.get("nodes", {}).get(name, {})
        if row.get("chart_safe") is not True:
            failures.append(f"{name}: chart is not safe")
        if not float(row.get("cayley_composition_denominator_lower", -1.0)) > 0.0:
            failures.append(f"{name}: denominator is not positive")
        if float(row.get("reset_inverse_operator_norm_upper", math.inf)) != 1.0:
            failures.append(f"{name}: reset inverse norm is not exact one")
        if not (isinstance(row.get("reset_attitude_defect_norm_upper"), (int, float))
                and math.isfinite(float(row["reset_attitude_defect_norm_upper"]))
                and float(row["reset_attitude_defect_norm_upper"]) >= 0.0):
            failures.append(f"{name}: reset defect bound invalid")
    if not failures and d.get("P5_EXACT_CORRECTION_TRANSPORT_ALGEBRA_CERTIFICATE") != "PASS":
        failures.append("exact correction transport algebra did not pass")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build(args.domain.resolve())
    vf = validate(out)
    out["validation_pass"] = not vf
    out["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": out["P5_EXACT_CORRECTION_TRANSPORT_ALGEBRA_CERTIFICATE"],
        "first_S_delta": out["first_due_S_correction_norm_upper_rad"],
        "nodes": out["nodes"],
        "next": out["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
