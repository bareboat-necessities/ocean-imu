#!/usr/bin/env python3
"""Source-correlated signed first-accelerometer bound for the 45 deg P5 entrance.

The preceding q<8 bridge deliberately forgets the direction of the first
accelerometer correction and permits every vector with ||d||<=dmax.  That is a
safe chart certificate, but it turns the actual gravity correction into an
adversarial rotation and expands q from about 0.83 to about 7.67.

At the deployed first Live packet the source is much more structured:

* the goLive yaw-covariance axis is the body image of world down;
* the first predicted specific-force axis is that same world-down axis;
* J_aw is orthogonal and the first theta/a_w covariance is exactly zero;
* after rotation gauge, H_theta acts only on the two gravity-tangent axes;
* the ideal K_theta H_theta is therefore a positive scalar on the tangent
  plane and exactly zero on yaw.

Write c=(c_t,c_y), t=||c_t||, q=||c||.  The physical correction is

    d = -K_theta H_theta c_t - K_theta n,

where n contains the rotated latent-acceleration error, accelerometer-bias
error, and the exact finite-angle attitude residual remainder.  The latter is
bounded by

    ||eta_R|| <= q/sqrt(4+q^2) * g * t,

which is the exact Cayley vector identity specialized to gravity; importantly
it vanishes for pure yaw.  The small one-step PSD covariance remainder is kept
as a resolvent perturbation of K, not discarded.

For each source/phase cell and a validated scalar subdivision of t in [0,q],
this producer bounds both ||d|| and the signed dot d^T c.  The deployed
correction Cayley vector is a=s(||d||)d with the exact positive source scale,
so a^T c keeps the same sign information.  The exact Cayley composition obeys

  ||c+||^2 = (||a||^2+||c||^2+2 a^T c
              + 1/4 (||a||^2||c||^2-(a^T c)^2))
             /(1-1/4 a^T c)^2.

Dropping only the favorable negative square term gives an outward upper bound.
No favorable correction direction is assumed, no replay sample is used, and no
filter/helper range is changed.  This is still a first-packet transient bound;
it does not by itself establish 30 deg recapture or the complete P4 word.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p4_candidate_first_accel_exact_source as FIRST
import ou3_p4_p5_entrance_search_domain as ENTRANCE
import ou3_p5_45deg_first_accel_q8_bridge as Q8
import ou3_p5_first_accel_rotation_gauge as RG
import ou3_p5_first_accel_structured_gain as SG
import ou3_p5_full_h_prefix_cells as FULL
import ou3_p5_signed_cayley_cell as SIGNED
import ou3_vector_uco_certificate as VECTOR

DEFAULT_DOMAIN = FIRST.DEFAULT_DOMAIN
SCHEMA = 1
DEFAULT_TANGENT_CELLS = 96


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def _tangent_cells(q: float, n: int) -> list[tuple[float, float]]:
    if not (math.isfinite(q) and q > 0.0 and n >= 1):
        raise ValueError("positive q and tangent cell count required")
    edges = [q * i / n for i in range(n + 1)]
    return [
        (0.0 if i == 0 else down(edges[i]), up(edges[i + 1]))
        for i in range(n)
    ]


def _sqrt_lower(x: float) -> float:
    if not x > 0.0:
        raise ValueError("positive square-root input required")
    return down(math.sqrt(x))


def _accepted_q_upper(q: float, d_norm_upper: float, d_dot_c_upper: float) -> dict:
    if not (0.0 <= d_norm_upper < math.pi):
        raise RuntimeError("signed first correction left the monotone pre-pi range")
    scale = SIGNED.correction_cayley_scale_interval(d_norm_upper)
    a_norm = up(scale.hi * d_norm_upper)
    if d_dot_c_upper >= 0.0:
        adot = up(scale.hi * d_dot_c_upper)
    else:
        # Multiplying a negative upper bound by the smallest positive scale is
        # the conservative (least negative) upper bound.
        adot = up(scale.lo * d_dot_c_upper)

    den = down(1.0 - up(0.25 * adot))
    if not den > 0.0:
        raise RuntimeError("signed first-accelerometer Cayley denominator can cross zero")

    a2 = up(a_norm * a_norm)
    q2 = up(q * q)
    # Exact numerator has an additional -0.25*(a^T c)^2 <= 0, which is dropped.
    num = up(up(a2 + q2) + up(2.0 * adot) + up(0.25 * up(a2 * q2)))
    num = max(0.0, num)
    q2plus = up(num / down(den * den))
    qplus = up(math.sqrt(q2plus))
    return {
        "correction_cayley_scale": scale.as_list(),
        "correction_cayley_norm_upper": a_norm,
        "a_dot_c_upper": adot,
        "composition_denominator_lower": den,
        "post_update_q_upper": qplus,
    }


def _row_bound(
    row: dict,
    *,
    q: float,
    gravity: float,
    ba: float,
    tilt: float,
    yaw: float,
    eps: float,
    racc_var: Interval,
    tangent_cells: list[tuple[float, float]],
) -> dict:
    paw = Interval.outward_bounds(*map(float, row["P_aw_variance_interval"]))
    m = Interval.outward_bounds(gravity, gravity)
    x0 = Interval.outward_bounds(0.0, 0.0)
    k, _kh, detail = SG._structured_gain_bounds(
        tilt=tilt,
        yaw=yaw,
        eps=eps,
        x=x0,
        m=m,
        paw=paw,
        racc_var=racc_var,
    )

    # E=0 tangent channel.  With x=0 both tangent directions have this same
    # gain and there is exactly no ideal yaw injection.
    lam_hi = up(paw.hi + racc_var.hi)
    m2_lo = down(gravity * gravity)
    m2_hi = up(gravity * gravity)
    kh0_lo = down(down(m2_lo * tilt) / up(up(m2_hi * tilt) + lam_hi))
    kh0_hi = float(detail["KH0_norm_upper"])
    k0_hi = float(detail["K0_norm_upper"])
    dk = float(detail["PSD_remainder_K_perturbation_upper"])
    dkh = up(dk * gravity)
    aw = float(row["predicted_aw_error_norm_upper_mps2"])

    den_eta = _sqrt_lower(down(4.0 + down(q * q)))
    eta_per_t = up(up(q * gravity) / den_eta)

    worst = None
    max_qplus = 0.0
    max_d = 0.0
    max_adot = -math.inf
    min_den = math.inf
    for ci, (tlo, thi) in enumerate(tangent_cells):
        # Exact nuisance decomposition:
        # R^T e_aw has norm ||e_aw||; bias is additive; eta_R is proportional
        # to the gravity-tangent Cayley component.
        rho = up(aw + ba + up(eta_per_t * thi))

        # Ideal KH is tangent only.  The source-certified PSD remainder can
        # perturb both the tangent gain and yaw row, so its action pays q.
        corrective_norm = up(up(kh0_hi * thi) + up(dkh * q))
        nuisance_norm = up(k * rho)
        dnorm = up(corrective_norm + nuisance_norm)

        # d=-KH c-K n.  Lower c^T KH c by the ideal tangent quadratic minus
        # the perturbation norm.  For K n, the ideal gain is tangent (pay t),
        # while only the tiny resolvent perturbation can pay the full q.
        quad_lower = down(down(kh0_lo * down(tlo * tlo)) - up(dkh * up(q * q)))
        nuisance_dot = up(up(k0_hi * rho * thi) + up(dk * rho * q))
        ddot_upper = up(-quad_lower + nuisance_dot)

        comp = _accepted_q_upper(q, dnorm, ddot_upper)
        qplus = float(comp["post_update_q_upper"])
        max_qplus = max(max_qplus, qplus)
        max_d = max(max_d, dnorm)
        max_adot = max(max_adot, float(comp["a_dot_c_upper"]))
        min_den = min(min_den, float(comp["composition_denominator_lower"]))
        if worst is None or qplus > worst["post_update_q_upper"]:
            worst = {
                "tangent_cell": ci,
                "tangent_norm_interval": [tlo, thi],
                "nuisance_residual_norm_upper_mps2": rho,
                "correction_norm_upper_rad": dnorm,
                "d_dot_c_upper": ddot_upper,
                **comp,
            }

    return {
        "source_phase_cell": int(row["source_phase_cell"]),
        "pseudo_phase": row["pseudo_phase"],
        "tau_s": row["tau_s"],
        "sigma_aw_mps2": row["sigma_aw_mps2"],
        "R_S_filter_std": row["R_S_filter_std"],
        "P_aw_variance_interval": row["P_aw_variance_interval"],
        "predicted_aw_error_norm_upper_mps2": aw,
        "Ktheta_norm_upper": k,
        "ideal_tangent_KH_lower": kh0_lo,
        "ideal_tangent_KH_upper": kh0_hi,
        "resolvent_K_perturbation_upper": dk,
        "resolvent_KH_perturbation_upper": dkh,
        "eta_attitude_per_tangent_cayley_norm_upper_mps2": eta_per_t,
        "max_correction_norm_upper_rad": max_d,
        "max_a_dot_c_upper": max_adot,
        "minimum_composition_denominator_lower": min_den,
        "accepted_post_update_q_upper": max_qplus,
        "accepted_or_rejected_post_update_q_upper": max(q, max_qplus),
        "worst_tangent_child": worst,
    }


def build(
    domain_path: Path = DEFAULT_DOMAIN,
    *,
    source_pieces: int = 2,
    tangent_cells: int = DEFAULT_TANGENT_CELLS,
) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("signed first-accelerometer source bound must not be trajectory fitted")
    if domain.get("configured_runtime", {}).get("imu_lever_arm_enabled") is not False:
        raise RuntimeError("signed first-accelerometer source bound requires lever arm disabled")

    first = FIRST.build(path, source_pieces=source_pieces)
    entrance = ENTRANCE.build(path)
    q8 = Q8.build(path, source_pieces=source_pieces)
    vector = VECTOR.build()
    failures = [f"first: {x}" for x in FIRST.validate(first)]
    failures += [f"entrance: {x}" for x in ENTRANCE.validate(entrance)]
    failures += [f"q8: {x}" for x in Q8.validate(q8)]
    failures += [f"vector: {x}" for x in VECTOR.validate(vector)]

    if first.get("first_accel_yaw_covariance_axis_aligned_with_force_axis") is not True:
        failures.append("first source does not pin yaw covariance axis to force axis")
    if first.get("yaw_alignment_x_equals_zero") is not True:
        failures.append("first source yaw-alignment scalar is not x=0")

    q = float(q8["P5_45deg_entrance_first_accel"]["post_prediction_q_upper"])
    scalar_qplus = float(q8["P5_45deg_entrance_first_accel"]["product"]["post_update_q_upper_from_scalar"])
    gravity = float(domain["startup"]["gravity_mps2"])
    ba = float(domain["startup"]["physical_handoff_coordinate_bounds"]["accelerometer_bias_error_norm_upper_mps2"])
    h = float(FULL._source_cell()["dt_s"])
    tilt, yaw, eps = RG._attitude_covariance_epsilon(path, h)
    vc = vector["configured_measurement_bounds"]
    Racc = FULL._R_diag(float(vc["acc_measurement_std_mps2"]))
    racc_var = Racc[0][0]
    tcells = _tangent_cells(q, int(tangent_cells))

    candidate_rows = first.get("candidate_rows", [])
    source_rows = candidate_rows[0].get("source_rows", []) if candidate_rows else []
    if not source_rows:
        failures.append("first-accelerometer source rows missing")

    rows = []
    first_bad = None
    for row in source_rows:
        try:
            r = _row_bound(
                row,
                q=q,
                gravity=gravity,
                ba=ba,
                tilt=tilt,
                yaw=yaw,
                eps=eps,
                racc_var=racc_var,
                tangent_cells=tcells,
            )
            rows.append(r)
        except Exception as exc:
            first_bad = {
                "source_phase_cell": row.get("source_phase_cell"),
                "pseudo_phase": row.get("pseudo_phase"),
                "reason": f"{type(exc).__name__}: {exc}",
            }
            break

    complete = len(rows) == len(source_rows) and first_bad is None and bool(rows)
    qplus = max((float(r["accepted_or_rejected_post_update_q_upper"]) for r in rows), default=math.inf)
    dmax = max((float(r["max_correction_norm_upper_rad"]) for r in rows), default=math.inf)
    min_den = min((float(r["minimum_composition_denominator_lower"]) for r in rows), default=-math.inf)
    improved = complete and math.isfinite(qplus) and qplus < scalar_qplus
    inside_q8 = complete and qplus < 8.0
    if not complete:
        failures.append("signed source/tangent family did not close")
    if not inside_q8:
        failures.append("signed source/tangent first-packet family is not inside q8")
    if not improved:
        failures.append("signed source/tangent bound did not improve the sign-agnostic scalar bridge")

    passed = not failures
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_45DEG_FIRST_ACCEL_SIGNED_SOURCE_CORRELATED_BOUND",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "deployed_correction_limit_increased": False,
        "first_accel_yaw_covariance_axis_aligned_with_force_axis": True,
        "ideal_first_accel_gain_tangent_only": True,
        "ideal_first_accel_yaw_injection_exact_zero": True,
        "PSD_remainder_resolvent_perturbation_retained": True,
        "finite_angle_attitude_eta_scales_with_tangent_error": True,
        "latent_rotation_plus_linear_error_norm_combined_exactly": True,
        "signed_a_dot_c_retained": True,
        "favorable_correction_direction_assumed": False,
        "exact_Cayley_cross_product_norm_identity_used": True,
        "negative_dot_square_term_dropped_only_for_upper_bound": True,
        "pre_update_q_upper": q,
        "tangent_subdivision_cell_count": len(tcells),
        "source_phase_cell_count": len(source_rows),
        "evaluated_source_phase_cells": len(rows),
        "sign_agnostic_scalar_post_update_q_upper": scalar_qplus,
        "signed_source_correlated_post_update_q_upper": qplus,
        "q_upper_improvement_factor": scalar_qplus / qplus if improved else 0.0,
        "max_signed_decomposition_correction_norm_upper_rad": dmax,
        "minimum_signed_composition_denominator_lower": min_den,
        "inside_q8": inside_q8,
        "strictly_improves_sign_agnostic_bridge": improved,
        "source_rows": rows,
        "first_unclosed_child": first_bad,
        "returned_to_30deg_P4_sector_here": qplus <= float(entrance["P4_complete_word_candidate_sectors"][0]["cayley_norm_upper"]),
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE": False,
        "P5_FINITE_CAPTURE_TO_P4_ESTABLISHED_HERE": False,
        "P5_45DEG_FIRST_ACCEL_SIGNED_SOURCE_BOUND": "PASS" if passed else "NOT_ESTABLISHED",
        "next_obligation": (
            "propagate this signed/source-correlated first-packet attitude child together with the physical H group-norm and Joseph/reset covariance children through sample1 S/accelerometer/magnetometer operations; test finite recapture into the 30deg P4 sector"
        ),
        "failures": list(dict.fromkeys(failures)),
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for k in (
        "source_generated_not_trajectory_fit",
        "first_accel_yaw_covariance_axis_aligned_with_force_axis",
        "ideal_first_accel_gain_tangent_only",
        "ideal_first_accel_yaw_injection_exact_zero",
        "PSD_remainder_resolvent_perturbation_retained",
        "finite_angle_attitude_eta_scales_with_tangent_error",
        "latent_rotation_plus_linear_error_norm_combined_exactly",
        "signed_a_dot_c_retained",
        "exact_Cayley_cross_product_norm_identity_used",
        "negative_dot_square_term_dropped_only_for_upper_bound",
        "inside_q8",
        "strictly_improves_sign_agnostic_bridge",
    ):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in (
        "source_replay_used",
        "filter_changed",
        "deployed_correction_limit_increased",
        "favorable_correction_direction_assumed",
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE",
        "P5_FINITE_CAPTURE_TO_P4_ESTABLISHED_HERE",
    ):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    if int(d.get("evaluated_source_phase_cells", 0)) != int(d.get("source_phase_cell_count", -1)):
        f.append("source phase family incomplete")
    if d.get("first_unclosed_child") is not None:
        f.append("signed source family retains an unclosed child")
    q0 = float(d.get("sign_agnostic_scalar_post_update_q_upper", math.inf))
    q1 = float(d.get("signed_source_correlated_post_update_q_upper", math.inf))
    if not (math.isfinite(q1) and 0.0 < q1 < q0 < 8.0):
        f.append("signed q bound is not a strict finite improvement inside q8")
    if not float(d.get("minimum_signed_composition_denominator_lower", -math.inf)) > 0.0:
        f.append("signed composition denominator not positive")
    st = d.get("P5_45DEG_FIRST_ACCEL_SIGNED_SOURCE_BOUND")
    if st not in ("PASS", "NOT_ESTABLISHED"):
        f.append("invalid status")
    if st == "PASS" and f:
        f.append("PASS carries validation failures")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces", type=int, default=2)
    ap.add_argument("--tangent-cells", type=int, default=DEFAULT_TANGENT_CELLS)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build(args.domain.resolve(), source_pieces=args.source_pieces, tangent_cells=args.tangent_cells)
    vf = validate(out)
    out["validation_pass"] = not vf
    out["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": out["P5_45DEG_FIRST_ACCEL_SIGNED_SOURCE_BOUND"],
        "q_pre": out["pre_update_q_upper"],
        "q_scalar": out["sign_agnostic_scalar_post_update_q_upper"],
        "q_signed": out["signed_source_correlated_post_update_q_upper"],
        "improvement_factor": out["q_upper_improvement_factor"],
        "max_d": out["max_signed_decomposition_correction_norm_upper_rad"],
        "min_den": out["minimum_signed_composition_denominator_lower"],
        "returned_to_30deg": out["returned_to_30deg_P4_sector_here"],
        "first_unclosed": out["first_unclosed_child"],
        "validation_failures": vf,
        "next": out["next_obligation"],
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
