#!/usr/bin/env python3
"""Dependency-preserving signed first-accelerometer check for 30 deg P4.

This stage repairs three conservative losses in the first signed finite-angle
accelerometer enclosure without changing the filter, the declared 30 deg P4
candidate, the 0.3 g startup a_w set, or the 3 rad correction-Cayley chart.

1. Ratios X/(X+lambda) are bounded as one positive rational expression instead
   of interval-dividing a numerator by a denominator that contains the same X.
2. The ideal K_theta H_theta correction is charged only against the force-
   tangent attitude component.  The accelerometer is not credited with axial
   yaw contraction; the source-varying axial cross gain is retained explicitly.
3. The latent acceleration and its finite-rotation cross term are combined by

       e_aw + (R^T-I)e_aw = R^T e_aw,

   so their joint norm is exactly ||e_aw|| rather than two adverse independent
   disturbances.  The separate finite-angle force-vector remainder and accel
   bias are still retained.

The output remains a first-operation sector lemma.  It does not yet propagate
the Joseph posterior/reset covariance or promote P4.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p4_30deg_signed_first_accel_sector as V1
import ou3_p4_candidate_aw_capture_budget as AWB
import ou3_p4_candidate_first_accel_range_v3 as RANGE
import ou3_p4_p5_entrance_search_domain as ENTRANCE
import ou3_p4_operation_matched_sector_certificate as SECTOR
import ou3_p5_effective_vector_input as VEFF
import ou3_p5_first_accel_rotation_gauge as RG
import ou3_p5_first_accel_rotation_gauge_v3 as RG3
import ou3_p5_first_accel_structured_gain as SG
import ou3_p5_full_h_prefix_cells as FULL
import ou3_p5_full_h_prefix_cells_v3 as FULL3
import ou3_vector_uco_certificate as VECTOR

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 2
CANDIDATE_INDEX = 0
CORRECTION_CAYLEY_MONOTONE_LIMIT_RAD = V1.CORRECTION_CAYLEY_MONOTONE_LIMIT_RAD


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def _sqrt_up(x: float) -> float:
    return up(math.sqrt(max(0.0, float(x))))


def _positive_ratio_same_x(x: Interval, lam: Interval) -> Interval:
    """Outward range of x/(x+lam), x>=0, lam>0, preserving shared x."""
    if x.lo < 0.0 or lam.lo <= 0.0:
        raise RuntimeError("positive-ratio domain lost positivity")
    lo_den = up(x.lo + lam.hi)
    hi_den = down(x.hi + lam.lo)
    lo = 0.0 if x.lo == 0.0 else down(x.lo / lo_den)
    hi = up(x.hi / hi_den)
    return Interval.outward_bounds(max(0.0, lo), min(1.0, hi))


def _m2_over_m2p_plus_lam(m2: Interval, p: Interval, lam: Interval) -> Interval:
    """Outward range M/(M p+lam), monotone in M,p,lam on positive domain."""
    if m2.lo < 0.0 or p.lo <= 0.0 or lam.lo <= 0.0:
        raise RuntimeError("canonical positive resolvent domain lost positivity")
    dlo = down(down(m2.lo * p.hi) + lam.hi)
    dhi = up(up(m2.hi * p.lo) + lam.lo)
    lo = 0.0 if m2.lo == 0.0 else down(m2.lo / up(down(m2.lo * p.hi) + lam.hi))
    hi = up(m2.hi / down(up(m2.hi * p.lo) + lam.lo))
    if dlo <= 0.0 or dhi <= 0.0:
        raise RuntimeError("canonical resolvent denominator lost positivity")
    return Interval.outward_bounds(max(0.0, lo), hi)


def _canonical_kh_bounds(*, tilt: float, yaw: float, x: Interval, m: Interval,
                         paw: Interval, racc_var: Interval, gain_detail: dict) -> dict:
    t = Interval.outward_bounds(tilt, tilt)
    delta = Interval.outward_bounds(yaw - tilt, yaw - tilt)
    pu = t + delta * x
    m2 = m.square()
    lam = paw + racc_var
    if t.lo <= 0.0 or pu.lo <= 0.0 or lam.lo <= 0.0:
        raise RuntimeError("canonical innovation domain lost positivity")

    xt = m2 * t
    xu = m2 * pu
    a = _positive_ratio_same_x(xu, lam)
    b = _positive_ratio_same_x(xt, lam)

    geom_hi = V1._sqrt_x1mx_upper(x)
    resolv = _m2_over_m2p_plus_lam(m2, pu, lam)
    delta_abs = max(abs(delta.lo), abs(delta.hi))
    z_hi = up(up(delta_abs * geom_hi) * resolv.hi)

    dk = float(gain_detail["PSD_remainder_K_perturbation_upper"])
    dkh = up(dk * m.hi)
    kh_tangent = _sqrt_up(max(
        up(a.hi * a.hi + z_hi * z_hi),
        up(b.hi * b.hi),
    ))
    return {
        "tangent_quadratic_lower": min(a.lo, b.lo),
        "axial_cross_upper": z_hi,
        "ideal_KH_on_force_tangent_norm_upper": kh_tangent,
        "KH_remainder_operator_upper": dkh,
        "tangent_K_n_gain_upper": up(max(
            float(gain_detail["g_perp_upper"]),
            float(gain_detail["g_u_upper"]),
        ) + dk),
        "axial_K_n_gain_upper": up(float(gain_detail["g_z_upper"]) + dk),
        "a_interval": a.as_list(),
        "b_interval": b.as_list(),
        "z_interval": [0.0, z_hi],
        "dependency_preserving_positive_ratios": True,
    }


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 2,
          alignment_pieces: int = 16, force_magnitude_pieces: int = 4,
          tangent_pieces: int = 32) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("signed P4 candidate domain must not be trajectory fitted")

    RG3._install_backend(path, source_pieces)
    FULL3._install_backend()
    entrance = ENTRANCE.build(path)
    sector = SECTOR.build(path)
    vector = VECTOR.build()
    failures = [f"entrance: {x}" for x in ENTRANCE.validate(entrance)]
    failures += [f"sector: {x}" for x in SECTOR.validate(sector)]
    failures += [f"vector: {x}" for x in VECTOR.validate(vector)]

    crow = entrance["P4_complete_word_search"]["candidate_rows"][CANDIDATE_INDEX]
    if float(crow["angle_deg"]) != 30.0:
        failures.append("candidate index 0 is no longer 30 deg")
    q0 = float(crow["cayley_norm_upper"])
    h = float(FULL._source_cell()["dt_s"])
    q = RG._q_after_first_prediction(q0, domain, h)
    outer_q = float(sector["design_cayley_norm_upper"])
    tilt, yaw, eps = RG._attitude_covariance_epsilon(path, h)
    vc = vector["configured_measurement_bounds"]
    racc_var = FULL._R_diag(float(vc["acc_measurement_std_mps2"]))[0][0]
    startup = domain["startup"]
    handoff = startup["physical_handoff_coordinate_bounds"]
    aw0 = float(handoff["latent_acceleration_error_norm_upper_mps2"])
    ba = float(handoff["accelerometer_bias_error_norm_upper_mps2"])
    pnorm = AWB._p5_position_norm_upper(domain)

    src_phases = RG._source_phase_children(source_pieces)
    xcells = SG._linear_cells(alignment_pieces)
    live = domain["normal_live"]
    force_cells = RG._geom_ranges(
        float(live["specific_force_norm_lower_mps2"]),
        float(live["specific_force_norm_upper_mps2"]),
        force_magnitude_pieces,
    )
    tcells = V1._tangent_cells(q, tangent_pieces)

    worst = None
    max_qplus = q
    max_d = 0.0
    min_den = math.inf
    evaluated = 0
    first_failure = None

    for si, (src, phase) in enumerate(src_phases):
        P0 = FULL._initial_covariance(src, path)
        F, Q, _ = FULL._transition_and_Q(src, domain)
        Pp = FULL._psd_tighten(FULL.matrix_add(
            FULL.matrix_mul(FULL.matrix_mul(F, P0), FULL.matrix_transpose(F)), Q))
        _pss, _psa, paw_pred = RG._scalar_axis_structure(Pp)
        paw = RG._due_paw_and_error_norm(Pp, src, 0.0, 0.0)[0] if phase == "due" else paw_pred
        intercept, slope, _sdetail = AWB._s_phase_affine_aw_bound(src, phase, Pp, domain, pnorm)
        aw = up(intercept + up(slope * aw0))

        for xi, x in enumerate(xcells):
            for mi, m in enumerate(force_cells):
                k, _kh_generic, gdetail = RANGE._tangent_structured_gain_bounds(
                    tilt=tilt, yaw=yaw, eps=eps, x=x, m=m,
                    paw=paw, racc_var=racc_var)
                canon = _canonical_kh_bounds(
                    tilt=tilt, yaw=yaw, x=x, m=m, paw=paw,
                    racc_var=racc_var, gain_detail=gdetail)

                # Exact latent rotation combination:
                # e_aw + (R^T-I)e_aw = R^T e_aw, hence norm exactly ||e_aw||.
                eta_force = up(VEFF.accel_attitude_eta_per_vector_norm_upper(q) * m.hi)
                rho = up(aw + up(eta_force + ba))

                for ti, (rlo, rhi) in enumerate(tcells):
                    s_hi = V1._axial_upper(q, rlo)
                    quad = down(
                        down(canon["tangent_quadratic_lower"] * down(rlo * rlo))
                        - up(canon["axial_cross_upper"] * up(rhi * s_hi))
                        - up(canon["KH_remainder_operator_upper"] * up(q * q))
                    )
                    nuisance_dot = up(rho * up(
                        up(canon["tangent_K_n_gain_upper"] * rhi)
                        + up(canon["axial_K_n_gain_upper"] * s_hi)
                    ))
                    ddot = up(-quad + nuisance_dot)

                    ideal_khc = up(canon["ideal_KH_on_force_tangent_norm_upper"] * rhi)
                    remainder_khc = up(canon["KH_remainder_operator_upper"] * q)
                    nuisance_corr = up(k * rho)
                    dnorm = up(up(ideal_khc + remainder_khc) + nuisance_corr)
                    try:
                        comp = V1._compose_q_upper(q, dnorm, ddot)
                    except Exception as exc:
                        first_failure = {
                            "source_phase_cell": si,
                            "pseudo_phase": phase,
                            "alignment_cell": xi,
                            "force_cell": mi,
                            "tangent_cell": ti,
                            "q": q,
                            "r_interval": [rlo, rhi],
                            "s_upper": s_hi,
                            "aw_after_prefix_upper_mps2": aw,
                            "force_attitude_remainder_upper_mps2": eta_force,
                            "nuisance_effective_aw_upper_mps2": rho,
                            "Ktheta_norm_upper": k,
                            "ideal_KHc_norm_upper_rad": ideal_khc,
                            "PSD_remainder_KHc_norm_upper_rad": remainder_khc,
                            "nuisance_correction_norm_upper_rad": nuisance_corr,
                            "correction_norm_upper_rad": dnorm,
                            "d_dot_c_upper": ddot,
                            "canonical": canon,
                            "reason": f"{type(exc).__name__}: {exc}",
                        }
                        break
                    evaluated += 1
                    qp = float(comp["post_update_q_upper"])
                    max_qplus = max(max_qplus, qp)
                    max_d = max(max_d, dnorm)
                    min_den = min(min_den, float(comp["denominator_lower"]))
                    row = {
                        "source_phase_cell": si,
                        "pseudo_phase": phase,
                        "alignment_cell": xi,
                        "force_cell": mi,
                        "tangent_cell": ti,
                        "force_magnitude_mps2": m.as_list(),
                        "alignment_x": x.as_list(),
                        "r_interval": [rlo, rhi],
                        "s_upper": s_hi,
                        "aw_after_prefix_upper_mps2": aw,
                        "force_attitude_remainder_upper_mps2": eta_force,
                        "nuisance_effective_aw_upper_mps2": rho,
                        "Ktheta_norm_upper": k,
                        "ideal_KHc_norm_upper_rad": ideal_khc,
                        "PSD_remainder_KHc_norm_upper_rad": remainder_khc,
                        "nuisance_correction_norm_upper_rad": nuisance_corr,
                        "correction_norm_upper_rad": dnorm,
                        "d_dot_c_upper": ddot,
                        "canonical": canon,
                        **comp,
                    }
                    if worst is None or qp > float(worst["post_update_q_upper"]):
                        worst = row
                if first_failure is not None:
                    break
            if first_failure is not None:
                break
        if first_failure is not None:
            break

    complete = first_failure is None and evaluated > 0
    inside_outer = complete and max_qplus < outer_q
    below_correction_chart = complete and max_d < CORRECTION_CAYLEY_MONOTONE_LIMIT_RAD
    if not complete:
        failures.append("dependency-preserving signed 30deg first-accelerometer family did not complete")
    if complete and not below_correction_chart:
        failures.append("signed 30deg correction family reached the 3 rad Cayley monotone limit")
    if complete and not inside_outer:
        failures.append("signed 30deg first-accelerometer family leaves operation-matched outer sector")
    passed = not failures

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_30DEG_SIGNED_FIRST_ACCEL_DEPENDENCY_PRESERVING",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "candidate_angle_deg": 30.0,
        "candidate_q_upper": q0,
        "post_prediction_q_upper": q,
        "operation_matched_outer_q_upper": outer_q,
        "declared_startup_aw_error_norm_upper_mps2": aw0,
        "declared_startup_aw_error_fraction_g": float(startup["latent_acceleration_error_fraction_g"]),
        "rotation_gauge_used": True,
        "accelerometer_yaw_contraction_assumed": False,
        "axial_yaw_cross_gain_retained": True,
        "PSD_gain_remainder_retained": True,
        "dependency_preserving_positive_ratios": True,
        "ideal_KH_charged_on_force_tangent_only": True,
        "latent_aw_rotation_combined_exactly": True,
        "effective_aw_finite_angle_force_remainder_retained": True,
        "accepted_and_rejected_branches_covered": True,
        "correction_cayley_monotone_limit_rad": CORRECTION_CAYLEY_MONOTONE_LIMIT_RAD,
        "evaluated_children": evaluated,
        "max_correction_norm_upper_rad": max_d,
        "minimum_signed_composition_denominator_lower": min_den,
        "max_accepted_or_rejected_post_update_q_upper": max_qplus,
        "all_children_below_correction_cayley_monotone_limit": below_correction_chart,
        "all_children_inside_operation_matched_outer_sector": inside_outer,
        "worst_child": worst,
        "first_failure": first_failure,
        "P4_30DEG_SIGNED_FIRST_ACCEL_SECTOR_CERTIFICATE": "PASS" if passed else "NOT_ESTABLISHED",
        "JOSEPH_RESET_COVARIANCE_PROPAGATED_HERE": False,
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE": False,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "next_obligation": (
            "if this dependency-preserving signed sector closes, propagate the same source/alignment/force/tangent children through the exact Joseph posterior and reset congruence before sample1; otherwise subdivide the first reported directional child without shrinking the declared 30deg/0.3g entrance"
        ),
        "failures": list(dict.fromkeys(failures)),
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for k in (
        "source_generated_not_trajectory_fit", "rotation_gauge_used",
        "axial_yaw_cross_gain_retained", "PSD_gain_remainder_retained",
        "dependency_preserving_positive_ratios", "ideal_KH_charged_on_force_tangent_only",
        "latent_aw_rotation_combined_exactly", "effective_aw_finite_angle_force_remainder_retained",
        "accepted_and_rejected_branches_covered",
    ):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in (
        "source_replay_used", "filter_changed", "accelerometer_yaw_contraction_assumed",
        "JOSEPH_RESET_COVARIANCE_PROPAGATED_HERE", "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE",
        "P4_USABLE_CERTIFICATE_PROMOTED",
    ):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    if float(d.get("candidate_angle_deg", -1.0)) != 30.0:
        f.append("candidate angle changed")
    if float(d.get("declared_startup_aw_error_fraction_g", -1.0)) != 0.3:
        f.append("startup a_w domain changed")
    if d.get("P4_30DEG_SIGNED_FIRST_ACCEL_SECTOR_CERTIFICATE") == "PASS":
        if d.get("all_children_below_correction_cayley_monotone_limit") is not True:
            f.append("PASS without correction-Cayley chart closure")
        if d.get("all_children_inside_operation_matched_outer_sector") is not True:
            f.append("PASS without outer-sector closure")
        if d.get("first_failure") is not None:
            f.append("PASS retains first failure")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces", type=int, default=2)
    ap.add_argument("--alignment-pieces", type=int, default=16)
    ap.add_argument("--force-magnitude-pieces", type=int, default=4)
    ap.add_argument("--tangent-pieces", type=int, default=32)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(a.domain.resolve(), source_pieces=a.source_pieces,
              alignment_pieces=a.alignment_pieces,
              force_magnitude_pieces=a.force_magnitude_pieces,
              tangent_pieces=a.tangent_pieces)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P4_30DEG_SIGNED_FIRST_ACCEL_SECTOR_CERTIFICATE"],
        "q_pre": d["post_prediction_q_upper"],
        "q_outer": d["operation_matched_outer_q_upper"],
        "max_d": d["max_correction_norm_upper_rad"],
        "max_qplus": d["max_accepted_or_rejected_post_update_q_upper"],
        "min_den": d["minimum_signed_composition_denominator_lower"],
        "evaluated": d["evaluated_children"],
        "first_failure": d["first_failure"],
        "validation_failures": vf,
        "next": d["next_obligation"],
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
