#!/usr/bin/env python3
"""Signed source-complete first-accelerometer sector check for the 30 deg P4 candidate.

The scalar candidate range certificate proves the first normal-Live
accelerometer correction is below the shipping 6 rad helper limit on the
30-degree P4 entrance, but a norm-only correction cannot prove finite-angle
sector invariance.  This producer retains the correction/current-attitude sign
needed by exact Cayley composition.

In a rotation gauge with predicted specific force ``f=m e3`` and yaw covariance
axis in the x-z plane, the nominal attitude covariance is

    P = t I + delta v v' ,  v=(sqrt(x),0,sqrt(1-x)).

For isotropic effective accelerometer noise ``lambda=P_aw+R_acc``, the ideal
``K_theta H_theta`` map has the exact canonical form

    [ a  0  0 ]
    [ 0  b  0 ]
    [ z  0  0 ]

with positive a,b and source-varying axial coupling z.  Thus the accelerometer
is not falsely credited with yaw contraction.  For a force-tangent Cayley norm
r and axial magnitude s inside ``||c||<=q``, we lower-bound

    c' K H c >= min(a,b) r^2 - |z| r s - ||dKH|| q^2,

where the last term retains the already-certified small PSD remainder.  The
nuisance correction uses the exact effective-a_w representation, including the
finite-angle attitude remainder and latent-acceleration rotation cross term.
Its tangent and axial gain channels are charged separately in ``c' K n``.

The resulting signed ``d'c`` and correction norm are composed through the
shipping correction Cayley scale only when the whole child stays below the
existing monotone correction radius (3 rad < pi).  Accepted and rejected
branches are both covered.  This is a first-operation sector lemma only: it
does not propagate Joseph/reset covariance, later packets, the full word, or
promote P4.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p4_candidate_aw_capture_budget as AWB
import ou3_p4_candidate_first_accel_range_v3 as V3
import ou3_p4_p5_entrance_search_domain as ENTRANCE
import ou3_p4_operation_matched_sector_certificate as SECTOR
import ou3_p5_effective_vector_input as VEFF
import ou3_p5_first_accel_rotation_gauge as RG
import ou3_p5_first_accel_rotation_gauge_v3 as RG3
import ou3_p5_first_accel_structured_gain as SG
import ou3_p5_full_h_prefix_cells as FULL
import ou3_p5_full_h_prefix_cells_v3 as FULL3
import ou3_p5_signed_cayley_cell as SIGNED
import ou3_vector_uco_certificate as VECTOR

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
CANDIDATE_INDEX = 0
CORRECTION_CAYLEY_MONOTONE_LIMIT_RAD = 3.0


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def _sqrt_up(x: float) -> float:
    return up(math.sqrt(max(0.0, float(x))))


def _sqrt_x1mx_upper(x: Interval) -> float:
    lo = max(0.0, x.lo)
    hi = min(1.0, x.hi)
    vals = [lo * (1.0 - lo), hi * (1.0 - hi)]
    if lo <= 0.5 <= hi:
        vals.append(0.25)
    return _sqrt_up(max(vals))


def _tangent_cells(q: float, pieces: int):
    if pieces < 1:
        raise ValueError("positive tangent piece count required")
    edges = [q * i / pieces for i in range(pieces + 1)]
    return [(0.0 if i == 0 else down(edges[i]), up(edges[i + 1])) for i in range(pieces)]


def _axial_upper(q: float, tangent_lower: float) -> float:
    return _sqrt_up(max(0.0, up(q * q) - down(tangent_lower * tangent_lower)))


def _canonical_kh_bounds(*, tilt: float, yaw: float, x: Interval, m: Interval,
                         paw: Interval, racc_var: Interval, gain_detail: dict):
    t = Interval.outward_bounds(tilt, tilt)
    delta = Interval.outward_bounds(yaw - tilt, yaw - tilt)
    pu = t + delta * x
    m2 = m.square()
    lam = paw + racc_var
    den_t = m2 * t + lam
    den_u = m2 * pu + lam
    if den_t.lo <= 0.0 or den_u.lo <= 0.0:
        raise RuntimeError("canonical innovation denominator lost positivity")
    a = m2 * pu / den_u
    b = m2 * t / den_t
    geom = Interval(0.0, _sqrt_x1mx_upper(x))
    z = m2 * delta * geom / den_u
    dk = float(gain_detail["PSD_remainder_K_perturbation_upper"])
    dkh = up(dk * m.hi)
    return {
        "tangent_quadratic_lower": min(a.lo, b.lo),
        "axial_cross_upper": z.hi,
        "KH_remainder_operator_upper": dkh,
        "tangent_K_n_gain_upper": up(max(
            float(gain_detail["g_perp_upper"]),
            float(gain_detail["g_u_upper"]),
        ) + dk),
        "axial_K_n_gain_upper": up(float(gain_detail["g_z_upper"]) + dk),
        "a_interval": a.as_list(),
        "b_interval": b.as_list(),
        "z_interval": z.as_list(),
    }


def _compose_q_upper(q: float, dnorm: float, ddot_upper: float) -> dict:
    if not (0.0 <= dnorm < CORRECTION_CAYLEY_MONOTONE_LIMIT_RAD):
        raise RuntimeError(
            f"correction norm {dnorm} does not stay below {CORRECTION_CAYLEY_MONOTONE_LIMIT_RAD} rad monotone Cayley limit")
    scale = SIGNED.correction_cayley_scale_interval(dnorm)
    anorm = up(scale.hi * dnorm)
    adot = up((scale.hi if ddot_upper >= 0.0 else scale.lo) * ddot_upper)
    den = down(1.0 - up(0.25 * adot))
    if den <= 0.0:
        raise RuntimeError("signed candidate correction can cross Cayley denominator")
    a2 = up(anorm * anorm)
    q2 = up(q * q)
    num = up(up(a2 + q2) + up(2.0 * adot) + up(0.25 * up(a2 * q2)))
    qplus = _sqrt_up(up(max(0.0, num) / down(den * den)))
    return {
        "correction_cayley_scale": scale.as_list(),
        "correction_cayley_norm_upper": anorm,
        "a_dot_c_upper": adot,
        "denominator_lower": den,
        "post_update_q_upper": qplus,
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
    tcells = _tangent_cells(q, tangent_pieces)

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
        intercept, slope, sdetail = AWB._s_phase_affine_aw_bound(src, phase, Pp, domain, pnorm)
        aw = up(intercept + up(slope * aw0))

        for xi, x in enumerate(xcells):
            for mi, m in enumerate(force_cells):
                k, kh, gdetail = V3._tangent_structured_gain_bounds(
                    tilt=tilt, yaw=yaw, eps=eps, x=x, m=m,
                    paw=paw, racc_var=racc_var)
                canon = _canonical_kh_bounds(
                    tilt=tilt, yaw=yaw, x=x, m=m, paw=paw,
                    racc_var=racc_var, gain_detail=gdetail)
                eta = up(
                    up(VEFF.accel_attitude_eta_per_vector_norm_upper(q) * m.hi)
                    + up(VEFF.accel_latent_cross_gain_upper(q) * aw)
                )
                rho = up(aw + up(eta + ba))
                dnorm = up(up(kh * q) + up(k * rho))

                for ti, (rlo, rhi) in enumerate(tcells):
                    s_hi = _axial_upper(q, rlo)
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
                    try:
                        comp = _compose_q_upper(q, dnorm, ddot)
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
                            "eta_effective_aw_upper_mps2": eta,
                            "nuisance_effective_aw_upper_mps2": rho,
                            "Ktheta_norm_upper": k,
                            "KHtheta_norm_upper": kh,
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
                        "eta_effective_aw_upper_mps2": eta,
                        "nuisance_effective_aw_upper_mps2": rho,
                        "Ktheta_norm_upper": k,
                        "KHtheta_norm_upper": kh,
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
        failures.append("signed 30deg first-accelerometer family did not complete")
    if complete and not below_correction_chart:
        failures.append("signed 30deg correction family reached the 3 rad Cayley monotone limit")
    if complete and not inside_outer:
        failures.append("signed 30deg first-accelerometer family leaves operation-matched outer sector")
    passed = not failures

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_30DEG_SIGNED_FIRST_ACCEL_FINITE_ANGLE_SECTOR",
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
        "effective_aw_finite_angle_remainder_retained": True,
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
            "if this signed sector closes, propagate the same source/alignment/force children through the exact Joseph posterior and reset congruence before sample1; otherwise subdivide the first reported directional child without shrinking the declared 30deg/0.3g entrance"
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
        "effective_aw_finite_angle_remainder_retained", "accepted_and_rejected_branches_covered",
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
