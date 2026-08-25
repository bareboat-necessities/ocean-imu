#!/usr/bin/env python3
"""Analytic structured-gain enclosure for the first P5 H accelerometer map.

The rotation gauge reduces the first accepted accelerometer correction to a
rank-two canonical problem.  In the gauged body frame J_aw=I and the predicted
specific-force direction is e3.  The goLive attitude covariance before the
first vector correction has the source form

    P_theta = t I + delta v v^T + E,
    0 <= E <= eps I,

where t is the tilt variance, t+delta the gauged-yaw variance, and v is the
source-varying yaw axis.  Rotational symmetry about e3 means the gain norm
depends on v only through x=||v_tangent||^2 in [0,1]; no 2-D cube-face cover is
needed.

For E=0, with lambda=P_aw+R_acc, m=||f|| and
p_u=t+delta x, the two singular channels of the attitude Kalman gain are
available without a matrix inverse:

    g_perp = m t / (m^2 t + lambda),
    g_u    = m p_u / (m^2 p_u + lambda),
    g_z    = m delta sqrt(x(1-x)) / (m^2 p_u + lambda).

Thus ||K_theta|| is bounded by max(g_perp,sqrt(g_u^2+g_z^2)); multiplying the
same expressions by m gives ||K_theta H_theta||.  The small PSD remainder E is
added with a resolvent perturbation bound; it is never silently discarded.
All denominators are positive scalar intervals, so the loose 3x3 S>=R inverse
fallback is absent by construction.

This is still only the first accelerometer range gate.  It does not promote the
complete q<=8 word or set N_H_words.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p5_effective_vector_input as VEFF
import ou3_p5_first_accel_rotation_gauge as G1
import ou3_p5_first_accel_rotation_gauge_v3 as G3
import ou3_p5_full_h_prefix_cells as FULL1
import ou3_p5_full_h_prefix_cells_v3 as FULL3
import ou3_vector_uco_certificate as VECTOR

DEFAULT_DOMAIN = G1.DEFAULT_DOMAIN
SCHEMA = 1
DEPLOYED_CORRECTION_LIMIT_RAD = 6.0


def I(x: float) -> Interval:
    return Interval.point(float(x))


def _linear_cells(pieces: int) -> list[Interval]:
    if pieces < 1:
        raise ValueError("alignment pieces must be positive")
    edges = [i / pieces for i in range(pieces + 1)]
    return [Interval.outward_bounds(edges[i], edges[i + 1]) for i in range(pieces)]


def _sqrt_x1mx_upper(x: Interval) -> float:
    lo = max(0.0, x.lo)
    hi = min(1.0, x.hi)
    if lo > hi:
        raise ValueError("alignment cell outside [0,1]")
    candidates = [lo * (1.0 - lo), hi * (1.0 - hi)]
    if lo <= 0.5 <= hi:
        candidates.append(0.25)
    return FULL1.up(math.sqrt(max(candidates)))


def _structured_gain_bounds(
    *,
    tilt: float,
    yaw: float,
    eps: float,
    x: Interval,
    m: Interval,
    paw: Interval,
    racc_var: Interval,
) -> tuple[float, float, dict]:
    t = I(tilt)
    delta = Interval.outward_bounds(yaw - tilt, yaw - tilt)
    pu = t + delta * x
    lam = paw + racc_var
    if lam.lo <= 0.0:
        raise RuntimeError("first accelerometer lambda floor is nonpositive")
    m2 = m.square()
    den_perp = m2 * t + lam
    den_u = m2 * pu + lam
    if den_perp.lo <= 0.0 or den_u.lo <= 0.0:
        raise RuntimeError("structured gain denominator is nonpositive")

    geom = Interval(0.0, _sqrt_x1mx_upper(x))
    g_perp = m * t / den_perp
    g_u = m * pu / den_u
    g_z = m * delta * geom / den_u
    k0 = FULL1.up(max(
        g_perp.hi,
        math.sqrt(FULL1.up(g_u.hi * g_u.hi + g_z.hi * g_z.hi)),
    ))

    kh_perp = m2 * t / den_perp
    kh_u = m2 * pu / den_u
    kh_z = m2 * delta * geom / den_u
    kh0 = FULL1.up(max(
        kh_perp.hi,
        math.sqrt(FULL1.up(kh_u.hi * kh_u.hi + kh_z.hi * kh_z.hi)),
    ))

    # Resolvent perturbation for 0<=E<=eps I:
    # K(P0+E)-K(P0)
    # = E H' S^-1 - K0 H E H' S^-1.
    # ||H||=m, ||S^-1||<=1/lambda, ||E||<=eps.
    mhi = m.hi
    first = FULL1.up(eps * mhi / lam.lo)
    second = FULL1.up(k0 * mhi * mhi * eps / lam.lo)
    dk = FULL1.up(first + second)
    k = FULL1.up(k0 + dk)
    kh = FULL1.up(kh0 + FULL1.up(dk * mhi))
    return k, kh, {
        "lambda_lower": lam.lo,
        "p_u": pu.as_list(),
        "g_perp_upper": g_perp.hi,
        "g_u_upper": g_u.hi,
        "g_z_upper": g_z.hi,
        "K0_norm_upper": k0,
        "KH0_norm_upper": kh0,
        "PSD_remainder_K_perturbation_upper": dk,
    }


def build(
    domain_path: Path = DEFAULT_DOMAIN,
    *,
    source_pieces: int = 2,
    alignment_pieces: int = 16,
    force_magnitude_pieces: int = 4,
) -> dict:
    domain_path = Path(domain_path).resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("structured first-accel gain domain must not be trajectory fitted")
    if domain.get("configured_runtime", {}).get("imu_lever_arm_enabled") is not False:
        raise RuntimeError("structured first-accel gain requires lever arm disabled")

    G3._install_backend(domain_path, source_pieces)
    # G1.build normally installs the active full V3 backend; do it explicitly
    # because this producer evaluates only the first covariance prediction.
    FULL3._install_backend()

    heading = G1.HEADING.build(domain_path)
    go = G1.GOLIVE.build(domain_path)
    veff = VEFF.build(domain_path)
    vector = VECTOR.build()
    failures = [f"heading: {x}" for x in G1.HEADING.validate(heading)]
    failures += [f"goLive: {x}" for x in G1.GOLIVE.validate(go)]
    failures += [f"effective-vector: {x}" for x in VEFF.validate(veff)]
    failures += [f"vector: {x}" for x in VECTOR.validate(vector)]

    src_phases = G1._source_phase_children(source_pieces)
    xcells = _linear_cells(alignment_pieces)
    live = domain["normal_live"]
    force_cells = G1._geom_ranges(
        float(live["specific_force_norm_lower_mps2"]),
        float(live["specific_force_norm_upper_mps2"]),
        force_magnitude_pieces,
    )
    h = float(FULL1._source_cell()["dt_s"])
    tilt, yaw, eps = G1._attitude_covariance_epsilon(domain_path, h)
    q0 = float(heading["gauged_timeout_subbranch"]["full_attitude_cayley_norm_upper"])
    qpred = G1._q_after_first_prediction(q0, domain, h)
    vc = vector["configured_measurement_bounds"]
    Racc = FULL1._R_diag(float(vc["acc_measurement_std_mps2"]))
    racc_var = Racc[0][0]
    ba = float(domain["startup"]["physical_handoff_coordinate_bounds"]["accelerometer_bias_error_norm_upper_mps2"])

    total = 0
    over = 0
    max_d = 0.0
    min_margin = math.inf
    max_k = 0.0
    max_kh = 0.0
    first_over = None

    for si, (src, phase) in enumerate(src_phases):
        P0 = FULL1._initial_covariance(src, domain_path)
        F, Q, _ = FULL1._transition_and_Q(src, domain)
        Pp = FULL1._psd_tighten(FULL1.matrix_add(FULL1.matrix_mul(FULL1.matrix_mul(F, P0), FULL1.matrix_transpose(F)), Q))
        _pss, _psa, paw_pred = G1._scalar_axis_structure(Pp)
        aw_pred, eS_pred = G1._prediction_norms(src, domain)
        if phase == "due":
            paw, aw_norm = G1._due_paw_and_error_norm(Pp, src, aw_pred, eS_pred)
        else:
            paw, aw_norm = paw_pred, aw_pred

        for xi, x in enumerate(xcells):
            for mi, m in enumerate(force_cells):
                k, kh, detail = _structured_gain_bounds(
                    tilt=tilt,
                    yaw=yaw,
                    eps=eps,
                    x=x,
                    m=m,
                    paw=paw,
                    racc_var=racc_var,
                )
                eta = FULL1.up(
                    VEFF.accel_attitude_eta_per_vector_norm_upper(qpred) * m.hi
                    + VEFF.accel_latent_cross_gain_upper(qpred) * aw_norm
                )
                rho = FULL1.up(aw_norm + eta + ba)
                d = FULL1.up(FULL1.up(kh * qpred) + FULL1.up(k * rho))
                total += 1
                max_k = max(max_k, k)
                max_kh = max(max_kh, kh)
                max_d = max(max_d, d)
                margin = DEPLOYED_CORRECTION_LIMIT_RAD - d
                min_margin = min(min_margin, margin)
                if not math.isfinite(d) or d > DEPLOYED_CORRECTION_LIMIT_RAD:
                    over += 1
                    if first_over is None:
                        first_over = {
                            "source_phase_cell": si,
                            "pseudo_phase": phase,
                            "tau_s": src["tau_s"].as_list(),
                            "sigma_aw_mps2": src["sigma_aw_mps2"].as_list(),
                            "R_S_filter_std": src["R_S_filter_std"].as_list(),
                            "alignment_x_tangent_yaw_fraction": x.as_list(),
                            "force_magnitude_cell": mi,
                            "force_magnitude_mps2": m.as_list(),
                            "predicted_aw_error_norm_upper_mps2": aw_norm,
                            "effective_aw_eta_norm_upper_mps2": eta,
                            "nuisance_residual_norm_upper_mps2": rho,
                            "Ktheta_norm_upper": k,
                            "Ktheta_Htheta_norm_upper": kh,
                            "correction_norm_upper_rad": d,
                            "gain_detail": detail,
                        }

    closed = total > 0 and over == 0 and not failures
    next_obligation = (
        "PROPAGATE_STRUCTURED_FIRST_ACCEL_CHILDREN_THROUGH_JOSEPH_RESET_AND_LATER_PREFIXES"
        if closed else
        "REFINE_ACCEPTED_ACCEL_EFFECTIVE_AW_RESIDUAL_DIRECTION_COUPLING"
    )
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_FIRST_ACCEL_ANALYTIC_STRUCTURED_GAIN_ENCLOSURE",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "uses_v3_dependency_preserving_covariance_backend": True,
        "first_prefix_source_sparsity_certified": True,
        "rotation_gauge_sets_J_aw_to_identity": True,
        "specific_force_direction_gauged_to_e3": True,
        "yaw_axis_reduced_by_axial_symmetry_to_scalar_alignment": True,
        "alignment_coordinate": "x=||v_tangent||^2",
        "alignment_cell_count": len(xcells),
        "matrix_inverse_used_for_first_accel_gain": False,
        "loose_spectral_inverse_fallback_used": False,
        "analytic_rank_two_gain_channels_used": True,
        "attitude_PSD_remainder_retained_by_resolvent_bound": True,
        "deployed_correction_limit_rad": DEPLOYED_CORRECTION_LIMIT_RAD,
        "deployed_correction_limit_increased": False,
        "source_phase_cell_count": len(src_phases),
        "force_magnitude_cell_count": len(force_cells),
        "evaluated_child_count": total,
        "predicted_cayley_norm_upper": qpred,
        "attitude_covariance_tilt_variance": tilt,
        "attitude_covariance_gauged_yaw_variance": yaw,
        "attitude_covariance_PSD_remainder_upper": eps,
        "max_Ktheta_norm_upper": max_k,
        "max_Ktheta_Htheta_norm_upper": max_kh,
        "children_above_validated_correction_limit": over,
        "max_first_accelerometer_correction_norm_upper_rad": max_d,
        "minimum_correction_range_margin_rad": min_margin,
        "all_first_accelerometer_children_inside_validated_correction_range": closed,
        "first_unclosed_child": first_over,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_FIRST_ACCEL_STRUCTURED_GAIN_CERTIFICATE": "PASS" if closed else "NOT_ESTABLISHED",
        "next_obligation": next_obligation,
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    for k in (
        "source_generated_not_trajectory_fit",
        "uses_v3_dependency_preserving_covariance_backend",
        "first_prefix_source_sparsity_certified",
        "rotation_gauge_sets_J_aw_to_identity",
        "specific_force_direction_gauged_to_e3",
        "yaw_axis_reduced_by_axial_symmetry_to_scalar_alignment",
        "analytic_rank_two_gain_channels_used",
        "attitude_PSD_remainder_retained_by_resolvent_bound",
    ):
        if d.get(k) is not True:
            failures.append(f"{k} is not true")
    for k in (
        "source_replay_used", "filter_changed", "matrix_inverse_used_for_first_accel_gain",
        "loose_spectral_inverse_fallback_used", "deployed_correction_limit_increased",
        "whole_word_promoted_here", "N_H_words_set_here",
    ):
        if d.get(k) is not False:
            failures.append(f"{k} is not false")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        failures.append("deployed correction range changed")
    if int(d.get("evaluated_child_count", 0)) <= 0:
        failures.append("no structured gain children evaluated")
    status = d.get("P5_FIRST_ACCEL_STRUCTURED_GAIN_CERTIFICATE")
    if status == "PASS":
        if d.get("all_first_accelerometer_children_inside_validated_correction_range") is not True:
            failures.append("PASS without complete first-accelerometer closure")
        if d.get("first_unclosed_child") is not None:
            failures.append("PASS retains an unclosed child")
    elif status == "NOT_ESTABLISHED":
        if d.get("first_unclosed_child") is None:
            failures.append("nonclosure missing source child witness")
    else:
        failures.append("invalid structured-gain status")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces", type=int, default=2)
    ap.add_argument("--alignment-pieces", type=int, default=16)
    ap.add_argument("--force-magnitude-pieces", type=int, default=4)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build(
        args.domain.resolve(),
        source_pieces=args.source_pieces,
        alignment_pieces=args.alignment_pieces,
        force_magnitude_pieces=args.force_magnitude_pieces,
    )
    vf = validate(out)
    out["validation_pass"] = not vf
    out["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": out["P5_FIRST_ACCEL_STRUCTURED_GAIN_CERTIFICATE"],
        "children": out["evaluated_child_count"],
        "over_limit": out["children_above_validated_correction_limit"],
        "max_d": out["max_first_accelerometer_correction_norm_upper_rad"],
        "margin": out["minimum_correction_range_margin_rad"],
        "max_K": out["max_Ktheta_norm_upper"],
        "max_KH": out["max_Ktheta_Htheta_norm_upper"],
        "first_unclosed": out["first_unclosed_child"],
        "next": out["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
