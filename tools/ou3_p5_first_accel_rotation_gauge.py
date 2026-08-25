#!/usr/bin/env python3
"""Rotation-gauged source subdivision for the first P5 H accelerometer map.

The V3 full-matrix prefix backend is structurally complete but its source-wide
accelerometer Jacobian encloses each entry of J_aw=R_wb independently in
[-1,1].  That box contains matrices which are not rotations, destroys the
isotropic first-prefix a_w covariance, and forces the 3x3 innovation inverse
onto a deliberately loose S>=R fallback.  The resulting correction bound is a
numerical dependency artifact, not a reachable filter correction.

At the first normal-Live accelerometer packet the source has a stronger exact
symmetry:

* the goLive attitude/a_w cross covariance is zero and one prediction plus an
  optional isotropic S=0 pseudo update preserves that zero;
* the a_w marginal and all first-prefix S/a_w axis blocks are isotropic;
* the shipping accelerometer Jacobian is J_aw=R_wb, an orthogonal matrix;
* the physical a_w error is supplied as a Euclidean ball.

Therefore a simultaneous orthogonal change of the world linear coordinates to
the current body coordinates sends J_aw to I without changing any of those
sets.  A second body-coordinate rotation sends the nonzero predicted specific
force to +e3.  The only attitude-covariance direction that remains is the
source-varying yaw axis.  It is covered by normalized cube-face cells, while
force magnitude, tau, sigma_aw, R_S and the first pseudo due/not-due phase are
kept source-correlated.

This stage only certifies whether the *first* accelerometer attitude correction
lies in the already validated deployed-quaternion range [0,6] rad.  It does not
promote the complete q<=8 word and does not set N_H_words.  If it closes, the
next obligation is to carry the same structured/gauged child through Joseph,
reset and later prefixes.  No replay samples or fitted extrema are used.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import (
    Interval,
    matrix_abs_col_sum_upper,
    matrix_abs_row_sum_upper,
    matrix_add,
    matrix_mul,
    matrix_transpose,
)
import ou3_full_process_ucc as PROCESS
import ou3_p5_effective_vector_input as VEFF
import ou3_p5_first_accel_subdivision as SUB
import ou3_p5_full_h_prefix_cells as V1
import ou3_p5_full_h_prefix_cells_v2 as V2
import ou3_p5_full_h_prefix_cells_v3 as V3
import ou3_p5_go_live_covariance_stage as GOLIVE
import ou3_p5_heading_handoff_contract as HEADING
import ou3_source_reachable_matrix_p3 as P3CELL
import ou3_vector_uco_certificate as VECTOR
import ou3_validated_transcendentals as VT

DEFAULT_DOMAIN = V1.DEFAULT_DOMAIN
SCHEMA = 1
DEPLOYED_CORRECTION_LIMIT_RAD = 6.0
FLOAT_EPS = 2.0 ** -23


def I(x: float) -> Interval:
    return Interval.point(float(x))


def _sqrt_interval(x: Interval) -> Interval:
    if x.lo < 0.0 or not math.isfinite(x.hi):
        raise ValueError("finite nonnegative square-root interval required")
    return Interval(V1.down(math.sqrt(x.lo)), V1.up(math.sqrt(x.hi)))


def _op2_upper(A) -> float:
    """Rigorous ||A||_2 upper from sqrt(||A||_1 ||A||_inf)."""
    n1 = matrix_abs_col_sum_upper(A)
    ni = matrix_abs_row_sum_upper(A)
    return V1.up(math.sqrt(V1.up(n1 * ni)))


def _unit_face_cells(pieces: int) -> list[list[Interval]]:
    """Finite cover of S^2 by normalized dominant-coordinate cube faces."""
    if pieces < 1:
        raise ValueError("face pieces must be positive")
    edges = [-1.0 + 2.0 * i / pieces for i in range(pieces + 1)]
    uv = [Interval.outward_bounds(edges[i], edges[i + 1]) for i in range(pieces)]
    out: list[list[Interval]] = []
    for axis in range(3):
        other = [i for i in range(3) if i != axis]
        for sign in (-1.0, 1.0):
            for u in uv:
                for v in uv:
                    raw = [I(0.0), I(0.0), I(0.0)]
                    raw[axis] = I(sign)
                    raw[other[0]] = u
                    raw[other[1]] = v
                    sq = I(1.0) + u.square() + v.square()
                    inv = _sqrt_interval(sq).reciprocal()
                    out.append([x * inv for x in raw])
    return out


def _geom_ranges(lo: float, hi: float, pieces: int) -> list[Interval]:
    if not (0.0 < lo <= hi and pieces >= 1):
        raise ValueError("positive geometric range required")
    if pieces == 1 or lo == hi:
        return [Interval.outward_bounds(lo, hi)]
    r = (hi / lo) ** (1.0 / pieces)
    edges = [lo]
    for _ in range(pieces - 1):
        edges.append(edges[-1] * r)
    edges.append(hi)
    return [Interval.outward_bounds(edges[i], edges[i + 1]) for i in range(pieces)]


def _split_interval_at(x: Interval, cut: float) -> list[Interval]:
    if cut <= x.lo or cut >= x.hi:
        return [x]
    return [Interval.outward_bounds(x.lo, cut), Interval.outward_bounds(cut, x.hi)]


def _source_phase_children(source_pieces: int) -> list[tuple[dict, str]]:
    """Split tuner source and retain only first-sample-reachable pseudo phases."""
    sched = P3CELL.source_schedule()
    src0 = V1._source_cell()
    h = float(src0["dt_s"])
    ratio = float(sched["pseudo_ratio"])
    # periodic_update_due<float> uses total + 16 eps max(1,period) < period.
    # All deployed first-prefix periods are below one second, so max=1.
    tol = 16.0 * FLOAT_EPS
    threshold_period = h + tol
    tau_cut = threshold_period / ratio

    tau_base = SUB._geom_split(src0["tau_s"], source_pieces)
    taus: list[Interval] = []
    for t in tau_base:
        taus.extend(_split_interval_at(t, tau_cut))
    sigmas = SUB._geom_split(src0["sigma_aw_mps2"], source_pieces)
    rss = SUB._geom_split(src0["R_S_filter_std"], source_pieces)

    out: list[tuple[dict, str]] = []
    for tau in taus:
        plo, phi = P3CELL.cadence_bounds(tau, sched)
        period = Interval.outward_bounds(plo, phi)
        due_possible = period.lo <= threshold_period
        not_due_possible = period.hi > threshold_period
        for sigma in sigmas:
            for rs in rss:
                src = dict(src0)
                src["tau_s"] = tau
                src["sigma_aw_mps2"] = sigma
                src["R_S_filter_std"] = rs
                src["pseudo_period_s"] = period
                if due_possible:
                    out.append((src, "due"))
                if not_due_possible:
                    out.append((src, "not_due"))
    return out


def _attitude_covariance_epsilon(domain_path: Path, h: float) -> tuple[float, float, float]:
    """Return tilt, yaw and isotropic one-step PSD remainder upper."""
    go = GOLIVE.build(domain_path)["goLive_H_covariance_seed"]
    tilt = float(go["attitude_covariance_seed"]["tilt_variance"])
    yaw = float(go["attitude_covariance_seed"]["gauged_yaw_variance"])
    proc = PROCESS.build()["source_constants"]
    qg = float(proc["gyro_noise_density_rad_sqrt_s"]) ** 2
    qb = float(proc["gyro_bias_rw_variance_density"])
    pb0 = V1._source_pb0()
    pbg = V1.up(pb0 + qb * V2._startup_timeout_s())
    bb = V1.up(qb * h * h * h / 3.0)
    # ||B||_2<=h for B=int_0^h R(s)ds.  V2's interval Q_theta has diagonal
    # upper qg*h+bb and two off-diagonal absolute uppers bb, hence Gershgorin
    # gives qg*h+3bb as a source-uniform spectral upper.
    eps = V1.up(V1.up(pbg * h * h) + V1.up(qg * h + 3.0 * bb))
    if not (0.0 < tilt <= yaw and eps >= 0.0):
        raise RuntimeError("invalid source attitude covariance seed")
    return tilt, yaw, eps


def _ptheta_cell(v: list[Interval], tilt: float, yaw: float, eps: float):
    delta = Interval.outward_bounds(yaw - tilt, yaw - tilt)
    P = V1._zero(3, 3)
    for i in range(3):
        for j in range(3):
            base = delta * v[i] * v[j]
            if i == j:
                P[i][j] = I(tilt) + base + Interval(0.0, V1.up(eps))
            else:
                # Any 0<=E<=eps I has |E_ij|<=eps/2.  eps itself is already
                # tiny; using eps here is a simpler conservative enclosure.
                P[i][j] = base + Interval(-V1.up(eps), V1.up(eps))
    return P


def _canonical_Htheta(force_mag: Interval):
    # Gauge f=(0,0,+m).  J_att=-[f]_x.
    H = V1._zero(3, 3)
    H[0][1] = force_mag
    H[1][0] = -force_mag
    return H


def _scalar_axis_structure(Pm) -> tuple[Interval, Interval, Interval]:
    """Extract first-prefix P_SS, P_Saw, P_aw and require exact axis symmetry."""
    pss = Pm[12][12]
    psa = Pm[12][15]
    paw = Pm[15][15]
    for ax in range(3):
        if Pm[12 + ax][12 + ax] != pss or Pm[12 + ax][15 + ax] != psa or Pm[15 + ax][15 + ax] != paw:
            raise RuntimeError("first-prefix linear covariance lost axis symmetry")
    for ai in range(3):
        for aj in range(3):
            if ai != aj:
                for a, b in ((12 + ai, 12 + aj), (12 + ai, 15 + aj), (15 + ai, 15 + aj)):
                    if Pm[a][b].lo != 0.0 or Pm[a][b].hi != 0.0:
                        raise RuntimeError("first-prefix linear covariance gained cross-axis terms")
    for ti in range(3):
        for aj in range(3):
            z = Pm[ti][15 + aj]
            if z.lo != 0.0 or z.hi != 0.0:
                raise RuntimeError("first-prefix theta/a_w covariance is not exactly zero")
    return pss, psa, paw


def _due_paw_and_error_norm(Pp, src: dict, aw_pred_norm: float, eS_pred_norm: float) -> tuple[Interval, float]:
    pss, psa, paw = _scalar_axis_structure(Pp)
    rs2 = src["R_S_filter_std"].square()
    den = pss + rs2
    if den.lo <= 0.0:
        raise RuntimeError("S pseudo innovation lost positive floor")
    k = psa / den
    reduction = psa.square() / den
    raw = paw - reduction
    paw_due = V1._intersect(raw, Interval(0.0, paw.hi))
    aw_due_norm = V1.up(aw_pred_norm + V1.up(k.abs_upper() * eS_pred_norm))
    return paw_due, aw_due_norm


def _prediction_norms(src: dict, domain: dict) -> tuple[float, float]:
    h = float(src["dt_s"])
    alpha, _pv, _pp, phiS = V2._monotone_coeff_hull(src["tau_s"], h)
    b = domain["startup"]["physical_handoff_coordinate_bounds"]
    aw0 = float(b["latent_acceleration_error_norm_upper_mps2"])
    v0 = float(b["velocity_error_norm_upper_mps"])
    p0 = float(b["position_error_norm_upper_m"])
    S0 = float(b["integral_displacement_error_norm_upper_m_s"])
    aw = V1.up(alpha.hi * aw0)
    s = V1.up(S0 + V1.up(h * p0) + V1.up(0.5 * h * h * v0) + V1.up(phiS.hi * aw0))
    return aw, s


def _q_after_first_prediction(q0: float, domain: dict, h: float) -> float:
    b = domain["startup"]["physical_handoff_coordinate_bounds"]
    bg = float(b["gyro_bias_error_norm_upper_rad_s"])
    wd = float(domain["startup"]["effective_deterministic_gyro_transport_disturbance_upper_rad_s"])
    theta = V1.up(h * (bg + wd))
    half = V1.up(0.5 * theta)
    sn = VT.sin_point(half)
    co = VT.cos_point(half)
    if co.lo <= 0.0:
        raise RuntimeError("one-step transport input reaches Cayley antipode")
    a = V1.up(2.0 * sn.hi / co.lo)
    den = V1.down(1.0 - V1.up(0.25 * a * q0))
    if den <= 0.0:
        raise RuntimeError("one-step Cayley norm denominator is nonpositive")
    num = V1.up(q0 + a + V1.up(0.5 * a * q0))
    return V1.up(num / den)


def _gain_bounds(Ptheta, paw: Interval, force_mag: Interval, Racc):
    H = _canonical_Htheta(force_mag)
    PHt = matrix_mul(Ptheta, matrix_transpose(H))
    lam = paw + Racc[0][0]
    Reff = V1._zero(3, 3)
    for i in range(3):
        Reff[i][i] = lam
    S = matrix_add(matrix_mul(H, PHt), Reff)
    Sinv, backend = V1._spd_inverse_enclosure(S, Reff)
    K = matrix_mul(PHt, Sinv)
    KH = matrix_mul(K, H)
    return _op2_upper(K), _op2_upper(KH), backend, S


def build(
    domain_path: Path = DEFAULT_DOMAIN,
    *,
    source_pieces: int = 2,
    yaw_axis_face_pieces: int = 4,
    force_magnitude_pieces: int = 4,
) -> dict:
    V3._install_backend()
    domain_path = Path(domain_path).resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("rotation-gauged P5 diagnostic must not be trajectory fitted")
    if domain.get("configured_runtime", {}).get("imu_lever_arm_enabled") is not False:
        raise RuntimeError("rotation-gauged first accel stage requires lever arm disabled")

    heading = HEADING.build(domain_path)
    go = GOLIVE.build(domain_path)
    veff = VEFF.build(domain_path)
    vector = VECTOR.build()
    failures = [f"heading: {x}" for x in HEADING.validate(heading)]
    failures += [f"goLive: {x}" for x in GOLIVE.validate(go)]
    failures += [f"effective-vector: {x}" for x in VEFF.validate(veff)]
    failures += [f"vector: {x}" for x in VECTOR.validate(vector)]

    src_phases = _source_phase_children(source_pieces)
    yaw_cells = _unit_face_cells(yaw_axis_face_pieces)
    live = domain["normal_live"]
    force_cells = _geom_ranges(
        float(live["specific_force_norm_lower_mps2"]),
        float(live["specific_force_norm_upper_mps2"]),
        force_magnitude_pieces,
    )
    h = float(V1._source_cell()["dt_s"])
    tilt, yaw, peps = _attitude_covariance_epsilon(domain_path, h)
    q0 = float(heading["gauged_timeout_subbranch"]["full_attitude_cayley_norm_upper"])
    qpred = _q_after_first_prediction(q0, domain, h)
    vc = vector["configured_measurement_bounds"]
    Racc = V1._R_diag(float(vc["acc_measurement_std_mps2"]))
    ba = float(domain["startup"]["physical_handoff_coordinate_bounds"]["accelerometer_bias_error_norm_upper_mps2"])

    max_d = 0.0
    min_margin = math.inf
    first_over = None
    over = 0
    fixed = 0
    fallback = 0
    source_structure_checked = True

    for si, (src, phase) in enumerate(src_phases):
        P0 = V1._initial_covariance(src, domain_path)
        e0 = V1._initial_error(domain)
        F, Q, _Rstep = V1._transition_and_Q(src, domain)
        Pp = V1._psd_tighten(matrix_add(matrix_mul(matrix_mul(F, P0), matrix_transpose(F)), Q))
        _pss, _psa, paw_pred = _scalar_axis_structure(Pp)
        aw_pred, eS_pred = _prediction_norms(src, domain)
        if phase == "due":
            paw, aw_norm = _due_paw_and_error_norm(Pp, src, aw_pred, eS_pred)
        else:
            paw, aw_norm = paw_pred, aw_pred

        for vi, v in enumerate(yaw_cells):
            Ptheta = _ptheta_cell(v, tilt, yaw, peps)
            for mi, m in enumerate(force_cells):
                try:
                    knorm, khnorm, backend, S = _gain_bounds(Ptheta, paw, m, Racc)
                except Exception as exc:
                    dupper = math.inf
                    backend = f"EXCEPTION:{type(exc).__name__}"
                    eta = math.inf
                    rho = math.inf
                else:
                    if backend == "FIXED_PIVOT_INTERVAL_GAUSS_JORDAN":
                        fixed += 1
                    else:
                        fallback += 1
                    eta = V1.up(
                        VEFF.accel_attitude_eta_per_vector_norm_upper(qpred) * m.hi
                        + VEFF.accel_latent_cross_gain_upper(qpred) * aw_norm
                    )
                    rho = V1.up(aw_norm + eta + ba)
                    dupper = V1.up(V1.up(khnorm * qpred) + V1.up(knorm * rho))
                max_d = max(max_d, dupper)
                margin = DEPLOYED_CORRECTION_LIMIT_RAD - dupper
                min_margin = min(min_margin, margin)
                if not math.isfinite(dupper) or dupper > DEPLOYED_CORRECTION_LIMIT_RAD:
                    over += 1
                    if first_over is None:
                        first_over = {
                            "source_phase_cell": si,
                            "pseudo_phase": phase,
                            "tau_s": src["tau_s"].as_list(),
                            "sigma_aw_mps2": src["sigma_aw_mps2"].as_list(),
                            "R_S_filter_std": src["R_S_filter_std"].as_list(),
                            "pseudo_period_s": src["pseudo_period_s"].as_list(),
                            "yaw_axis_cell": vi,
                            "yaw_axis_box": [x.as_list() for x in v],
                            "force_magnitude_cell": mi,
                            "force_magnitude_mps2": m.as_list(),
                            "predicted_cayley_norm_upper": qpred,
                            "predicted_aw_error_norm_upper_mps2": aw_norm,
                            "effective_aw_eta_norm_upper_mps2": eta,
                            "nuisance_residual_norm_upper_mps2": rho,
                            "correction_norm_upper_rad": dupper,
                            "inverse_backend": backend,
                        }

    total = len(src_phases) * len(yaw_cells) * len(force_cells)
    closed = total > 0 and over == 0 and not failures
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_FIRST_ACCEL_ROTATION_GAUGED_FULL_MATRIX_SUBDIVISION",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "uses_v3_dependency_preserving_backend": True,
        "deployed_correction_limit_rad": DEPLOYED_CORRECTION_LIMIT_RAD,
        "deployed_correction_limit_increased": False,
        "rotation_gauge_sets_J_aw_to_identity": True,
        "rotation_gauge_sets_specific_force_direction_to_e3": True,
        "gauge_requires_first_prefix_theta_aw_cross_zero": True,
        "gauge_requires_first_prefix_aw_axis_isotropy": True,
        "first_prefix_source_structure_checked": source_structure_checked,
        "attitude_seed_rank_one_yaw_axis_retained": True,
        "yaw_axis_cube_face_cover_complete": True,
        "pseudo_phase_coupled_to_tau_before_branching": True,
        "source_parameter_subdivision": ["tau", "sigma_aw", "R_S", "pseudo_phase"],
        "source_phase_cell_count": len(src_phases),
        "yaw_axis_cell_count": len(yaw_cells),
        "force_magnitude_cell_count": len(force_cells),
        "evaluated_child_count": total,
        "predicted_cayley_norm_upper": qpred,
        "attitude_covariance_tilt_variance": tilt,
        "attitude_covariance_gauged_yaw_variance": yaw,
        "attitude_covariance_one_step_psd_remainder_spectral_upper": peps,
        "fixed_pivot_inverse_count": fixed,
        "spectral_fallback_inverse_count": fallback,
        "children_above_validated_correction_limit": over,
        "max_first_accelerometer_correction_norm_upper_rad": max_d,
        "minimum_correction_range_margin_rad": min_margin,
        "all_first_accelerometer_children_inside_validated_correction_range": closed,
        "first_unclosed_child": first_over,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_FIRST_ACCEL_ROTATION_GAUGED_CERTIFICATE": "PASS" if closed else "NOT_ESTABLISHED",
        "next_obligation": (
            "PROPAGATE_ROTATION_GAUGED_CHILDREN_THROUGH_ACCEL_JOSEPH_RESET_AND_LATER_PREFIXES"
            if closed else
            "REFINE_FIRST_ACCEL_ATTITUDE_COVARIANCE_AND_EFFECTIVE_AW_DIRECTION_COUPLING"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    for k in (
        "source_generated_not_trajectory_fit",
        "uses_v3_dependency_preserving_backend",
        "rotation_gauge_sets_J_aw_to_identity",
        "rotation_gauge_sets_specific_force_direction_to_e3",
        "gauge_requires_first_prefix_theta_aw_cross_zero",
        "gauge_requires_first_prefix_aw_axis_isotropy",
        "first_prefix_source_structure_checked",
        "attitude_seed_rank_one_yaw_axis_retained",
        "yaw_axis_cube_face_cover_complete",
        "pseudo_phase_coupled_to_tau_before_branching",
    ):
        if d.get(k) is not True:
            failures.append(f"{k} is not true")
    for k in (
        "source_replay_used", "filter_changed", "deployed_correction_limit_increased",
        "whole_word_promoted_here", "N_H_words_set_here",
    ):
        if d.get(k) is not False:
            failures.append(f"{k} is not false")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        failures.append("deployed correction range changed")
    if int(d.get("evaluated_child_count", 0)) <= 0:
        failures.append("no rotation-gauged children evaluated")
    if int(d.get("spectral_fallback_inverse_count", 0)) != 0:
        failures.append("rotation-gauged innovation still required loose spectral inverse fallback")
    status = d.get("P5_FIRST_ACCEL_ROTATION_GAUGED_CERTIFICATE")
    if status == "PASS":
        if d.get("all_first_accelerometer_children_inside_validated_correction_range") is not True:
            failures.append("PASS without complete first-accelerometer correction-range closure")
        if d.get("first_unclosed_child") is not None:
            failures.append("PASS retains an unclosed child")
    elif status == "NOT_ESTABLISHED":
        if d.get("first_unclosed_child") is None:
            failures.append("nonclosure missing source child witness")
    else:
        failures.append("invalid rotation-gauged status")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces", type=int, default=2)
    ap.add_argument("--yaw-axis-face-pieces", type=int, default=4)
    ap.add_argument("--force-magnitude-pieces", type=int, default=4)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build(
        args.domain.resolve(),
        source_pieces=args.source_pieces,
        yaw_axis_face_pieces=args.yaw_axis_face_pieces,
        force_magnitude_pieces=args.force_magnitude_pieces,
    )
    vf = validate(out)
    out["validation_pass"] = not vf
    out["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": out["P5_FIRST_ACCEL_ROTATION_GAUGED_CERTIFICATE"],
        "children": out["evaluated_child_count"],
        "over_limit": out["children_above_validated_correction_limit"],
        "max_d": out["max_first_accelerometer_correction_norm_upper_rad"],
        "margin": out["minimum_correction_range_margin_rad"],
        "fixed_inverse": out["fixed_pivot_inverse_count"],
        "fallback_inverse": out["spectral_fallback_inverse_count"],
        "first_unclosed": out["first_unclosed_child"],
        "next": out["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
