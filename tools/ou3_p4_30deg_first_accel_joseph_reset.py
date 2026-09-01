#!/usr/bin/env python3
"""Joseph/reset lift for the dependency-preserving 30 deg P4 first accelerometer.

This producer is deliberately narrower than a full-word covariance proof.  It
lifts the *same source/alignment/force family* used by
``ou3_p4_30deg_signed_first_accel_sector_v2`` through the part of the shipping
Joseph/reset map that previously destroyed the interval enclosure: the small
one-step attitude PSD remainder.

The canonical gravity gauge has

    P_theta = t I + delta v v' + E,   0 <= E <= eps I,

with ``v=(sqrt(x),0,sqrt(1-x))``.  The diagonal of E is absorbed into the two
force-tangent variances.  Matrix order then gives every off-diagonal entry of
the remaining zero-diagonal O an absolute bound eps/2 and ||O||_2 <= eps, the
same lemma already certified by P5 V36/V40.

For the nominal accelerometer Joseph factor B=I-KH, the three columns acting on
O are orthogonal in the rotation gauge.  Their norms are bounded directly from
the two tangent innovation channels; the force/yaw covariance cross term is
retained in the first tangent column.  Therefore B O B' is bounded by the
maximum row sum of the scaled zero-diagonal component matrix instead of by a
full entrywise matrix box.  Gain-perturbation cross and quadratic terms are
then added exactly as in V40.

Finally the shipping MEKF reset congruence is applied with

    ||G_reset(d)||_2 <= sqrt(1 + ||d||^2/4),

using the signed producer's certified correction norm.  No correction range is
increased, no source set is changed, and no P4/full-word claim is made here.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p4_30deg_signed_first_accel_sector_v2 as SIGNED
import ou3_p4_candidate_first_accel_range_v3 as RANGE
import ou3_p5_first_accel_rotation_gauge as RG
import ou3_p5_first_accel_structured_gain as SG
import ou3_p5_full_h_prefix_cells as FULL
import ou3_p5_full_h_prefix_cells_v3 as FULL3
import ou3_vector_uco_certificate as VECTOR

DEFAULT_DOMAIN = SIGNED.DEFAULT_DOMAIN
SCHEMA = 1


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def _sum_up(*xs: float) -> float:
    y = 0.0
    for x in xs:
        y = up(y + float(x))
    return y


def _mul_up(*xs: float) -> float:
    y = 1.0
    for x in xs:
        y = up(y * float(x))
    return y


def _ratio_lambda_over_x_plus_lambda(x: Interval, lam: Interval) -> float:
    """Upper bound lambda/(x+lambda), preserving the shared positive lambda."""
    if x.lo < 0.0 or lam.lo <= 0.0:
        raise RuntimeError("positive Joseph ratio domain lost")
    den = down(x.lo + lam.hi)
    if den <= 0.0:
        raise RuntimeError("positive Joseph ratio denominator lost")
    return min(1.0, up(lam.hi / den))


def _component_transport(*, tilt: float, yaw: float, eps: float,
                         x: Interval, m: Interval, paw: Interval,
                         racc: Interval, dnorm: float) -> dict:
    """Transport the absorbed-diagonal PSD remainder through Joseph + reset."""
    I = SG.I
    t0 = I(tilt)
    delta = Interval.outward_bounds(yaw - tilt, yaw - tilt)
    geom = Interval(0.0, SG._sqrt_x1mx_upper(x))
    m2 = m.square()
    lam = paw + racc
    if lam.lo <= 0.0:
        raise RuntimeError("Joseph lambda floor is nonpositive")

    # Absorb diag(E), 0<=E<=eps I, into the tangent variances.  Widen the two
    # nominal tangent channels independently; keep the source yaw-axis cross
    # covariance delta*sqrt(x(1-x)) explicit.
    te = Interval.outward_bounds(tilt, up(tilt + eps))
    pu0 = t0 + delta * x
    pue = Interval.outward_bounds(pu0.lo, up(pu0.hi + eps))
    den_perp = m2 * te + lam
    den_u = m2 * pue + lam
    if den_perp.lo <= 0.0 or den_u.lo <= 0.0:
        raise RuntimeError("Joseph tangent innovation floor lost")

    beta_perp = _ratio_lambda_over_x_plus_lambda(m2 * te, lam)
    beta_u = _ratio_lambda_over_x_plus_lambda(m2 * pue, lam)

    # aw rows of B=I-KH: |m P_aw / (m^2 P_theta + lambda)|.
    gamma_perp = (m * paw / den_perp).abs_upper()
    gamma_u = (m * paw / den_u).abs_upper()

    # Base yaw-axis covariance produces the exact axial component in the u
    # tangent column.  It is not an E perturbation and therefore must stay in
    # the nominal column norm.
    z_u = (m2 * delta * geom / den_u).abs_upper()
    s_perp = up(math.sqrt(max(0.0, _sum_up(beta_perp * beta_perp,
                                           gamma_perp * gamma_perp))))
    s_u = up(math.sqrt(max(0.0, _sum_up(beta_u * beta_u,
                                       gamma_u * gamma_u,
                                       z_u * z_u))))
    b0 = max(1.0, s_perp, s_u)

    # V36/V40 matrix-order lemma after diagonal absorption:
    # each offdiag <= eps/2 and ||O|| <= eps.  With unequal tangent column
    # norms the transformed zero-diagonal matrix has these three row sums.
    e = up(0.5 * eps)
    eop = up(eps)
    row1 = _mul_up(e, _sum_up(up(s_perp * s_u), s_perp))
    row2 = _mul_up(e, _sum_up(up(s_perp * s_u), s_u))
    row3 = _mul_up(e, _sum_up(s_perp, s_u))
    nominal = max(row1, row2, row3)

    # Gain perturbation from O.  Use the same tangent-only resolvent as the P4
    # signed/range stage, but include the aw gain rows as V40 does.
    ktheta, _kh, detail = RANGE._tangent_structured_gain_bounds(
        tilt=tilt, yaw=yaw, eps=eps, x=x, m=m, paw=paw,
        racc_var=racc)
    del ktheta
    dS = float(detail["tangent_innovation_perturbation_upper"])
    inv = float(detail["perturbed_tangent_inverse_operator_upper"])
    dkth = float(detail["PSD_remainder_K_perturbation_upper"])
    kaw0 = max((paw / (m2 * t0 + lam)).abs_upper(),
               (paw / (m2 * pu0 + lam)).abs_upper())
    dkaw = up(up(kaw0 * dS) * inv)
    dk = up(math.sqrt(max(0.0, _sum_up(dkth * dkth, dkaw * dkaw))))
    delta_B = _mul_up(m.hi, dk)
    cross = _mul_up(2.0, delta_B, eop, b0)
    quadratic = _mul_up(delta_B, delta_B, eop)
    posterior = _sum_up(nominal, cross, quadratic)

    if not (math.isfinite(dnorm) and 0.0 <= dnorm < math.pi):
        raise RuntimeError("signed correction is outside the reset chart")
    reset_op = up(math.sqrt(up(1.0 + up(0.25 * up(dnorm * dnorm)))))
    after_reset = _mul_up(reset_op, reset_op, posterior)
    return {
        "tangent_column_norm_perp_upper": s_perp,
        "tangent_column_norm_u_upper": s_u,
        "nominal_B_operator_upper": b0,
        "PSD_offdiagonal_entry_abs_upper": e,
        "PSD_remainder_operator_upper": eop,
        "nominal_component_transport_upper": nominal,
        "deltaK_theta_operator_upper": dkth,
        "deltaK_aw_operator_upper": dkaw,
        "deltaK_full_operator_upper": dk,
        "deltaB_operator_upper": delta_B,
        "cross_transport_upper": cross,
        "quadratic_transport_upper": quadratic,
        "posterior_PSD_remainder_operator_upper": posterior,
        "reset_operator_upper": reset_op,
        "post_reset_PSD_remainder_operator_upper": after_reset,
        "yaw_axis_cross_component_retained": True,
        "Joseph_zero_diagonal_component_matrix_used": True,
        "reset_congruence_used": True,
    }


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 2,
          alignment_pieces: int = 16, force_magnitude_pieces: int = 4,
          tangent_pieces: int = 32) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    signed = SIGNED.build(
        path, source_pieces=source_pieces, alignment_pieces=alignment_pieces,
        force_magnitude_pieces=force_magnitude_pieces,
        tangent_pieces=tangent_pieces)
    sf = SIGNED.validate(signed)
    failures = [f"signed: {x}" for x in sf]
    if signed.get("P4_30DEG_SIGNED_FIRST_ACCEL_SECTOR_CERTIFICATE") != "PASS":
        failures.append("signed 30deg prerequisite did not close")

    FULL3._install_backend()
    vec = VECTOR.build()
    failures += [f"vector: {x}" for x in VECTOR.validate(vec)]
    h = float(FULL._source_cell()["dt_s"])
    tilt, yaw, eps = RG._attitude_covariance_epsilon(path, h)
    racc = FULL._R_diag(float(
        vec["configured_measurement_bounds"]["acc_measurement_std_mps2"]))[0][0]
    live = domain["normal_live"]
    xcells = SG._linear_cells(alignment_pieces)
    mcells = RG._geom_ranges(
        float(live["specific_force_norm_lower_mps2"]),
        float(live["specific_force_norm_upper_mps2"]),
        force_magnitude_pieces)
    src_phases = RG._source_phase_children(source_pieces)
    dnorm = float(signed.get("max_correction_norm_upper_rad", math.inf))

    rows = []
    worst = None
    for si, (src, phase) in enumerate(src_phases):
        P0 = FULL._initial_covariance(src, path)
        F, Q, _ = FULL._transition_and_Q(src, domain)
        Pp = FULL._psd_tighten(FULL.matrix_add(
            FULL.matrix_mul(FULL.matrix_mul(F, P0), FULL.matrix_transpose(F)), Q))
        _pss, _psa, paw_pred = RG._scalar_axis_structure(Pp)
        paw = RG._due_paw_and_error_norm(Pp, src, 0.0, 0.0)[0] if phase == "due" else paw_pred
        for xi, x in enumerate(xcells):
            for mi, m in enumerate(mcells):
                try:
                    d = _component_transport(
                        tilt=tilt, yaw=yaw, eps=eps, x=x, m=m,
                        paw=paw, racc=racc, dnorm=dnorm)
                except Exception as exc:
                    failures.append(
                        f"source={si} phase={phase} alignment={xi} force={mi}: "
                        f"{type(exc).__name__}: {exc}")
                    continue
                row = {"source_phase_cell": si, "pseudo_phase": phase,
                       "alignment_cell": xi, "force_cell": mi,
                       "alignment_x": x.as_list(),
                       "force_magnitude_mps2": m.as_list(), **d}
                rows.append(row)
                if worst is None or d["post_reset_PSD_remainder_operator_upper"] > \
                        worst["post_reset_PSD_remainder_operator_upper"]:
                    worst = row

    complete = bool(rows) and not failures
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_30DEG_FIRST_ACCEL_JOSEPH_RESET_PSD_REMAINDER",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "candidate_angle_deg": 30.0,
        "declared_startup_aw_error_fraction_g": float(
            domain["startup"]["latent_acceleration_error_fraction_g"]),
        "signed_30deg_prerequisite_status": signed.get(
            "P4_30DEG_SIGNED_FIRST_ACCEL_SECTOR_CERTIFICATE"),
        "signed_max_correction_norm_upper_rad": dnorm,
        "same_source_alignment_force_family_as_signed_stage": True,
        "PSD_matrix_order_diagonal_absorption_used": True,
        "PSD_offdiagonal_abs_le_eps_over_2_used": True,
        "PSD_zero_diagonal_operator_le_eps_used": True,
        "unequal_canonical_tangent_column_norms_used": True,
        "yaw_axis_cross_component_retained": True,
        "aw_gain_rows_retained_in_deltaK": True,
        "Joseph_cross_terms_retained": True,
        "Joseph_quadratic_term_retained": True,
        "shipping_reset_congruence_used": True,
        "evaluated_source_alignment_force_cells": len(rows),
        "worst_cell": worst,
        "max_post_reset_PSD_remainder_operator_upper": (
            None if worst is None else worst["post_reset_PSD_remainder_operator_upper"]),
        "P4_30DEG_FIRST_ACCEL_JOSEPH_RESET_PSD_REMAINDER": (
            "PASS" if complete else "NOT_ESTABLISHED"),
        "FULL_18_STATE_COVARIANCE_PROPAGATED_HERE": False,
        "A21_MODE_PROPAGATED_HERE": False,
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE": False,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "next_obligation": (
            "lift the nominal and PSD-remainder Joseph/reset enclosure into the full H=18 source-correlated covariance cell, then replace the candidate full-word backend's entrywise J_aw accelerometer operation; only after H closes extend the same operation to A=21 with the active-bias projection guard"
        ),
        "failures": list(dict.fromkeys(failures)),
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for k in (
        "source_generated_not_trajectory_fit",
        "same_source_alignment_force_family_as_signed_stage",
        "PSD_matrix_order_diagonal_absorption_used",
        "PSD_offdiagonal_abs_le_eps_over_2_used",
        "PSD_zero_diagonal_operator_le_eps_used",
        "unequal_canonical_tangent_column_norms_used",
        "yaw_axis_cross_component_retained",
        "aw_gain_rows_retained_in_deltaK",
        "Joseph_cross_terms_retained", "Joseph_quadratic_term_retained",
        "shipping_reset_congruence_used",
    ):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in (
        "source_replay_used", "filter_changed", "FULL_18_STATE_COVARIANCE_PROPAGATED_HERE",
        "A21_MODE_PROPAGATED_HERE", "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE",
        "P4_USABLE_CERTIFICATE_PROMOTED",
    ):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    if float(d.get("candidate_angle_deg", -1.0)) != 30.0:
        f.append("candidate angle changed")
    if float(d.get("declared_startup_aw_error_fraction_g", -1.0)) != 0.3:
        f.append("startup a_w domain changed")
    if d.get("P4_30DEG_FIRST_ACCEL_JOSEPH_RESET_PSD_REMAINDER") == "PASS":
        if d.get("signed_30deg_prerequisite_status") != "PASS":
            f.append("Joseph/reset PASS without signed prerequisite")
        m = d.get("max_post_reset_PSD_remainder_operator_upper")
        if not (isinstance(m, (int, float)) and math.isfinite(float(m)) and float(m) >= 0.0):
            f.append("Joseph/reset PASS without finite remainder bound")
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
        "status": d["P4_30DEG_FIRST_ACCEL_JOSEPH_RESET_PSD_REMAINDER"],
        "signed": d["signed_30deg_prerequisite_status"],
        "dmax": d["signed_max_correction_norm_upper_rad"],
        "cells": d["evaluated_source_alignment_force_cells"],
        "max_post_reset_PSD": d["max_post_reset_PSD_remainder_operator_upper"],
        "worst": d["worst_cell"],
        "validation_failures": vf,
        "next": d["next_obligation"],
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
