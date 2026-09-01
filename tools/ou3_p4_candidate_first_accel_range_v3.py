#!/usr/bin/env python3
"""V3 candidate first-accelerometer range with tangent-only PSD resolvent.

V2 fixes the A-mode structural proof, but the underlying rank-two gain helper
still perturbs the small attitude-covariance remainder E using the scalar floor
P_aw+R_acc.  That is the axial innovation channel.  Since

    H_theta = -[m e3]_x,

H_theta E H_theta^T has an exact zero axial row/column.  The perturbation acts
only on the two force-tangent innovation channels, whose nominal denominators
already contain m^2 P_theta.  This is the same algebraic refinement certified
by the existing P5 V12D producer, generalized here to the candidate force and
yaw-alignment cells.

No source set, filter gain, 6-rad helper range, or P4 angle is changed.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p4_candidate_first_accel_range_v2 as V2
import ou3_p5_first_accel_structured_gain as SG

DEFAULT_DOMAIN = V2.DEFAULT_DOMAIN
SCHEMA = 3


def _tangent_structured_gain_bounds(*, tilt: float, yaw: float, eps: float,
                                    x: Interval, m: Interval, paw: Interval,
                                    racc_var: Interval):
    F = SG.FULL1
    I = SG.I
    t = I(tilt)
    delta = Interval.outward_bounds(yaw - tilt, yaw - tilt)
    pu = t + delta * x
    lam = paw + racc_var
    if lam.lo <= 0.0:
        raise RuntimeError("candidate first-accel lambda floor is nonpositive")
    m2 = m.square()
    den_perp = m2 * t + lam
    den_u = m2 * pu + lam
    if den_perp.lo <= 0.0 or den_u.lo <= 0.0:
        raise RuntimeError("candidate tangent gain denominator is nonpositive")

    geom = Interval(0.0, SG._sqrt_x1mx_upper(x))
    g_perp = m * t / den_perp
    g_u = m * pu / den_u
    g_z = m * delta * geom / den_u
    k0 = F.up(max(
        g_perp.hi,
        math.sqrt(F.up(g_u.hi * g_u.hi + g_z.hi * g_z.hi)),
    ))

    kh_perp = m2 * t / den_perp
    kh_u = m2 * pu / den_u
    kh_z = m2 * delta * geom / den_u
    kh0 = F.up(max(
        kh_perp.hi,
        math.sqrt(F.up(kh_u.hi * kh_u.hi + kh_z.hi * kh_z.hi)),
    ))

    # V12D tangent-channel resolvent.  Use 2*eps, matching the established
    # off-diagonal/remainder operator enclosure in that proof stage.
    eoff = F.up(2.0 * eps)
    mhi = m.hi
    dS = F.up(mhi * mhi * eoff)
    tangent_nominal_floor = min(den_perp.lo, den_u.lo)
    tangent_floor = F.down(tangent_nominal_floor - dS)
    if tangent_floor <= 0.0:
        raise RuntimeError("candidate perturbed tangent innovation floor lost positivity")
    inv_tangent = F.up(1.0 / tangent_floor)
    dCtheta = F.up(mhi * eoff)
    dk = F.up(F.up(dCtheta + F.up(k0 * dS)) * inv_tangent)
    k = F.up(k0 + dk)
    kh = F.up(kh0 + F.up(dk * mhi))
    return k, kh, {
        "lambda_lower": lam.lo,
        "p_u": pu.as_list(),
        "g_perp_upper": g_perp.hi,
        "g_u_upper": g_u.hi,
        "g_z_upper": g_z.hi,
        "K0_norm_upper": k0,
        "KH0_norm_upper": kh0,
        "attitude_PSD_remainder_operator_upper": eoff,
        "PSD_innovation_perturbation_tangent_only": True,
        "PSD_innovation_axial_row_column_exact_zero": True,
        "nominal_tangent_innovation_lower": tangent_nominal_floor,
        "tangent_innovation_perturbation_upper": dS,
        "perturbed_tangent_innovation_lower": tangent_floor,
        "perturbed_tangent_inverse_operator_upper": inv_tangent,
        "PSD_remainder_K_perturbation_upper": dk,
        "retired_axial_lambda_inverse_for_PSD_remainder": True,
    }


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 2,
          alignment_pieces: int = 16, force_magnitude_pieces: int = 4) -> dict:
    old = SG._structured_gain_bounds
    try:
        SG._structured_gain_bounds = _tangent_structured_gain_bounds
        out = dict(V2.build(
            Path(domain_path).resolve(),
            source_pieces=source_pieces,
            alignment_pieces=alignment_pieces,
            force_magnitude_pieces=force_magnitude_pieces,
        ))
    finally:
        SG._structured_gain_bounds = old
    out["schema"] = SCHEMA
    out["qualification"] = "OU3_P4_CANDIDATE_FIRST_ACCEL_RANGE_TANGENT_PSD_RESOLVENT"
    out["V12D_tangent_PSD_resolvent_used"] = True
    out["PSD_remainder_axial_noise_floor_inverse_used"] = False
    return out


def validate(d: dict) -> list[str]:
    base = dict(d)
    base["schema"] = V2.SCHEMA
    failures = V2.validate(base)
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("V12D_tangent_PSD_resolvent_used") is not True:
        failures.append("V12D tangent PSD resolvent is not active")
    if d.get("PSD_remainder_axial_noise_floor_inverse_used") is not False:
        failures.append("retired axial noise-floor inverse is still active")
    return list(dict.fromkeys(failures))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces", type=int, default=2)
    ap.add_argument("--alignment-pieces", type=int, default=16)
    ap.add_argument("--force-magnitude-pieces", type=int, default=4)
    ap.add_argument("--output", type=Path, required=True)
    x = ap.parse_args()
    d = build(x.domain.resolve(), source_pieces=x.source_pieces,
              alignment_pieces=x.alignment_pieces,
              force_magnitude_pieces=x.force_magnitude_pieces)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    x.output.parent.mkdir(parents=True, exist_ok=True)
    x.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P4_CANDIDATE_FIRST_ACCEL_RANGE_CERTIFICATE"],
        "widest_safe_deg": d["widest_candidate_first_accel_range_safe_deg"],
        "rows": [{
            "angle_deg": r["angle_deg"],
            "q_pred": r["post_prediction_q_upper"],
            "max_Ktheta": r["max_Ktheta_norm_upper"],
            "max_residual": r["max_combined_residual_norm_upper_mps2"],
            "max_d": r["max_first_accelerometer_correction_norm_upper_rad"],
            "margin": r["minimum_correction_range_margin_rad"],
            "safe": r["first_accelerometer_range_safe"],
            "first_unclosed": r["first_unclosed_child"],
        } for r in d["candidate_rows"]],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
