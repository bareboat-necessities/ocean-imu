#!/usr/bin/env python3
"""V13 signed radial subcells for the OU-III P5 sample-1 correction.

V10/V12 certify the sample-1 correction norm but do not retain a signed
three-vector suitable for the winding-aware deployed-quaternion/Cayley V2
primitive.  This producer restores only the directional information that is
actually present in the analytic one-plus-two innovation structure.

In V7's un-Rx proof gauge the nominal attitude gain is exactly block diagonal:

    r_x  -> (d_y,d_z),       r_yz -> d_x.

For one V10 source child, split the signed nominal residual coordinate r_x and
cover the scalar nominal d_x coordinate by signed intervals.  The exact
(d_y,d_z) gain components are reconstructed from the same positive 1+2 block
identities used by V7.  The y/z residual allocation is bounded by the residual
ball after fixing an r_x subcell, so the scalar d_x cover remains source
complete without counting a second independent residual ball.

V12 bounds all omitted PSD/S effects by an additive correction perturbation
norm.  V13 adds that perturbation as a symmetric component box, while retaining
V12's sharper global correction-norm upper bound.  Therefore every real
correction belongs to at least one component box and also obeys the independent
V12 radial upper bound.

The winding-aware V2 primitive needs a positive radial lower bound only when the
correction norm upper exceeds the deployed six-radian proof range.  V13 uses the
component-box guaranteed norm lower for that purpose.  A >6-rad subcell whose
component box still reaches the origin fails closed and requests finer signed
subdivision; no correction limit is widened.

This stage is conditional on V12 PASS.  It does not compose the Cayley state,
claim q<8, promote sample 1 or a complete word, or set N_H_words.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p5_sample1_structured_full_gain_v7 as V7
import ou3_p5_sample1_structured_full_gain_v12 as V12
import ou3_p5_deployed_quaternion_cayley_cell_v2 as CAYLEY2

DEFAULT_DOMAIN = V12.DEFAULT_DOMAIN
SCHEMA = 13
OLD_RANGE = 6.0
EXTENDED_RANGE = 9.0
FULL = V12.FULL


def _norm_upper(v) -> float:
    s = 0.0
    for x in v:
        a = x.abs_upper()
        s = FULL.up(s + FULL.up(a * a))
    return FULL.up(math.sqrt(max(0.0, s)))


def _minimum_abs(x: Interval) -> float:
    if x.lo <= 0.0 <= x.hi:
        return 0.0
    return min(abs(x.lo), abs(x.hi))


def _signed_gain_components(*, a:Interval, Y:Interval, c0:Interval,
                            alpha:Interval, qaw:Interval, b:Interval,
                            bz:Interval, det_first:Interval, d:Interval,
                            fy:Interval, fz:Interval, r:Interval):
    """Exact interval components of the V7 one-plus-two attitude gain."""
    h = FULL.I(0.5) * d
    Bt = alpha.square() * b + qaw
    Bz = alpha.square() * bz + qaw
    delta = alpha.square() * det_first + a * qaw
    if not (a.lo > 0.0 and Y.lo > 0.0 and Bz.lo > 0.0 and delta.lo > 0.0 and r.lo > 0.0):
        raise RuntimeError("positive V13 covariance/noise floors required")

    U = fz - h * fy
    V = -(h * fz) - fy
    cu = alpha * c0
    Nu = a * U + cu
    Sx = Nu.square() / a + delta / a + Y * V.square() + r
    if not Sx.lo > 0.0:
        raise RuntimeError("V13 scalar innovation lost positive floor")

    # Cov(theta_yz, y_x) = [Nu-hYV, hNu+YV].
    gy = (Nu - h * Y * V) / Sx
    gz = (h * Nu + Y * V) / Sx

    cx = -(alpha * c0)
    q = cx - a * fz
    A = delta + a * r
    det = fy.square() * A + (q.square() + A) * (Bz + r) / a
    if not (A.lo > 0.0 and det.lo > 0.0):
        raise RuntimeError("V13 two-by-two innovation lost positivity")
    kxy = q * (Bz + r) / det
    kxz = fy * A / det
    return (gy, gz), (kxy, kxz), {
        "scalar_Sx": Sx.as_list(),
        "two_by_two_det": det.as_list(),
        "perpendicular_gain_components": [gy.as_list(), gz.as_list()],
        "parallel_gain_components": [kxy.as_list(), kxz.as_list()],
    }


def _prerequisite_failure(v12:dict, failures:list[str]) -> dict:
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_SIGNED_RADIAL_SUBCELLS_V13",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "V12_actual_innovation_PSD_S_prerequisite_required": True,
        "V12_prerequisite_passed": False,
        "signed_one_plus_two_correction_cover_used": False,
        "V12_correction_perturbation_ball_retained": False,
        "radial_lower_bound_required_above_6_rad": True,
        "deployed_correction_limit_rad": 6.0,
        "deployed_correction_limit_increased": False,
        "signed_cayley_q8_composed_here": False,
        "complete_sample1_branch_closed_here": False,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "evaluated_source_rows": 0,
        "evaluated_signed_subcells": 0,
        "above_6rad_subcells": 0,
        "unclosed_radial_subcells": 0,
        "max_radial_upper": 0.0,
        "minimum_radial_lower_above_6": None,
        "first_unclosed_radial_subcell": None,
        "P5_SAMPLE1_SIGNED_RADIAL_SUBCELLS_V13": "NOT_ESTABLISHED",
        "next_obligation": "ESTABLISH_V12_ACTUAL_INNOVATION_PSD_S_CLOSURE",
        "failures": failures,
    }


def build(domain_path:Path=DEFAULT_DOMAIN, *, source_pieces:int=4,
          source_cell_index:int=0, p_pieces:int=24,
          tangent_pieces:int=24, axial_pieces:int=24,
          residual_x_pieces:int=12, parallel_pieces:int=12) -> dict:
    if residual_x_pieces < 2 or parallel_pieces < 2:
        raise ValueError("signed subdivision counts must be >=2")
    path = Path(domain_path).resolve()
    v12 = V12.build(path, source_pieces=source_pieces,
                    source_cell_index=source_cell_index,
                    p_pieces=p_pieces, tangent_pieces=tangent_pieces,
                    axial_pieces=axial_pieces)
    failures = [f"V12: {x}" for x in V12.validate(v12)]
    if v12.get("P5_SAMPLE1_PSD_S_ACTUAL_INNOVATION_V12") != "PASS":
        failures.append("V12 prerequisite did not pass")
        return _prerequisite_failure(v12, failures)

    # Rebuild V10 once so V13 can recover the exact source coordinates and the
    # nominal two-block residual allocation.  V12 already certified that every
    # one of these source rows survives the PSD/S perturbation stage.
    core = V12.V11.V10.build(path, source_pieces=source_pieces,
                             source_cell_index=source_cell_index,
                             p_pieces=p_pieces, tangent_pieces=tangent_pieces,
                             axial_pieces=axial_pieces)
    failures += [f"V10: {x}" for x in V12.V11.V10.validate(core)]
    if core.get("P5_SAMPLE1_COMBINED_PERPENDICULAR_RESIDUAL_V10") != "PASS":
        failures.append("V10 prerequisite did not pass")

    dom = json.loads(path.read_text(encoding="utf-8"))
    first = V12.V11.FIRST.build(path, source_pieces=source_pieces)
    vec = V12.V11.VECTOR.build()
    failures += [f"first: {x}" for x in V12.V11.FIRST.validate(first)]
    failures += [f"vector: {x}" for x in V12.V11.VECTOR.validate(vec)]
    src, phase = V12.V11.RG._source_phase_children(source_pieces)[source_cell_index]
    if phase != "due": failures.append("V13 focused family requires first due source cell")

    fr = first["source_cells"][source_cell_index]
    p_all = Interval.outward_bounds(*map(float, fr["P_aw_variance_interval"]))
    pcells = V12.V11.SUB.parts(p_all.lo, p_all.hi, p_pieces)
    hstep = float(src["dt_s"]); g = float(dom["startup"]["gravity_mps2"])
    tilt, yaw, eps = V12.V11.RG._attitude_covariance_epsilon(path, hstep)
    t = Interval.outward_bounds(tilt, FULL.up(tilt + eps))
    Y = Interval.outward_bounds(yaw, FULL.up(yaw + eps))
    r = FULL._R_diag(float(vec["configured_measurement_bounds"]["acc_measurement_std_mps2"]))[0][0]
    F, Q, _ = FULL._transition_and_Q(src, dom)
    alpha = F[15][15]; qaw = Q[15][15]

    if len(core["rows"]) != len(v12["rows"]):
        failures.append("V10/V12 row counts differ")

    total = high = unclosed = 0
    max_hi = 0.0; min_high_lo = math.inf
    first_bad = None; worst = None
    samples = []

    for base, vr in zip(core["rows"], v12["rows"]):
        ids = ("p_cell", "tangent_residual_cell", "axial_residual_cell")
        if any(int(base[k]) != int(vr[k]) for k in ids):
            failures.append("V10/V12 row ordering mismatch")
            break
        if "V12_correction_norm_upper_rad" not in vr:
            failures.append("V12 PASS row lacks correction bound")
            break

        pi = int(base["p_cell"]); p = pcells[pi]
        rt = Interval.outward_bounds(*map(float, base["first_tangent_residual_magnitude_mps2"]))
        rz = Interval.outward_bounds(*map(float, base["first_axial_residual_mps2"]))
        d0 = Interval.outward_bounds(*map(float, base["first_attitude_correction_rad"]))
        D = FULL.I(g * g) * t + p + r
        a = t * (p + r) / D
        c0 = -(FULL.I(g) * t * p / D)
        b = p * (FULL.I(g * g) * t + r) / D
        bz = p * r / (p + r)
        det_first = t * p * r / D
        fy = -(alpha * (p / D) * rt)
        fz = FULL.I(g) + alpha * (p / (p + r)) * rz
        (gy, gz), (_kxy, _kxz), gain_detail = _signed_gain_components(
            a=a, Y=Y, c0=c0, alpha=alpha, qaw=qaw, b=b, bz=bz,
            det_first=det_first, d=d0, fy=fy, fz=fz, r=r)

        kperp = float(base["Ktheta_perpendicular_block_upper"])
        kpar = float(base["Ktheta_parallel_block_upper"])
        if _norm_upper((gy, gz)) > FULL.up(kperp + 1e-10):
            failures.append("signed perpendicular gain exceeds V10 norm parent")
            break
        if _norm_upper((_kxy, _kxz)) > FULL.up(kpar + 1e-10):
            failures.append("signed parallel gain exceeds V10 norm parent")
            break

        drho = float(vr["total_residual_perturbation_upper_mps2"])
        dk = float(vr["sample1_gain_operator_perturbation_upper"])
        rho = float(base["sample1_full_residual_norm_upper_mps2"])
        rho_x = min(rho, float(base["sample1_combined_source_x_residual_upper_mps2"]))
        k0 = max(kperp, kpar)
        eta = FULL.up(FULL.up(k0 * drho) + FULL.up(dk * FULL.up(rho + drho)))
        global_hi = float(vr["V12_correction_norm_upper_rad"])
        if FULL.up(float(base["combined_directional_correction_norm_upper_rad"]) + eta) > FULL.up(global_hi + 1e-9):
            failures.append("V13 perturbation reconstruction exceeds V12 parent")
            break

        rx_cells = V12.V11.SUB.parts(-rho_x, rho_x, residual_x_pieces)
        for rxc in rx_cells:
            rx_min = _minimum_abs(rxc)
            rem2 = max(0.0, FULL.up(rho * rho) - FULL.down(rx_min * rx_min))
            ryz_hi = FULL.up(math.sqrt(rem2))
            ux_hi = FULL.up(kpar * ryz_hi)
            u_cells = V12.V11.SUB.parts(-ux_hi, ux_hi, parallel_pieces)
            for uc in u_cells:
                e = Interval.outward_bounds(-eta, eta)
                dbox = [uc + e, gy * rxc + e, gz * rxc + e]
                component_box_hi = CAYLEY2.V1._norm_upper(dbox)
                component_box_lo = CAYLEY2._norm_lower(dbox)
                rx_abs = rxc.abs_upper(); u_abs = uc.abs_upper()
                nominal_hi = FULL.up(math.sqrt(FULL.up(FULL.up(u_abs * u_abs) + FULL.up((kperp * rx_abs) ** 2))))
                radial_hi = min(global_hi, FULL.up(nominal_hi + eta))
                radial_lo = min(component_box_lo, radial_hi)
                total += 1; max_hi = max(max_hi, radial_hi)
                ready = radial_hi <= OLD_RANGE
                if radial_hi > OLD_RANGE:
                    high += 1
                    min_high_lo = min(min_high_lo, radial_lo)
                    ready = radial_lo > 0.0 and radial_hi <= EXTENDED_RANGE
                row = {
                    "p_cell": pi,
                    "tangent_residual_cell": int(base["tangent_residual_cell"]),
                    "axial_residual_cell": int(base["axial_residual_cell"]),
                    "nominal_rx_subcell_mps2": rxc.as_list(),
                    "nominal_parallel_correction_subcell_rad": uc.as_list(),
                    "correction_component_box_rad": [x.as_list() for x in dbox],
                    "correction_perturbation_norm_upper_rad": eta,
                    "component_box_norm_lower_rad": component_box_lo,
                    "component_box_norm_upper_rad": component_box_hi,
                    "certified_radial_lower_rad": radial_lo,
                    "certified_radial_upper_rad": radial_hi,
                    "above_6rad_requires_positive_radial_lower": radial_hi > OLD_RANGE,
                    "radial_ready_for_deployed_cayley_v2": ready,
                    "gain_detail": gain_detail,
                }
                if len(samples) < 24 and radial_hi > OLD_RANGE:
                    samples.append(row)
                if worst is None or radial_hi > worst["certified_radial_upper_rad"]:
                    worst = row
                if not ready:
                    unclosed += 1
                    if first_bad is None: first_bad = row

    ok = (not failures and total > 0 and unclosed == 0 and first_bad is None)
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_SIGNED_RADIAL_SUBCELLS_V13",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "V12_actual_innovation_PSD_S_prerequisite_required": True,
        "V12_prerequisite_passed": v12.get("P5_SAMPLE1_PSD_S_ACTUAL_INNOVATION_V12") == "PASS",
        "signed_one_plus_two_correction_cover_used": True,
        "residual_ball_not_double_counted_across_blocks": True,
        "V12_correction_perturbation_ball_retained": True,
        "V12_global_radial_upper_intersected_with_component_cover": True,
        "radial_lower_bound_required_above_6_rad": True,
        "deployed_correction_limit_rad": 6.0,
        "deployed_correction_limit_increased": False,
        "extended_cayley_proof_primitive_limit_rad": 9.0,
        "signed_cayley_q8_composed_here": False,
        "complete_sample1_branch_closed_here": False,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "residual_x_pieces": residual_x_pieces,
        "parallel_pieces": parallel_pieces,
        "evaluated_source_rows": len(core["rows"]),
        "evaluated_signed_subcells": total,
        "above_6rad_subcells": high,
        "unclosed_radial_subcells": unclosed,
        "max_radial_upper": max_hi,
        "minimum_radial_lower_above_6": None if high == 0 else min_high_lo,
        "first_unclosed_radial_subcell": first_bad,
        "worst_radial_subcell": worst,
        "sample_above_6rad_subcells": samples,
        "P5_SAMPLE1_SIGNED_RADIAL_SUBCELLS_V13": "PASS" if ok else "NOT_ESTABLISHED",
        "next_obligation": (
            "COMPOSE_SIGNED_SAMPLE1_CORRECTION_WITH_CURRENT_CAYLEY_AND_REQUIRE_Q_LT_8"
            if ok else "REFINE_SIGNED_RX_PARALLEL_SUBDIVISION_OR_ESTABLISH_V12"
        ),
        "failures": failures,
    }


def validate(d:dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA: f.append("schema mismatch")
    for k in ("source_generated_not_trajectory_fit", "V12_actual_innovation_PSD_S_prerequisite_required",
              "radial_lower_bound_required_above_6_rad"):
        if d.get(k) is not True: f.append(f"{k} is not true")
    for k in ("source_replay_used", "filter_changed", "deployed_correction_limit_increased",
              "signed_cayley_q8_composed_here", "complete_sample1_branch_closed_here",
              "q8_word_promoted_here", "whole_word_promoted_here", "N_H_words_set_here"):
        if d.get(k) is not False: f.append(f"{k} is not false")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        f.append("deployed correction limit changed")
    prereq = d.get("V12_prerequisite_passed") is True
    st = d.get("P5_SAMPLE1_SIGNED_RADIAL_SUBCELLS_V13")
    if prereq:
        for k in ("signed_one_plus_two_correction_cover_used", "residual_ball_not_double_counted_across_blocks",
                  "V12_correction_perturbation_ball_retained", "V12_global_radial_upper_intersected_with_component_cover"):
            if d.get(k) is not True: f.append(f"{k} is not true")
        if int(d.get("evaluated_signed_subcells", 0)) <= 0: f.append("no signed V13 subcells")
        if st == "PASS" and int(d.get("unclosed_radial_subcells", -1)) != 0:
            f.append("V13 PASS retains unclosed radial subcell")
        if st == "PASS" and float(d.get("max_radial_upper", math.inf)) > EXTENDED_RANGE:
            f.append("V13 PASS exceeds extended Cayley primitive range")
        if st == "NOT_ESTABLISHED" and d.get("first_unclosed_radial_subcell") is None and not f:
            f.append("missing V13 radial witness")
    else:
        if st != "NOT_ESTABLISHED": f.append("V13 passed without V12 prerequisite")
        if not any("V12 prerequisite did not pass" in x for x in f):
            f.append("missing V12 prerequisite failure")
    return f


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces", type=int, default=4)
    ap.add_argument("--source-cell-index", type=int, default=0)
    ap.add_argument("--p-pieces", type=int, default=24)
    ap.add_argument("--tangent-pieces", type=int, default=24)
    ap.add_argument("--axial-pieces", type=int, default=24)
    ap.add_argument("--residual-x-pieces", type=int, default=12)
    ap.add_argument("--parallel-pieces", type=int, default=12)
    ap.add_argument("--output", type=Path, required=True)
    x = ap.parse_args()
    d = build(x.domain, source_pieces=x.source_pieces, source_cell_index=x.source_cell_index,
              p_pieces=x.p_pieces, tangent_pieces=x.tangent_pieces,
              axial_pieces=x.axial_pieces, residual_x_pieces=x.residual_x_pieces,
              parallel_pieces=x.parallel_pieces)
    vf = validate(d); d["validation_failures"] = vf
    x.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P5_SAMPLE1_SIGNED_RADIAL_SUBCELLS_V13"],
        "V12_passed": d["V12_prerequisite_passed"],
        "source_rows": d["evaluated_source_rows"],
        "signed_subcells": d["evaluated_signed_subcells"],
        "above_6": d["above_6rad_subcells"],
        "unclosed": d["unclosed_radial_subcells"],
        "max_radial_upper": d["max_radial_upper"],
        "min_radial_lower_above_6": d["minimum_radial_lower_above_6"],
        "first_unclosed": d["first_unclosed_radial_subcell"],
        "worst": d.get("worst_radial_subcell"),
        "next": d["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2

if __name__ == "__main__":
    raise SystemExit(main())
