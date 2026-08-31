#!/usr/bin/env python3
"""Source-faithful finite-angle P4 candidate full-word prefix backend.

This is the first heavy numerical consumer of the entrance/search domain added
by PR #445.  The older P5 full-H diagnostic starts from one Cartesian cube and
uses q<=8.  That is useful for finding large-angle numerical obstructions, but
it is not the P4 search requested by the theorem architecture.

This producer instead:

* consumes the declared P4 candidate ladder (30,25,20,15 deg);
* covers the actual Cayley norm ball for each candidate with an outward union
  of boxes instead of replacing ||c||<=q by the much larger [-q,q]^3 cube;
* uses the operation-matched outer sector only as a prefix-safety chart;
* replaces the old |delta p_i|<=20 m Cartesian initialization by the declared
  P5 entrance |delta p_i|<=0.5 Hs with the declared Hs proof ceiling;
* carries the shipping 18-state H mode and a source-faithful 21-state A mode;
* retains full S->attitude cross gain, Joseph covariance update, immediate
  reset congruence, exact deployed quaternion composition, accelerometer
  effective-a_w reduction, magnetometer radial annihilation, and all
  accepted/rejected identity hulls used by the active H backend.

The A mode adds the shipping residual accelerometer-bias OU block and J_ba=I.
It fails closed if the deterministic A error enclosure reaches the 0.5 m/s^2
projection surface; this producer does not linearize through the projection.

This stage proves finite prefix safety/closure of the numerical word.  It does
NOT relabel prefix safety as the final P4 Lyapunov-dissipation certificate.  A
candidate can be promoted only after the complete H/A word dissipation and the
normalized translation/nontranslation cross-block obligation also close.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

from ou3_interval import Interval
import ou3_full_process_ucc as PROCESS
import ou3_p4_operation_matched_sector_certificate as SECTOR
import ou3_p4_p5_entrance_search_domain as ENTRANCE
import ou3_p5_full_h_prefix_cells as H
import ou3_p5_full_h_prefix_cells_v2 as H2
import ou3_p5_full_h_prefix_cells_v3 as H3
import ou3_validated_transcendentals as VT

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def I(x: float) -> Interval:
    return Interval.point(float(x))


def _box(a: float) -> Interval:
    a = abs(float(a))
    return Interval(down(-a), up(a))


def _norm_bounds_box(box) -> tuple[float, float]:
    lo2 = 0.0
    hi2 = 0.0
    for a, b in box:
        if a <= 0.0 <= b:
            mn = 0.0
        else:
            mn = min(abs(a), abs(b))
        mx = max(abs(a), abs(b))
        lo2 = down(lo2 + down(mn * mn))
        hi2 = up(hi2 + up(mx * mx))
    return down(math.sqrt(max(0.0, lo2))), up(math.sqrt(max(0.0, hi2)))


def _ball_box_cover(q: float, max_box_norm_factor: float = 1.5) -> list[list[tuple[float, float]]]:
    """Outward box cover of ||c||<=q with bounded box-corner inflation.

    Boxes wholly outside the ball are discarded.  Boundary boxes are bisected
    until their farthest corner is <= factor*q.  The union therefore contains
    the whole ball, while avoiding the sqrt(3) inflation of one global cube.
    """
    q = float(q)
    factor = float(max_box_norm_factor)
    if not (q > 0.0 and 1.0 < factor < math.sqrt(3.0)):
        raise ValueError("invalid Cayley ball-cover parameters")
    todo = [([(-q, q), (-q, q), (-q, q)], 0)]
    out = []
    while todo:
        box, depth = todo.pop()
        qmin, qmax = _norm_bounds_box(box)
        if qmin > q:
            continue
        if qmax <= factor * q:
            out.append([(down(a), up(b)) for a, b in box])
            continue
        if depth >= 18:
            raise RuntimeError("Cayley ball cover failed to reach requested inflation")
        widths = [b - a for a, b in box]
        axis = max(range(3), key=lambda i: widths[i])
        a, b = box[axis]
        m = 0.5 * (a + b)
        left = list(box); right = list(box)
        left[axis] = (a, m); right[axis] = (m, b)
        todo.append((right, depth + 1))
        todo.append((left, depth + 1))
    out.sort(key=lambda b: _norm_bounds_box(b)[1], reverse=True)
    return out


def _configure_mode(mode: str) -> None:
    if mode not in ("H", "A"):
        raise ValueError("mode must be H or A")
    H.N = 18 if mode == "H" else 21
    H.TH = range(0, 3)
    H.BG = range(3, 6)
    H.V = range(6, 9)
    H.P = range(9, 12)
    H.SS = range(12, 15)
    H.AW = range(15, 18)
    H.BA = range(18, 21) if mode == "A" else range(18, 18)


def _source_sigma_bacc0() -> float:
    text = H.MEKF.read_text(encoding="utf-8")
    m = re.search(r"sigma_bacc0_\s*=\s*T\(([0-9.eE+-]+)\)", text)
    if not m:
        raise RuntimeError("cannot extract source sigma_bacc0")
    return float(m.group(1))


def _transition_and_Q(mode: str, src: dict, domain: dict):
    F, Q, Rstep = H2._tight_transition_and_Q(src, domain)
    if mode == "H":
        return F, Q, Rstep, None

    h = float(src["dt_s"])
    pc = PROCESS.build()["source_constants"]
    tau_b = float(pc["accel_bias_tau_s"])
    q_b = float(pc["accel_bias_process_variance_density"])
    if not (tau_b > 0.0 and q_b > 0.0):
        raise RuntimeError("active accelerometer-bias source constants lost positivity")

    x1 = Interval.outward_bounds(h / tau_b, h / tau_b)
    phi = VT.exp_interval(-x1)
    x2 = Interval.outward_bounds(2.0 * h / tau_b, 2.0 * h / tau_b)
    em1 = VT.expm1_interval(-x2)
    qd_scale = Interval.outward_bounds(-0.5 * tau_b, -0.5 * tau_b) * em1
    qd = Interval.outward_bounds(q_b, q_b) * qd_scale
    for i in H.BA:
        F[i][i] = Interval(phi.lo, min(1.0, phi.hi))
        Q[i][i] = qd
    return F, H._psd_tighten(Q), Rstep, {
        "phi_interval": phi.as_list(),
        "Qd_variance_interval": qd.as_list(),
        "tau_bacc_s": tau_b,
        "process_variance_density": q_b,
    }


def _initial_covariance(mode: str, src: dict, domain_path: Path):
    Pm = H2._corrected_initial_covariance(src, domain_path)
    if mode == "A":
        s = _source_sigma_bacc0()
        v = Interval.outward_bounds(s * s, s * s)
        for i in H.BA:
            Pm[i][i] = v
            for j in range(H.N):
                if j != i:
                    Pm[i][j] = I(0.0)
                    Pm[j][i] = I(0.0)
    return H._psd_tighten(Pm)


def _initial_error(mode: str, domain: dict) -> tuple[list[Interval], float | None, dict]:
    old = domain["startup"]["physical_handoff_coordinate_bounds"]
    pos = domain["initial_filter_entrance"]["position"]
    hs = float(pos["significant_wave_height_Hs_upper_m"])
    p_factor = float(pos["component_abs_error_upper_Hs_factor"])
    if not (hs > 0.0 and p_factor == 0.5):
        raise RuntimeError("finite 0.5 Hs entrance box not declared")
    p_component = up(p_factor * hs)

    e = [I(0.0) for _ in range(H.N)]
    for idxs, key in (
        (H.BG, "gyro_bias_error_norm_upper_rad_s"),
        (H.V, "velocity_error_norm_upper_mps"),
        (H.SS, "integral_displacement_error_norm_upper_m_s"),
        (H.AW, "latent_acceleration_error_norm_upper_mps2"),
    ):
        a = float(old[key])
        for i in idxs:
            e[i] = _box(a)
    for i in H.P:
        e[i] = _box(p_component)

    ba_cap = None
    if mode == "A":
        ba_cap = float(domain["normal_live"]["active_accelerometer_bias_state_norm_upper_mps2"])
        if not 0.0 < ba_cap < float(domain["normal_live"]["active_accelerometer_bias_projection_limit_mps2"]):
            raise RuntimeError("A entrance bias ball is not strictly inside projection ball")
        # Component intervals are a conservative projection of the norm ball.
        # A separate norm cap is propagated below to avoid turning the 0.45 ball
        # into a fictitious sqrt(3)*0.45 projection-boundary violation.
        for i in H.BA:
            e[i] = _box(ba_cap)

    return e, ba_cap, {
        "Hs_upper_m": hs,
        "position_component_abs_upper_m": p_component,
        "legacy_P1_position_norm_upper_m_not_used_as_P4_entry": float(old["position_error_norm_upper_m"]),
    }


def _predict_error(mode: str, e, F, ba_cap):
    out = H._predict_error(e, F)
    if mode == "A":
        phi_hi = max(F[i][i].abs_upper() for i in H.BA)
        old_ba = list(e[i] for i in H.BA)
        for row, i in enumerate(H.BA):
            out[i] = F[i][i] * old_ba[row]
        ba_cap = up(float(ba_cap) * phi_hi)
        for i in H.BA:
            out[i] = H._intersect(out[i], Interval(-ba_cap, ba_cap))
    return out, ba_cap


def _H_acc(mode: str, domain: dict):
    M = H._H_acc(domain)
    if mode == "A":
        for ax, i in enumerate(H.BA):
            M[ax][i] = I(1.0)
    return M


def _acc_residual(mode: str, e, c, domain: dict, q_hi: float):
    if mode == "H":
        return H._acc_residual(e, c, domain, q_hi)
    M = _H_acc(mode, domain)
    fhi = float(domain["normal_live"]["specific_force_norm_upper_mps2"])
    aw_hi = max(e[i].abs_upper() for i in H.AW)
    eta = up(
        H.VEFF.accel_attitude_eta_per_vector_norm_upper(q_hi) * fhi
        + H.VEFF.accel_latent_cross_gain_upper(q_hi) * aw_hi
    )
    z = [I(0.0) for _ in range(H.N)]
    for i in range(3):
        z[i] = c[i]
    for i in H.AW:
        z[i] = e[i] + _box(eta)
    for i in H.BA:
        z[i] = e[i]
    return M, H._mat_vec(M, z), eta


def _update_ba_cap(mode: str, e, cell, ba_cap: float | None, projection_limit: float):
    if mode == "H":
        return e, ba_cap
    dx = [cell["dx"][i] for i in H.BA]
    dnorm = H._norm_upper(dx)
    ba_cap = up(float(ba_cap) + dnorm)
    if not ba_cap < projection_limit:
        raise RuntimeError(
            f"active accelerometer-bias enclosure reaches projection surface: {ba_cap} >= {projection_limit}"
        )
    for i in H.BA:
        e[i] = H._intersect(e[i], Interval(-ba_cap, ba_cap))
    return e, ba_cap


def _run_cell(mode: str, domain_path: Path, domain: dict, cbox, candidate_q: float,
              outer_q: float, max_samples: int | None) -> dict:
    _configure_mode(mode)
    src = H._source_cell()
    F, Q, Rstep, ba_process = _transition_and_Q(mode, src, domain)
    Pm = _initial_covariance(mode, src, domain_path)
    e, ba_cap, position = _initial_error(mode, domain)
    c = [Interval(a, b) for a, b in cbox]
    h = float(src["dt_s"])
    Tword = float(domain["normal_live"]["vector_pe_recurrence_window_s"])
    full_samples = int(math.ceil(Tword / h)) + 2
    samples = full_samples if max_samples is None else min(full_samples, max(0, int(max_samples)))
    vc = H.VECTOR.build()["configured_measurement_bounds"]
    Racc = H._R_diag(float(vc["acc_measurement_std_mps2"]))
    Rmag = H._R_diag(float(vc["mag_measurement_std_uT"]))
    RS = H._R_S(src)
    projection_limit = float(domain["normal_live"]["active_accelerometer_bias_projection_limit_mps2"])

    max_q = H._norm_upper(c)
    first_failure = None
    inverse_counts = {"FIXED_PIVOT_INTERVAL_GAUSS_JORDAN": 0, "SPD_S_GE_R_SPECTRAL_ENTRY_ENCLOSURE": 0}
    first_S_done = False
    last_cells = {}

    for k in range(samples):
        try:
            Pm = H._psd_tighten(H.matrix_add(H.matrix_mul(H.matrix_mul(F, Pm), H.matrix_transpose(F)), Q))
            e, ba_cap = _predict_error(mode, e, F, ba_cap)
            c = H._predict_c(c, Rstep, domain, h)
            qnow = H._norm_upper(c)
            max_q = max(max_q, qnow)
            if not qnow < outer_q:
                raise RuntimeError(f"prediction prefix leaves operation-matched outer sector: q={qnow} >= {outer_q}")

            HS = H._H_S()
            rS = [-e[12 + i] for i in range(3)]
            Pm, e, c, Scell, Ssigned = H._measurement_branch_hull(Pm, e, c, HS, RS, rS, allow_rejected=True)
            inverse_counts[Scell["inverse_backend"]] += 1
            e, ba_cap = _update_ba_cap(mode, e, Scell, ba_cap, projection_limit)
            first_S_done = True
            last_cells["S"] = {
                "sample": k,
                "d_norm_upper": H._norm_upper(H._vec_neg(Scell["dx"][0:3])),
                "signed_denominator": Ssigned["denominator"].as_list(),
            }
            qnow = H._norm_upper(c)
            max_q = max(max_q, qnow)
            if not qnow < outer_q:
                raise RuntimeError(f"S prefix leaves operation-matched outer sector: q={qnow} >= {outer_q}")

            Hacc, racc, eta = _acc_residual(mode, e, c, domain, qnow)
            Pm, e, c, Acell, Asigned = H._measurement_branch_hull(Pm, e, c, Hacc, Racc, racc, allow_rejected=True)
            inverse_counts[Acell["inverse_backend"]] += 1
            e, ba_cap = _update_ba_cap(mode, e, Acell, ba_cap, projection_limit)
            last_cells["accelerometer"] = {
                "sample": k,
                "d_norm_upper": H._norm_upper(H._vec_neg(Acell["dx"][0:3])),
                "effective_aw_eta_norm_upper": eta,
                "signed_denominator": Asigned["denominator"].as_list(),
            }
            qnow = H._norm_upper(c)
            max_q = max(max_q, qnow)
            if not qnow < outer_q:
                raise RuntimeError(f"accelerometer prefix leaves operation-matched outer sector: q={qnow} >= {outer_q}")

            Hmag, rmag, deff = H._mag_residual(c, domain, 0.0, qnow)
            Pm, e, c, Mcell, Msigned = H._measurement_branch_hull(Pm, e, c, Hmag, Rmag, rmag, allow_rejected=True)
            inverse_counts[Mcell["inverse_backend"]] += 1
            e, ba_cap = _update_ba_cap(mode, e, Mcell, ba_cap, projection_limit)
            last_cells["magnetometer"] = {
                "sample": k,
                "d_norm_upper": H._norm_upper(H._vec_neg(Mcell["dx"][0:3])),
                "signed_denominator": Msigned["denominator"].as_list(),
            }
            qnow = H._norm_upper(c)
            max_q = max(max_q, qnow)
            if not qnow < outer_q:
                raise RuntimeError(f"magnetometer prefix leaves operation-matched outer sector: q={qnow} >= {outer_q}")

            if (k + 1) % max(1, int(round(0.1 / h))) == 0:
                u = max(src["sigma_aw_mps2"].square().hi, max(Pm[i][i].hi for i in H.AW))
                for i in H.AW:
                    Pm[i][i] = Interval(0.0, up(u))
                    for j in H.AW:
                        if i != j:
                            Pm[i][j] = H._intersect(Pm[i][j], Interval(-up(u), up(u)))
                Pm = H._psd_tighten(Pm)
        except Exception as exc:
            first_failure = {
                "sample": k,
                "reason": f"{type(exc).__name__}: {exc}",
                "cayley_norm_upper_before_failure": H._norm_upper(c),
                "ba_error_norm_cap_before_failure": ba_cap,
                "covariance": H._matrix_summary(Pm),
                "last_cells": last_cells,
            }
            break

    q_final = H._norm_upper(c)
    prefix_safe = first_failure is None and (samples == 0 or first_S_done) and max_q < outer_q
    returned_to_candidate = prefix_safe and q_final <= candidate_q
    return {
        "mode": mode,
        "dimension": H.N,
        "samples_run": samples,
        "full_word_samples": full_samples,
        "full_word_completed": samples == full_samples and prefix_safe,
        "prefix_safe_in_outer_sector": prefix_safe,
        "returned_to_candidate_cayley_ball": returned_to_candidate,
        "candidate_q_upper": candidate_q,
        "outer_q_upper": outer_q,
        "initial_box_q_range": list(_norm_bounds_box(cbox)),
        "max_prefix_q_upper": max_q,
        "final_q_upper": q_final,
        "position_entrance": position,
        "active_bias_process": ba_process,
        "active_bias_final_norm_cap": ba_cap,
        "inverse_backend_counts": inverse_counts,
        "first_failure": first_failure,
    }


def build(domain_path: Path = DEFAULT_DOMAIN, max_samples: int | None = None,
          candidate_index: int | None = None, cover_factor: float = 1.5) -> dict:
    domain_path = Path(domain_path).resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("candidate full-word domain must not be trajectory fitted")

    H3._install_backend()
    entrance = ENTRANCE.build(domain_path)
    ef = ENTRANCE.validate(entrance)
    sector = SECTOR.build(domain_path)
    sf = SECTOR.validate(sector)
    failures = [f"entrance: {x}" for x in ef] + [f"sector: {x}" for x in sf]
    outer_q = float(sector["design_cayley_norm_upper"])
    candidates = entrance["P4_complete_word_search"]["candidate_rows"]
    if candidate_index is not None:
        if not 0 <= int(candidate_index) < len(candidates):
            raise ValueError("candidate index out of range")
        candidates = [candidates[int(candidate_index)]]

    rows = []
    widest_prefix_safe = None
    for crow in candidates:
        q = float(crow["cayley_norm_upper"])
        cover = _ball_box_cover(q, cover_factor)
        cover_q = max(_norm_bounds_box(b)[1] for b in cover)
        if not cover_q < outer_q:
            failures.append(
                f"candidate {crow['angle_deg']} deg ball cover exceeds outer sector: {cover_q} >= {outer_q}"
            )
            rows.append({
                "angle_deg": crow["angle_deg"],
                "candidate_q_upper": q,
                "cover_cells": len(cover),
                "cover_q_upper": cover_q,
                "H": None,
                "A": None,
                "both_modes_prefix_safe": False,
                "first_unclosed_cell": "cover exceeds outer sector",
            })
            continue

        mode_results = {}
        first_unclosed = None
        for mode in ("H", "A"):
            cell_results = []
            for ci, cbox in enumerate(cover):
                r = _run_cell(mode, domain_path, domain, cbox, q, outer_q, max_samples)
                r["cell_index"] = ci
                cell_results.append(r)
                if not r["prefix_safe_in_outer_sector"]:
                    first_unclosed = {"mode": mode, "cell_index": ci, "result": r}
                    break
            mode_results[mode] = {
                "dimension": 18 if mode == "H" else 21,
                "evaluated_cells": len(cell_results),
                "all_evaluated_cells_prefix_safe": bool(cell_results) and all(x["prefix_safe_in_outer_sector"] for x in cell_results),
                "all_cover_cells_evaluated": len(cell_results) == len(cover),
                "all_cover_cells_prefix_safe": len(cell_results) == len(cover) and all(x["prefix_safe_in_outer_sector"] for x in cell_results),
                "all_cover_cells_return_to_candidate": len(cell_results) == len(cover) and all(x["returned_to_candidate_cayley_ball"] for x in cell_results),
                "worst_prefix_q_upper": max((x["max_prefix_q_upper"] for x in cell_results), default=None),
                "first_failure": next((x["first_failure"] for x in cell_results if x["first_failure"] is not None), None),
            }
            if first_unclosed is not None:
                break

        both = (
            mode_results.get("H", {}).get("all_cover_cells_prefix_safe") is True
            and mode_results.get("A", {}).get("all_cover_cells_prefix_safe") is True
        )
        row = {
            "angle_deg": crow["angle_deg"],
            "candidate_q_upper": q,
            "cover_factor": cover_factor,
            "cover_cells": len(cover),
            "cover_q_upper": cover_q,
            "H": mode_results.get("H"),
            "A": mode_results.get("A"),
            "both_modes_prefix_safe": both,
            "first_unclosed_cell": first_unclosed,
        }
        rows.append(row)
        if both and widest_prefix_safe is None:
            widest_prefix_safe = float(crow["angle_deg"])
            if max_samples is None:
                break

    full_requested = max_samples is None
    status = "PASS" if not failures and widest_prefix_safe is not None else "NOT_ESTABLISHED"
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_FINITE_ANGLE_CANDIDATE_FULL_HA_WORD_PREFIX",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "full_S_to_attitude_cross_gain_retained": True,
        "shipping_Joseph_and_immediate_reset_retained": True,
        "deployed_quaternion_composition_retained": True,
        "H_dimension": 18,
        "A_dimension": 21,
        "active_accelerometer_bias_J_ba_identity": True,
        "active_accelerometer_bias_projection_linearized_through_boundary": False,
        "p4_candidate_ball_covered_not_single_cube": True,
        "cover_is_outward_union": True,
        "cover_factor": cover_factor,
        "old_q8_chart_used": False,
        "operation_matched_outer_q_upper": outer_q,
        "position_entrance_uses_half_Hs": True,
        "Hs_upper_m": float(domain["initial_filter_entrance"]["position"]["significant_wave_height_Hs_upper_m"]),
        "max_samples_override": max_samples,
        "full_word_requested": full_requested,
        "candidate_rows": rows,
        "widest_candidate_prefix_safe_deg": widest_prefix_safe,
        "P4_CANDIDATE_FULL_HA_PREFIX_CERTIFICATE": status,
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE": False,
        "P4_NORMALIZED_CROSS_BLOCK_ESTABLISHED_HERE": False,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "next_obligation": (
            "if a full-word H/A candidate prefix is safe, compose the operation-matched nonlinear defect with the P3 word decrease and validate the normalized translation/nontranslation cross block on that same candidate; otherwise subdivide the first reported source/vector/covariance obstruction"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    for k in (
        "source_generated_not_trajectory_fit", "full_S_to_attitude_cross_gain_retained",
        "shipping_Joseph_and_immediate_reset_retained", "deployed_quaternion_composition_retained",
        "active_accelerometer_bias_J_ba_identity", "p4_candidate_ball_covered_not_single_cube",
        "cover_is_outward_union", "position_entrance_uses_half_Hs",
    ):
        if d.get(k) is not True:
            failures.append(f"{k} is not true")
    for k in (
        "source_replay_used", "filter_changed", "old_q8_chart_used",
        "active_accelerometer_bias_projection_linearized_through_boundary",
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE",
        "P4_NORMALIZED_CROSS_BLOCK_ESTABLISHED_HERE",
        "P4_USABLE_CERTIFICATE_PROMOTED",
    ):
        if d.get(k) is not False:
            failures.append(f"{k} is not false")
    if d.get("H_dimension") != 18 or d.get("A_dimension") != 21:
        failures.append("H/A dimensions changed")
    if not 0.0 < float(d.get("Hs_upper_m", 0.0)) <= 8.5:
        failures.append("declared Hs proof envelope is missing or widened without review")
    if not d.get("candidate_rows"):
        failures.append("candidate ladder was not evaluated")
    return list(dict.fromkeys(failures))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--max-samples", type=int, default=None,
                    help="diagnostic prefix length; omit for complete recurrence word")
    ap.add_argument("--candidate-index", type=int, default=None,
                    help="0=30deg,1=25deg,2=20deg,3=15deg")
    ap.add_argument("--cover-factor", type=float, default=1.5)
    args = ap.parse_args()
    d = build(args.domain.resolve(), args.max_samples, args.candidate_index, args.cover_factor)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P4_CANDIDATE_FULL_HA_PREFIX_CERTIFICATE"],
        "full_word_requested": d["full_word_requested"],
        "widest_candidate_prefix_safe_deg": d["widest_candidate_prefix_safe_deg"],
        "candidates": [{
            "angle_deg": x["angle_deg"],
            "cover_cells": x["cover_cells"],
            "cover_q_upper": x["cover_q_upper"],
            "both_modes_prefix_safe": x["both_modes_prefix_safe"],
            "first_unclosed_cell": x["first_unclosed_cell"],
        } for x in d["candidate_rows"]],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
