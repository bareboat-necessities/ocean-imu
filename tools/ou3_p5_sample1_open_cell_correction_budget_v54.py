#!/usr/bin/env python3
"""V54: correction budget at the cells V53 leaves open.

V53 re-ran V41's signed-chart `q<8` cover with V51's exact monotone block and
retired the survivor V42 through V50 could not move, but left roughly 48 percent
of the cover open.  Its next obligation names the nominal geometry at the
remaining cells without saying how far away they are.  V54 measures that.

For each named open cell, V54 reconstructs the refined V8/V10 chain from source
with V51's exact monotone enclosure installed, decomposes

    corr^2 = k_perp^2 rho_x^2 + k_par^2 (rho^2 - rho_x^2),
    rho    = |f| chord + eaw1 + ba,

and then reports, for each of the four factors, the reduction that would bring
the cell's composed correction to the largest value the geodesic branch admits.

The admissible correction follows from the SO(3) triangle alone.  A cell closes
if its composed q falls below the target, and `q_+ = 2 tan((phi_c + phi_d)/2)`
is increasing in the correction principal angle `phi_d`, so

    phi_d <= 2 atan(Q_TARGET/2) - phi_c,   phi_c = 2 atan(q_current/2)

is exactly the geodesic-branch admissible correction for a cell whose archived
current chart is `q_current`.

Scope, stated plainly.  The decomposition and the reconstructed correction are
outward-rounded enclosures; the retired witness ties them to V51 by reproducing
its refined correction bit-for-bit.  V41's own recorded radius is reported
beside each reconstruction but is not an ordering constraint: it is the composed
radius after V41's signed/component intersections, which can sit below V10's raw
correction magnitude as well as above it.

The *required reduction* figures are arithmetic on those enclosures, not
enclosures themselves: they say where to aim, not what is proved.  They also
cover the geodesic branch only.  V41 closes a cell on `min(geodesic, product)`,
so a cell whose geodesic branch reaches the target closes, but a cell that
misses it may still close through the product branch.  V54 therefore never
reports a cell as closed or as unreachable; it reports a distance.

Nothing here composes `q<8`, promotes sample 1 or P5, or sets `N_H_words`, and
no filter setting, source domain, six-radian correction limit, `q<8` target, or
source language changes.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p5_sample1_structured_full_gain_v8 as V8
import ou3_p5_sample1_structured_full_gain_v9 as V9
import ou3_p5_sample1_structured_full_gain_v10 as V10
import ou3_p5_sample1_exact_monotone_source_gain_v51 as V51

DEFAULT_DOMAIN = V51.DEFAULT_DOMAIN
SCHEMA = 5400
Q_TARGET = V51.Q_TARGET

V4 = V8.V4
V5 = V8.V5
V6 = V8.V6
FULL = V8.FULL
SUB = V8.SUB

#: Cells recorded by V53, with the run each record came from.  The first two are
#: the refined cover's first-open and worst cells; the third is the survivor the
#: refinement retired, whose record therefore only exists in the parent cover.
#: Pairing that parent chart with a refined correction is the same conservative
#: pairing V51 makes: the parent chart is an upper bound computed with the
#: unrefined gains.
#:
#: `correction_radial_upper_rad` is reported for comparison only.  It is V41's
#: composed radius after its signed/component intersections, which can fall
#: below or rise above V10's raw correction magnitude, so V54 does not require
#: an ordering between the two.
V53_OPEN_CELLS = {
    "first_open": {
        "cell": (0, 2, 23),
        "record_source": "V53 refined cover",
        "sample1_current_cayley_norm_upper": 0.7444250125377312,
        "correction_radial_upper_rad": 1.9654578067229402,
        "post_sample1_cayley_norm_upper": 8.475205389989586,
    },
    "worst": {
        "cell": (23, 18, 4),
        "record_source": "V53 refined cover",
        "sample1_current_cayley_norm_upper": 1.674977608400414,
        "correction_radial_upper_rad": 3.1593913566884573,
        "post_sample1_cayley_norm_upper": 499303.8238549043,
    },
    "retired_witness": {
        "cell": (0, 0, 23),
        "record_source": "V53 parent cover",
        "sample1_current_cayley_norm_upper": V51.V41_Q_CURRENT,
        "correction_radial_upper_rad": 2.050326092645528,
        "post_sample1_cayley_norm_upper": V51.V41_Q_POST,
    },
}

#: V51's refined correction at the retired witness.  V54 reconstructs the same
#: cell through its own code path and must reproduce it bit-for-bit; that is the
#: provenance tie between the two producers.
V51_RETIRED_WITNESS_CORRECTION_RAD = 1.7313776836494923


def _admissible_correction(q_current: float) -> float:
    """Largest correction principal angle the geodesic branch admits.

    Rounded down, so a correction at or below this value is genuinely inside the
    target on that branch.
    """
    q = float(q_current)
    if not (math.isfinite(q) and q >= 0.0):
        raise ValueError("finite nonnegative current chart radius required")
    phi_c = FULL.up(2.0 * math.atan(FULL.up(0.5 * q)))
    phi_target = FULL.down(2.0 * math.atan(FULL.down(0.5 * Q_TARGET)))
    return FULL.down(phi_target - phi_c)


def _required_reductions(*, k_perp: float, k_par: float, rho: float,
                         rho_x: float, target: float) -> dict:
    """Reduction in each factor that reaches ``target`` with the others fixed.

    Solves `k_perp^2 rho_x^2 + k_par^2 (rho^2 - rho_x^2) = target^2` for one
    factor at a time.  A factor that cannot reach the target on its own, because
    the other three already exceed it, is reported as unreachable.
    """
    vals = (k_perp, k_par, rho, rho_x, target)
    if not all(math.isfinite(float(x)) and float(x) >= 0.0 for x in vals):
        raise ValueError("finite nonnegative correction-budget inputs required")
    if rho_x > rho:
        raise ValueError("perpendicular residual exceeds the full residual")

    t2 = target * target
    perp = k_perp * k_perp * rho_x * rho_x
    par_span = rho * rho - rho_x * rho_x
    par = k_par * k_par * par_span
    out = {
        "corr_squared_perpendicular_term": perp,
        "corr_squared_parallel_term": par,
        "corr_squared_total": perp + par,
        "target_corr_squared": t2,
    }

    def _factor(name, current, solve):
        try:
            need = solve()
        except (ValueError, ZeroDivisionError):
            need = None
        if need is None or not math.isfinite(need) or need < 0.0:
            out[name] = {"reachable_alone": False, "current": current}
            return
        out[name] = {
            "reachable_alone": True,
            "current": current,
            "required": need,
            "absolute_reduction": current - need,
            "fractional_reduction": (
                0.0 if current == 0.0 else (current - need) / current),
        }

    def _rho():
        rem = t2 - perp
        if rem < 0.0 or k_par == 0.0:
            return None
        return math.sqrt(rem / (k_par * k_par) + rho_x * rho_x)

    def _rho_x():
        denom = k_perp * k_perp - k_par * k_par
        if denom <= 0.0:
            return None
        rem = t2 - k_par * k_par * rho * rho
        if rem < 0.0:
            return None
        return math.sqrt(rem / denom)

    def _k_par():
        rem = t2 - perp
        if rem < 0.0 or par_span <= 0.0:
            return None
        return math.sqrt(rem / par_span)

    def _k_perp():
        rem = t2 - par
        if rem < 0.0 or rho_x == 0.0:
            return None
        return math.sqrt(rem) / rho_x

    _factor("rho", rho, _rho)
    _factor("rho_x", rho_x, _rho_x)
    _factor("k_parallel", k_par, _k_par)
    _factor("k_perpendicular", k_perp, _k_perp)
    return out


def _evaluate_cell(ctx: dict, pi: int, ti: int, zi: int, *, p_pieces: int,
                   tangent_pieces: int, axial_pieces: int) -> dict | None:
    """Reconstruct one V8/V10 cell with V51's exact monotone block installed."""
    p = SUB.parts(ctx["p_all"].lo, ctx["p_all"].hi, p_pieces)[pi]
    rt0 = SUB.parts(0.0, ctx["rt_cap"], tangent_pieces)[ti]
    rz0 = SUB.parts(-ctx["rz_cap"], ctx["rz_cap"], axial_pieces)[zi]
    child = V6._residual_child(rt0, rz0, ctx["rho0"])
    if child is None:
        return None
    rt, rz = child

    g = ctx["g"]
    t = ctx["t"]
    Y = ctx["Y"]
    r = ctx["r"]
    alpha = ctx["alpha"]
    qaw = ctx["qaw"]
    ba = ctx["ba"]
    blk = V51._first_block(t=t, p=p, r=r, g=g, exact=True)

    d = blk["ktheta"] * rt
    awt = blk["kaw_t"] * rt
    az = blk["kz"] * rz
    fy = -(alpha * awt)
    fz = FULL.I(g) + alpha * az
    fn = V5._norm2_upper(fy.abs_upper(), fz.abs_upper())
    _kn, det = V8._block_gain_bounds(
        blk["a"], Y, blk["c0"], alpha, qaw, blk["b"], blk["bz"],
        blk["det_first"], d, fy, fz, r)

    dhi = max(0.0, d.hi)
    chord = min(2.0, FULL.up(
        ctx["chord0"] + FULL.up(V4._correction_chord_upper(dhi) + ctx["pred_chord"])))
    left_t = ((FULL.I(1.0) - blk["kaw_t"]) * rt).abs_upper()
    left_z = ((FULL.I(1.0) - blk["kz"]) * rz).abs_upper()
    post_r = FULL.up(V5._norm2_upper(
        FULL.up(left_t + ctx["yRt"]), FULL.up(left_z + ctx["yRz"])) + ba)
    post_tri = FULL.up(ctx["aw_pred"] + V5._norm2_upper(
        awt.abs_upper(), az.abs_upper()))
    post = min(post_r, post_tri)
    eaw1 = FULL.up(alpha.hi * post)
    force_term = FULL.up(fn * chord)
    rho = FULL.up(force_term + FULL.up(eaw1 + ba))

    comb = V10._combined_x_residual_upper(
        alpha_lo=float(alpha.lo), alpha_hi=float(alpha.hi),
        first_rot_x_upper=ctx["yRt"], bias_upper=ba,
        error_transport_rotation_norm_upper=ctx["transport_chord"],
        series_rotation_mismatch_upper=V10._series_vs_axis_rotation_mismatch_upper(
            float(d.lo), float(d.hi)),
        pre_first_aw_error_norm_upper=ctx["aw_pred"], gravity=g)
    rho_x = min(rho, float(comb["combined_x_residual_upper_mps2"]))
    k_perp = float(det["scalar_gain_yz_norm_upper"])
    k_par = float(det["two_by_two_theta_x_gain_norm_upper"])
    corr = V9._two_block_correction_upper(k_perp, k_par, rho_x, rho)
    return {
        "cell": [pi, ti, zi],
        "sample1_force_norm_upper_mps2": fn,
        "attitude_rotation_chord_upper": chord,
        "force_term_upper_mps2": force_term,
        "post_prediction_aw_error_norm_upper_mps2": eaw1,
        "accelerometer_bias_bound_mps2": ba,
        "post_residual_identity_upper_mps2": post_r,
        "post_triangle_upper_mps2": post_tri,
        "post_branch_used": "residual_identity" if post == post_r else "triangle",
        "post_left_tangent_upper_mps2": left_t,
        "post_left_axial_upper_mps2": left_z,
        "sample1_residual_norm_upper_mps2": rho,
        "sample1_combined_source_x_residual_upper_mps2": rho_x,
        "Ktheta_perpendicular_block_upper": k_perp,
        "Ktheta_parallel_block_upper": k_par,
        "combined_directional_correction_norm_upper_rad": corr,
        "rho_share_force": 0.0 if rho == 0.0 else force_term / rho,
        "rho_share_post_aw": 0.0 if rho == 0.0 else eaw1 / rho,
        "rho_share_bias": 0.0 if rho == 0.0 else ba / rho,
    }


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24) -> dict:
    """Measure the correction budget at the cells V53 leaves open."""
    path = Path(domain_path).resolve()
    failures: list[str] = []
    cells: dict = {}
    try:
        ctx = V51._source_context(path, source_pieces=source_pieces,
                                  source_cell_index=source_cell_index)
        for label, spec in V53_OPEN_CELLS.items():
            pi, ti, zi = spec["cell"]
            got = _evaluate_cell(ctx, pi, ti, zi, p_pieces=p_pieces,
                                 tangent_pieces=tangent_pieces,
                                 axial_pieces=axial_pieces)
            if got is None:
                failures.append(f"{label} cell {spec['cell']} was pruned")
                continue
            corr = float(got["combined_directional_correction_norm_upper_rad"])
            recorded = float(spec["correction_radial_upper_rad"])
            if label == "retired_witness" and corr != V51_RETIRED_WITNESS_CORRECTION_RAD:
                failures.append(
                    f"retired witness reconstruction {corr} does not reproduce "
                    f"V51's {V51_RETIRED_WITNESS_CORRECTION_RAD}")
            admissible = _admissible_correction(
                spec["sample1_current_cayley_norm_upper"])
            got.update({
                "V53_record_source": spec["record_source"],
                "V53_recorded_correction_radial_upper_rad": recorded,
                "reconstruction_over_V53_recorded_ratio": (
                    math.inf if recorded == 0.0 else corr / recorded),
                "V53_recorded_current_chart_upper": spec[
                    "sample1_current_cayley_norm_upper"],
                "V53_recorded_post_q_upper": spec["post_sample1_cayley_norm_upper"],
                "geodesic_admissible_correction_rad": admissible,
                "correction_gap_to_geodesic_target_rad": corr - admissible,
                "correction_fractional_gap": (
                    math.inf if corr == 0.0 else (corr - admissible) / corr),
                "already_inside_geodesic_target": corr <= admissible,
                "required_reductions": _required_reductions(
                    k_perp=float(got["Ktheta_perpendicular_block_upper"]),
                    k_par=float(got["Ktheta_parallel_block_upper"]),
                    rho=float(got["sample1_residual_norm_upper_mps2"]),
                    rho_x=float(got["sample1_combined_source_x_residual_upper_mps2"]),
                    target=max(0.0, admissible)),
            })
            cells[label] = got
    except Exception as exc:
        failures.append(f"V54 open-cell correction budget: {exc}")

    ranked = sorted(
        ((k, float(v["correction_fractional_gap"])) for k, v in cells.items()
         if not v["already_inside_geodesic_target"]),
        key=lambda kv: kv[1])
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_OPEN_CELL_CORRECTION_BUDGET_V54",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "exact_monotone_corner_enclosure_used": True,
        "required_reductions_are_diagnostics_not_enclosures": True,
        "geodesic_branch_only": True,
        "cell_reported_closed_or_unreachable_here": False,
        "deployed_correction_limit_rad": 6.0,
        "deployed_correction_limit_increased": False,
        "q_target": Q_TARGET,
        "q8_composed_here": False,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_established_here": False,
        "open_cells": cells,
        "nearest_open_cell": ranked[0][0] if ranked else None,
        "nearest_open_cell_fractional_gap": ranked[0][1] if ranked else None,
        "P5_SAMPLE1_OPEN_CELL_CORRECTION_BUDGET_V54": (
            "PASS" if not failures else "NOT_ESTABLISHED"),
        "next_obligation": (
            "REFINE_NOMINAL_SAMPLE1_RESIDUAL_AT_THE_NEAREST_OPEN_Q8_CELL"
            if ranked else "RECHECK_OPEN_Q8_CELLS_AGAINST_A_FRESH_V53_COVER"),
        "failures": list(dict.fromkeys(failures)),
    }


def validate(d: dict) -> list[str]:
    """Validate the fail-closed V54 proof-artifact contract."""
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P5_SAMPLE1_OPEN_CELL_CORRECTION_BUDGET_V54":
        f.append("qualification mismatch")
    for k in ("source_generated_not_trajectory_fit",
              "exact_monotone_corner_enclosure_used",
              "required_reductions_are_diagnostics_not_enclosures",
              "geodesic_branch_only"):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in ("source_replay_used", "filter_changed",
              "cell_reported_closed_or_unreachable_here",
              "deployed_correction_limit_increased", "q8_composed_here",
              "q8_word_promoted_here", "whole_word_promoted_here",
              "N_H_words_set_here", "P5_established_here"):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        f.append("deployed correction limit changed")
    if float(d.get("q_target", 0.0)) != Q_TARGET:
        f.append("q target changed")

    cells = d.get("open_cells") or {}
    if set(cells) != set(V53_OPEN_CELLS):
        f.append("open-cell set changed")
    for label, cell in cells.items():
        spec = V53_OPEN_CELLS.get(label) or {}
        if tuple(cell.get("cell", ())) != tuple(spec.get("cell", ())):
            f.append(f"{label} cell index changed")
        corr = float(cell.get("combined_directional_correction_norm_upper_rad",
                              math.inf))
        if not (math.isfinite(corr) and corr >= 0.0):
            f.append(f"{label} reconstruction is not a finite correction")
        if label == "retired_witness" and corr != V51_RETIRED_WITNESS_CORRECTION_RAD:
            f.append("retired witness does not reproduce V51's correction")
        rho = float(cell.get("sample1_residual_norm_upper_mps2", -1.0))
        rho_x = float(cell.get("sample1_combined_source_x_residual_upper_mps2",
                               math.inf))
        if not (0.0 <= rho_x <= rho):
            f.append(f"{label} residual decomposition is inconsistent")
    if d.get("P5_SAMPLE1_OPEN_CELL_CORRECTION_BUDGET_V54") not in (
            "PASS", "NOT_ESTABLISHED"):
        f.append("invalid V54 status")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces", type=int, default=4)
    ap.add_argument("--source-cell-index", type=int, default=0)
    ap.add_argument("--p-pieces", type=int, default=24)
    ap.add_argument("--tangent-pieces", type=int, default=24)
    ap.add_argument("--axial-pieces", type=int, default=24)
    ap.add_argument("--output", type=Path, required=True)
    x = ap.parse_args()
    d = build(x.domain, source_pieces=x.source_pieces,
              source_cell_index=x.source_cell_index, p_pieces=x.p_pieces,
              tangent_pieces=x.tangent_pieces, axial_pieces=x.axial_pieces)
    vf = validate(d)
    d["validation_failures"] = vf
    x.output.parent.mkdir(parents=True, exist_ok=True)
    x.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P5_SAMPLE1_OPEN_CELL_CORRECTION_BUDGET_V54"],
        "nearest": d.get("nearest_open_cell"),
        "nearest_fractional_gap": d.get("nearest_open_cell_fractional_gap"),
        "open_cells": d.get("open_cells"),
        "next": d.get("next_obligation"),
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
