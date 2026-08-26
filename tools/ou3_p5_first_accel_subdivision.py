#!/usr/bin/env python3
"""Source-complete subdivision diagnostic for the first P5 H accelerometer prefix.

The V3 full-matrix producer now fails at sample zero only because its single
source cell gives the first accelerometer correction an interval norm above the
validated deployed-quaternion range.  This producer does not widen that range
and does not promote P5.  It asks the narrower question needed before building
an expensive whole-word branch tree: which dependencies must be retained to
make that first correction finite enough?

It uses exactly the active V3 primitives (dependency-preserving OU transition,
18x18 covariance prediction, Joseph/reset S branch and deployed-quaternion
backend), but keeps the pseudo due/not-due branches separate instead of hulling
them.  The continuous tuner source cell is split geometrically in tau,
sigma_aw and R_S.  The initial Cayley ball is covered by signed Cartesian boxes
whose interiors intersect the certified Euclidean ball.  Each child is then
propagated independently through prediction and, when selected, the due S map.

The accelerometer Jacobian is intentionally still the V3 source-uniform box in
this diagnostic.  Therefore failure after the source/Cayley split is useful:
it proves that the next mandatory split is the physical vector/orientation
cell, rather than encouraging arbitrary refinement of tuner intervals.  No
sample/replay values enter this calculation.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, matrix_add, matrix_mul, matrix_transpose
import ou3_p5_full_h_prefix_cells as V1
import ou3_p5_full_h_prefix_cells_v3 as V3
import ou3_p5_heading_handoff_contract as HEADING

DEFAULT_DOMAIN = V1.DEFAULT_DOMAIN
SCHEMA = 1


def _geom_split(x: Interval, pieces: int) -> list[Interval]:
    if pieces < 1 or not (0.0 < x.lo <= x.hi):
        raise ValueError("positive interval and pieces>=1 required")
    if pieces == 1 or x.lo == x.hi:
        return [x]
    ratio = (x.hi / x.lo) ** (1.0 / pieces)
    edges = [x.lo]
    for _ in range(pieces - 1):
        edges.append(edges[-1] * ratio)
    edges.append(x.hi)
    return [Interval.outward_bounds(edges[i], edges[i + 1]) for i in range(pieces)]


def _linear_edges(radius: float, pieces_per_axis: int) -> list[float]:
    if not (radius > 0.0 and pieces_per_axis >= 2):
        raise ValueError("positive radius and at least two axis pieces required")
    return [-radius + (2.0 * radius * i) / pieces_per_axis for i in range(pieces_per_axis + 1)]


def _min_abs(x: Interval) -> float:
    if x.lo <= 0.0 <= x.hi:
        return 0.0
    return min(abs(x.lo), abs(x.hi))


def _ball_cover(radius: float, pieces_per_axis: int) -> list[list[Interval]]:
    """Finite Cartesian cover of ||x||<=radius; exterior-only boxes are dropped."""
    edges = _linear_edges(radius, pieces_per_axis)
    one = [Interval.outward_bounds(edges[i], edges[i + 1]) for i in range(pieces_per_axis)]
    out: list[list[Interval]] = []
    r2 = radius * radius
    for x in one:
        for y in one:
            for z in one:
                mn2 = _min_abs(x) ** 2 + _min_abs(y) ** 2 + _min_abs(z) ** 2
                if mn2 <= math.nextafter(r2, math.inf):
                    out.append([x, y, z])
    return out


def _source_children(src: dict, pieces: int) -> list[dict]:
    sched = V1.P3CELL.source_schedule()
    out = []
    for tau in _geom_split(src["tau_s"], pieces):
        cadence = V1.P3CELL.cadence_bounds(tau, sched)
        for sigma in _geom_split(src["sigma_aw_mps2"], pieces):
            for rs in _geom_split(src["R_S_filter_std"], pieces):
                lo = max(src["pseudo_period_s"].lo, cadence[0])
                hi = min(src["pseudo_period_s"].hi, cadence[1])
                child = dict(src)
                child["tau_s"] = tau
                child["sigma_aw_mps2"] = sigma
                child["R_S_filter_std"] = rs
                child["pseudo_period_s"] = Interval.outward_bounds(lo, hi)
                out.append(child)
    return out


def _predict_state(Pm, e, src: dict, domain: dict):
    F, Q, Rstep = V1._transition_and_Q(src, domain)
    Pp = V1._psd_tighten(matrix_add(matrix_mul(matrix_mul(F, Pm), matrix_transpose(F)), Q))
    ep = V1._predict_error(e, F)
    return Pp, ep, Rstep


def _due_s_state(Pm, e, src: dict):
    H = V1._H_S()
    R = V1._R_S(src)
    r = [-e[12 + i] for i in range(3)]
    cell = V1._measurement_cell(Pm, H, R, r)
    d = V1._vec_neg(cell["dx"][0:3])
    e2 = list(e)
    for i in range(3, V1.N):
        e2[i] = e[i] - cell["dx"][i]
    return cell["P_accepted"], e2, d, cell


def _acc_gain(Pm, domain: dict, Racc):
    H = V1._H_acc(domain)
    PHt, S = V1._innovation(Pm, H, Racc)
    Sinv, backend = V1._spd_inverse_enclosure(S, Racc)
    return H, matrix_mul(PHt, Sinv), backend


def _acc_correction_from_gain(K, H, e, c, domain: dict):
    q = min(8.0, V1._norm_upper(c))
    fhi = float(domain["normal_live"]["specific_force_norm_upper_mps2"])
    aw_hi = max(e[i].abs_upper() for i in V1.AW)
    eta = V1.up(
        V1.VEFF.accel_attitude_eta_per_vector_norm_upper(q) * fhi
        + V1.VEFF.accel_latent_cross_gain_upper(q) * aw_hi
    )
    z = [V1.I(0.0) for _ in range(V1.N)]
    for i in range(3):
        z[i] = c[i]
    for i in V1.AW:
        z[i] = e[i] + V1._box(eta)
    r = V1._mat_vec(H, z)
    ba = float(domain["startup"]["physical_handoff_coordinate_bounds"]["accelerometer_bias_error_norm_upper_mps2"])
    r = [x + V1._box(ba) for x in r]
    dx = V1._mat_vec(K, r)
    return V1._norm_upper(V1._vec_neg(dx[0:3])), eta


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 2, cayley_axis_pieces: int = 4) -> dict:
    V3._install_backend()
    domain_path = Path(domain_path).resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("P5 subdivision diagnostic must not be trajectory fitted")

    heading = HEADING.build(domain_path)
    failures = [f"heading: {x}" for x in HEADING.validate(heading)]
    src0 = V1._source_cell()
    source_cells = _source_children(src0, source_pieces)
    q0 = float(heading["gauged_timeout_subbranch"]["full_attitude_cayley_norm_upper"])
    c_cells = _ball_cover(q0, cayley_axis_pieces)
    vc = V1.VECTOR.build()["configured_measurement_bounds"]
    Racc = V1._R_diag(float(vc["acc_measurement_std_mps2"]))

    max_d = 0.0
    min_d = math.inf
    over_limit = 0
    backend_counts = {"FIXED_PIVOT_INTERVAL_GAUSS_JORDAN": 0, "SPD_S_GE_R_SPECTRAL_ENTRY_ENCLOSURE": 0}
    first_over = None

    for si, src in enumerate(source_cells):
        P0 = V1._initial_covariance(src, domain_path)
        e0 = V1._initial_error(domain)
        Pp, ep, Rstep = _predict_state(P0, e0, src, domain)
        Pdue, edue, dS, Scell = _due_s_state(Pp, ep, src)
        phases = {
            "not_due": (Pp, ep, None, None),
            "due": (Pdue, edue, dS, Scell["inverse_backend"]),
        }
        phase_gain = {}
        for phase, (P1, _e1, _dS, _sbackend) in phases.items():
            H, K, backend = _acc_gain(P1, domain, Racc)
            backend_counts[backend] += 1
            phase_gain[phase] = (H, K, backend)

        for ci, c0 in enumerate(c_cells):
            cp = V1._predict_c(c0, Rstep, domain, float(src["dt_s"]))
            for phase, (_P1, e1, ds, s_backend) in phases.items():
                c1 = cp if ds is None else V1.SIGNED.compose_cell(cp, ds)["c_plus"]
                try:
                    H, K, backend = phase_gain[phase]
                    dnorm, eta = _acc_correction_from_gain(K, H, e1, c1, domain)
                except Exception as exc:
                    row = {
                        "source_cell": si,
                        "cayley_cell": ci,
                        "pseudo_phase": phase,
                        "exception": f"{type(exc).__name__}: {exc}",
                    }
                    if first_over is None:
                        first_over = row
                    over_limit += 1
                    continue
                max_d = max(max_d, dnorm)
                min_d = min(min_d, dnorm)
                if dnorm > 6.0:
                    over_limit += 1
                    if first_over is None:
                        first_over = {
                            "source_cell": si,
                            "cayley_cell": ci,
                            "pseudo_phase": phase,
                            "correction_norm_upper_rad": dnorm,
                            "acc_effective_aw_eta_norm_upper": eta,
                            "inverse_backend": backend,
                            "S_inverse_backend": s_backend,
                            "tau_s": src["tau_s"].as_list(),
                            "sigma_aw_mps2": src["sigma_aw_mps2"].as_list(),
                            "R_S_filter_std": src["R_S_filter_std"].as_list(),
                            "cayley_cell_box": [x.as_list() for x in c0],
                        }

    total = len(source_cells) * len(c_cells) * 2
    all_inside = total > 0 and over_limit == 0 and not failures
    next_obligation = (
        "PROPAGATE_SUBDIVIDED_CHILDREN_TO_NEXT_PREFIX"
        if all_inside else
        "ACCEL_VECTOR_ORIENTATION_AND_STATE_DIRECTION_SUBDIVISION_REQUIRED"
    )
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_FIRST_ACCEL_FULL_MATRIX_SOURCE_CAYLEY_SUBDIVISION_DIAGNOSTIC",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "uses_v3_dependency_preserving_full_matrix_backend": True,
        "deployed_correction_limit_rad": 6.0,
        "deployed_correction_limit_increased": False,
        "pseudo_due_and_not_due_branches_kept_separate": True,
        "source_parameter_subdivision": ["tau", "sigma_aw", "R_S"],
        "source_pieces_per_parameter": source_pieces,
        "cayley_ball_direction_subdivided": True,
        "cayley_axis_pieces": cayley_axis_pieces,
        "cayley_initial_norm_upper": q0,
        "source_cell_count": len(source_cells),
        "cayley_cell_count": len(c_cells),
        "evaluated_child_count": total,
        "children_above_validated_correction_limit": over_limit,
        "all_first_accelerometer_children_inside_validated_correction_range": all_inside,
        "min_first_accelerometer_correction_norm_upper_rad": None if math.isinf(min_d) else min_d,
        "max_first_accelerometer_correction_norm_upper_rad": max_d,
        "first_unclosed_child": first_over,
        "inverse_backend_counts": backend_counts,
        "accelerometer_vector_orientation_still_source_uniform_box": True,
        "whole_word_promoted_here": False,
        "P5_FIRST_ACCEL_SOURCE_CAYLEY_SUBDIVISION": "PASS" if all_inside else "NOT_ESTABLISHED",
        "next_obligation": next_obligation,
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    for k in (
        "source_generated_not_trajectory_fit",
        "uses_v3_dependency_preserving_full_matrix_backend",
        "pseudo_due_and_not_due_branches_kept_separate",
        "cayley_ball_direction_subdivided",
        "accelerometer_vector_orientation_still_source_uniform_box",
    ):
        if d.get(k) is not True:
            failures.append(f"{k} is not true")
    for k in ("source_replay_used", "filter_changed", "deployed_correction_limit_increased", "whole_word_promoted_here"):
        if d.get(k) is not False:
            failures.append(f"{k} is not false")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        failures.append("deployed correction limit changed")
    if int(d.get("evaluated_child_count", 0)) <= 0:
        failures.append("no subdivision children evaluated")
    status = d.get("P5_FIRST_ACCEL_SOURCE_CAYLEY_SUBDIVISION")
    if status == "PASS":
        if d.get("all_first_accelerometer_children_inside_validated_correction_range") is not True:
            failures.append("PASS without complete first-accelerometer closure")
    elif status == "NOT_ESTABLISHED":
        if d.get("first_unclosed_child") is None:
            failures.append("nonclosure missing child witness")
    else:
        failures.append("invalid subdivision status")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces", type=int, default=2)
    ap.add_argument("--cayley-axis-pieces", type=int, default=4)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build(args.domain.resolve(), source_pieces=args.source_pieces, cayley_axis_pieces=args.cayley_axis_pieces)
    vf = validate(out)
    out["validation_pass"] = not vf
    out["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": out["P5_FIRST_ACCEL_SOURCE_CAYLEY_SUBDIVISION"],
        "source_cells": out["source_cell_count"],
        "cayley_cells": out["cayley_cell_count"],
        "children": out["evaluated_child_count"],
        "over_limit": out["children_above_validated_correction_limit"],
        "min_d": out["min_first_accelerometer_correction_norm_upper_rad"],
        "max_d": out["max_first_accelerometer_correction_norm_upper_rad"],
        "first_unclosed": out["first_unclosed_child"],
        "next": out["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
