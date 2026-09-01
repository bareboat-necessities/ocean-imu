#!/usr/bin/env python3
"""H=18 interval-AD word backend for the parallel path-dependent P4 route.

This is the first numerical differential backend behind PR #450.  It propagates
an 18-dimensional interval first derivative through the finite-angle error map,
source-bounded H-mode linear chain, Joseph measurement gains, and every
sequential deployed quaternion correction used by one theorem-qualified word.

The backend is deliberately split into a **screen** and a later promotion step.
The screen already differentiates the nonlinear state map; it is not a scalar
Lipschitz estimate.  However, three pieces are still deliberately broader or
less complete than the final P4 theorem:

* the vector pair is put in one canonical rotation gauge with source norm and
  separation intervals instead of exhausting every source orientation cell;
* the word uses the theorem's mandatory four-S/two-vector events and omits
  optional accepted corrections between them instead of hulling the complete
  branch family;
* whitening uses the P3 computational congruence as a conditioning screen until
  actual per-node Sigma_KF(g), Sigma_KF(h) factors are materialized.

Consequently ``H18_SCREEN_GAMMA_LT_ONE`` is useful numerical evidence but can
never set ``P4_USABLE_CERTIFICATE_PROMOTED``.  The route remains fail-closed.

The attitude prediction preserves the source tangent blocks ``R`` and ``B``.
The finite-angle transport uses exact Cayley composition for the transported
attitude and sends ``B*delta_b_g`` through the deployed quaternion/Cayley AD
primitive.  The deterministic gyro-transport disturbance enlarges only the AD
value interval, not the state derivative.

Accepted vector residuals are differentiated as exact rational rotation maps in
a canonical common-rotation gauge:

    r_a = (R(c)-I) f + R(c) delta_a_w,
    r_m = (R(c)-I) m,
    r_S = -delta_S.

At c=0 these have the shipping Jacobians ``-[f]_x``, ``J_aw=I`` and
``-[m]_x``.  The common-rotation gauge is exact for the vector geometry; the
remaining source-orientation/covariance correlation is why this stage is still
called a screen rather than the final edge certificate.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Sequence

from ou3_interval import Interval, hull, matrix_add, matrix_mul, matrix_transpose
import ou3_explicit_information_word_certificate as P3
import ou3_implementation_word_language as WORDS
import ou3_interval_ad as AD
import ou3_p4_candidate_full_word as CAND
import ou3_p4_operation_matched_sector_certificate as SECTOR
import ou3_p5_full_h_prefix_cells as H
import ou3_p5_full_h_prefix_cells_v2 as H2
import ou3_vector_uco_certificate as VECTOR

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
N = 18


def I(x: float) -> Interval:
    return Interval.point(float(x))


def _box(a: float) -> Interval:
    a = abs(float(a))
    return Interval(math.nextafter(-a, -math.inf), math.nextafter(a, math.inf))


def _ad_constant_interval(x: Interval, n: int = N) -> AD.AD:
    return AD.constant(x, n)


def _ad_matvec_interval(A, x: Sequence[AD.AD]) -> list[AD.AD]:
    out = []
    for row in A:
        y = AD.constant(0.0, N)
        for a, b in zip(row, x):
            y = y + AD.constant(a, N) * b
        out.append(y)
    return out


def _ad_hull_vectors(a: Sequence[AD.AD], b: Sequence[AD.AD]) -> list[AD.AD]:
    if len(a) != len(b):
        raise ValueError("AD branch hull length mismatch")
    return [AD.hull_ad(x, y) for x, y in zip(a, b)]


def _state_values(z: Sequence[AD.AD]) -> list[Interval]:
    return [x.val for x in z]


def _initial_state(domain: dict, cbox) -> list[AD.AD]:
    CAND._configure_mode("H")
    e, _ba, _position = CAND._initial_error("H", domain)
    vals = [Interval(a, b) for a, b in cbox] + list(e[3:])
    if len(vals) != N:
        raise RuntimeError("H18 initial state dimension mismatch")
    return [AD.independent(vals[i], i, N) for i in range(N)]


def _prediction(z: Sequence[AD.AD], F, Rstep, domain: dict, h: float) -> list[AD.AD]:
    """Finite-angle attitude prediction plus exact H linear-chain transport."""
    c = list(z[:3])
    bg = list(z[3:6])
    transported = _ad_matvec_interval(Rstep, c)
    Bstep = [[F[i][3 + j] for j in range(3)] for i in range(3)]
    db = _ad_matvec_interval(Bstep, bg)
    wdist = float(domain["startup"]["effective_deterministic_gyro_transport_disturbance_upper_rad_s"])
    extra = _box(h * wdist)
    db = [AD.AD(x.val + extra, x.der) for x in db]
    cp = AD.deployed_correct_cayley(transported, db)

    out = list(z)
    out[:3] = cp
    # H mode: gyro bias is held; v,p,S,a_w use the exact integrated-OU chain.
    for i in range(3, N):
        y = AD.constant(0.0, N)
        for j in range(3, N):
            y = y + AD.constant(F[i][j], N) * z[j]
        out[i] = y
    return out


def _rotation_residual_acc(z: Sequence[AD.AD], force: Sequence[Interval]) -> list[AD.AD]:
    R = AD.rotation_from_cayley(z[:3])
    f = [_ad_constant_interval(x) for x in force]
    aw = list(z[15:18])
    Rf = AD.matvec(R, f)
    Raw = AD.matvec(R, aw)
    return [Rf[i] - f[i] + Raw[i] for i in range(3)]


def _rotation_residual_mag(z: Sequence[AD.AD], mag: Sequence[Interval]) -> list[AD.AD]:
    R = AD.rotation_from_cayley(z[:3])
    m = [_ad_constant_interval(x) for x in mag]
    Rm = AD.matvec(R, m)
    return [Rm[i] - m[i] for i in range(3)]


def _residual_S(z: Sequence[AD.AD]) -> list[AD.AD]:
    return [-z[12 + i] for i in range(3)]


def _zero_matrix(rows: int, cols: int):
    return [[I(0.0) for _ in range(cols)] for _ in range(rows)]


def _H_acc_canonical(force: Sequence[Interval]):
    # Local derivative of (R(c)-I)f + R(c)aw at c=aw=0.
    fx, fy, fz = force
    Hm = _zero_matrix(3, N)
    # -[f]_x
    Hm[0][1] = fz; Hm[0][2] = -fy
    Hm[1][0] = -fz; Hm[1][2] = fx
    Hm[2][0] = fy; Hm[2][1] = -fx
    for i in range(3):
        Hm[i][15 + i] = I(1.0)
    return Hm


def _H_mag_canonical(mag: Sequence[Interval]):
    mx, my, mz = mag
    Hm = _zero_matrix(3, N)
    Hm[0][1] = mz; Hm[0][2] = -my
    Hm[1][0] = -mz; Hm[1][2] = mx
    Hm[2][0] = my; Hm[2][1] = -mx
    return Hm


def _accepted_update(Pm, z: Sequence[AD.AD], Hm, Rm, r: Sequence[AD.AD]):
    cell = H._measurement_cell(Pm, Hm, Rm, [x.val for x in r])
    K = cell["K"]
    dx = _ad_matvec_interval(K, r)
    d = [-x for x in dx[:3]]
    cp = AD.deployed_correct_cayley(z[:3], d)
    out = list(z)
    out[:3] = cp
    for i in range(3, N):
        out[i] = z[i] - dx[i]
    return cell["P_accepted"], out, cell


def _conditioned_jacobian(J, scale2: Sequence[float]):
    if len(scale2) != N:
        raise RuntimeError("P3 conditioning scale dimension mismatch")
    s = [math.sqrt(float(x)) for x in scale2]
    if not all(math.isfinite(x) and x > 0.0 for x in s):
        raise RuntimeError("P3 conditioning scale is not positive")
    out = []
    for i in range(N):
        row = []
        for j in range(N):
            q = Interval.outward_bounds(s[i] / s[j], s[i] / s[j])
            row.append(q * J[i][j])
        out.append(row)
    return out


def _block_norms(J) -> dict:
    blocks = {
        "theta": range(0, 3), "bg": range(3, 6), "v": range(6, 9),
        "p": range(9, 12), "S": range(12, 15), "aw": range(15, 18),
    }
    out = {}
    for rname, rr in blocks.items():
        for cname, cc in blocks.items():
            A = [[J[i][j] for j in cc] for i in rr]
            out[f"{rname}<-{cname}"] = AD.interval_matrix_op2_upper(A)
    return out


def _canonical_vector_cells(domain: dict) -> tuple[list[Interval], list[Interval], dict]:
    live = domain["normal_live"]
    f = Interval.outward_bounds(
        float(live["specific_force_norm_lower_mps2"]),
        float(live["specific_force_norm_upper_mps2"]),
    )
    m = Interval.outward_bounds(
        float(live["magnetic_vector_norm_lower_uT"]),
        float(live["magnetic_vector_norm_upper_uT"]),
    )
    s = float(live["vector_sine_separation_lower"])
    if not 0.0 < s < 1.0:
        raise RuntimeError("invalid vector separation lower")
    c_hi = math.nextafter(math.sqrt(max(0.0, 1.0 - s * s)), math.inf)
    # Common rotation sends f to +e3; rotate about e3 so m lies in x-z plane.
    force = [I(0.0), I(0.0), f]
    mag = [m * Interval.outward_bounds(s, 1.0), I(0.0), m * Interval.outward_bounds(-c_hi, c_hi)]
    return force, mag, {"sine_lower": s, "cos_abs_upper": c_hi}


def _mandatory_schedule(words: dict, samples: int, h: float) -> dict:
    wc = words["word_contract"]
    trans = wc["translation_recurrence"]
    q = int(trans["spread_index_q_W"])
    pseudo_gap_lo = float(trans["pseudo_gap_min_s"])
    spacing_steps = max(1, int(math.floor(q * pseudo_gap_lo / h)))
    # Place the validated spread family as late as possible while retaining four firings.
    last = samples - 1
    s_steps = [last - 3 * spacing_steps, last - 2 * spacing_steps, last - spacing_steps, last]
    if min(s_steps) < 0:
        # Fallback to an evenly spread four-event schedule for a truncated smoke run.
        s_steps = sorted(set([0, max(0, last // 3), max(0, 2 * last // 3), last]))
    gap = VECTOR.build()["operating_envelope"]["packet_gap_s"]
    mag_gap_steps = max(1, int(math.ceil(float(gap[1]) / h)))
    v_steps = [max(0, last - mag_gap_steps), last]
    return {
        "S_steps": sorted(set(s_steps)),
        "vector_steps": sorted(set(v_steps)),
        "spread_index_q": q,
        "spread_spacing_steps_from_lower_gap": spacing_steps,
        "vector_gap_steps": mag_gap_steps,
        "source_complete_timing_realization_proved": False,
    }


def _aw_sync_covariance_overapprox(Pm, sigma_hi: float):
    """H-mode source screen: allow the filter a_w covariance floor at every prefix."""
    cap = math.nextafter(float(sigma_hi) * float(sigma_hi), math.inf)
    out = [[Pm[i][j] for j in range(N)] for i in range(N)]
    for i in range(15, 18):
        out[i][i] = hull(out[i][i], Interval(0.0, cap))
    return H._psd_tighten(out)


def _run_cell(domain_path: Path, domain: dict, cbox, samples: int, scale2: Sequence[float]) -> dict:
    src = H._source_cell()
    F, Q, Rstep = H2._tight_transition_and_Q(src, domain)
    Pm = H2._corrected_initial_covariance(src, domain_path)
    z = _initial_state(domain, cbox)
    force, mag, geometry = _canonical_vector_cells(domain)
    Hs = H._H_S()
    Ha = _H_acc_canonical(force)
    Hm = _H_mag_canonical(mag)
    vc = VECTOR.build()["configured_measurement_bounds"]
    Racc = H._R_diag(float(vc["acc_measurement_std_mps2"]))
    Rmag = H._R_diag(float(vc["mag_measurement_std_uT"]))
    RS = H._R_S(src)
    h = float(src["dt_s"])
    words = WORDS.build(domain_path)
    schedule = _mandatory_schedule(words, samples, h)
    sigma_hi = src["sigma_aw_mps2"].hi

    max_raw = 0.0
    max_conditioned = 0.0
    worst = None
    inverse_counts = {}
    operation_count = 0

    def record(label: str, step: int):
        nonlocal max_raw, max_conditioned, worst, operation_count
        J = AD.jacobian(z)
        raw = AD.interval_matrix_op2_upper(J)
        Jc = _conditioned_jacobian(J, scale2)
        cond = AD.interval_matrix_op2_upper(Jc)
        operation_count += 1
        if cond >= max_conditioned:
            max_conditioned = cond
            max_raw = max(max_raw, raw)
            worst = {
                "step": step,
                "operation": label,
                "raw_jacobian_norm_upper": raw,
                "P3_congruence_conditioned_norm_upper": cond,
                "block_norms_conditioned": _block_norms(Jc),
                "state_value_abs_upper": [x.abs_upper() for x in _state_values(z)],
            }
        else:
            max_raw = max(max_raw, raw)

    record("entry", -1)
    for k in range(samples):
        Pm = H._psd_tighten(matrix_add(matrix_mul(matrix_mul(F, Pm), matrix_transpose(F)), Q))
        Pm = _aw_sync_covariance_overapprox(Pm, sigma_hi)
        z = _prediction(z, F, Rstep, domain, h)
        record("prediction", k)

        if k in schedule["S_steps"]:
            Pm, z, cell = _accepted_update(Pm, z, Hs, RS, _residual_S(z))
            inverse_counts[cell["inverse_backend"]] = inverse_counts.get(cell["inverse_backend"], 0) + 1
            record("mandatory_S", k)

        if k in schedule["vector_steps"]:
            Pm, z, cell = _accepted_update(Pm, z, Ha, Racc, _rotation_residual_acc(z, force))
            inverse_counts[cell["inverse_backend"]] = inverse_counts.get(cell["inverse_backend"], 0) + 1
            record("mandatory_accelerometer", k)
            Pm, z, cell = _accepted_update(Pm, z, Hm, Rmag, _rotation_residual_mag(z, mag))
            inverse_counts[cell["inverse_backend"]] = inverse_counts.get(cell["inverse_backend"], 0) + 1
            record("mandatory_magnetometer", k)

    J = AD.jacobian(z)
    Jc = _conditioned_jacobian(J, scale2)
    return {
        "screen_completed": True,
        "samples": samples,
        "schedule": schedule,
        "canonical_vector_geometry": geometry,
        "entry_cayley_box": cbox,
        "endpoint_raw_jacobian_norm_upper": AD.interval_matrix_op2_upper(J),
        "endpoint_P3_congruence_conditioned_norm_upper": AD.interval_matrix_op2_upper(Jc),
        "endpoint_block_norms_conditioned": _block_norms(Jc),
        "maximum_prefix_raw_jacobian_norm_upper": max_raw,
        "maximum_prefix_P3_congruence_conditioned_norm_upper": max_conditioned,
        "worst_prefix": worst,
        "inverse_backend_counts": inverse_counts,
        "differentiated_operation_count": operation_count,
        "endpoint_state_value_abs_upper": [x.abs_upper() for x in _state_values(z)],
    }


def build(domain_path: Path = DEFAULT_DOMAIN, *, samples: int | None = None,
          cell_limit: int = 1, ball_inflation: float = 1.5) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("H18 interval-AD screen must not be trajectory fitted")

    words = WORDS.build(path)
    wf = WORDS.validate(words)
    sector = SECTOR.build(path)
    sf = SECTOR.validate(sector)
    p3 = P3.build(path)
    pf = P3.validate(p3)
    failures = [f"word-language: {x}" for x in wf] + [f"sector: {x}" for x in sf] + [f"P3: {x}" for x in pf]

    full_samples = int(words["word_contract"]["conditional_word_language"]["word_samples_upper_at_configured_dt"])
    use_samples = full_samples if samples is None else max(1, min(full_samples, int(samples)))
    q = float(sector["design_cayley_norm_upper"])
    covers = CAND._ball_box_cover(q, max_box_norm_factor=float(ball_inflation))
    total_cells = len(covers)
    limit = total_cells if int(cell_limit) <= 0 else min(total_cells, int(cell_limit))
    selected = covers[:limit]

    comparison = p3["modes"]["H"]["matrix_comparison"]
    scale2 = [float(x) for x in comparison["comparison_scale_diagonal_squared"]]
    cells = []
    first_failure = None
    for i, cbox in enumerate(selected):
        try:
            row = _run_cell(path, domain, cbox, use_samples, scale2)
            row["index"] = i
            cells.append(row)
        except Exception as exc:
            first_failure = {"index": i, "entry_cayley_box": cbox, "error": f"{type(exc).__name__}: {exc}"}
            break

    max_endpoint = max((float(x["endpoint_P3_congruence_conditioned_norm_upper"]) for x in cells), default=math.inf)
    max_prefix = max((float(x["maximum_prefix_P3_congruence_conditioned_norm_upper"]) for x in cells), default=math.inf)
    screen_complete = first_failure is None and len(cells) == len(selected)
    whole_word = use_samples == full_samples
    whole_ball = limit == total_cells and screen_complete
    gamma_lt_one = bool(whole_word and whole_ball and max_endpoint < 1.0)

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_H18_INTERVAL_AD_WORD_SCREEN",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "dimension": N,
        "outer_angle_rad": float(sector["design_full_attitude_angle_rad"]),
        "outer_cayley_norm_upper": q,
        "outer_ball_box_cover_total": total_cells,
        "outer_ball_box_cells_requested": limit,
        "outer_ball_box_cells_completed": len(cells),
        "all_outer_ball_cells_checked": whole_ball,
        "full_word_samples": full_samples,
        "samples_checked": use_samples,
        "full_word_horizon_checked": whole_word,
        "interval_AD_used_for_state_return_map": True,
        "finite_difference_used": False,
        "deployed_quaternion_generalized_jacobian_used": True,
        "full_18_state_cross_derivatives_retained": True,
        "canonical_common_rotation_PE_gauge_used": True,
        "all_source_vector_orientation_covariance_correlations_checked": False,
        "mandatory_four_S_events_used": True,
        "mandatory_two_packet_vector_PE_used": True,
        "optional_accepted_branch_family_between_required_events_checked": False,
        "aw_covariance_sync_overapproximated_at_every_prefix": True,
        "actual_per_node_Sigma_KF_whitening_used": False,
        "P3_computational_congruence_used_for_screening_only": True,
        "P3_delta_used_as_nonlinear_radius": False,
        "source_graph_all_reachable_edges_checked": False,
        "screen_cells": cells,
        "first_failure": first_failure,
        "max_endpoint_P3_congruence_conditioned_norm_upper": max_endpoint,
        "max_prefix_P3_congruence_conditioned_norm_upper": max_prefix,
        "H18_SCREEN_GAMMA_LT_ONE": gamma_lt_one,
        "H18_COMPLETE_SOURCE_EDGE_CONTRACTION_ESTABLISHED_HERE": False,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "validation_scope": "NUMERICAL_H18_DIFFERENTIAL_SCREEN_NOT_THEOREM_PROMOTION",
        "next_obligation": (
            "materialize actual source-node H covariance/information factors; cover every source vector orientation/covariance cell and every optional accepted/rejected/not-due branch; then rerun the full outer-ball H18 word on every reachable g->h edge and replace the P3 congruence screen by the exact generalized metric test"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P4_H18_INTERVAL_AD_WORD_SCREEN":
        f.append("wrong H18 interval-AD qualification")
    for key in ("source_generated_not_trajectory_fit", "interval_AD_used_for_state_return_map",
                "deployed_quaternion_generalized_jacobian_used", "full_18_state_cross_derivatives_retained",
                "canonical_common_rotation_PE_gauge_used", "mandatory_four_S_events_used",
                "mandatory_two_packet_vector_PE_used", "aw_covariance_sync_overapproximated_at_every_prefix",
                "P3_computational_congruence_used_for_screening_only"):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in ("source_replay_used", "filter_changed", "finite_difference_used",
                "all_source_vector_orientation_covariance_correlations_checked",
                "optional_accepted_branch_family_between_required_events_checked",
                "actual_per_node_Sigma_KF_whitening_used", "P3_delta_used_as_nonlinear_radius",
                "source_graph_all_reachable_edges_checked", "H18_COMPLETE_SOURCE_EDGE_CONTRACTION_ESTABLISHED_HERE",
                "P4_USABLE_CERTIFICATE_PROMOTED"):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if int(d.get("dimension", 0)) != N:
        f.append("H18 screen dimension is not 18")
    if float(d.get("outer_angle_rad", math.nan)) != 0.80:
        f.append("H18 screen outer angle is not exactly 0.80 rad")
    if int(d.get("outer_ball_box_cover_total", 0)) <= 0:
        f.append("H18 screen outer-ball cover is empty")
    for key in ("max_endpoint_P3_congruence_conditioned_norm_upper", "max_prefix_P3_congruence_conditioned_norm_upper"):
        x = d.get(key)
        if not isinstance(x, (int, float)) or not math.isfinite(float(x)) or float(x) < 0.0:
            if d.get("first_failure") is None:
                f.append(f"{key} is invalid without a numerical failure witness")
    # A sub-unity screen is never enough to promote while the source/metric flags above are false.
    if d.get("H18_SCREEN_GAMMA_LT_ONE") is True and d.get("P4_USABLE_CERTIFICATE_PROMOTED") is not False:
        f.append("H18 screening norm improperly promoted P4")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--samples", type=int, default=None)
    ap.add_argument("--cell-limit", type=int, default=1, help="0 means all Cayley ball-cover cells")
    ap.add_argument("--ball-inflation", type=float, default=1.5)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build(args.domain.resolve(), samples=args.samples, cell_limit=args.cell_limit,
                ball_inflation=args.ball_inflation)
    vf = validate(out)
    out["validation_pass"] = not vf
    out["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "validation_pass": not vf,
        "samples": [out["samples_checked"], out["full_word_samples"]],
        "cells": [out["outer_ball_box_cells_completed"], out["outer_ball_box_cells_requested"], out["outer_ball_box_cover_total"]],
        "max_endpoint_conditioned_norm": out["max_endpoint_P3_congruence_conditioned_norm_upper"],
        "max_prefix_conditioned_norm": out["max_prefix_P3_congruence_conditioned_norm_upper"],
        "screen_gamma_lt_one": out["H18_SCREEN_GAMMA_LT_ONE"],
        "first_failure": out["first_failure"],
        "next": out["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
