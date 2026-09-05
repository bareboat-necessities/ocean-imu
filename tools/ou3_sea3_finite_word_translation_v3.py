#!/usr/bin/env python3
"""Commit-aligned source-history-free translation word for canonical OU-III P3.

The active MEKF OU schedule changes only when the staged online tuner candidate
is applied at an IMU boundary.  The shipping staging test is on a 0.1 s clock;
the dynamic-source certificate already proves a 22-sample upper gap.  For this
proof we conservatively admit 20, 21, or 22 samples between active commits.
Ten commit segments therefore span at least 1.0 s, enough to contain the
Normal-Live vector-PE recurrence window.

Within one commit segment tau is fixed, but it may be any value in the complete
SEA3 invariant.  Sigma and R_S are replaced by their global lower endpoints.
At every segment boundary the instantaneous tau cell is forgotten and a single
common Loewner lower is propagated.  Thus the construction has history depth
zero even though successive segments may select unrelated tau values.

The key numerical point is also theorem-level: a 5 ms integrated-OU process
covariance is almost rank deficient in the [v,p,S,a_w] chain.  We therefore do
NOT require intermediate per-sample lower matrices to be SPD.  Scalar Kalman
covariance updates are Loewner monotone for every symmetric lower with positive
innovation denominator.  Adaptive x subdivision is continued until each whole
commit-segment endpoint is strict SPD; only then are all current-source cells
collapsed to the common boundary lower.  This is the retained finite-segment
mechanism, stripped of the retired P2 history/correlation interface.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, symmetric_positive_definite_ldlt
import ou3_sea3_dynamic_source_certificate as DYNAMIC
import ou3_sea3_finite_word_translation as BASE
import ou3_sea3_riccati_tube as TUBE_BASE
import ou3_sea3_riccati_tube_factored as TUBE
import ou3_vector_uco_certificate as VECTOR

SCHEMA = 3
QUALIFICATION = "OU3_SEA3_COMMIT_ALIGNED_FINITE_WORD_TRANSLATION_LOWER"
DEFAULT_DOMAIN = BASE.DEFAULT_DOMAIN
USEFUL_GATE = BASE.USEFUL_GATE
COMMIT_GAPS = (20, 21, 22)
WORD_SEGMENTS = 10
INITIAL_X_CELLS = 96
MAX_X_DEPTH = 10


def _split_x(x: Interval):
    for cut in (TUBE.BRANCH_X,):
        if x.lo < cut < x.hi:
            return (
                Interval(x.lo, math.nextafter(cut, -math.inf)),
                Interval(cut, x.hi),
            )
    mid = math.sqrt(x.lo * x.hi)
    if not (x.lo < mid < x.hi):
        return None
    return Interval.outward_bounds(x.lo, mid), Interval.outward_bounds(mid, x.hi)


def _geom_cells(lo: float, hi: float, count: int):
    if not (0.0 < lo < hi and count > 0):
        raise ValueError("invalid geometric x cover")
    ratio = (hi / lo) ** (1.0 / count)
    edges = [lo]
    for k in range(1, count):
        edges.append(lo * ratio ** k)
    edges.append(hi)
    out = []
    for a, b in zip(edges[:-1], edges[1:]):
        out.append(Interval.outward_bounds(a, b))
    return out


def _one_step(P, x: Interval, sigma2: float, R_aw: float, R_S_z: float):
    # Collapse the incoming lower once.  It need not be SPD; the scalar update
    # below only requires a positive innovation denominator.
    L0, _eps0, _route0 = BASE._common_point_lower(P)
    F = BASE._transition(x)
    Q = BASE._scale(TUBE.step_scaled_q(x), sigma2)
    pred_interval = BASE._sym(
        BASE._add(BASE._mul(BASE._mul(F, L0), BASE._transpose(F)), Q)
    )
    pred, _ep, _rp = BASE._common_point_lower(pred_interval)
    post_aw, _ea, _ra = BASE._measurement_lower(pred, 3, R_aw)
    post_s, _es, _rs = BASE._measurement_lower(post_aw, 2, R_S_z)
    return post_s


def _propagate(P, gap: int, x: Interval, sigma2: float, R_aw: float, R_S_z: float):
    out = [[P[i][j] for j in range(4)] for i in range(4)]
    for _ in range(int(gap)):
        out = _one_step(out, x, sigma2, R_aw, R_S_z)
    return out


def _segment_images(P, gap: int, x: Interval, sigma2: float,
                    R_aw: float, R_S_z: float, depth: int = 0):
    try:
        out = _propagate(P, gap, x, sigma2, R_aw, R_S_z)
        if not symmetric_positive_definite_ldlt(out)[0]:
            raise RuntimeError("segment endpoint lower is not SPD")
        return [(x, out)]
    except (RuntimeError, ValueError) as exc:
        if depth >= MAX_X_DEPTH:
            raise RuntimeError(
                f"cannot certify commit-segment endpoint on x={x.as_list()} gap={gap}: {exc}"
            ) from exc
        halves = _split_x(x)
        if halves is None:
            raise
        left, right = halves
        return (
            _segment_images(P, gap, left, sigma2, R_aw, R_S_z, depth + 1)
            + _segment_images(P, gap, right, sigma2, R_aw, R_S_z, depth + 1)
        )


def _minus(A, alpha: float, B):
    q = BASE.I(float(alpha))
    return BASE._sym([
        [A[i][j] - q * B[i][j] for j in range(4)]
        for i in range(4)
    ])


def _relative_floor(A, B) -> float:
    if not symmetric_positive_definite_ldlt(A)[0]:
        return 0.0
    if not symmetric_positive_definite_ldlt(B)[0]:
        raise RuntimeError("common reference is not SPD")
    if symmetric_positive_definite_ldlt(_minus(A, 1.0, B))[0]:
        return 1.0
    lo, hi = 0.0, 1.0
    for _ in range(54):
        mid = 0.5 * (lo + hi)
        if symmetric_positive_definite_ldlt(_minus(A, mid, B))[0]:
            lo = mid
        else:
            hi = mid
    return BASE.down(lo)


def _common_boundary(images):
    if not images:
        raise RuntimeError("empty segment image family")
    # The lowest-x / shortest-gap image is only a reference shape.  Every image
    # is explicitly compared, so no monotonicity in tau or gap is assumed.
    ref = images[0][2]
    if not symmetric_positive_definite_ldlt(ref)[0]:
        raise RuntimeError("reference segment endpoint is not SPD")
    alpha = 1.0
    limiting = (images[0][0], images[0][1])
    for gap, x, A in images[1:]:
        a = _relative_floor(A, ref)
        if a < alpha:
            alpha = a
            limiting = (gap, x)
    if not alpha > 0.0:
        raise RuntimeError("no positive common Loewner lower across current-source segment images")
    common = BASE._scale(ref, BASE.down(alpha))
    if not symmetric_positive_definite_ldlt(common)[0]:
        raise RuntimeError("common boundary lower lost SPD")
    return common, alpha, limiting


def build(domain_path: Path = DEFAULT_DOMAIN, tube_path: Path | None = None) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("finite-word translation proof may not be trajectory fitted")

    dynamic = DYNAMIC.build(path)
    df = DYNAMIC.validate(dynamic)
    if df:
        raise RuntimeError(f"SEA3 dynamic source prerequisite failed: {df}")
    tube = BASE._load_tube(path, tube_path)

    inv = dynamic["dynamic_invariant"]
    rates = dynamic["validated_rate_and_jump_bounds"]
    h = float(rates["dt_s"])
    tau_lo, tau_hi = map(float, inv["tau_applied_s"])
    sigma_lo = float(inv["sigma_aw_filter_mps2"][0])
    if not (h > 0.0 and tau_lo > 0.0 and tau_hi >= tau_lo and sigma_lo > 0.0):
        raise RuntimeError("SEA3 dynamic invariant lost positive translation source")

    adapt_gap_upper = int(rates["active_commit_gap_samples_upper"])
    if adapt_gap_upper != 22:
        raise RuntimeError("shipping commit-gap upper changed; re-audit commit alphabet")
    # 20 is deliberately admitted although the strict >0.1 s scheduler normally
    # makes 21 the first firing.  This absorbs binary32 equality/clock details
    # without depending on a lower-gap floating-point proof.
    if max(COMMIT_GAPS) != adapt_gap_upper or min(COMMIT_GAPS) * h < 0.099:
        raise RuntimeError("commit gap alphabet no longer covers the shipping timing safely")

    global_x = Interval.outward_bounds(BASE.down(h / tau_hi), BASE.up(h / tau_lo))
    xcells = _geom_cells(global_x.lo, global_x.hi, INITIAL_X_CELLS)

    vector = VECTOR.build()
    vf = VECTOR.validate(vector)
    if vf:
        raise RuntimeError(f"vector measurement prerequisite failed: {vf}")
    R_aw = BASE.down(float(vector["configured_measurement_bounds"]["acc_measurement_std_mps2"]) ** 2)
    axis_factor = min(map(float, TUBE_BASE._axis_factors()))
    rs_lo = float(inv["R_S_applied"][0])
    rS = BASE.down(rs_lo * axis_factor)
    R_S_z = BASE.down(BASE.down(rS * rS) / BASE.up(h ** 6))
    sigma2 = BASE.down(sigma_lo * sigma_lo)
    if not (R_aw > 0.0 and R_S_z > 0.0 and sigma2 > 0.0):
        raise RuntimeError("translation lower constants lost positivity")

    P = BASE._zero()
    segment_rows = []
    max_images = 0
    min_alpha = 1.0
    for segment in range(1, WORD_SEGMENTS + 1):
        images = []
        for gap in COMMIT_GAPS:
            for cell in xcells:
                for xx, endpoint in _segment_images(P, gap, cell, sigma2, R_aw, R_S_z):
                    images.append((gap, xx, endpoint))
        P, alpha, limiting = _common_boundary(images)
        max_images = max(max_images, len(images))
        min_alpha = min(min_alpha, alpha)
        segment_rows.append({
            "segment": segment,
            "current_source_images": len(images),
            "common_scale": alpha,
            "limiting_gap_samples": int(limiting[0]),
            "limiting_x_h_over_tau": limiting[1].as_list(),
            "common_diagonal_z": [P[i][i].lo for i in range(4)],
        })

    pdiag = list(map(float, tube["modes"]["H"]["Pbar_diagonal_variance_upper"]))
    upper = [pdiag[6], pdiag[9], pdiag[12], pdiag[15]]
    dscale = [h, h * h, h * h * h, 1.0]
    upper_z = [BASE.up(upper[i] / BASE.down(dscale[i] * dscale[i])) for i in range(4)]
    rho = BASE._generalized_rho(P, upper_z)
    passed = rho >= USEFUL_GATE

    word_min_s = BASE.down(WORD_SEGMENTS * min(COMMIT_GAPS) * h)
    word_max_s = BASE.up(WORD_SEGMENTS * max(COMMIT_GAPS) * h)
    pe_window = float(domain["normal_live"]["vector_pe_recurrence_window_s"])

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_shrunk": False,
        "source_history_graph_consumed": False,
        "predecessor_path_enumeration_consumed": False,
        "old_P2_800_state_graph_consumed": False,
        "old_P2_history_frontier_consumed": False,
        "P2_correlation_interface_consumed": False,
        "SEA3_dynamic_source_consumed": True,
        "commit_aligned_word": True,
        "commit_gap_samples_admitted": list(COMMIT_GAPS),
        "commit_segments_per_word": WORD_SEGMENTS,
        "word_duration_s_lower": word_min_s,
        "word_duration_s_upper": word_max_s,
        "vector_PE_recurrence_window_s": pe_window,
        "word_covers_PE_recurrence": word_min_s >= BASE.down(pe_window),
        "instantaneous_source_cell_partition_consumed": True,
        "history_depth": 0,
        "cell_index_forgotten_after_each_segment": True,
        "tau_constant_only_inside_shipping_commit_segment": True,
        "tau_path_correlation_assumed_across_segments": False,
        "sigma_path_correlation_assumed": False,
        "R_S_path_correlation_assumed": False,
        "full_4x4_translation_matrix_retained": True,
        "intermediate_per_sample_SPD_required": False,
        "segment_endpoint_SPD_required": True,
        "strongest_accelerometer_measurement_each_sample_for_lower": True,
        "strongest_S_zero_measurement_each_sample_for_lower": True,
        "nuisance_states_conditioned_known_for_translation_lower": True,
        "fixed_physical_scaling": "z=[v/h,p/h^2,S/h^3,a_w]",
        "global_x_h_over_tau": global_x.as_list(),
        "initial_x_cells": INITIAL_X_CELLS,
        "maximum_adaptive_x_depth": MAX_X_DEPTH,
        "maximum_segment_image_count": max_images,
        "minimum_cross_cell_common_scale": min_alpha,
        "segment_rows": segment_rows,
        "translation_covariance_upper": upper,
        "translation_covariance_upper_z": upper_z,
        "endpoint_common_lower_diagonal_z": [P[i][i].lo for i in range(4)],
        "relative_word_injection_floor_lower": rho,
        "useful_gate": USEFUL_GATE,
        "useful_margin_pass": passed,
        "theorem_identity": {
            "finite_word_concavity": "R_W(P)-D R_W(P)[P] >= R_W(0)",
            "selected_process_lower": "R_W(0) >= L_translation_commit_word",
            "comparison": "L_translation_commit_word >= delta_translation * Pbar_translation",
        },
        "pass": passed,
    }


def validate(d: dict) -> list[str]:
    f = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    for key in (
        "source_generated_not_trajectory_fit", "SEA3_dynamic_source_consumed",
        "commit_aligned_word", "word_covers_PE_recurrence",
        "instantaneous_source_cell_partition_consumed",
        "cell_index_forgotten_after_each_segment",
        "tau_constant_only_inside_shipping_commit_segment",
        "full_4x4_translation_matrix_retained", "segment_endpoint_SPD_required",
        "strongest_accelerometer_measurement_each_sample_for_lower",
        "strongest_S_zero_measurement_each_sample_for_lower",
        "nuisance_states_conditioned_known_for_translation_lower",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_shrunk",
        "source_history_graph_consumed", "predecessor_path_enumeration_consumed",
        "old_P2_800_state_graph_consumed", "old_P2_history_frontier_consumed",
        "P2_correlation_interface_consumed",
        "tau_path_correlation_assumed_across_segments",
        "sigma_path_correlation_assumed", "R_S_path_correlation_assumed",
        "intermediate_per_sample_SPD_required",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("history_depth") != 0:
        f.append("commit word acquired source-history depth")
    if d.get("commit_gap_samples_admitted") != list(COMMIT_GAPS):
        f.append("commit gap alphabet changed")
    if d.get("commit_segments_per_word") != WORD_SEGMENTS:
        f.append("commit word segment count changed")
    if float(d.get("useful_gate", math.nan)) != USEFUL_GATE:
        f.append("translation useful gate changed")
    rho = d.get("relative_word_injection_floor_lower")
    if not isinstance(rho, (int, float)) or not math.isfinite(float(rho)) or float(rho) < 0.0:
        f.append("translation finite-word floor is not finite nonnegative")
    expected = isinstance(rho, (int, float)) and float(rho) >= USEFUL_GATE
    if d.get("useful_margin_pass") is not expected or d.get("pass") is not expected:
        f.append("translation pass flag disagrees with finite-word margin")
    rows = d.get("segment_rows", [])
    if len(rows) != WORD_SEGMENTS or any(float(r.get("common_scale", 0.0)) <= 0.0 for r in rows):
        f.append("commit-segment common lower did not remain strict")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--tube", type=Path, default=None)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain, args.tube)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "word_duration_s": [d["word_duration_s_lower"], d["word_duration_s_upper"]],
        "max_segment_images": d["maximum_segment_image_count"],
        "min_common_scale": d["minimum_cross_cell_common_scale"],
        "translation_delta": d["relative_word_injection_floor_lower"],
        "useful_gate": d["useful_gate"],
        "pass": d["pass"],
        "endpoint_diag": d["endpoint_common_lower_diagonal_z"],
        "segments": d["segment_rows"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
