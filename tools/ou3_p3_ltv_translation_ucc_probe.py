#!/usr/bin/env python3
"""Source-uniform time-varying translation process-floor probe for OU-III P3.

This module asks one narrow question before canonical P3 is changed: does the
translation process inject enough *directional* covariance over a finite suffix
when tau/sigma may vary arbitrarily inside the source word?

The one-axis normal-Live translation chain, ordered [S,p,v,a_w], is

    S' = p,  p' = v,  v' = a_w,
    a_w' = -lambda(t) a_w + sqrt(q_c(t)) w,

with lambda(t)=1/tau(t).  No frozen-parameter assumption is used below.
For four process-input times s0<s1<s2<s3, multiply each endpoint response
column by exp(-Lambda(si)), Lambda(t)=int_0^t lambda.  After consecutive column
subtractions the determinant is an integral of the positive polynomial
Vandermonde against exp(-Lambda).  Restoring the four column factors cancels
all but at most one endpoint-to-start decay interval, giving the same robust
bound as the constant-lambda calculation,

    |det K(s0..s3)| >= V(s0..s3) exp(-lambda_max H) / 12.

Thus the existing four-subinterval Andreief construction remains valid for
arbitrary measurable lambda(t) in the declared source box:

    det G_unit >= (2025/144) (H/7)^16 exp(-2 lambda_max H).

Because q_c(t)>=2 sigma_min^2/tau_max, the actual process Gramian dominates
q_c,min G_unit.  The calculation is normalized directly by the endpoint P3
translation covariance dominator diag(Uv,Up,US,Ua), avoiding a global covariance
condition-number conversion.

This is intentionally a *pre-measurement* process-floor probe.  Interleaved
accelerometer/S updates can reduce the floor and must be enclosed separately
before P3 promotion.  A positive/useful result here therefore establishes only
that time-varying controllability itself is not the remaining bottleneck.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p3_scaled_process as SCALED
import ou3_p4_source_node_cells as NODES
import ou3_source_reachable_matrix_p3 as BASE

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def ltv_relative_process_floor(upper: list[float], horizon_s: float,
                               tau_min: float, tau_max: float,
                               sigma_min: float) -> dict:
    """Return rho with W_process >= rho*diag(upper) on [v,p,S,a_w]."""
    if len(upper) != 4 or any(not (math.isfinite(float(v)) and float(v) > 0.0) for v in upper):
        raise ValueError("four finite positive translation covariance uppers required")
    H = float(horizon_s)
    if not (H > 0.0 and 0.0 < tau_min <= tau_max and sigma_min > 0.0):
        raise ValueError("positive LTV source bounds required")

    lam_max = up(1.0 / tau_min)
    qc_min = down(2.0 * sigma_min * sigma_min / tau_max)
    width = down(H / 7.0)
    det_unit = down(
        (2025.0 / 144.0) * (width ** 16) * math.exp(-up(2.0 * lam_max * H))
    )
    if not det_unit > 0.0:
        return {"relative_process_floor_lower": 0.0}

    Uv, Up, US, Ua = map(float, upper)
    det_normalized = down(det_unit / up(up(Uv * Up) * up(US * Ua)))

    # Integral response bounds valid for arbitrary nonnegative lambda(t):
    # |a|<=1, |v|<=r, |p|<=r^2/2, |S|<=r^3/6.
    trace_normalized = up(
        H / Ua
        + (H ** 3) / (3.0 * Uv)
        + (H ** 5) / (20.0 * Up)
        + (H ** 7) / (252.0 * US)
    )
    if not (det_normalized > 0.0 and trace_normalized > 0.0):
        return {"relative_process_floor_lower": 0.0}
    gram_lambda = down(det_normalized / up(trace_normalized ** 3))
    rho = down(qc_min * gram_lambda)
    return {
        "horizon_s": H,
        "lambda_max_per_s_upper": lam_max,
        "q_c_min_lower": qc_min,
        "unit_gramian_det_lower": det_unit,
        "Sigma_normalized_gramian_det_lower": det_normalized,
        "Sigma_normalized_gramian_trace_upper": trace_normalized,
        "Sigma_normalized_unit_gramian_lambda_min_lower": gram_lambda,
        "relative_process_floor_lower": rho,
    }


def _mode_node(mode: str, node: dict, domain: dict, candidates: list[float]) -> dict:
    live = domain["normal_live"]
    vector = BASE.VECTOR.build()
    process = BASE.PROCESS.build()
    sched = BASE.source_schedule()
    h = float(sched["dt_s"])
    tau_lo, tau_hi = map(float, node["tau_s"])
    x = Interval.outward_bounds(BASE.down(h / tau_hi), BASE.up(h / tau_lo))
    sigma = Interval(*map(float, node["sigma_filter_committed_mps2"]))
    rs = Interval(*map(float, node["R_S_filter_std"]))
    alpha6 = BASE.vector_alpha6(live, vector)

    global_tau_lo, global_tau_hi = map(float, sched["tau_applied_invariant_s"])
    sigma_min = float(NODES.build()["filter_sigma_floor_mps2"])
    rows = []
    for xcell, rho_x in SCALED.split_x_cell(x):
        raw = BASE.mode_cell(mode, xcell, rho_x, sigma, rs, live, vector, process, sched, alpha6)
        upper = [
            float(raw["Sigma_diagonal_upper"][6]),
            float(raw["Sigma_diagonal_upper"][9]),
            float(raw["Sigma_diagonal_upper"][12]),
            float(raw["Sigma_diagonal_upper"][15]),
        ]
        word_lo = float(raw["word_horizon_s_lower"])
        hs = [H for H in candidates if H <= word_lo]
        if not hs:
            raise RuntimeError("no LTV probe horizon fits the certified P3 word")
        probes = [
            ltv_relative_process_floor(upper, H, global_tau_lo, global_tau_hi, sigma_min)
            for H in hs
        ]
        best = max(probes, key=lambda d: d["relative_process_floor_lower"])
        rows.append({
            "x_h_over_tau": xcell.as_list(),
            "word_horizon_s_lower": word_lo,
            "Sigma_translation_diagonal_upper": upper,
            "best": best,
            "candidates": probes,
        })
    worst = min(rows, key=lambda r: r["best"]["relative_process_floor_lower"])
    rho = float(worst["best"]["relative_process_floor_lower"])
    return {
        "mode": mode,
        "relative_premeasurement_process_floor_lower": rho,
        "useful_gate": BASE.MIN_USEFUL_DELTA,
        "premeasurement_floor_above_useful_gate": rho >= BASE.MIN_USEFUL_DELTA,
        "worst_x_subcell": worst,
        "x_subcell_count": len(rows),
    }


def build(domain_path: Path = DEFAULT_DOMAIN, source_node_indices=(0, 729)) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("LTV UCC probe must not be trajectory fitted")
    nodes = NODES.build()
    nf = NODES.validate(nodes)
    if nf:
        raise RuntimeError(f"P2 source nodes invalid: {nf}")

    candidates = [0.25, 0.5, 0.75, 1.0, 1.25, 1.5]
    results = {}
    for index in source_node_indices:
        node = NODES.node(int(index), nodes)
        results[str(index)] = {
            "source_node": node,
            "H": _mode_node("H", node, domain, candidates),
            "A": _mode_node("A", node, domain, candidates),
        }

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P3_SOURCE_UNIFORM_LTV_TRANSLATION_PROCESS_FLOOR_PROBE",
        "source_generated_not_trajectory_fit": True,
        "arbitrary_time_varying_tau_inside_window_covered": True,
        "arbitrary_time_varying_sigma_inside_window_covered_by_qc_min": True,
        "frozen_parameter_Q_Nh_identity_used": False,
        "endpoint_covariance_condition_number_conversion_used": False,
        "interleaved_measurement_attenuation_enclosed_here": False,
        "P3_PROMOTED": False,
        "candidate_horizons_s": candidates,
        "nodes": results,
        "next_obligation": (
            "enclose the interleaved normal-Live measurement attenuation of this LTV process floor; "
            "only the resulting post-measurement source-uniform H/A comparison may replace canonical P3"
        ),
    }


def validate(d: dict) -> list[str]:
    f = []
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for key in (
        "source_generated_not_trajectory_fit",
        "arbitrary_time_varying_tau_inside_window_covered",
        "arbitrary_time_varying_sigma_inside_window_covered_by_qc_min",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "frozen_parameter_Q_Nh_identity_used",
        "endpoint_covariance_condition_number_conversion_used",
        "interleaved_measurement_attenuation_enclosed_here",
        "P3_PROMOTED",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    for node, row in d.get("nodes", {}).items():
        for mode in ("H", "A"):
            rho = row.get(mode, {}).get("relative_premeasurement_process_floor_lower")
            if not isinstance(rho, (int, float)) or not (math.isfinite(float(rho)) and float(rho) > 0.0):
                f.append(f"node {node} {mode}: LTV floor is not strict")
    return f


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--nodes", type=int, nargs="*", default=[0, 729])
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain, args.nodes)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    summary = {}
    for node, row in d["nodes"].items():
        summary[node] = {
            mode: {
                "rho_premeasurement": row[mode]["relative_premeasurement_process_floor_lower"],
                "useful": row[mode]["premeasurement_floor_above_useful_gate"],
                "best_horizon_s": row[mode]["worst_x_subcell"]["best"]["horizon_s"],
            }
            for mode in ("H", "A")
        }
    print(json.dumps({"validation_failures": vf, "nodes": summary}, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
