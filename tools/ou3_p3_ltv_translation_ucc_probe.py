#!/usr/bin/env python3
"""Source-uniform time-varying translation process-floor probe for OU-III P3.

For one translation axis ordered [S,p,v,a_w], normal-Live dynamics are
S'=p, p'=v, v'=a_w, a_w'=-lambda(t)a_w+sqrt(q_c(t))w.  No frozen tuner
parameter is assumed.  Applying the integrating factor to four endpoint
response columns and taking consecutive column differences turns their
determinant into a positive polynomial-Vandermonde integral.  Restoring the
column factors leaves one source decay integral, hence

    |det K(s0..s3)| >= V(s0..s3) exp(-int lambda)/12,

and the four-subinterval Andreief construction gives

    det G_unit >= (2025/144)(H/7)^16 exp(-2 int lambda).

The active decay bound comes from the retained 13..26-sample staged tuner path,
not from the Cartesian H/tau_min shortcut.  Process intensity is bounded by the
source-uniform q_c(t)>=2 sigma_min^2/tau_max.

Two comparisons are emitted:

* exact-node endpoint-static covariance uppers for nodes 0 and 729 remain only
  diagnostics, to preserve continuity with the earlier whole-word experiment;
* the active premeasurement comparison uses
  :mod:`ou3_p3_source_uniform_translation_covariance`, a changing-parameter
  finite-memory covariance upper valid on the complete invariant source box,
  and scans all ten endpoint tau cells through the clock-phase decay budget.

That second result is source complete for the covariance upper and process
floor.  It still stops before interleaved accelerometer/S measurement
attenuation, so it cannot promote P3 by itself.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p3_scaled_process as SCALED
import ou3_p3_source_uniform_translation_covariance as COV
import ou3_p3_tau_decay_budget as DECAY
import ou3_p4_source_node_cells as NODES
import ou3_source_reachable_matrix_p3 as BASE
import ou3_translational_uco_ucc as TRANS

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 4


def point(x: float) -> Interval:
    return Interval.outward_bounds(float(x), float(x))


def ltv_relative_process_floor(upper: list[float], horizon_s: float,
                               tau_min: float, tau_max: float,
                               sigma_min: float, *,
                               decay_exponent_upper: float | None = None) -> dict:
    """Return rho such that W_process >= rho*diag(upper) on [v,p,S,a_w]."""
    if len(upper) != 4 or any(not (math.isfinite(float(v)) and float(v) > 0.0) for v in upper):
        raise ValueError("four finite positive translation covariance uppers required")
    H = float(horizon_s)
    if not (H > 0.0 and 0.0 < tau_min <= tau_max and sigma_min > 0.0):
        raise ValueError("positive LTV source bounds required")

    if decay_exponent_upper is None:
        decay_exponent_upper = math.nextafter(H / tau_min, math.inf)
        decay_route = "GLOBAL_H_OVER_TAU_MIN_FALLBACK"
    else:
        decay_exponent_upper = float(decay_exponent_upper)
        decay_route = "CLOCK_PHASE_TAU_PATH_BUDGET"
    if not (math.isfinite(decay_exponent_upper) and decay_exponent_upper > 0.0):
        raise ValueError("positive finite decay exponent required")

    Hiv = point(H)
    decay = TRANS._exp_negative_wide(point(decay_exponent_upper))
    decay2_lower = point(decay.lo).square().lo
    width = Hiv / point(7.0)
    det_unit = (
        point(2025.0 / 144.0)
        * TRANS._pow_nonnegative(width, 16)
        * point(decay2_lower)
    ).lo
    qc_min = (point(2.0) * point(sigma_min).square() / point(tau_max)).lo
    if not (det_unit > 0.0 and qc_min > 0.0):
        return {"relative_process_floor_lower": 0.0}

    Uv, Up, US, Ua = map(float, upper)
    product_U = point(Uv) * point(Up) * point(US) * point(Ua)
    det_normalized = (point(det_unit) / product_U).lo

    # Arbitrary nonnegative lambda only shortens these endpoint responses.
    H3 = TRANS._pow_nonnegative(Hiv, 3)
    H5 = TRANS._pow_nonnegative(Hiv, 5)
    H7 = TRANS._pow_nonnegative(Hiv, 7)
    trace_normalized = (
        Hiv / point(Ua)
        + H3 / point(3.0 * Uv)
        + H5 / point(20.0 * Up)
        + H7 / point(252.0 * US)
    ).hi
    trace3 = TRANS._pow_nonnegative(point(trace_normalized), 3)
    if not (det_normalized > 0.0 and trace_normalized > 0.0):
        return {"relative_process_floor_lower": 0.0}
    gram_lambda = (point(det_normalized) / trace3).lo
    rho = (point(qc_min) * point(gram_lambda)).lo
    return {
        "horizon_s": H,
        "decay_route": decay_route,
        "decay_exponent_upper": decay_exponent_upper,
        "exp_minus_decay_lower": decay.lo,
        "q_c_min_lower": qc_min,
        "unit_gramian_det_lower": det_unit,
        "Sigma_normalized_gramian_det_lower": det_normalized,
        "Sigma_normalized_gramian_trace_upper": trace_normalized,
        "Sigma_normalized_unit_gramian_lambda_min_lower": gram_lambda,
        "relative_process_floor_lower": rho,
        "validated_interval_arithmetic": True,
        "validated_exponential_arithmetic": True,
    }


def _mode_node(mode: str, node: dict, domain: dict, candidates: list[float],
               sigma_floor: float, domain_path: Path) -> dict:
    """Historical exact-node diagnostic; not the active source-complete route."""
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
    endpoint_tau_index = int(node["tau_index"])

    rows = []
    for xcell, rho_x in SCALED.split_x_cell(x):
        raw = BASE.mode_cell(mode, xcell, rho_x, sigma, rs, live, vector, process, sched, alpha6)
        upper = [float(raw["Sigma_diagonal_upper"][i]) for i in (6, 9, 12, 15)]
        word_lo = float(raw["word_horizon_s_lower"])
        hs = [H for H in candidates if H <= word_lo]
        if not hs:
            raise RuntimeError("no LTV probe horizon fits the certified P3 word")
        probes = []
        for H in hs:
            samples = max(1, int(math.ceil(H / h)))
            budget = DECAY.decay_budget(endpoint_tau_index, samples, domain_path)
            row = ltv_relative_process_floor(
                upper, H, global_tau_lo, global_tau_hi, sigma_floor,
                decay_exponent_upper=float(budget["decay_exponent_upper"]),
            )
            row["clock_phase_decay_budget"] = budget
            probes.append(row)
        best = max(probes, key=lambda d: d["relative_process_floor_lower"])
        rows.append({
            "x_h_over_tau": xcell.as_list(),
            "word_horizon_s_lower": word_lo,
            "Sigma_translation_diagonal_upper": upper,
            "endpoint_tau_index": endpoint_tau_index,
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


def _source_complete_translation(path: Path, candidates: list[float], sigma_floor: float) -> dict:
    cov = COV.build(path)
    cf = COV.validate(cov)
    if cf:
        raise RuntimeError(f"source-uniform translation covariance upper failed: {cf}")
    upper = [float(x) for x in cov["Sigma_translation_diagonal_upper"]]
    word_lo = float(cov["timing"]["word_horizon_s_lower"])
    sched = BASE.source_schedule()
    h = float(sched["dt_s"])
    tau_lo, tau_hi = map(float, sched["tau_applied_invariant_s"])
    hs = [H for H in candidates if H <= word_lo]
    if not hs:
        raise RuntimeError("no LTV horizon fits source-uniform covariance word")

    endpoint_rows = []
    for tau_index in range(10):
        probes = []
        for H in hs:
            samples = max(1, int(math.ceil(H / h)))
            budget = DECAY.decay_budget(tau_index, samples, path)
            row = ltv_relative_process_floor(
                upper, H, tau_lo, tau_hi, sigma_floor,
                decay_exponent_upper=float(budget["decay_exponent_upper"]),
            )
            row["clock_phase_decay_budget"] = budget
            probes.append(row)
        best = max(probes, key=lambda d: d["relative_process_floor_lower"])
        endpoint_rows.append({
            "endpoint_tau_index": tau_index,
            "best": best,
            "candidates": probes,
        })
    worst = min(endpoint_rows, key=lambda r: r["best"]["relative_process_floor_lower"])
    rho = float(worst["best"]["relative_process_floor_lower"])
    return {
        "covariance_upper": cov,
        "endpoint_tau_cells_scanned": 10,
        "endpoint_rows": endpoint_rows,
        "worst_endpoint": worst,
        "relative_premeasurement_process_floor_lower": rho,
        "useful_gate": BASE.MIN_USEFUL_DELTA,
        "premeasurement_floor_above_useful_gate": rho >= BASE.MIN_USEFUL_DELTA,
        "source_complete_covariance_upper_consumed": True,
        "source_complete_process_floor_before_measurements": True,
    }


def build(domain_path: Path = DEFAULT_DOMAIN, source_node_indices=(0, 729)) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("LTV UCC probe must not be trajectory fitted")
    nodes = NODES.build()
    failures = NODES.validate(nodes)
    if failures:
        raise RuntimeError(f"P2 source nodes invalid: {failures}")
    sigma_floor = float(nodes["filter_sigma_floor_mps2"])
    candidates = [0.25, 0.5, 0.75, 1.0, 1.25, 1.5]

    source_complete = _source_complete_translation(path, candidates, sigma_floor)
    diagnostics = {}
    for index in source_node_indices:
        node = NODES.node(int(index), nodes)
        diagnostics[str(index)] = {
            "source_node": node,
            "H": _mode_node("H", node, domain, candidates, sigma_floor, path),
            "A": _mode_node("A", node, domain, candidates, sigma_floor, path),
        }
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P3_CLOCK_PHASE_LTV_TRANSLATION_PROCESS_FLOOR_PROBE",
        "source_generated_not_trajectory_fit": True,
        "validated_interval_arithmetic": True,
        "validated_exponential_arithmetic": True,
        "arbitrary_time_varying_tau_inside_window_covered": True,
        "tau_variation_uses_clock_phase_path_budget": True,
        "arbitrary_time_varying_sigma_inside_window_covered_by_qc_min": True,
        "source_complete_covariance_upper_consumed": True,
        "source_complete_process_floor_before_measurements": True,
        "frozen_parameter_Q_Nh_identity_used": False,
        "endpoint_covariance_condition_number_conversion_used": False,
        "interleaved_measurement_attenuation_enclosed_here": False,
        "P3_PROMOTED": False,
        "candidate_horizons_s": candidates,
        "source_complete_translation": source_complete,
        "diagnostic_exact_nodes": diagnostics,
        "next_obligation": (
            "enclose the interleaved accelerometer and S=0 measurement attenuation of the source-complete LTV process floor; "
            "if that post-measurement floor remains above the unchanged useful gate, combine it with the existing attitude/gyro-bias "
            "and active accelerometer-bias blocks to promote canonical P3"
        ),
    }


def validate(d: dict) -> list[str]:
    f = []
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for key in (
        "source_generated_not_trajectory_fit", "validated_interval_arithmetic",
        "validated_exponential_arithmetic", "arbitrary_time_varying_tau_inside_window_covered",
        "tau_variation_uses_clock_phase_path_budget",
        "arbitrary_time_varying_sigma_inside_window_covered_by_qc_min",
        "source_complete_covariance_upper_consumed",
        "source_complete_process_floor_before_measurements",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "frozen_parameter_Q_Nh_identity_used", "endpoint_covariance_condition_number_conversion_used",
        "interleaved_measurement_attenuation_enclosed_here", "P3_PROMOTED",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")

    sc = d.get("source_complete_translation", {})
    rho = sc.get("relative_premeasurement_process_floor_lower")
    if not isinstance(rho, (int, float)) or not (math.isfinite(float(rho)) and float(rho) > 0.0):
        f.append("source-complete translation LTV floor is not strict")
    if sc.get("source_complete_covariance_upper_consumed") is not True:
        f.append("source-complete covariance upper was not consumed")
    if sc.get("source_complete_process_floor_before_measurements") is not True:
        f.append("source-complete premeasurement process floor missing")
    if int(sc.get("endpoint_tau_cells_scanned", 0)) != 10:
        f.append("source-complete LTV scan did not cover ten tau cells")
    best = sc.get("worst_endpoint", {}).get("best", {})
    if best.get("decay_route") != "CLOCK_PHASE_TAU_PATH_BUDGET":
        f.append("source-complete worst endpoint did not use clock-phase decay budget")

    for node, row in d.get("diagnostic_exact_nodes", {}).items():
        for mode in ("H", "A"):
            x = row.get(mode, {}).get("relative_premeasurement_process_floor_lower")
            if not isinstance(x, (int, float)) or not (math.isfinite(float(x)) and float(x) > 0.0):
                f.append(f"diagnostic node {node} {mode}: LTV floor is not strict")
    return list(dict.fromkeys(f))


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
    sc = d["source_complete_translation"]
    summary = {
        "source_complete_translation": {
            "rho_premeasurement": sc["relative_premeasurement_process_floor_lower"],
            "useful": sc["premeasurement_floor_above_useful_gate"],
            "worst_tau_index": sc["worst_endpoint"]["endpoint_tau_index"],
            "best_horizon_s_at_worst": sc["worst_endpoint"]["best"]["horizon_s"],
            "decay_exponent_upper": sc["worst_endpoint"]["best"]["decay_exponent_upper"],
            "Sigma_translation_diagonal_upper": sc["covariance_upper"]["Sigma_translation_diagonal_upper"],
        },
        "diagnostic_nodes": {},
        "validation_failures": vf,
    }
    for node, row in d["diagnostic_exact_nodes"].items():
        summary["diagnostic_nodes"][node] = {
            mode: {
                "rho_premeasurement": row[mode]["relative_premeasurement_process_floor_lower"],
                "useful": row[mode]["premeasurement_floor_above_useful_gate"],
                "best_horizon_s": row[mode]["worst_x_subcell"]["best"]["horizon_s"],
            }
            for mode in ("H", "A")
        }
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
