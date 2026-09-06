#!/usr/bin/env python3
"""Run finite nonlinear same-history P4 feasibility checks along linear worst directions.

This is a non-promoting diagnostic.  It consumes the worst H18/A21 directions
from ``ou3_p4_complete_sea3_word_feasibility.py`` and injects scaled versions
into the retained host-only ``ou3-neighborhood-sim``.  Both estimators receive
the same genuine coupled measurements and the observer rejects a case from the
same-history set as soon as Live/mode/acceptance/tuner/pseudo-period histories
differ.  In particular, the nominal and perturbed filters must retain the same
actual applied ``R_S`` sequence.

The scales are restricted by the already-declared startup/error envelope and
the widest P4 attitude candidate; this script does not widen the theorem domain.
A rho >= 1 finding is reported, never hidden by changing source language,
retuning the filter, or replacing due S updates.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import subprocess
from pathlib import Path

GROUPS = {
    "theta": slice(0, 3),
    "bg": slice(3, 6),
    "v": slice(6, 9),
    "p": slice(9, 12),
    "S": slice(12, 15),
    "aw": slice(15, 18),
    "ba": slice(18, 21),
}


def norm(v: list[float]) -> float:
    return math.sqrt(sum(x * x for x in v))


def group_norm(x: list[float], name: str) -> float:
    s = GROUPS[name]
    return norm(x[s])


def parse_hs(path: Path) -> float | None:
    m = re.search(r"_H([0-9]+(?:\.[0-9]+)?)_", path.name)
    return float(m.group(1)) if m else None


def max_scale_in_declared_envelope(direction: list[float], mode: str, domain: dict, input_path: Path) -> tuple[float, dict]:
    handoff = domain["startup"]["physical_handoff_coordinate_bounds"]
    angle_deg = max(float(x) for x in domain["certificate_search"]["p4_complete_word_full_attitude_candidate_deg"])
    caps = {
        "theta": math.radians(angle_deg),
        "bg": float(handoff["gyro_bias_error_norm_upper_rad_s"]),
        "v": float(handoff["velocity_error_norm_upper_mps"]),
        "p": float(handoff["position_error_norm_upper_m"]),
        "S": float(handoff["integral_displacement_error_norm_upper_m_s"]),
        "aw": float(handoff["latent_acceleration_error_norm_upper_mps2"]),
        "ba": float(handoff["accelerometer_bias_error_norm_upper_mps2"]),
    }
    if mode == "H18":
        caps.pop("ba")
    candidates: dict[str, float] = {}
    for name, cap in caps.items():
        g = group_norm(direction, name)
        if g > 0.0:
            candidates[name] = cap / g

    # P5 entrance additionally declares componentwise |delta p_i| <= 0.5 Hs.
    # Preserve that restriction in this diagnostic when Hs is encoded in the
    # genuine source filename.  This is stricter than the 20 m handoff norm.
    hs = parse_hs(input_path)
    if hs is not None:
        p_cap = 0.5 * hs
        for j, value in enumerate(direction[9:12]):
            if abs(value) > 0.0:
                candidates[f"p_component_{j}"] = p_cap / abs(value)

    if not candidates:
        raise RuntimeError(f"zero maximizing direction for {mode}")
    limiting = min(candidates, key=candidates.get)
    return float(candidates[limiting]), {
        "limiting_constraint": limiting,
        "all_scale_limits": candidates,
        "widest_attitude_candidate_deg": angle_deg,
        "Hs_m": hs,
    }


def choose_scales(max_scale: float) -> list[float]:
    base = [0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0]
    out = [x for x in base if x <= max_scale * (1.0 + 1e-12)]
    if not out or max_scale > out[-1] * (1.0 + 1e-6):
        out.append(max_scale)
    uniq: list[float] = []
    for x in out:
        if not uniq or abs(x - uniq[-1]) > 1e-6 * max(1.0, abs(x)):
            uniq.append(x)
    return uniq


def read_trace(path: Path) -> dict:
    with path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    if len(rows) < 2:
        raise RuntimeError(f"nonlinear trace has fewer than two rows: {path}")
    first = rows[0]
    endpoints = [r for r in rows if int(r["endpoint"]) == 1]
    if not endpoints:
        raise RuntimeError(f"nonlinear trace has no endpoint: {path}")
    last = endpoints[-1]
    injection_cov = float(first["covariance_rel_fro"])
    covariance_identical_at_injection = math.isfinite(injection_cov) and injection_cov <= 2.0e-7
    source_match = (
        covariance_identical_at_injection
        and all(int(r["source_match"]) == 1 for r in rows)
    )
    w0 = float(first["W_nominal"])
    w1 = float(last["W_nominal"])
    rho = w1 / w0 if w0 > 0.0 and math.isfinite(w0) and math.isfinite(w1) else None
    return {
        "rows": len(rows),
        "source_match": source_match,
        "covariance_identical_at_injection": covariance_identical_at_injection,
        "injection_covariance_rel_fro": injection_cov,
        "W0_nominal_metric": w0,
        "W1_nominal_metric": w1,
        "rho_nonlinear": rho,
        "distance_to_one": 1.0 - rho if rho is not None else None,
        "initial_theta_rad": float(first["theta_rad"]),
        "endpoint_theta_rad": float(last["theta_rad"]),
        "endpoint_error_norm": float(last["error_norm"]),
        "endpoint_covariance_rel_fro": float(last["covariance_rel_fro"]),
        "nominal_RS_start": float(first["nom_rs"]),
        "nominal_RS_end": float(last["nom_rs"]),
        "perturbed_RS_start": float(first["pert_rs"]),
        "perturbed_RS_end": float(last["pert_rs"]),
        "acc_accept_match_all": all(int(r["acc_accept_match"]) == 1 for r in rows),
        "mag_accept_match_all": all(int(r["mag_accept_match"]) == 1 for r in rows),
    }


def run_case(sim: Path, input_path: Path, mode: str, t0: float, direction: list[float], scale: float, out_dir: Path) -> dict:
    delta = [scale * x for x in direction]
    if mode == "H18":
        delta[18:21] = [0.0, 0.0, 0.0]
    trace = out_dir / f"nonlinear_{mode}_scale_{scale:.12g}.csv"
    env = os.environ.copy()
    env.update({
        "OU3_NEIGHBOR_TRACE": str(trace),
        "OU3_NEIGHBOR_DELTA": ",".join(f"{x:.17g}" for x in delta),
        "OU3_NEIGHBOR_INJECT_TIME_S": f"{t0:.17g}",
        "OU3_NEIGHBOR_HORIZON_S": "3.0",
        "OU3_NEIGHBOR_MODE": "H" if mode == "H18" else "A",
        "OU3_NEIGHBOR_TRACE_STRIDE": "10",
        "W3D_WRITE_TIMESERIES": "0",
        "W3D_VALIDATION_WINDOW_SEC": "0",
    })
    cp = subprocess.run(
        [str(sim), "--input", str(input_path)],
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    result = {
        "mode": mode,
        "scale": scale,
        "requested_injection_s": t0,
        "delta": delta,
        "delta_group_norms": {name: group_norm(delta, name) for name in GROUPS},
        "returncode": cp.returncode,
        "stdout_tail": "\n".join(cp.stdout.splitlines()[-12:]),
        "trace": str(trace),
    }
    if cp.returncode == 0 and trace.exists():
        result.update(read_trace(trace))
    else:
        result.update({
            "source_match": False,
            "covariance_identical_at_injection": False,
            "rho_nonlinear": None,
            "distance_to_one": None,
        })
    return result


def finite_case(c: dict) -> bool:
    rho = c.get("rho_nonlinear")
    return c.get("source_match") is True and isinstance(rho, (int, float)) and math.isfinite(float(rho))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--linear-json", type=Path, required=True)
    ap.add_argument("--domain", type=Path, required=True)
    ap.add_argument("--sim", type=Path, required=True)
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    linear = json.loads(args.linear_json.read_text(encoding="utf-8"))
    domain = json.loads(args.domain.read_text(encoding="utf-8"))
    args.output_dir.mkdir(parents=True, exist_ok=True)

    report = {
        "qualification": "NON_PROMOTING_COMPLETE_SEA3_P4_NONLINEAR_SAME_HISTORY_FEASIBILITY",
        "canonical_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "point_same_history_diagnostic_only": True,
        "P4_promoted": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "same_shipping_covariance_at_injection_required": True,
        "same_actual_applied_RS_history_required": True,
        "same_accelerometer_acceptance_history_required": True,
        "same_vector_acceptance_history_required": True,
        "packet_count_remainder_budget_used": False,
        "selected_S_subset_used": False,
        "input": str(args.input),
        "modes": {},
    }

    for mode in ("H18", "A21"):
        worst = linear["modes"][mode]["worst"]
        direction = [float(x) for x in worst["maximizing_direction"]["components"]]
        max_scale, limit = max_scale_in_declared_envelope(direction, mode, domain, args.input)
        scales = choose_scales(max_scale)
        cases = [
            run_case(args.sim.resolve(), args.input.resolve(), mode, float(worst["t0"]), direction, s, args.output_dir.resolve())
            for s in scales
        ]
        valid = [c for c in cases if finite_case(c)]
        report["modes"][mode] = {
            "rho_linear_worst_word": float(worst["rho_linear"]),
            "linear_distance_to_one": float(worst["distance_to_one"]),
            "word_t0": float(worst["t0"]),
            "word_t1": float(worst["t1"]),
            "S_update_count": int(worst["S_update_count"]),
            "acc_count": int(worst["acc_count"]),
            "mag_count": int(worst["mag_count"]),
            "RS_scalar_start": float(worst["RS_scalar_start"]),
            "RS_scalar_end": float(worst["RS_scalar_end"]),
            "declared_scale_limit": max_scale,
            "declared_scale_limit_detail": limit,
            "cases": cases,
            "same_history_cases": len(valid),
            "worst_same_history_nonlinear_rho": max((float(c["rho_nonlinear"]) for c in valid), default=None),
            "largest_same_history_scale": max((float(c["scale"]) for c in valid), default=None),
            "strict_contraction_on_all_same_history_cases": bool(valid) and all(float(c["rho_nonlinear"]) < 1.0 for c in valid),
        }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        mode: {
            "linear_rho": report["modes"][mode]["rho_linear_worst_word"],
            "scale_limit": report["modes"][mode]["declared_scale_limit"],
            "limiter": report["modes"][mode]["declared_scale_limit_detail"]["limiting_constraint"],
            "same_history_cases": report["modes"][mode]["same_history_cases"],
            "largest_same_history_scale": report["modes"][mode]["largest_same_history_scale"],
            "worst_nonlinear_rho": report["modes"][mode]["worst_same_history_nonlinear_rho"],
            "all_contract": report["modes"][mode]["strict_contraction_on_all_same_history_cases"],
        }
        for mode in ("H18", "A21")
    }, indent=2, sort_keys=True))
    return 0 if all(report["modes"][m]["same_history_cases"] > 0 for m in ("H18", "A21")) else 2


if __name__ == "__main__":
    raise SystemExit(main())