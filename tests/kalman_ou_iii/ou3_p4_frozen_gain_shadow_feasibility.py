#!/usr/bin/env python3
"""Point nonlinear complete-word feasibility on one frozen shipping word.

The shipping estimator alone generates the complete SEA3 source, covariance,
branch decisions, gains, actual R_S sequence and endpoint metric.  The host-only
shadow has state but no covariance and therefore cannot create a second Riccati
history.  It recomputes only the nonlinear residuals under those frozen shipping
operations.  This is a falsification/feasibility experiment, never P4 promotion.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import re
import subprocess
from pathlib import Path

GROUPS = {
    "theta": slice(0, 3), "bg": slice(3, 6), "v": slice(6, 9),
    "p": slice(9, 12), "S": slice(12, 15), "aw": slice(15, 18),
    "ba": slice(18, 21),
}
DONE_RE = re.compile(
    r"OU3_FROZEN_SHADOW_DONE mode=(?P<mode>H18|A21) scale=(?P<scale>[-+0-9.eE]+) "
    r"t0=(?P<t0>[-+0-9.eE]+) t1=(?P<t1>[-+0-9.eE]+) "
    r"V0=(?P<V0>[-+0-9.eE]+) V1=(?P<V1>[-+0-9.eE]+) rho=(?P<rho>[-+0-9.eE]+) "
    r"reconstruction_max=(?P<recon>[-+0-9.eE]+) prediction_count=(?P<pred>\d+) "
    r"S_count=(?P<S>\d+) accel_count=(?P<acc>\d+) vector_count=(?P<vector>\d+)"
)


def norm(v: list[float]) -> float:
    return math.sqrt(sum(x * x for x in v))


def group_norm(v: list[float], name: str) -> float:
    s = GROUPS[name]
    return norm(v[s])


def parse_hs(path: Path) -> float | None:
    m = re.search(r"_H([0-9]+(?:\.[0-9]+)?)_", path.name)
    return float(m.group(1)) if m else None


def max_scale(direction: list[float], mode: str, domain: dict, input_path: Path) -> tuple[float, dict]:
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
    if mode == "H18": caps.pop("ba")
    limits: dict[str, float] = {}
    for name, cap in caps.items():
        g = group_norm(direction, name)
        if g > 0.0: limits[name] = cap / g
    hs = parse_hs(input_path)
    if hs is not None:
        for j, x in enumerate(direction[9:12]):
            if abs(x) > 0.0: limits[f"p_component_{j}"] = 0.5 * hs / abs(x)
    if not limits: raise RuntimeError(f"zero direction for {mode}")
    limiter = min(limits, key=limits.get)
    return float(limits[limiter]), {
        "limiting_constraint": limiter,
        "all_scale_limits": limits,
        "widest_attitude_candidate_deg": angle_deg,
        "Hs_m": hs,
    }


def choose_magnitudes(limit: float) -> list[float]:
    base = [
        0.001, 0.002, 0.004, 0.008, 0.015625, 0.03125, 0.0625,
        0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0,
    ]
    out = [x for x in base if x <= limit * (1.0 + 1e-12)]
    if not out or limit > out[-1] * (1.0 + 1e-6): out.append(limit)
    return out


def run_case(sim: Path, input_path: Path, out_dir: Path, mode: str,
             worst: dict, direction: list[float], scale: float) -> dict:
    tag = f"{scale:+.12g}".replace("+", "p").replace("-", "m")
    trace = out_dir / f"frozen_shadow_{mode}_{tag}.csv"
    env = os.environ.copy()
    env.update({
        "OU3_SHADOW_TRACE": str(trace),
        "OU3_SHADOW_T0": f"{float(worst['t0']):.17g}",
        "OU3_SHADOW_T1": f"{float(worst['t1']):.17g}",
        "OU3_SHADOW_MODE": mode,
        "OU3_SHADOW_DIRECTION": ",".join(f"{x:.17g}" for x in direction),
        "OU3_SHADOW_SCALE": f"{scale:.17g}",
        "W3D_WRITE_TIMESERIES": "0",
        "W3D_VALIDATION_WINDOW_SEC": "0",
    })
    cp = subprocess.run(
        [str(sim), "--input", str(input_path)], env=env, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
    )
    matches = list(DONE_RE.finditer(cp.stdout))
    result = {
        "mode": mode, "scale": scale, "absolute_scale": abs(scale),
        "returncode": cp.returncode, "trace": str(trace),
        "stdout_tail": "\n".join(cp.stdout.splitlines()[-10:]),
    }
    if cp.returncode != 0 or not matches:
        result["valid"] = False
        return result
    m = matches[-1].groupdict()
    result.update({
        "valid": True,
        "t0": float(m["t0"]), "t1": float(m["t1"]),
        "V0": float(m["V0"]), "V1": float(m["V1"]), "rho": float(m["rho"]),
        "reconstruction_max": float(m["recon"]),
        "prediction_count": int(m["pred"]), "S_count": int(m["S"]),
        "acc_count": int(m["acc"]), "vector_count": int(m["vector"]),
    })
    result["endpoint_matches"] = abs(result["t1"] - float(worst["t1"])) <= 2e-6
    result["event_counts_match"] = (
        result["prediction_count"] == 600
        and result["acc_count"] == int(worst["acc_count"])
        and result["S_count"] == int(worst["S_update_count"])
        and result["vector_count"] == int(worst["mag_count"])
    )
    result["valid"] = bool(result["endpoint_matches"] and result["event_counts_match"])
    return result


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
        "qualification": "NON_PROMOTING_COMPLETE_SEA3_P4_FROZEN_GAIN_NONLINEAR_SHADOW",
        "canonical_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "point_same_history_diagnostic_only": True,
        "P4_promoted": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "single_shipping_estimator_owns_covariance_and_gains": True,
        "shadow_covariance_exists": False,
        "shipping_acceptance_branches_frozen": True,
        "same_actual_applied_RS_word_retained": True,
        "all_due_S_updates_retained": True,
        "packet_count_remainder_budget_used": False,
        "selected_S_subset_used": False,
        "both_signs_tested": True,
        "modes": {},
    }

    for mode in ("H18", "A21"):
        worst = linear["modes"][mode]["worst"]
        direction = [float(x) for x in worst["maximizing_direction"]["components"]]
        limit, limit_detail = max_scale(direction, mode, domain, args.input)
        magnitudes = choose_magnitudes(limit)
        cases = []
        for mag in magnitudes:
            cases.append(run_case(args.sim.resolve(), args.input.resolve(), args.output_dir.resolve(),
                                  mode, worst, direction, -mag))
            cases.append(run_case(args.sim.resolve(), args.input.resolve(), args.output_dir.resolve(),
                                  mode, worst, direction, +mag))
        valid = [c for c in cases if c.get("valid")]
        pairs = []
        for mag in magnitudes:
            neg = next((c for c in valid if math.isclose(c["scale"], -mag, rel_tol=0, abs_tol=1e-12)), None)
            pos = next((c for c in valid if math.isclose(c["scale"], +mag, rel_tol=0, abs_tol=1e-12)), None)
            if neg and pos:
                central = 0.5 * (float(neg["rho"]) + float(pos["rho"]))
                pairs.append({
                    "absolute_scale": mag,
                    "rho_negative": float(neg["rho"]),
                    "rho_positive": float(pos["rho"]),
                    "rho_central_average": central,
                    "central_minus_linear": central - float(worst["rho_linear"]),
                })
        smallest_pair = min(pairs, key=lambda p: p["absolute_scale"]) if pairs else None
        crossings = [c for c in valid if float(c["rho"]) >= 1.0]
        report["modes"][mode] = {
            "rho_linear_worst_word": float(worst["rho_linear"]),
            "linear_distance_to_one": float(worst["distance_to_one"]),
            "word_t0": float(worst["t0"]), "word_t1": float(worst["t1"]),
            "acc_count": int(worst["acc_count"]), "S_update_count": int(worst["S_update_count"]),
            "mag_count": int(worst["mag_count"]),
            "RS_scalar_start": float(worst["RS_scalar_start"]),
            "RS_scalar_end": float(worst["RS_scalar_end"]),
            "declared_scale_limit": limit,
            "declared_scale_limit_detail": limit_detail,
            "cases": cases, "valid_cases": len(valid), "paired_scales": pairs,
            "smallest_pair_absolute_scale": smallest_pair["absolute_scale"] if smallest_pair else None,
            "smallest_pair_central_rho": smallest_pair["rho_central_average"] if smallest_pair else None,
            "smallest_pair_central_minus_linear": smallest_pair["central_minus_linear"] if smallest_pair else None,
            "worst_rho": max((float(c["rho"]) for c in valid), default=None),
            "strict_contraction_on_all_tested_frozen_word_cases": bool(valid) and not crossings,
            "rho_ge_one_cases": crossings,
            "first_rho_ge_one_absolute_scale": min((abs(float(c["scale"])) for c in crossings), default=None),
            "max_nominal_reconstruction_error": max((float(c["reconstruction_max"]) for c in valid), default=None),
        }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({m: {
        "linear_rho": report["modes"][m]["rho_linear_worst_word"],
        "smallest_pair_central_rho": report["modes"][m]["smallest_pair_central_rho"],
        "central_minus_linear": report["modes"][m]["smallest_pair_central_minus_linear"],
        "worst_rho": report["modes"][m]["worst_rho"],
        "first_rho_ge_one_absolute_scale": report["modes"][m]["first_rho_ge_one_absolute_scale"],
        "declared_scale_limit": report["modes"][m]["declared_scale_limit"],
        "max_reconstruction_error": report["modes"][m]["max_nominal_reconstruction_error"],
    } for m in ("H18", "A21")}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
