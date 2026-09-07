#!/usr/bin/env python3
"""Non-promoting 6/9 s physical complete-SEA3 P4 feasibility probe.

This is deliberately different from the retired multiword estimator-pair
shadow experiment.  One shipping observer scans ONE source-contiguous window of
the requested length, recomputes the worst linear H18/A21 direction for that
same horizon, replays only the selected nominal shipping source history through
the same observer to emit the physical P/H/R/event payload, and then evaluates
the paper's true-minus-estimated finite Cayley map.

No source word is generated here.  No independent scheduler, tuner, R_S path,
Riccati recursion or second estimator is introduced.  Every valid accelerometer
update and every due S=0 operation are required, and each S event keeps the
actual applied anisotropic SpectralMSE R_S from the same shipping history.

The result is a falsification/feasibility diagnostic only.  Endpoint rho is
reported but is never a CI pass/fail condition and cannot promote P4.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess

import ou3_p4_complete_sea3_word_feasibility as WORD
import ou3_p4_physical_finite_map_feasibility as BASE
import ou3_p4_physical_finite_map_feasibility_fast as FAST
import ou3_p4_canonicalize_physical_payload as CANON

DT = 1.0 / 200.0


def analyze_linear_windows(map_path: Path, cov_path: Path, stride: int) -> dict:
    map_stride, maps = WORD.read_maps(map_path)
    cov_stride, covs = WORD.read_covariances(cov_path)
    if map_stride != cov_stride or map_stride != stride:
        raise RuntimeError(
            f"long-window map/cov stride mismatch map={map_stride} cov={cov_stride} expected={stride}"
        )
    expected_horizon = stride * DT
    n = min(len(maps), len(covs))
    if n == 0:
        raise RuntimeError("long-window scan produced no records")

    rows: list[dict] = []
    rejected: dict[str, int] = {}
    min_s_updates = max(1, int(round(20.0 * stride / 600.0)))
    for i in range(n):
        m, c = maps[i], covs[i]
        if abs(m["t0"] - c["t0"]) > 2.0e-4 or abs(m["t1"] - c["t1"]) > 2.0e-4:
            raise RuntimeError(f"long-window map/cov block alignment lost at record {i}")
        flags = int(m["flags"])
        reasons: list[str] = []
        if not (flags & 1):
            reasons.append("map_invalid")
        if flags & (1 << 7):
            reasons.append("hybrid_jump")
        if not (flags & (1 << 1)) or not (flags & (1 << 2)):
            reasons.append("not_live")
        mode = WORD._mode(flags)
        if mode == "HYBRID":
            reasons.append("mode_change")
        horizon = float(m["t1"] - m["t0"])
        if abs(horizon - expected_horizon) > 0.03:
            reasons.append("wrong_horizon")
        if int(m["acc_count"]) != stride:
            reasons.append("missing_accelerometer_update")
        if int(m["pseudo_count"]) < min_s_updates:
            reasons.append("missing_due_S_regularization")
        if reasons:
            for r in reasons:
                rejected[r] = rejected.get(r, 0) + 1
            continue

        dim = 18 if mode == "H18" else 21
        rho, direction = WORD._linear_ratio(
            m["M"][:dim, :dim], c["P0"][:dim, :dim], c["P1"][:dim, :dim]
        )
        rows.append(
            {
                "record": i,
                "mode": mode,
                "t0": float(m["t0"]),
                "t1": float(m["t1"]),
                "rho_linear": rho,
                "distance_to_one": 1.0 - rho,
                "acc_count": int(m["acc_count"]),
                "mag_count": int(m["mag_count"]),
                "S_update_count": int(m["pseudo_count"]),
                "tau_start": float(m["tau0"]),
                "tau_end": float(m["tau1"]),
                "sigma_start": float(m["sigma0"]),
                "sigma_end": float(m["sigma1"]),
                "RS_scalar_start": float(m["rs0"]),
                "RS_scalar_end": float(m["rs1"]),
                "exact_map_linearization_recovery_residual": float(m["linearization_residual"]),
                "maximizing_direction": WORD._direction_json(direction, dim),
            }
        )

    modes: dict[str, dict] = {}
    for mode in ("H18", "A21"):
        candidates = [x for x in rows if x["mode"] == mode]
        if not candidates:
            raise RuntimeError(f"no legal {mode} window at {expected_horizon:g} s")
        worst = max(candidates, key=lambda x: x["rho_linear"])
        modes[mode] = {
            "legal_blocks": len(candidates),
            "rho_min": min(x["rho_linear"] for x in candidates),
            "rho_max": worst["rho_linear"],
            "worst": worst,
        }

    return {
        "qualification": "NON_PROMOTING_COMPLETE_SEA3_P4_LONG_PHYSICAL_WINDOW_LINEAR_SCAN",
        "canonical_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "same_history_point_diagnostic_only": True,
        "stride_samples": stride,
        "horizon_s": expected_horizon,
        "all_valid_accelerometer_updates_required": True,
        "all_due_S_updates_required": True,
        "actual_applied_RS_used_inside_each_S_gain": True,
        "selected_S_subset_used": False,
        "packet_count_remainder_budget_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "P4_promoted": False,
        "records_seen": n,
        "records_rejected": rejected,
        "modes": modes,
    }


def emit_payload(sim: Path, input_path: Path, mode: str, worst: dict, out_dir: Path) -> Path:
    ledger = out_dir / f"event_ledger_{mode}.csv"
    payload = out_dir / f"payload_{mode}.bin"
    env = os.environ.copy()
    env.pop("OU3_LEDGER_MAP_TRACE", None)
    env.pop("OU3_LEDGER_COV_TRACE", None)
    env.update(
        {
            "OU3_LEDGER_TRACE": str(ledger),
            "OU3_PHYSICAL_PAYLOAD_TRACE": str(payload),
            "OU3_LEDGER_T0": format(float(worst["t0"]), ".17g"),
            "OU3_LEDGER_T1": format(float(worst["t1"]), ".17g"),
            "OU3_LEDGER_MODE": mode,
            "OU3_LEDGER_DIRECTION": worst["maximizing_direction"]["csv"],
            "W3D_WRITE_TIMESERIES": "0",
            "W3D_VALIDATION_WINDOW_SEC": "0",
        }
    )
    cp = subprocess.run(
        [str(sim.resolve()), "--input", str(input_path.resolve())],
        env=env,
        cwd=str(sim.resolve().parent),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    print(cp.stdout)
    if cp.returncode != 0:
        raise RuntimeError(f"{mode} long-window physical payload replay failed rc={cp.returncode}")
    if not payload.exists() or payload.stat().st_size == 0:
        raise RuntimeError(f"{mode} long-window physical payload was not emitted")
    return payload


def physical_report(
    payload_path: Path,
    domain_path: Path,
    input_path: Path,
    expected_stride: int,
    out_dir: Path,
    mode: str,
) -> dict:
    canonical = out_dir / f"payload_canonical_{mode}.bin"
    projection_json = out_dir / f"payload_projection_{mode}.json"
    projection = CANON.canonicalize(payload_path, canonical)
    projection_json.write_text(json.dumps(projection, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    payload = BASE.read_payload(canonical)
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    projection_limit = float(domain["normal_live"]["active_accelerometer_bias_projection_limit_mps2"])
    linear = BASE.build_reset_normalized_linear_path(payload, projection_limit)
    parity = BASE.zero_state_parity(payload, linear, projection_limit)
    if not parity["pass"]:
        raise RuntimeError(f"{mode} long-window zero-state parity failed: {parity}")
    if linear["minimum_measurement_R_eigenvalue"] <= 0.0:
        raise RuntimeError(f"{mode} long-window payload has non-positive measurement R")
    if linear["counts"]["prediction"] != expected_stride:
        raise RuntimeError(
            f"{mode} prediction count {linear['counts']['prediction']} != {expected_stride}"
        )
    if linear["counts"]["accelerometer"] != expected_stride:
        raise RuntimeError(
            f"{mode} accelerometer count {linear['counts']['accelerometer']} != {expected_stride}"
        )
    if linear["counts"]["S_zero"] <= 0:
        raise RuntimeError(f"{mode} long-window payload lost all due S updates")
    if linear["actual_RS_std_ratio_max_error"] > 2.0e-5:
        raise RuntimeError(f"{mode} long-window payload lost actual R_S anisotropy")

    kernels = FAST.prepare_kernels(payload, linear)
    limit, limit_detail = BASE.scale_limit(
        linear["direction"], payload["mode_dim"], domain, input_path
    )
    scales = BASE.choose_scales(limit)
    # Add the known 3-second transition region so 6/9 s results are directly
    # comparable even when the generic scale grid skips 6.5.
    for x in (6.5, 7.0, 7.5):
        if x <= limit * (1.0 + 1.0e-12):
            scales.append(float(x))
    scales = sorted(set(scales))

    cases: list[dict] = []
    for s in scales:
        cases.append(
            FAST.propagate_case(
                payload, linear, kernels, domain, input_path, projection_limit, -s
            )
        )
        cases.append(
            FAST.propagate_case(
                payload, linear, kernels, domain, input_path, projection_limit, +s
            )
        )

    valid = [c for c in cases if c["domain_retained"]]
    crossings = [c for c in valid if c["rho_endpoint"] >= 1.0]
    return {
        "qualification": "NON_PROMOTING_COMPLETE_SEA3_PHYSICAL_LONG_WINDOW_POINT_V1",
        "canonical_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "point_same_history_diagnostic_only": True,
        "physical_true_minus_estimated_map": True,
        "reset_normalized_metric_diagnostic_only": True,
        "same_single_shipping_observer_source_payload": True,
        "second_estimator_or_Riccati_history_used": False,
        "all_due_S_updates_with_actual_RS_retained": True,
        "packet_count_remainder_budget_used": False,
        "state_elimination_used": False,
        "longer_estimator_pair_shadow_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "P4_promoted": False,
        "P5_may_start": False,
        "mode": mode,
        "dimension": payload["mode_dim"],
        "word_t0": payload["t0"],
        "word_t1": payload["t1"],
        "word_horizon_s": payload["t1"] - payload["t0"],
        "event_counts": linear["counts"],
        "reset_normalized_linear_rho": linear["rho_linear"],
        "reset_normalized_linear_distance_to_one": 1.0 - linear["rho_linear"],
        "zero_state_parity": parity,
        "payload_projection": projection,
        "actual_RS_std_ratio_max_error": linear["actual_RS_std_ratio_max_error"],
        "minimum_measurement_R_eigenvalue": linear["minimum_measurement_R_eigenvalue"],
        "direction": linear["direction"].tolist(),
        "declared_scale_limit": limit,
        "declared_scale_limit_detail": limit_detail,
        "tested_absolute_scales": scales,
        "cases": cases,
        "valid_domain_cases": len(valid),
        "strict_endpoint_contraction_on_all_domain_retained_tested_cases": bool(valid)
        and not crossings,
        "first_endpoint_rho_ge_one_absolute_scale": min(
            (c["absolute_scale"] for c in crossings), default=None
        ),
        "worst_endpoint_rho_on_domain_retained_cases": max(
            (c["rho_endpoint"] for c in valid), default=None
        ),
        "worst_prefix_ratio_on_domain_retained_cases": max(
            (c["max_prefix_ratio"] for c in valid), default=None
        ),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--map", type=Path, required=True)
    ap.add_argument("--cov", type=Path, required=True)
    ap.add_argument("--sim", type=Path, required=True)
    ap.add_argument("--input", type=Path, required=True)
    ap.add_argument("--domain", type=Path, required=True)
    ap.add_argument("--stride", type=int, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    if args.stride not in (1200, 1800):
        raise ValueError("long-window diagnostic accepts only 1200 (6 s) or 1800 (9 s)")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    linear_scan = analyze_linear_windows(args.map, args.cov, args.stride)
    (args.output_dir / "linear_scan.json").write_text(
        json.dumps(linear_scan, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    modes: dict[str, dict] = {}
    for mode in ("H18", "A21"):
        worst = linear_scan["modes"][mode]["worst"]
        payload = emit_payload(args.sim, args.input, mode, worst, args.output_dir)
        modes[mode] = physical_report(
            payload, args.domain, args.input, args.stride, args.output_dir, mode
        )
        (args.output_dir / f"physical_{mode}.json").write_text(
            json.dumps(modes[mode], indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )

    result = {
        "qualification": "NON_PROMOTING_COMPLETE_SEA3_P4_PHYSICAL_LONG_WINDOW_FEASIBILITY_V1",
        "canonical_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "point_same_history_diagnostic_only": True,
        "window_stride_samples": args.stride,
        "window_horizon_s": args.stride * DT,
        "paper_permits_longer_finite_observability_window": True,
        "retired_estimator_pair_multiword_shadow_used": False,
        "single_shipping_observer_used": True,
        "actual_applied_RS_retained": True,
        "filter_changed": False,
        "declared_domain_changed": False,
        "P4_promoted": False,
        "P5_may_start": False,
        "linear_scan": linear_scan,
        "modes": modes,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8"
    )

    print(
        json.dumps(
            {
                "horizon_s": result["window_horizon_s"],
                "H18_linear_rho": linear_scan["modes"]["H18"]["worst"]["rho_linear"],
                "A21_linear_rho": linear_scan["modes"]["A21"]["worst"]["rho_linear"],
                "H18_worst_physical_rho": modes["H18"]["worst_endpoint_rho_on_domain_retained_cases"],
                "A21_worst_physical_rho": modes["A21"]["worst_endpoint_rho_on_domain_retained_cases"],
                "A21_first_crossing": modes["A21"]["first_endpoint_rho_ge_one_absolute_scale"],
                "A21_worst_prefix_ratio": modes["A21"]["worst_prefix_ratio_on_domain_retained_cases"],
                "P4_promoted": False,
            },
            indent=2,
            sort_keys=True,
        )
    )
    # Contraction is a result, never an infrastructure pass condition.
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
