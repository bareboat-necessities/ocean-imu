#!/usr/bin/env python3
"""Generate the eight OU-III exact-map/covariance replays without solving a metric.

The old numerical driver coupled expensive replay generation to the coarse
fixed-node SDP.  The latter is now a diagnostic because it merges distinct
Riccati source points into one node.  This driver performs only the pieces that
are common to every valid certificate route:

* run the unchanged adaptive OU-III filter on the exact eight reference seas;
* preserve the normal RMS/quality regression gate;
* record exact closed-loop map blocks and (with the information observer) full
  covariance endpoints;
* reconstruct truth errors for attitude-domain diagnostics;
* enforce the exact-map reconstruction-integrity gate.

No transition identification and no Lyapunov metric are computed here. Raw
simulator stdout is diagnostic rather than scientific result data and is kept
outside ``reports/results`` so result-tree fingerprints cannot depend on runner
workspace paths or other execution-environment text.
"""
from __future__ import annotations

import argparse
import json
import math
import subprocess
from pathlib import Path

import numpy as np

import ou3_numerical_certificate as BASE

DEFAULT_DIAGNOSTIC_DIR = (
    BASE.REPO / "reports" / "diagnostics" / "ou3_numerical_certificate" / "logs"
)


def markdown(report: dict) -> str:
    out = [
        "# OU-III exact replay inputs", "",
        f"Filter regression: **{report['filter_regression']}**",
        f"Exact-map integrity: **{'PASS' if report['map_integrity']['pass'] else 'FAIL'}**",
        f"Maximum map reconstruction residual: {report['map_integrity']['max_linearization_residual']:.3e}",
        "",
        "The coarse fixed-node SDP is deliberately not run in this stage. Linear certification is performed by the source-varying information-metric stage.",
        "",
        "| Sea | RMS | max theta | map blocks | invalid | hybrid |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for r in report["records"]:
        out.append(
            f"| {r['family']} {r['Hs_m']:.2f} | {'PASS' if r['rms_regression_pass'] else 'FAIL'} | "
            f"{r['theta_max_deg']:.2f} deg | {r['map_blocks']} | {r['invalid_map_blocks']} | {r['hybrid_map_blocks']} |"
        )
    return "\n".join(out)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", type=Path, default=BASE.DEFAULT_DATA_DIR)
    ap.add_argument("--output-dir", type=Path, default=BASE.DEFAULT_OUT)
    ap.add_argument("--diagnostic-log-dir", type=Path, default=DEFAULT_DIAGNOSTIC_DIR)
    ap.add_argument("--sim", type=Path,
                    default=BASE.TEST_DIR / "ou3-information-certificate-sim")
    ap.add_argument("--no-build", action="store_true")
    args = ap.parse_args()

    data_dir = args.data_dir.resolve()
    out = args.output_dir.resolve()
    log_dir = args.diagnostic_log_dir.resolve()
    exe = args.sim.resolve()
    out.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    if not args.no_build:
        subprocess.run(["make", "-C", str(BASE.TEST_DIR), exe.name], check=True)

    records = []
    meta = {}
    for family, hs, name in BASE.RECORDS:
        data = (data_dir / name).resolve()
        if not data.exists():
            raise FileNotFoundError(data)
        trace_path, map_path, timeseries, metrics, ok, log = BASE.run_record(exe, data, out)
        slug = f"{family.lower().replace('-', '_')}_{hs:.2f}".replace(".", "_")
        (log_dir / f"{slug}.log").write_text(log)
        trace = np.genfromtxt(trace_path, delimiter=",", names=True, dtype=None, encoding=None)
        E, theta = BASE.build_error_states(trace, timeseries)
        blocks, m = BASE.load_exact_maps(map_path, slug)
        meta[slug] = m
        residual = max((b.linearization_residual for b in blocks), default=math.nan)
        records.append({
            "family": family,
            "Hs_m": hs,
            "slug": slug,
            "rms_regression_pass": bool(ok),
            "rms_metrics": metrics,
            "theta_max_deg": math.degrees(float(np.max(theta))),
            "group_energy_max": float(max(BASE.group_energy(x[:3]) for x in E)),
            "all_attitude_inside_pi": bool(np.max(theta) < math.pi),
            "map_blocks": len(blocks),
            "invalid_map_blocks": sum(not b.valid for b in blocks),
            "hybrid_map_blocks": sum(b.hybrid_jump for b in blocks),
            "max_linearization_residual": residual,
            "handoff_hybrid": BASE.handoff_hybrid(trace, E),
        })

    rms = all(r["rms_regression_pass"] for r in records)
    maxres = max(r["max_linearization_residual"] for r in records)
    integrity = math.isfinite(maxres) and maxres < 5e-3
    report = {
        "schema": 5,
        "scope": "eight_noisy_reference_replays_exact_filter_maps_and_covariance_inputs",
        "record_count": len(records),
        "filter_regression": "PASS" if rms else "FAIL",
        "map_integrity": {
            "max_linearization_residual": maxres,
            "pass": integrity,
            "per_record": meta,
        },
        "coarse_fixed_node_metric": "SKIPPED_DIAGNOSTIC_ONLY",
        "primary_linear_certificate": "SOURCE_VARYING_INFORMATION_METRIC_PENDING",
        "deployment_theorem_certificate": "NOT_ESTABLISHED",
        "records": records,
    }
    (out / "certificate.json").write_text(json.dumps(report, indent=2, sort_keys=True))
    text = markdown(report)
    (out / "certificate.md").write_text(text)
    print(text)
    return 1 if (not rms or not integrity) else 0


if __name__ == "__main__":
    raise SystemExit(main())
