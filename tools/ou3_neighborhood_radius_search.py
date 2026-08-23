#!/usr/bin/env python3
"""Locate sampled nonlinear-neighborhood boundaries for OU-III H and A modes.

For each requested mode this tool repeatedly runs the full information-normalized
coordinate sweep from ``ou3_neighborhood_diagnostic.py`` at a common information
energy W.  It expands W geometrically until that mode's full eight-sea basis
ceases to pass, then bisects the last-pass/first-fail bracket in log W.

The result is deliberately a *sampled numerical basin estimate*.  It cannot
promote the nonlinear or deployment theorem; promotion still requires validated
outward-rounded bounds for the complete continuous source/neighborhood family.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import subprocess
import sys

import ou3_numerical_certificate as BASE

REPO = Path(__file__).resolve().parents[1]
DRIVER = REPO / "tools" / "ou3_neighborhood_diagnostic.py"
DEFAULT_OUT = BASE.DEFAULT_OUT / "neighborhood_radius_search"
DEFAULT_DIAG = REPO / "reports" / "diagnostics" / "ou3_neighborhood_radius"


def classify_case_failure(case: dict) -> str | None:
    """Classify the first observable sampled safety obstruction for one case."""
    if case.get("pass_sampled"):
        return None
    status = case.get("status")
    if status == "SOURCE_TIME_MISALIGNMENT":
        return "SOURCE_TIME_ALIGNMENT"
    if status not in ("PASS_SAMPLED", "FAIL_SAMPLED"):
        return "TRACE_OR_TOOL_FAILURE"

    # The simulator's historical source_match field also drops when an
    # acceptance decision differs.  Use the dedicated gate fields first so a
    # gate bifurcation is not mislabeled as a source-parameter divergence.
    if case.get("measurement_acceptance_match_all") is False:
        return "MEASUREMENT_GATING_CONSISTENCY"
    if case.get("source_match_all") is False:
        return "SOURCE_WORD_IDENTITY"

    theta = case.get("theta_max_rad")
    try:
        if theta is None or not math.isfinite(float(theta)):
            return "PREFIX_FINITE_SAFETY"
        if float(theta) >= math.pi:
            return "SO3_CHART_SAFETY"
    except (TypeError, ValueError):
        return "PREFIX_FINITE_SAFETY"
    try:
        W0 = float(case.get("W0"))
        W1 = float(case.get("W1"))
        if not (math.isfinite(W0) and math.isfinite(W1) and W0 > 0.0):
            return "PREFIX_FINITE_SAFETY"
        if not W1 < W0:
            return "POSITIVE_W_DECREASE"
    except (TypeError, ValueError):
        return "PREFIX_FINITE_SAFETY"
    return "PREFIX_OR_UNCLASSIFIED_SAFETY"


def summarize_round(report: dict, target_W: float, mode: str | None = None) -> dict:
    cases = list(report.get("cases") or [])
    failures = []
    for case in cases:
        reason = classify_case_failure(case)
        if reason is not None:
            failures.append({
                "case": case.get("case"),
                "reason": reason,
                "mode": case.get("mode"),
                "direction": case.get("direction"),
                "sign": case.get("sign"),
                "family": case.get("family"),
                "Hs_m": case.get("Hs_m"),
                "relative_decrement": case.get("relative_decrement"),
                "theta_max_rad": case.get("theta_max_rad"),
                "prefix_W_gain_max": case.get("prefix_W_gain_max"),
                "W0": case.get("W0"),
                "W1": case.get("W1"),
            })
    passed = report.get("status") == "PASS_SAMPLED" and not failures
    reason_counts: dict[str, int] = {}
    for failure in failures:
        reason = failure["reason"]
        reason_counts[reason] = reason_counts.get(reason, 0) + 1
    return {
        "mode": mode,
        "target_W": float(target_W),
        "pass_all_sampled": bool(passed),
        "case_count": report.get("case_count"),
        "valid_endpoint_case_count": report.get("valid_endpoint_case_count"),
        "relative_decrement_min": report.get("relative_decrement_min"),
        "theta_max_rad": report.get("theta_max_rad"),
        "prefix_W_gain_max": report.get("prefix_W_gain_max"),
        "failure_reason_counts": reason_counts,
        "first_failures": failures[:32],
    }


def run_round(args, mode: str, target_W: float, round_index: int) -> dict:
    round_root = (
        args.diagnostic_dir.resolve() / mode /
        f"round_{round_index:02d}_W{target_W:.12g}"
    )
    result_dir = round_root / "result"
    trace_dir = round_root / "traces"
    result_dir.mkdir(parents=True, exist_ok=True)
    trace_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable, str(DRIVER),
        "--no-build",
        "--certificate-dir", str(args.certificate_dir.resolve()),
        "--data-dir", str(args.data_dir.resolve()),
        "--sim", str(args.sim.resolve()),
        "--output-dir", str(result_dir),
        "--diagnostic-dir", str(trace_dir),
        "--modes", mode,
        "--full-basis",
        "--target-W", f"{target_W:.17g}",
        "--jobs", str(args.jobs),
    ]
    p = subprocess.run(
        cmd, cwd=REPO, text=True, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT
    )
    (round_root / "driver.log").write_text(p.stdout)
    result_path = result_dir / "neighborhood_diagnostic.json"
    if p.returncode != 0 or not result_path.exists():
        return {
            "mode": mode,
            "target_W": float(target_W),
            "pass_all_sampled": False,
            "case_count": 0,
            "valid_endpoint_case_count": 0,
            "failure_reason_counts": {"TRACE_OR_TOOL_FAILURE": 1},
            "first_failures": [{
                "case": None, "reason": "TRACE_OR_TOOL_FAILURE", "mode": mode
            }],
            "driver_returncode": int(p.returncode),
        }
    report = json.loads(result_path.read_text())
    ans = summarize_round(report, target_W, mode)
    ans["driver_returncode"] = int(p.returncode)
    return ans


def search_mode(args, mode: str) -> dict:
    rounds = []
    round_index = 0
    w = args.start_W
    last_pass = None
    first_fail = None

    while w <= args.max_W * (1.0 + 1e-12):
        r = run_round(args, mode, w, round_index)
        rounds.append(r)
        round_index += 1
        if r["pass_all_sampled"]:
            last_pass = w
            w *= args.growth
        else:
            first_fail = w
            break

    if first_fail is not None and last_pass is not None:
        lo = last_pass
        hi = first_fail
        for _ in range(args.bisection_steps):
            mid = math.sqrt(lo * hi)
            r = run_round(args, mode, mid, round_index)
            rounds.append(r)
            round_index += 1
            if r["pass_all_sampled"]:
                lo = mid
            else:
                hi = mid
        last_pass = lo
        first_fail = hi

    if first_fail is not None and last_pass is not None:
        status = "BRACKETED_SAMPLED_BOUNDARY"
    elif first_fail is None and last_pass is not None:
        status = "NO_SAMPLED_FAILURE_THROUGH_MAX_W"
    else:
        status = "FAILED_AT_INITIAL_RADIUS"

    first_failure_summary = None
    if first_fail is not None:
        failing = min(
            (r for r in rounds if not r["pass_all_sampled"]),
            key=lambda r: abs(math.log(r["target_W"] / first_fail)),
        )
        first_failure_summary = {
            "target_W": failing["target_W"],
            "failure_reason_counts": failing.get("failure_reason_counts", {}),
            "limiting_cases": failing.get("first_failures", []),
        }

    return {
        "mode": mode,
        "status": status,
        "last_all_pass_W": last_pass,
        "first_fail_W": first_fail,
        "information_radius_last_pass_sqrt_W": (
            math.sqrt(last_pass) if last_pass is not None else None
        ),
        "information_radius_first_fail_sqrt_W": (
            math.sqrt(first_fail) if first_fail is not None else None
        ),
        "first_failure": first_failure_summary,
        "rounds": rounds,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--certificate-dir", type=Path, default=BASE.DEFAULT_OUT)
    ap.add_argument("--data-dir", type=Path, default=BASE.DEFAULT_DATA_DIR)
    ap.add_argument("--sim", type=Path, default=BASE.TEST_DIR / "ou3-neighborhood-sim")
    ap.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--diagnostic-dir", type=Path, default=DEFAULT_DIAG)
    ap.add_argument("--modes", default="H,A")
    ap.add_argument("--start-W", type=float, default=0.05)
    ap.add_argument("--growth", type=float, default=16.0)
    ap.add_argument("--max-W", type=float, default=838860.8)
    ap.add_argument("--bisection-steps", type=int, default=5)
    ap.add_argument("--jobs", type=int, default=4)
    args = ap.parse_args()

    if not (args.start_W > 0.0 and math.isfinite(args.start_W)):
        raise ValueError("--start-W must be finite positive")
    if not (args.growth > 1.0 and math.isfinite(args.growth)):
        raise ValueError("--growth must be finite and > 1")
    if not (args.max_W >= args.start_W and math.isfinite(args.max_W)):
        raise ValueError("--max-W must be finite and >= --start-W")
    if args.bisection_steps < 0 or args.jobs < 1:
        raise ValueError("--bisection-steps must be >= 0 and --jobs >= 1")
    if not args.sim.resolve().exists():
        raise FileNotFoundError(args.sim.resolve())
    modes = [x.strip() for x in args.modes.split(",") if x.strip()]
    if not modes or any(mode not in ("H", "A") for mode in modes):
        raise ValueError("--modes accepts H,A")
    modes = list(dict.fromkeys(modes))

    out = args.output_dir.resolve()
    diag = args.diagnostic_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    diag.mkdir(parents=True, exist_ok=True)

    mode_reports = {mode: search_mode(args, mode) for mode in modes}
    statuses = [mode_reports[m]["status"] for m in modes]
    if all(s == "BRACKETED_SAMPLED_BOUNDARY" for s in statuses):
        status = "BRACKETED_H_A_SAMPLED_BOUNDARIES" if set(modes) == {"H", "A"} else statuses[0]
    elif all(s == "NO_SAMPLED_FAILURE_THROUGH_MAX_W" for s in statuses):
        status = "NO_SAMPLED_FAILURE_THROUGH_MAX_W"
    else:
        status = "INCOMPLETE_OR_MIXED_SAMPLED_BOUNDARIES"

    report = {
        "schema": 2,
        "status": status,
        "qualification": "SAMPLED_NUMERICAL_BASIN_ESTIMATE_ONLY_NOT_A_VALIDATED_CERTIFICATE",
        "metric_radius_definition": "W=zeta^T Sigma_nominal^-1 zeta",
        "all_coordinates_both_signs": True,
        "all_eight_reference_seas": True,
        "search": {
            "start_W": float(args.start_W),
            "growth": float(args.growth),
            "max_W": float(args.max_W),
            "bisection_steps": int(args.bisection_steps),
        },
        "modes": mode_reports,
        "numerical_neighborhood_certificate": "NOT_ESTABLISHED",
        "deployment_theorem_certificate": "NOT_ESTABLISHED",
    }
    (out / "neighborhood_radius_search.json").write_text(
        json.dumps(report, indent=2, sort_keys=True)
    )
    print(json.dumps({
        "status": report["status"],
        "modes": {
            mode: {
                "status": mode_reports[mode]["status"],
                "last_all_pass_W": mode_reports[mode]["last_all_pass_W"],
                "first_fail_W": mode_reports[mode]["first_fail_W"],
                "first_failure": mode_reports[mode]["first_failure"],
            }
            for mode in modes
        },
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
