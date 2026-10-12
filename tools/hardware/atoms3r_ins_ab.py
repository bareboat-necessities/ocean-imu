#!/usr/bin/env python3
"""Summarize AtomS3R INS hardware A/B serial captures.

Input logs must use the non-NMEA, non-plotter debug stream and contain the
standardized [INS_AB] rows. [MAGTIME] rows are optional; when present their
last cumulative counters are included in the summary.
"""
from __future__ import annotations

import argparse
import csv
import math
import re
import statistics
import sys
from pathlib import Path

KV_RE = re.compile(r"([A-Za-z0-9_|]+)=([^\s]+)")


def kv(line: str) -> dict[str, str]:
    return {m.group(1): m.group(2) for m in KV_RE.finditer(line)}


def f(v: str | None) -> float:
    try:
        return float(v) if v is not None else math.nan
    except ValueError:
        return math.nan


def i(v: str | None) -> int:
    try:
        return int(v) if v is not None else 0
    except ValueError:
        return 0


def angle_error_deg(a: float, b: float) -> float:
    return (a - b + 180.0) % 360.0 - 180.0


def slope(xs: list[float], ys: list[float]) -> float:
    pts = [(x, y) for x, y in zip(xs, ys) if math.isfinite(x) and math.isfinite(y)]
    if len(pts) < 2:
        return math.nan
    mx = sum(x for x, _ in pts) / len(pts)
    my = sum(y for _, y in pts) / len(pts)
    den = sum((x - mx) ** 2 for x, _ in pts)
    if den == 0:
        return math.nan
    return sum((x - mx) * (y - my) for x, y in pts) / den


def summarize(label: str, path: Path, heading_ref: float | None, heave_ref: float | None) -> dict[str, object]:
    rows: list[dict[str, str]] = []
    last_mag: dict[str, str] = {}
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("[MAGTIME]"):
            last_mag = kv(line)
        elif line.startswith("[INS_AB]"):
            rows.append(kv(line))

    if not rows:
        raise SystemExit(f"{path}: no [INS_AB] rows found")

    ts = [f(r.get("t_ms")) * 1e-3 for r in rows]
    hdg = [f(r.get("heading_deg")) for r in rows]
    valid = [i(r.get("heading_valid")) != 0 for r in rows]
    raw = [f(r.get("heave_raw_m")) for r in rows]
    clean = [f(r.get("heave_clean_m")) for r in rows]

    valid_hdg = [x for x, ok in zip(hdg, valid) if ok and math.isfinite(x)]
    heading_steps = [
        abs(angle_error_deg(b, a))
        for a, b in zip(valid_hdg, valid_hdg[1:])
    ]
    flips = sum(x >= 150.0 for x in heading_steps)

    out: dict[str, object] = {
        "label": label,
        "filter": rows[-1].get("filter", ""),
        "samples": len(rows),
        "heading_valid_fraction": sum(valid) / len(valid),
        "heading_max_step_deg": max(heading_steps) if heading_steps else math.nan,
        "heading_flip_count": flips,
        "heave_raw_mean_m": statistics.fmean(x for x in raw if math.isfinite(x)),
        "heave_clean_mean_m": statistics.fmean(x for x in clean if math.isfinite(x)),
        "heave_clean_peak_to_peak_m": max(clean) - min(clean),
        "heave_clean_drift_m_per_min": slope(ts, clean) * 60.0,
    }

    if heading_ref is not None:
        errs = [angle_error_deg(x, heading_ref) for x, ok in zip(hdg, valid) if ok and math.isfinite(x)]
        out["heading_rms_error_deg"] = math.sqrt(statistics.fmean(e * e for e in errs)) if errs else math.nan
        out["heading_peak_error_deg"] = max((abs(e) for e in errs), default=math.nan)

    if heave_ref is not None:
        errs = [x - heave_ref for x in clean if math.isfinite(x)]
        out["heave_clean_rms_error_m"] = math.sqrt(statistics.fmean(e * e for e in errs)) if errs else math.nan
        out["heave_clean_peak_error_m"] = max((abs(e) for e in errs), default=math.nan)

    for name in (
        "source_total", "complete_total", "service_total", "rejected_total",
        "duplicates", "history_miss", "max_source_gap_ms", "max_service_gap_ms",
        "host_age_ms", "read_to_update_ms", "jitter_range_ms",
    ):
        if name in last_mag:
            out[f"mag_{name}"] = f(last_mag[name]) if "_ms" in name else i(last_mag[name])

    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="append", required=True,
                    help="LABEL=serial.log; repeat for v2.4.2/main and filters")
    ap.add_argument("--heading-reference-deg", type=float)
    ap.add_argument("--heave-reference-m", type=float)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()

    summaries = []
    for spec in args.run:
        if "=" not in spec:
            ap.error("--run must be LABEL=path")
        label, filename = spec.split("=", 1)
        summaries.append(summarize(label, Path(filename), args.heading_reference_deg,
                                   args.heave_reference_m))

    keys: list[str] = []
    for row in summaries:
        for k in row:
            if k not in keys:
                keys.append(k)
    fp = args.output.open("w", newline="") if args.output else sys.stdout
    try:
        w = csv.DictWriter(fp, fieldnames=keys)
        w.writeheader()
        w.writerows(summaries)
    finally:
        if args.output:
            fp.close()


if __name__ == "__main__":
    main()
