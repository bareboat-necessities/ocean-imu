#!/usr/bin/env python3
"""Re-derive OU-III's seven deterministic quality gates from the shipped filter.

The rule, from docs/quality-gate-regauge.md:

    the worst value the filter currently produces across the eight scored
    records, plus about half a percent, rounded up in the last digit the
    channel is quoted in.

The protocol is the deterministic one the gates are written against: default
seeds, the trailing 900 s of each 1200 s replay, the four JONSWAP and four
PM-Stokes records, run with W3D_COLLECT_ALL_GATES=1 so a breach scores the
remaining records instead of exiting at the first one.

Run it after any change to the filter that could move a gated quantity, and
paste the printed table into the FAIL_LIMITS comment block in
tests/kalman_ou_iii/kalman_ou_iii-sim.cpp.
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import math
import os
import re
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SIM = REPO_ROOT / "tests" / "kalman_ou_iii" / "kalman_ou_iii-sim"

RECORDS = (
    "wave_data_jonswap_H0.270_L14.047_A30.00_P60.00.csv",
    "wave_data_jonswap_H1.500_L50.710_A-30.00_P120.00.csv",
    "wave_data_jonswap_H4.000_L112.766_A30.00_P30.00.csv",
    "wave_data_jonswap_H8.500_L202.839_A-30.00_P72.00.csv",
    "wave_data_pmstokes_H0.270_L14.047_A30.00_P60.00.csv",
    "wave_data_pmstokes_H1.500_L50.710_A-30.00_P120.00.csv",
    "wave_data_pmstokes_H4.000_L112.766_A30.00_P30.00.csv",
    "wave_data_pmstokes_H8.500_L202.839_A-30.00_P72.00.csv",
)

# gate field -> (heading, which records it is taken over, current shipped value)
GATES = {
    "err_limit_percent_z_jonswap":   ("Z %Hs JONSWAP",     "jonswap",  4.72),
    "err_limit_percent_z_pmstokes":  ("Z %Hs PM-Stokes",   "pmstokes", 4.69),
    "err_limit_yaw_deg":             ("yaw deg",           "all",      1.068),
    "err_limit_percent_3d_jonswap":  ("3D % JONSWAP",      "jonswap",  21.05),
    "err_limit_percent_3d_pmstokes": ("3D % PM-Stokes",    "pmstokes", 20.83),
    "acc_z_bias_percent":            ("acc Z bias %",      "all",      4.93),
    "bias_3d_percent":               ("bias 3D %",         "all",      98.4),
}


def scrape(record: str, env_extra: dict[str, str]) -> dict[str, float]:
    """One deterministic replay; returns the gated quantities it produced."""
    env = os.environ.copy()
    env.update(env_extra)
    env.update({
        "W3D_COLLECT_ALL_GATES": "1",
        "W3D_WRITE_TIMESERIES": "0",
        "W3D_VALIDATION_WINDOW_SEC": "900",
    })
    proc = subprocess.run([str(SIM), "--input", record], cwd=SIM.parent, env=env,
                          capture_output=True, text=True, check=False)
    out = proc.stdout
    metrics = next((l for l in out.splitlines()
                    if l.startswith("VALIDATION_METRICS")), "")

    def metric(name: str) -> float:
        hit = re.search(rf"\b{name}=([-\d.eE+]+)", metrics)
        return float(hit.group(1)) if hit else math.nan

    # The two bias gates are quoted as a percentage of the largest true bias in
    # the window, which only the summary block carries.
    def bias_pct(kind: str, axis: str) -> float:
        hit = re.search(
            rf"Bias error RMS \(% of max TRUE bias\) \({kind}\):.*?"
            rf"{axis}=([-\d.eE+]+)%", out)
        return float(hit.group(1)) if hit else math.nan

    return {
        "record": record,
        "z_pct_hs": metric("disp_z_pct_hs"),
        "d3_pct": metric("disp_3d_pct_refmax"),
        "yaw_deg": metric("yaw_rms_deg"),
        "acc_z_bias_pct": bias_pct("acc", "Z"),
        "acc_3d_bias_pct": bias_pct("acc", r"\|3D\|"),
        "gyro_3d_bias_pct": bias_pct("gyro", r"\|3D\|"),
    }


def round_to_rule(worst: float) -> tuple[float, float]:
    """Half a percent above `worst`, rounded up to this channel's quantum.

    The quantum is chosen by magnitude rather than fixed at a tenth, which is
    the change docs/quality-gate-regauge.md records: a tenth is 0.03% of a
    400-valued gate and 3.5% of a 1-valued one.  A thousandth where the value
    is near 1, a hundredth for everything larger.  Reproduces all seven shipped
    gates from the shipped filter's worst values.
    """
    quantum = 0.001 if worst < 2.0 else 0.01
    limit = math.ceil(worst * 1.005 / quantum) * quantum
    return round(limit, 3), 100.0 * (limit - worst) / worst


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--env", action="append", default=[],
                    help="extra VAR=VALUE for the replays, e.g. OU_III_S_FACTOR=1.0")
    ap.add_argument("--jobs", type=int, default=4)
    args = ap.parse_args()

    if not SIM.is_file():
        raise SystemExit(f"{SIM} not built; run make -C tests/kalman_ou_iii build")
    env_extra = dict(kv.split("=", 1) for kv in args.env)
    if env_extra:
        print("overrides: " + " ".join(f"{k}={v}" for k, v in env_extra.items()))

    with cf.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        rows = list(pool.map(lambda r: scrape(r, env_extra), RECORDS))
    for row in rows:
        print(f"  {row['record'].split('_')[2]:9} "
              f"{row['record'].split('_')[3]:8} "
              f"z={row['z_pct_hs']:7.4f}%  3D={row['d3_pct']:8.4f}%  "
              f"yaw={row['yaw_deg']:7.4f}  accZbias={row['acc_z_bias_pct']:8.4f}%  "
              f"acc3Dbias={row['acc_3d_bias_pct']:8.4f}%  "
              f"gyr3Dbias={row['gyro_3d_bias_pct']:8.4f}%")

    def worst_of(key: str, family: str) -> tuple[float, str]:
        picked = [r for r in rows
                  if family == "all" or f"_{family}_" in r["record"]]
        best = max(picked, key=lambda r: r[key])
        return best[key], f"{best['record'].split('_')[2]} {best['record'].split('_')[3]}"

    print(f"\n{'gate':32} {'worst':>10}  {'record':22} {'shipped':>9} "
          f"{'rule':>9} {'margin':>8}")
    print("-" * 96)
    verdict_lines = []
    for field, (heading, family, shipped) in GATES.items():
        if field == "bias_3d_percent":
            # One limit covers both the accelerometer and the gyro 3D bias.
            worst_a, rec_a = worst_of("acc_3d_bias_pct", "all")
            worst_g, rec_g = worst_of("gyro_3d_bias_pct", "all")
            worst, rec = ((worst_a, rec_a + ", accel") if worst_a >= worst_g
                          else (worst_g, rec_g + ", gyro"))
        else:
            key = {"err_limit_percent_z_jonswap": "z_pct_hs",
                   "err_limit_percent_z_pmstokes": "z_pct_hs",
                   "err_limit_yaw_deg": "yaw_deg",
                   "err_limit_percent_3d_jonswap": "d3_pct",
                   "err_limit_percent_3d_pmstokes": "d3_pct",
                   "acc_z_bias_percent": "acc_z_bias_pct"}[field]
            worst, rec = worst_of(key, family)
        limit, margin = round_to_rule(worst)
        shipped_margin = 100.0 * (shipped - worst) / worst
        if shipped <= worst:
            flag = "  <-- SHIPPED GATE NOW FAILS"
        elif 0.45 <= shipped_margin <= 0.70:
            # Already at the rule; moving it would be churn, not a re-gauge.
            limit, margin, flag = shipped, shipped_margin, "  (shipped already at the rule)"
        else:
            flag = ""
        print(f"{heading:32} {worst:10.4f}  {rec:22} {shipped:9} "
              f"{limit:9} {margin:7.2f}%{flag}")
        verdict_lines.append(f"    .{field:30} = {limit}f,")
    print("\nFAIL_LIMITS body for the measured filter:\n")
    print("\n".join(verdict_lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
