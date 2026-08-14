#!/usr/bin/env python3
"""Re-derive the OU families' deterministic quality gates from the shipped filter.

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
tests/kalman_ou_ii/kalman_ou_ii-sim.cpp,
tests/kalman_ou_iii/kalman_ou_iii-sim.cpp, or the kRegressionBars block in
tests/kalman_tfg/kalman_tfg-sim.cpp.

`--json` writes every channel of every record, which is what the -march drift
measurement in docs/quality-gate-regauge.md compares between two builds.
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import math
import os
import re
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

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

# gate field -> (heading, which records it is taken over, scraped channel)
GATE_CHANNELS = {
    "err_limit_percent_z_jonswap":   ("Z %Hs JONSWAP",   "jonswap",  "z_pct_hs"),
    "err_limit_percent_z_pmstokes":  ("Z %Hs PM-Stokes", "pmstokes", "z_pct_hs"),
    "err_limit_yaw_deg":             ("yaw deg",         "all",      "yaw_deg"),
    "err_limit_roll_deg":            ("roll deg",        "all",      "roll_deg"),
    "err_limit_pitch_deg":           ("pitch deg",       "all",      "pitch_deg"),
    "err_limit_percent_3d_jonswap":  ("3D % JONSWAP",    "jonswap",  "d3_pct"),
    "err_limit_percent_3d_pmstokes": ("3D % PM-Stokes",  "pmstokes", "d3_pct"),
    "acc_z_bias_percent":            ("acc Z bias %",    "all",      "acc_z_bias_pct"),
    # One limit covers both the accelerometer and the gyro 3D bias; handled below.
    "bias_3d_percent":               ("bias 3D %",       "all",      None),
}

# The gates each family currently ships, in the order above.
FAMILIES = {
    "ou_ii": {
        "sim": REPO_ROOT / "tests" / "kalman_ou_ii" / "kalman_ou_ii-sim",
        "shipped": {
            "err_limit_percent_z_jonswap":   6.865,
            "err_limit_percent_z_pmstokes":  6.848,
            "err_limit_yaw_deg":             1.089,
            "err_limit_roll_deg":            0.4792,
            "err_limit_pitch_deg":           0.3657,
            "err_limit_percent_3d_jonswap":  20.92,
            "err_limit_percent_3d_pmstokes": 21.03,
            "acc_z_bias_percent":            5.324,
            "bias_3d_percent":               94.47,
        },
    },
    "ou_iii": {
        "sim": REPO_ROOT / "tests" / "kalman_ou_iii" / "kalman_ou_iii-sim",
        "shipped": {
            "err_limit_percent_z_jonswap":   4.735,
            "err_limit_percent_z_pmstokes":  4.682,
            "err_limit_yaw_deg":             1.297,
            "err_limit_roll_deg":            0.42,
            "err_limit_pitch_deg":           0.2211,
            "err_limit_percent_3d_jonswap":  20.95,
            "err_limit_percent_3d_pmstokes": 20.86,
            "acc_z_bias_percent":            4.624,
            "bias_3d_percent":               81.84,
        },
    },
    # TFG is gated by the same rule on the same protocol and had been
    # re-derived by hand every time, which is the arithmetic this script exists
    # to remove.  It carries seven bars rather than nine: roll and pitch were
    # gated for the two OU families and not for this one, so those two fields
    # are absent here and the loop below skips whatever a family does not
    # carry.
    "tfg": {
        "sim": REPO_ROOT / "tests" / "kalman_tfg" / "kalman_tfg-sim",
        "shipped": {
            "err_limit_percent_z_jonswap":   4.803,
            "err_limit_percent_z_pmstokes":  4.709,
            "err_limit_yaw_deg":             1.536,
            "err_limit_percent_3d_jonswap":  21.14,
            "err_limit_percent_3d_pmstokes": 20.71,
            "acc_z_bias_percent":            5.026,
            "bias_3d_percent":               167.6,
        },
    },
}


def scrape(sim: Path, record: str, env_extra: dict[str, str]) -> dict[str, float]:
    """One deterministic replay; returns the gated quantities it produced."""
    env = os.environ.copy()
    env.update(env_extra)
    env.update({
        "W3D_COLLECT_ALL_GATES": "1",
        "W3D_WRITE_TIMESERIES": "0",
        "W3D_VALIDATION_WINDOW_SEC": "900",
    })
    proc = subprocess.run([str(sim), "--input", record], cwd=sim.parent, env=env,
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
        "roll_deg": metric("roll_rms_deg"),
        "pitch_deg": metric("pitch_rms_deg"),
        "acc_z_bias_pct": bias_pct("acc", "Z"),
        "acc_3d_bias_pct": bias_pct("acc", r"\|3D\|"),
        "gyro_3d_bias_pct": bias_pct("gyro", r"\|3D\|"),
    }


def round_to_rule(worst: float) -> tuple[float, float]:
    """Half a percent above `worst`, rounded up to this channel's quantum.

    The quantum is a fixed *relative* precision -- four significant digits, so
    one part in a thousand of the value -- rather than a fixed absolute one.
    That is what lets a single half-percent rule land the same margin on a
    channel quoted near 0.2 and one quoted near 100: an absolute quantum is
    0.05% of a 94-valued gate and 45% of a 0.22-valued one, so it decides the
    margin instead of the rule doing it.  With this quantum the rounding costs
    at most a tenth of a percent on top of the half, which is the whole point
    of quoting the rounding rather than the rule.
    """
    quantum = 10.0 ** (math.floor(math.log10(worst)) - 3)
    limit = math.ceil(worst * 1.005 / quantum) * quantum
    # Renormalize: the ceil is exact in decimal but the product is not.
    limit = round(limit, max(0, 3 - math.floor(math.log10(worst))))
    return limit, 100.0 * (limit - worst) / worst


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--family", choices=sorted(FAMILIES), default="ou_iii",
                    help="which simulator to re-gauge (default: ou_iii)")
    ap.add_argument("--env", action="append", default=[],
                    help="extra VAR=VALUE for the replays, e.g. OU_III_S_FACTOR=1.0")
    ap.add_argument("--json", type=Path, default=None,
                    help="also write every channel of every record here")
    ap.add_argument("--jobs", type=int, default=4)
    args = ap.parse_args()

    family = FAMILIES[args.family]
    sim = family["sim"]
    if not sim.is_file():
        raise SystemExit(f"{sim} not built; run make -C {sim.parent} build")
    env_extra = dict(kv.split("=", 1) for kv in args.env)
    print(f"family: {args.family}")
    if env_extra:
        print("overrides: " + " ".join(f"{k}={v}" for k, v in env_extra.items()))

    with cf.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        rows = list(pool.map(lambda r: scrape(sim, r, env_extra), RECORDS))
    for row in rows:
        print(f"  {row['record'].split('_')[2]:9} "
              f"{row['record'].split('_')[3]:8} "
              f"z={row['z_pct_hs']:7.4f}%  3D={row['d3_pct']:8.4f}%  "
              f"roll={row['roll_deg']:6.4f}  pitch={row['pitch_deg']:6.4f}  "
              f"yaw={row['yaw_deg']:7.4f}  accZbias={row['acc_z_bias_pct']:8.4f}%  "
              f"acc3Dbias={row['acc_3d_bias_pct']:8.4f}%  "
              f"gyr3Dbias={row['gyro_3d_bias_pct']:8.4f}%")
    if args.json:
        args.json.write_text(json.dumps(rows, indent=1) + "\n")

    def worst_of(key: str, group: str) -> tuple[float, str]:
        picked = [r for r in rows
                  if group == "all" or f"_{group}_" in r["record"]]
        best = max(picked, key=lambda r: r[key])
        return best[key], f"{best['record'].split('_')[2]} {best['record'].split('_')[3]}"

    print(f"\n{'gate':32} {'worst':>10}  {'record':22} {'shipped':>9} "
          f"{'rule':>9} {'margin':>8}")
    print("-" * 96)
    verdict_lines = []
    for field, (heading, group, channel) in GATE_CHANNELS.items():
        if field not in family["shipped"]:
            continue  # a gate this family does not carry; see FAMILIES
        shipped = family["shipped"][field]
        if channel is None:
            # One limit covers both the accelerometer and the gyro 3D bias.
            worst_a, rec_a = worst_of("acc_3d_bias_pct", "all")
            worst_g, rec_g = worst_of("gyro_3d_bias_pct", "all")
            worst, rec = ((worst_a, rec_a + ", accel") if worst_a >= worst_g
                          else (worst_g, rec_g + ", gyro"))
        else:
            worst, rec = worst_of(channel, group)
        limit, margin = round_to_rule(worst)
        shipped_margin = 100.0 * (shipped - worst) / worst
        if shipped <= worst:
            flag = "  <-- SHIPPED GATE NOW FAILS"
        elif abs(shipped - limit) < 0.5 * 10.0 ** (math.floor(math.log10(worst)) - 3):
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
