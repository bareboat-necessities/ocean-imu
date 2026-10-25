#!/usr/bin/env python3
"""Monte Carlo benchmark of the experimental CARMA heave filters against OU-III.

Subcommands
-----------
run       Execute the paired synthetic grid with tests/kalman_carma/carma_bench
          (resumable: rows already in raw_runs.csv are skipped).
dataset   Replay the pinned vessel-RAO dataset with the unchanged
          kalman_ou_iii-sim and with carma-sim, over several W3D_SEED values.
matched   Level 1 oracle-shape diagnostic: CARMA(5,2) refitted to each case's
          true vessel-heave spectrum (needs mismatch_analysis.json).
timing    Serial full-INS runs on an otherwise idle host for the runtime table
          (the grid's per-update times are taken under parallel load).
traces    Write the example time series used by the figures.
truth     Realized motion levels of every scenario (truth_motion.csv).
report    Aggregate, compute paired statistics, write summary tables, the
          Markdown report and the SVG figures.

Every synthetic run feeds one truth realization and one corrupted sensor
stream to all estimators, so OU-III and CARMA results for a (scenario, seed)
are paired.  Statistics are paired over seeds: mean difference with a 95 %
bootstrap interval (10 000 resamples, fixed seed), relative change of the
mean, the fraction of seeds on which the candidate is better, and the
Wilcoxon signed-rank p-value.

Usage:
  python3 tools/carma_benchmark.py run [--jobs N] [--quick]
  python3 tools/carma_benchmark.py dataset [--seeds 10]
  python3 tools/carma_benchmark.py matched [--seeds 30]
  python3 tools/carma_benchmark.py timing
  python3 tools/carma_benchmark.py traces
  python3 tools/carma_benchmark.py truth
  python3 tools/carma_benchmark.py report
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BENCH = ROOT / "tests" / "kalman_carma" / "carma_bench"
CARMA_SIM = ROOT / "tests" / "kalman_carma" / "carma-sim"
OU3_SIM = ROOT / "tests" / "kalman_ou_iii" / "kalman_ou_iii-sim"
OUT = ROOT / "reports" / "results" / "carma_benchmark"
RAW = OUT / "raw_runs.csv"
DATASET_RAW = OUT / "dataset_runs.csv"
MATCHED_RAW = OUT / "matched_runs.csv"
TIMING_RAW = OUT / "timing_runs.csv"
TRUTH_STATS = OUT / "truth_motion.csv"
TRACE_DIR = OUT / "traces"

PRIMARY = ["jonswap_fp0.29_hs0.8_g3.3", "pm_fp0.25_hs1", "mixed_swell_wind", "changing_fp"]
SECONDARY = [
    "jonswap_fp0.29_hs0.8_g3.3_norao",
    "jonswap_fp0.20_hs1_g3.3", "jonswap_fp0.20_hs2_g3.3",
    "jonswap_fp0.45_hs0.2_g3.3", "jonswap_fp0.45_hs0.4_g3.3",
    "jonswap_fp0.29_hs0.4_g3.3", "jonswap_fp0.29_hs1.1_g3.3",
    "jonswap_fp0.29_hs0.8_g1.5", "jonswap_fp0.29_hs0.8_g7",
    "jonswap_fp0.10_hs4_g3.3",
    "pm_fp0.15_hs2.5", "pm_fp0.40_hs0.4",
    "broadband", "bimodal", "abrupt_change", "small_amplitude",
    "bias_drift",
]
QUIET = ["stationary", "stationary_calbias", "small_amplitude"]
FULL_ONLY = ["tilted", "large_rotation", "jonswap_fp0.29_hs0.8_g3.3_jit"]
LONG = ["jonswap_fp0.29_hs0.8_g3.3_T3600", "mixed_swell_wind_T3600"]
ORACLE = ["jonswap_fp0.29_hs0.8_g3.3", "pm_fp0.25_hs1", "mixed_swell_wind",
          "jonswap_fp0.20_hs1_g3.3", "jonswap_fp0.45_hs0.4_g3.3", "jonswap_fp0.29_hs0.8_g3.3_norao"]
MATCHED_CASES = ["jonswap_fp0.29_hs0.8_g3.3", "pm_fp0.25_hs1", "mixed_swell_wind", "bimodal"]
MATCHED_FILTERS = "OU3v,C41v,C52v,C52Mv"

FULL_FILTERS = "OU3,C30,C41,C52,C52S,C52F"
VERT_FILTERS = "OU3v,C30v,C41v,C52v,C52Sv,C52Fv"


def plan(quick: bool):
    """(level, params, scenario, seeds, filters) work items."""
    n_primary, n_secondary, n_long = (4, 2, 1) if quick else (100, 30, 20)
    items = []
    for sc in PRIMARY:
        items.append(("vertical", "tuner", sc, n_primary, VERT_FILTERS))
        items.append(("full", "tuner", sc, n_primary, FULL_FILTERS))
    for sc in SECONDARY:
        items.append(("vertical", "tuner", sc, n_secondary, VERT_FILTERS))
        items.append(("full", "tuner", sc, n_secondary, FULL_FILTERS))
    for sc in QUIET:
        items.append(("vertical", "tuner", sc, n_secondary, VERT_FILTERS + ",C52Nv"))
        items.append(("full", "tuner", sc, n_secondary, FULL_FILTERS + ",C52N"))
    for sc in FULL_ONLY:
        items.append(("full", "tuner", sc, n_secondary, FULL_FILTERS))
    for sc in ORACLE:
        items.append(("vertical", "oracle", sc, n_secondary, VERT_FILTERS))
    for sc in LONG:
        items.append(("full", "tuner", sc, n_long, "OU3,C41,C52,C52S"))
    jobs = []
    for level, params, sc, n, filters in items:
        for seed in range(1, n + 1):
            jobs.append((level, params, sc, seed, filters))
    return jobs


def read_rows(path: Path):
    if not path.exists():
        return []
    with path.open() as f:
        return list(csv.DictReader(f))


def run_bench(cmds, path: Path, jobs: int):
    """Run carma_bench argument lists in parallel, appending their CSV rows to path."""
    header = subprocess.run([str(BENCH), "--header"], capture_output=True, text=True, check=True).stdout.strip()
    new_file = not path.exists()
    with path.open("a") as out:
        if new_file:
            out.write(header + "\n")
        def work(argv):
            return argv, subprocess.run([str(BENCH)] + argv, capture_output=True, text=True)
        with ThreadPoolExecutor(max_workers=jobs) as ex:
            futs = [ex.submit(work, c) for c in cmds]
            for i, fut in enumerate(as_completed(futs), 1):
                argv, p = fut.result()
                if p.returncode != 0:
                    print(f"ERROR {' '.join(argv)}: {p.stderr.strip()[-400:]}", flush=True)
                    continue
                out.write(p.stdout)
                out.flush()
                if i % 50 == 0 or i == len(cmds):
                    print(f"  {i}/{len(cmds)}", flush=True)


def require_bench():
    if not BENCH.exists():
        sys.exit(f"{BENCH} not built: make -C tests/kalman_carma build")
    OUT.mkdir(parents=True, exist_ok=True)


def cmd_run(args):
    require_bench()
    done = {(r["level"], r["params"], r["scenario"], r["seed"]) for r in read_rows(RAW)}
    jobs = [j for j in plan(args.quick) if (j[0], j[1], j[2], str(j[3])) not in done]
    print(f"{len(jobs)} runs to do ({len(done)} done)", flush=True)
    cmds = []
    for level, params, sc, seed, filters in jobs:
        cmd = ["--level", level, "--scenario", sc, "--seed", str(seed), "--filters", filters]
        if level == "vertical":
            cmd += ["--params", params]
        cmds.append(cmd)
    run_bench(cmds, RAW, args.jobs)


def cmd_matched(args):
    """CARMA(5,2) with the shape refitted to each case's true spectrum (oracle shape and operating point)."""
    require_bench()
    mm = OUT / "mismatch_analysis.json"
    if not mm.exists():
        sys.exit(f"{mm} missing: python3 tools/carma_mismatch_analysis.py")
    shapes = {r["case"]: r["matched_c52_shape"] for r in json.loads(mm.read_text())}
    done = {(r["scenario"], r["seed"]) for r in read_rows(MATCHED_RAW)}
    cmds = []
    for sc in MATCHED_CASES:
        shape = ",".join(f"{shapes[sc][k]:.6f}" for k in ("f1", "z1", "f2", "z2", "fl", "fz", "zz"))
        for seed in range(1, args.seeds + 1):
            if (sc, str(seed)) not in done:
                cmds.append(["--level", "vertical", "--params", "oracle", "--scenario", sc, "--seed", str(seed),
                             "--filters", MATCHED_FILTERS, "--shape52", shape])
    print(f"{len(cmds)} matched-shape runs to do", flush=True)
    run_bench(cmds, MATCHED_RAW, args.jobs)


def parse_metrics(text: str):
    rows = []
    for line in text.splitlines():
        if not line.startswith("VALIDATION_METRICS "):     # not the _MEANS lines
            continue
        rows.append(dict(kv.split("=", 1) for kv in line.split()[1:] if "=" in kv))
    return rows


def cmd_timing(args):
    require_bench()
    if TIMING_RAW.exists():
        TIMING_RAW.unlink()
    cmds = [["--level", "full", "--scenario", sc, "--seed", str(seed), "--filters", FULL_FILTERS]
            for sc in ("jonswap_fp0.29_hs0.8_g3.3", "jonswap_fp0.29_hs0.8_g3.3_jit") for seed in range(1, args.seeds + 1)]
    run_bench(cmds, TIMING_RAW, 1)


def cmd_dataset(args):
    """Pinned vessel-RAO records: unchanged OU-III simulator vs carma-sim."""
    data_dir = ROOT / "tests" / "kalman_carma"
    files = sorted(str(p) for p in data_dir.glob("wave_data_*.csv"))
    if not files:
        sys.exit("no dataset in tests/kalman_carma: make ensure-sim-data")
    rows_done = read_rows(DATASET_RAW)
    done = {(r["family"], r["seed"]) for r in rows_done}
    work = []
    for seed in range(1, args.seeds + 1):
        work.append(("OU_III", seed, [str(OU3_SIM)], {}))
        for v in args.variants.split(","):
            work.append((f"CARMA_{v}", seed, [str(CARMA_SIM)], {"W3D_CARMA_VARIANT": v}))
    work = [w for w in work if (w[0], str(w[1])) not in done]
    print(f"{len(work)} dataset invocations to do", flush=True)

    def run(w):
        family, seed, cmd, extra = w
        # W3D_COLLECT_ALL_GATES: the OU-III simulator otherwise stops at its first
        # failed regression gate (as in CI, every record is processed).
        env = dict(os.environ, W3D_SEED=str(seed), W3D_VALIDATION_WINDOW_SEC="900",
                   W3D_WRITE_TIMESERIES="0", W3D_COLLECT_ALL_GATES="1", **extra)
        argv = cmd + sum((["--input", f] for f in files), [])
        p = subprocess.run(argv, capture_output=True, text=True, env=env, cwd=str(data_dir))
        return w, p

    fields = ["family", "seed", "input", "disp_z_rms_m", "disp_z_pct_hs", "disp_x_rms_m", "disp_y_rms_m",
              "roll_rms_deg", "pitch_rms_deg", "yaw_rms_deg", "accel_z_rms_mps2"]
    new_file = not DATASET_RAW.exists()
    with DATASET_RAW.open("a") as out:
        wr = csv.DictWriter(out, fieldnames=fields, extrasaction="ignore")
        if new_file:
            wr.writeheader()
        with ThreadPoolExecutor(max_workers=args.jobs) as ex:
            for fut in as_completed([ex.submit(run, w) for w in work]):
                (family, seed, _, _), p = fut.result()
                metrics = parse_metrics(p.stdout)
                if p.returncode not in (0, 1) or not metrics:
                    print(f"ERROR {family} seed {seed}: rc={p.returncode} {p.stderr[-300:]}", flush=True)
                    continue
                for m in metrics:
                    m["family"] = family
                    m["seed"] = str(seed)
                    m["input"] = Path(m.get("input", "")).name
                    wr.writerow(m)
                out.flush()
                print(f"  {family} seed {seed}: {len(metrics)} records", flush=True)


def cmd_traces(args):
    TRACE_DIR.mkdir(parents=True, exist_ok=True)
    for level, sc, seed, filters in [
        ("full", "jonswap_fp0.29_hs0.8_g3.3", 1, "OU3,C41,C52"),
        ("full", "mixed_swell_wind", 1, "OU3,C41,C52"),
        ("full", "stationary_calbias", 1, "OU3,C52,C52N"),
        ("vertical", "jonswap_fp0.29_hs0.8_g3.3", 1, "OU3v,C41v,C52v"),
    ]:
        path = TRACE_DIR / f"{level}_{sc}_s{seed}.csv"
        subprocess.run([str(BENCH), "--level", level, "--scenario", sc, "--seed", str(seed),
                        "--filters", filters, "--timeseries", str(path), "--decimate", "4"],
                       check=True, capture_output=True)
        print("wrote", path.relative_to(ROOT))


def cmd_truth(args):
    """Motion levels of each scenario's realization (seed 1), for the plausibility table."""
    require_bench()
    scenarios = list(dict.fromkeys(PRIMARY + SECONDARY + QUIET + FULL_ONLY))
    rows = []
    for sc in scenarios:
        out = subprocess.run([str(BENCH), "--level", "truth", "--scenario", sc, "--seed", "1"],
                             capture_output=True, text=True, check=True).stdout
        rows.append(dict(kv.split("=", 1) for kv in out.split()))
    with TRUTH_STATS.open("w") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0]))
        wr.writeheader()
        wr.writerows(rows)
    print("wrote", TRUTH_STATS.relative_to(ROOT))


# --------------------------------------------------------------------------
# Statistics
# --------------------------------------------------------------------------

def bootstrap_ci(x, stat=np.mean, n=10000, seed=20261025):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) < 2:
        return (np.nan, np.nan)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(x), size=(n, len(x)))
    s = stat(x[idx], axis=1)
    return (float(np.percentile(s, 2.5)), float(np.percentile(s, 97.5)))


def paired(a, b):
    """Candidate a vs baseline b on the same seeds (lower is better)."""
    from scipy.stats import wilcoxon

    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    m = np.isfinite(a) & np.isfinite(b)
    a, b = a[m], b[m]
    d = a - b
    rel = 1.0 - a.mean() / b.mean() if b.mean() > 0 else np.nan
    rng = np.random.default_rng(20261025)
    idx = rng.integers(0, len(a), size=(10000, len(a)))
    rel_bs = 1.0 - a[idx].mean(axis=1) / b[idx].mean(axis=1)
    try:
        p = float(wilcoxon(a, b).pvalue) if len(a) >= 6 and np.any(d != 0) else np.nan
    except ValueError:
        p = np.nan
    return {
        "n": int(len(a)),
        "mean_diff": float(d.mean()),
        "diff_ci": bootstrap_ci(d),
        "improvement": float(rel),
        "improvement_ci": (float(np.percentile(rel_bs, 2.5)), float(np.percentile(rel_bs, 97.5))),
        "win_rate": float(np.mean(d < 0)),
        "wilcoxon_p": p,
    }


def fnum(r, k):
    try:
        return float(r[k])
    except (KeyError, ValueError):
        return np.nan


METRICS = ["rms_h", "rms_h_pct_hs", "mean_h", "peak_h", "rms_v", "rms_a", "lf_rms_h", "trend_m_per_hr",
           "nees2", "nis", "nis_in95", "f60_rms_h", "f60_mean_h", "conv_s", "runtime_ns", "memory_bytes", "failed"]


def summarize(rows):
    """{(level, params, scenario): {filter: {metric: array over seeds}}}, seed-aligned."""
    groups = {}
    for r in rows:
        key = (r["level"], r["params"], r["scenario"])
        groups.setdefault(key, {}).setdefault(r["filter"], {})[int(r["seed"])] = r
    out = {}
    for key, by_filter in groups.items():
        seeds = sorted(set.intersection(*(set(v) for v in by_filter.values())))
        out[key] = {f: {m: np.array([fnum(v[s], m) for s in seeds]) for m in METRICS}
                    for f, v in by_filter.items()}
        out[key]["_seeds"] = seeds
    return out


def baseline_name(level):
    return "OU3" if level == "full" else "OU3v"


def fmt_ci(ci, scale=1.0, digits=1):
    return f"[{ci[0] * scale:.{digits}f}, {ci[1] * scale:.{digits}f}]"


def cmd_report(args):
    rows = read_rows(RAW)
    if not rows:
        sys.exit("no raw runs; run first")
    S = summarize(rows)
    stats = {}
    lines = []
    lines.append("# CARMA vs OU-III benchmark -- generated tables\n")
    lines.append("Generated by `tools/carma_benchmark.py report` from `raw_runs.csv`"
                 " and `dataset_runs.csv`. Heave errors are scored on t >= 300 s of"
                 " each record (final-60-s columns on the last minute). `RMS %Hs` is"
                 " relative to the incident significant wave height. Paired"
                 " statistics: improvement = 1 - mean(candidate)/mean(OU-III), 95 %"
                 " bootstrap interval over seeds; win = fraction of seeds on which"
                 " the candidate has lower RMS; p = Wilcoxon signed-rank.\n")
    # The requested headline table: complete INS, primary sea states.
    lines.append("\n## Headline -- complete INS, OU-III vs CARMA(5,2), primary sea states\n")
    lines.append("Mean over seeds; improvement = 1 - RMS(C52)/RMS(OU-III) with 95 % bootstrap interval;"
                 " drift = RMS of the 60 s moving average of the heave error and median |linear trend|.\n")
    lines.append("| Sea state | n | OU-III RMS | CARMA(5,2) RMS | Improvement | Win | Drift: LF RMS OU-III / C52 |"
                 " Drift: median abs trend OU-III / C52 (mm/h) | Final 60 s RMS OU-III / C52 |")
    lines.append("|" + "---|" * 9)
    for sc in PRIMARY:
        g = S.get(("full", "tuner", sc))
        if not g or "OU3" not in g or "C52" not in g:
            continue
        b, c = g["OU3"], g["C52"]
        ps = paired(c["rms_h"], b["rms_h"])
        lines.append(
            f"| {sc} | {len(g['_seeds'])} | {np.nanmean(b['rms_h']):.3f} m ({np.nanmean(b['rms_h_pct_hs']):.2f} % Hs) | "
            f"{np.nanmean(c['rms_h']):.3f} m ({np.nanmean(c['rms_h_pct_hs']):.2f} % Hs) | "
            f"{100 * ps['improvement']:+.1f} % {fmt_ci(ps['improvement_ci'], 100)} | {100 * ps['win_rate']:.0f} % | "
            f"{1000 * np.nanmean(b['lf_rms_h']):.1f} / {1000 * np.nanmean(c['lf_rms_h']):.1f} mm | "
            f"{1000 * np.nanmedian(np.abs(b['trend_m_per_hr'])):.1f} / {1000 * np.nanmedian(np.abs(c['trend_m_per_hr'])):.1f} | "
            f"{np.nanmean(b['f60_rms_h']):.3f} / {np.nanmean(c['f60_rms_h']):.3f} m |")

    order_levels = [("vertical", "tuner"), ("vertical", "oracle"), ("full", "tuner")]
    titles = {("vertical", "tuner"): "Level 1 -- vertical subproblem, prescribed attitude, shared measured operating point",
              ("vertical", "oracle"): "Level 1 -- vertical subproblem, prescribed attitude, oracle operating point",
              ("full", "tuner"): "Level 2 -- complete INS (shipping OU-III vs CARMA INS), identical sensors"}
    for level, params in order_levels:
        keys = [k for k in S if k[0] == level and k[1] == params]
        if not keys:
            continue
        base = baseline_name(level)
        lines.append(f"\n## {titles[(level, params)]}\n")
        cands = [f for f in ["C30", "C41", "C52", "C52S", "C52F", "C52N"]
                 if any((f if level == "full" else f + "v") in S[k] for k in keys)]
        hdr = "| Scenario | n | OU-III RMS %Hs | " + " | ".join(
            f"{c} RMS %Hs | {c} vs OU-III" for c in cands) + " |"
        lines.append(hdr)
        lines.append("|" + "---|" * (3 + 2 * len(cands)))
        scen_order = [s for s in PRIMARY + SECONDARY + QUIET + FULL_ONLY + LONG + ORACLE]
        seen = set()
        for sc in scen_order:
            if sc in seen:
                continue
            seen.add(sc)
            key = (level, params, sc)
            if key not in S or base not in S[key]:
                continue
            g = S[key]
            b = g[base]
            cells = [sc, str(len(g["_seeds"])),
                     f"{np.nanmean(b['rms_h_pct_hs']):.2f}" if np.isfinite(np.nanmean(b['rms_h_pct_hs']))
                     else f"{1000 * np.nanmean(b['rms_h']):.1f} mm"]
            for c in cands:
                name = c if level == "full" else c + "v"
                if name not in g:
                    cells += ["", ""]
                    continue
                a = g[name]
                metric = "rms_h_pct_hs" if np.isfinite(np.nanmean(b["rms_h_pct_hs"])) else "rms_h"
                ps = paired(a["rms_h"], b["rms_h"])
                stats.setdefault(f"{level}/{params}/{sc}", {})[name] = ps
                val = (f"{np.nanmean(a[metric]):.2f}" if metric == "rms_h_pct_hs"
                       else f"{1000 * np.nanmean(a['rms_h']):.1f} mm")
                cells += [val, f"{100 * ps['improvement']:+.1f} % {fmt_ci(ps['improvement_ci'], 100)}"
                               f" win {100 * ps['win_rate']:.0f} %"]
            lines.append("| " + " | ".join(cells) + " |")

        # Drift and consistency table for this level.
        lines.append(f"\n### Drift, final minute and consistency ({level}, {params})\n")
        lines.append("| Scenario | Filter | RMS h (m) | mean h (m) | peak h (m) | LF RMS h (m) | final-60 s RMS h (m) |"
                     " final-60 s mean h (m) | RMS v (m/s) | RMS a (m/s^2) | NEES [p,v] (2) | NIS | conv (s) | failures |")
        lines.append("|" + "---|" * 14)
        for sc in [s for s in PRIMARY + QUIET + ["bias_drift"] + LONG if (level, params, s) in S]:
            g = S[(level, params, sc)]
            for f in [base] + [c if level == "full" else c + "v" for c in cands]:
                if f not in g:
                    continue
                a = g[f]
                lines.append(
                    f"| {sc} | {f} | {np.nanmean(a['rms_h']):.4f} | {np.nanmean(a['mean_h']):+.4f} | "
                    f"{np.nanmean(a['peak_h']):.3f} | {np.nanmean(a['lf_rms_h']):.4f} | "
                    f"{np.nanmean(a['f60_rms_h']):.4f} | {np.nanmean(a['f60_mean_h']):+.4f} | "
                    f"{np.nanmean(a['rms_v']):.4f} | {np.nanmean(a['rms_a']):.4f} | "
                    f"{np.nanmean(a['nees2']):.2f} | {np.nanmean(a['nis']):.3f} | "
                    f"{np.nanmedian(a['conv_s']):.0f} | {int(np.nansum(a['failed']))} |")

    # Oracle-shape diagnostic.
    mrows = read_rows(MATCHED_RAW)
    if mrows:
        M = summarize(mrows)
        lines.append("\n## Level 1 -- oracle-shape diagnostic (not an independent result)\n")
        lines.append("C52M is CARMA(5,2) with its shape refitted to the true vessel-heave spectrum of the case"
                     " (`matched_c52_shape` in `mismatch_analysis.json`) and the oracle operating point: it uses"
                     " the truth and bounds what this model class can do with a perfect shape. C52 is the fixed"
                     " JONSWAP-fitted shape on the same oracle operating point.\n")
        lines.append("| Scenario | n | OU-III RMS %Hs | C41 RMS %Hs | C41 vs OU-III | C52 RMS %Hs | C52 vs OU-III |"
                     " C52M RMS %Hs | C52M vs OU-III | NEES OU-III / C41 / C52 / C52M |")
        lines.append("|" + "---|" * 10)
        for sc in MATCHED_CASES:
            g = M.get(("vertical", "oracle", sc))
            if not g or "OU3v" not in g:
                continue
            b = g["OU3v"]
            cells = [sc, str(len(g["_seeds"])), f"{np.nanmean(b['rms_h_pct_hs']):.2f}"]
            for c in ("C41v", "C52v", "C52Mv"):
                ps = paired(g[c]["rms_h"], b["rms_h"])
                stats.setdefault(f"matched/{sc}", {})[c] = ps
                cells += [f"{np.nanmean(g[c]['rms_h_pct_hs']):.2f}",
                          f"{100 * ps['improvement']:+.1f} % {fmt_ci(ps['improvement_ci'], 100)}"
                          f" win {100 * ps['win_rate']:.0f} %"]
            cells.append(" / ".join(f"{np.nanmean(g[c]['nees2']):.2f}" for c in ("OU3v", "C41v", "C52v", "C52Mv")))
            lines.append("| " + " | ".join(cells) + " |")

    trows = read_rows(TRUTH_STATS)
    if trows:
        lines.append("\n## Realized motion levels (seed 1)\n")
        lines.append("Truth only (`carma_bench --level truth`): vessel significant heave, vertical acceleration,"
                     " roll and pitch about the mean attitude, peak body rate and specific force, and the incident"
                     " steepness H_s / L_p of the steepest spectral part.\n")
        lines.append("| Scenario | H_s heave (m) | a RMS / max (m/s^2) | roll RMS / max (deg) | pitch RMS / max (deg) |"
                     " max rate (rad/s) | max |f| (m/s^2) | H_s/L_p |")
        lines.append("|" + "---|" * 8)
        for r in trows:
            lines.append(f"| {r['scenario']} | {float(r['hs_heave_m']):.2f} | {float(r['a_rms']):.2f} / {float(r['a_max']):.2f} | "
                         f"{float(r['roll_rms_deg']):.1f} / {float(r['roll_max_deg']):.1f} | "
                         f"{float(r['pitch_rms_deg']):.1f} / {float(r['pitch_max_deg']):.1f} | "
                         f"{float(r['rate_max_rad_s']):.2f} | {float(r['f_max_mps2']):.1f} | {float(r['steepness']):.3f} |")

    # Runtime / memory (host).
    trows = read_rows(TIMING_RAW)
    lines.append("\n## Host runtime and memory (full INS, x86-64, -O3)\n")
    lines.append(("Serial runs on an otherwise idle host (`timing`), JONSWAP 0.29 Hz, "
                  f"{len({r['seed'] for r in trows})} seeds; " if trows else "Grid runs under parallel load; ")
                 + "time per `update()` after Live, all filters interleaved sample by sample on one stream.\n")
    lines.append("| Filter | states | mean update (us) | object bytes |")
    lines.append("|---|---|---|---|")
    full_rows = [r for r in (trows or rows) if r["level"] == "full"]
    for f in ["OU3", "C30", "C41", "C52", "C52S", "C52F"]:
        rr = [r for r in full_rows if r["filter"] == f and "_jit" not in r["scenario"]]
        rj = [r for r in full_rows if r["filter"] == f and "_jit" in r["scenario"]]
        if not rr:
            continue
        lines.append(f"| {f} | {rr[0]['states']} | {np.mean([fnum(r, 'runtime_ns') for r in rr]) / 1000:.2f}"
                     f" (jittered dt: {np.mean([fnum(r, 'runtime_ns') for r in rj]) / 1000 if rj else float('nan'):.2f})"
                     f" | {rr[0]['memory_bytes']} |")

    # Dataset replay.
    drows = read_rows(DATASET_RAW)
    if drows:
        lines.append("\n## Pinned vessel-RAO dataset (oceanography-waves-lib v1.2.3), shared W3dSimulationRunner\n")
        lines.append("Trailing 900 s of each 1200 s record; `disp_z_pct_hs` relative to the record's incident Hs."
                     " OU-III is the unchanged `kalman_ou_iii-sim`.\n")
        fams = sorted({r["family"] for r in drows}, key=lambda s: (s != "OU_III", s))
        recs = sorted({r["input"] for r in drows})
        lines.append("| Record | " + " | ".join(fams) + " |")
        lines.append("|" + "---|" * (1 + len(fams)))
        table = {}
        for r in drows:
            table.setdefault((r["input"], r["family"]), []).append(fnum(r, "disp_z_pct_hs"))
        for rec in recs:
            cells = [rec.replace("wave_data_", "").replace(".csv", "")]
            for fam in fams:
                v = table.get((rec, fam), [])
                cells.append(f"{np.mean(v):.2f} (n={len(v)})" if v else "")
            lines.append("| " + " | ".join(cells) + " |")
        # Paired pooled comparison over (record, seed).
        def pooled(fam):
            d = {}
            for r in drows:
                if r["family"] == fam:
                    d[(r["input"], r["seed"])] = fnum(r, "disp_z_rms_m")
            return d
        base = pooled("OU_III")
        for fam in fams:
            if fam == "OU_III":
                continue
            c = pooled(fam)
            common = sorted(set(base) & set(c))
            if len(common) < 4:
                continue
            ps = paired([c[k] for k in common], [base[k] for k in common])
            stats.setdefault("dataset", {})[fam] = ps
            lines.append(f"\n{fam} vs OU-III over {ps['n']} (record, seed) pairs: improvement "
                         f"{100 * ps['improvement']:+.1f} % {fmt_ci(ps['improvement_ci'], 100)}, "
                         f"win {100 * ps['win_rate']:.0f} %, Wilcoxon p = {ps['wilcoxon_p']:.2g}")

    (OUT / "summary.md").write_text("\n".join(lines) + "\n")
    (OUT / "paired_stats.json").write_text(json.dumps(stats, indent=2, sort_keys=True, default=float) + "\n")
    print("wrote", (OUT / "summary.md").relative_to(ROOT))
    if not args.no_plots:
        figures(S, rows)


# --------------------------------------------------------------------------
# Figures
# --------------------------------------------------------------------------

COLORS = {"OU3": "#1f77b4", "C30": "#2ca02c", "C41": "#ff7f0e", "C52": "#d62728",
          "C52S": "#9467bd", "C52F": "#8c564b", "C52N": "#e377c2"}


def _plt():
    import importlib.util
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    spec = importlib.util.spec_from_file_location("svg_determinism", ROOT / "plots" / "svg_determinism.py")
    svgd = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(svgd)
    svgd.configure_svg(matplotlib)
    return plt, svgd


def figures(S, rows):
    plt, svgd = _plt()

    # RMS distributions over seeds, primary scenarios, both levels.
    for level, params in (("full", "tuner"), ("vertical", "tuner")):
        fig, axes = plt.subplots(1, len(PRIMARY), figsize=(11.0, 3.6), sharey=False)
        for ax, sc in zip(axes, PRIMARY):
            key = (level, params, sc)
            if key not in S:
                continue
            g = S[key]
            names = [n for n in ["OU3", "C30", "C41", "C52", "C52S", "C52F"]
                     if (n if level == "full" else n + "v") in g]
            data = [g[n if level == "full" else n + "v"]["rms_h_pct_hs"] for n in names]
            bp = ax.boxplot(data, widths=0.6, patch_artist=True, showfliers=True,
                            flierprops=dict(markersize=2))
            for patch, n in zip(bp["boxes"], names):
                patch.set_facecolor(COLORS[n]); patch.set_alpha(0.6)
            ax.set_xticks(range(1, len(names) + 1))
            ax.set_xticklabels([n.replace("OU3", "OU-III") for n in names], rotation=45, fontsize=8)
            # A few startup excursions would otherwise set the scale; count them instead.
            top = 1.25 * max(np.nanpercentile(x, 98) for x in data)
            n_above = sum(int(np.sum(x > top)) for x in data)
            ax.set_ylim(0.0, top)
            ax.set_title(f"{sc}\n(n={len(g['_seeds'])}" + (f", {n_above} points above axis)" if n_above else ")"),
                         fontsize=8)
            ax.grid(True, axis="y", color="#e5e5e5", lw=0.6)
        axes[0].set_ylabel("Heave RMS error (% of incident Hs)")
        fig.suptitle(("Complete INS" if level == "full" else "Vertical subproblem, prescribed attitude")
                     + ": heave RMS over Monte Carlo seeds", fontsize=10)
        fig.tight_layout()
        svgd.save_svg(fig, OUT / f"rms_distribution_{level}.svg", bbox_inches="tight")
        plt.close(fig)

    # Paired improvement across all scenarios (full INS) for C41 and C52.
    keys = [k for k in S if k[0] == "full" and k[1] == "tuner" and "OU3" in S[k]]
    keys = [k for k in keys if np.isfinite(np.nanmean(S[k]["OU3"]["rms_h_pct_hs"]))]
    keys.sort(key=lambda k: (k[2] not in PRIMARY, k[2]))
    fig, ax = plt.subplots(figsize=(7.5, 0.28 * len(keys) + 1.2))
    y = np.arange(len(keys))
    for off, c in ((-0.18, "C41"), (0.0, "C52"), (0.18, "C30")):
        vals, lo, hi = [], [], []
        for k in keys:
            if c not in S[k]:
                vals.append(np.nan); lo.append(np.nan); hi.append(np.nan); continue
            ps = paired(S[k][c]["rms_h"], S[k]["OU3"]["rms_h"])
            vals.append(100 * ps["improvement"]); lo.append(100 * ps["improvement_ci"][0]); hi.append(100 * ps["improvement_ci"][1])
        vals, lo, hi = map(np.array, (vals, lo, hi))
        ax.errorbar(vals, y + off, xerr=[vals - lo, hi - vals], fmt="o", ms=3, color=COLORS[c], label=c, lw=1)
    ax.axvline(0, color="#000000", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([k[2] for k in keys], fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("Heave RMS improvement vs OU-III (%), 95 % bootstrap CI; > 0 is better")
    ax.set_title("Complete INS: paired improvement over OU-III", fontsize=10)
    ax.grid(True, axis="x", color="#e5e5e5", lw=0.6)
    ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    svgd.save_svg(fig, OUT / "paired_improvement_full.svg", bbox_inches="tight")
    plt.close(fig)

    # Runtime / memory.
    full_rows = [r for r in (read_rows(TIMING_RAW) or rows) if r["level"] == "full" and "_jit" not in r["scenario"]]
    names = [n for n in ["OU3", "C30", "C41", "C52", "C52S", "C52F"] if any(r["filter"] == n for r in full_rows)]
    rt = [np.mean([fnum(r, "runtime_ns") for r in full_rows if r["filter"] == n]) / 1000 for n in names]
    mem = [fnum(next(r for r in full_rows if r["filter"] == n), "memory_bytes") / 1024 for n in names]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(8.0, 3.0))
    a1.bar(range(len(names)), rt, color=[COLORS[n] for n in names])
    a2.bar(range(len(names)), mem, color=[COLORS[n] for n in names])
    for a, lab in ((a1, "Mean update time (us, host x86-64)"), (a2, "Filter + wrapper object size (KiB)")):
        a.set_xticks(range(len(names)))
        a.set_xticklabels([n.replace("OU3", "OU-III") for n in names], fontsize=8)
        a.set_title(lab, fontsize=9)
        a.grid(True, axis="y", color="#e5e5e5", lw=0.6)
    fig.tight_layout()
    svgd.save_svg(fig, OUT / "runtime_memory.svg", bbox_inches="tight")
    plt.close(fig)

    trace_figures(plt, svgd)


def trace_figures(plt, svgd):
    from scipy.signal import welch

    path = TRACE_DIR / "full_jonswap_fp0.29_hs0.8_g3.3_s1.csv"
    if not path.exists():
        return
    d = np.genfromtxt(path, delimiter=",", names=True)
    t = d["t"]
    names = [n[2:] for n in d.dtype.names if n.startswith("h_") and n != "h_true"]
    # Time series: a late 60 s window.
    m = (t >= 1000) & (t <= 1060)
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(9.0, 5.0), sharex=True)
    a1.plot(t[m], d["h_true"][m], color="#000000", lw=1.4, label="truth (vessel heave)")
    a1.plot(t[m], d["eta_up"][m], color="#999999", lw=0.8, ls=":", label="incident elevation")
    for n in names:
        a1.plot(t[m], d["h_" + n][m], color=COLORS.get(n, "#333333"), lw=0.9, label=n.replace("OU3", "OU-III"))
        a2.plot(t[m], d["h_" + n][m] - d["h_true"][m], color=COLORS.get(n, "#333333"), lw=0.9)
    a1.set_ylabel("Heave up (m)")
    a1.legend(fontsize=7, ncol=5, loc="upper right")
    a2.set_ylabel("Estimate - truth (m)")
    a2.set_xlabel("Time (s)")
    a1.set_title("JONSWAP fp 0.29 Hz, Hs 0.8 m, seed 1, complete INS", fontsize=10)
    for a in (a1, a2):
        a.grid(True, color="#e5e5e5", lw=0.6)
    fig.tight_layout()
    svgd.save_svg(fig, OUT / "timeseries_jonswap.svg", bbox_inches="tight")
    plt.close(fig)

    # Error over the whole record, then PSDs.
    fig, ax = plt.subplots(figsize=(9.0, 3.0))
    for n in names:
        e = d["h_" + n] - d["h_true"]
        ax.plot(t, e, color=COLORS.get(n, "#333333"), lw=0.5, label=n.replace("OU3", "OU-III"))
    ax.set_ylim(-1.0, 1.0)
    ax.axvline(300, color="#000000", lw=0.6, ls="--")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Heave error (m)")
    ax.set_title("Heave error over the record (scoring starts at 300 s)", fontsize=10)
    ax.legend(fontsize=7, ncol=4)
    ax.grid(True, color="#e5e5e5", lw=0.6)
    fig.tight_layout()
    svgd.save_svg(fig, OUT / "error_timeseries_jonswap.svg", bbox_inches="tight")
    plt.close(fig)

    fs = 1.0 / np.median(np.diff(t))
    s = t >= 300
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.0, 3.6))
    f, Pt = welch(d["h_true"][s], fs=fs, nperseg=4096)
    a1.semilogy(f, Pt, color="#000000", lw=1.4, label="truth")
    for n in names:
        _, Pe = welch(d["h_" + n][s], fs=fs, nperseg=4096)
        _, Pr = welch((d["h_" + n] - d["h_true"])[s], fs=fs, nperseg=4096)
        a1.semilogy(f, Pe, color=COLORS.get(n, "#333333"), lw=0.9, label=n.replace("OU3", "OU-III"))
        a2.semilogy(f, Pr, color=COLORS.get(n, "#333333"), lw=0.9, label=n.replace("OU3", "OU-III"))
    for a, title in ((a1, "Heave PSD: truth and estimates"), (a2, "Heave error PSD")):
        a.set_xlim(0, 1.2)
        a.set_xlabel("Frequency (Hz)")
        a.set_ylabel("PSD (m$^2$/Hz)")
        a.set_title(title, fontsize=9)
        a.grid(True, which="both", color="#eeeeee", lw=0.5)
        a.legend(fontsize=7)
    a1.set_ylim(1e-6, 1e2)
    a2.set_ylim(1e-7, 1e0)
    fig.tight_layout()
    svgd.save_svg(fig, OUT / "psd_jonswap.svg", bbox_inches="tight")
    plt.close(fig)

    q = TRACE_DIR / "full_stationary_calbias_s1.csv"
    if q.exists():
        d = np.genfromtxt(q, delimiter=",", names=True)
        names = [n[2:] for n in d.dtype.names if n.startswith("h_") and n != "h_true"]
        fig, ax = plt.subplots(figsize=(9.0, 3.0))
        for n in names:
            ax.plot(d["t"], d["h_" + n], color=COLORS.get(n, "#333333"), lw=0.8, label=n.replace("OU3", "OU-III"))
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Heave estimate (m); truth = 0")
        ax.set_title("Still IMU, post-calibration accelerometer bias (+-0.01 m/s^2 range), seed 1", fontsize=10)
        ax.legend(fontsize=7)
        ax.grid(True, color="#e5e5e5", lw=0.6)
        fig.tight_layout()
        svgd.save_svg(fig, OUT / "stationary_timeseries.svg", bbox_inches="tight")
        plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--jobs", type=int, default=os.cpu_count() or 2)
    r.add_argument("--quick", action="store_true")
    d = sub.add_parser("dataset")
    d.add_argument("--jobs", type=int, default=os.cpu_count() or 2)
    d.add_argument("--seeds", type=int, default=10)
    d.add_argument("--variants", default="C52,C41,C30")
    mt = sub.add_parser("matched")
    mt.add_argument("--jobs", type=int, default=os.cpu_count() or 2)
    mt.add_argument("--seeds", type=int, default=30)
    tm = sub.add_parser("timing")
    tm.add_argument("--seeds", type=int, default=3)
    sub.add_parser("traces")
    sub.add_parser("truth")
    rp = sub.add_parser("report")
    rp.add_argument("--no-plots", action="store_true")
    args = ap.parse_args(argv)
    {"run": cmd_run, "dataset": cmd_dataset, "matched": cmd_matched, "timing": cmd_timing,
     "traces": cmd_traces, "truth": cmd_truth, "report": cmd_report}[args.cmd](args)


if __name__ == "__main__":
    main()
