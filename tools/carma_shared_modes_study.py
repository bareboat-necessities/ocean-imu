#!/usr/bin/env python3
"""Numerical study of the shared-modes QMEKF (EXPERIMENTAL).

Companion of docs/kalman-carma-shared-modes.md and
doc/kalman_carma_shared_modes/kalman_carma_shared_modes-validation.tex.
Drives tests/kalman_carma_shared_modes/kalman_carma_shared_modes-bench
(paired Monte Carlo, carma3_bench's engine: every filter of a row sees the
same sea and the same IMU stream) and kalman_carma_shared_modes-sim (dataset
replay), and writes reports/results/kalman_carma_shared_modes/.

Filters (CsmBenchAdapter.h):
    A41@C       CARMA3 baseline, Evidence excitation (21 states)
    CSM         shared-modes QMEKF, deployable configuration (27 states)
    CSM:noconf  ... with the coupling confidence disabled (robustness stage)
    CSM:off     ... with the coupling off: A41@C exactly (smoke stage)

Stages (resumable; each appends to its own CSV):
    smoke     one paired row with the reduction identity CSM:off = A41@C
    run       paired Monte Carlo: SCENARIOS x seeds 1..N (default 30)
    robust    wrong coupling priors (external, from the truth, deliberately
              wrong) with and without the confidence; heading changes
    traces    10 Hz time series around sea-state changes, seeds 1..10
    dataset   replay of the pinned vessel-RAO records (default seeds)
    gates     the regression-gate quantities of the default-seed replays: CSM
              and TFG-CARMA3 (T41; Legacy and Evidence), CSM's bars
    timing    host update time, one job (ranking only, not an ESP32-S3 time)
    size      ESP32-S3 size / stack probe (needs the Xtensa toolchain: --toolchain)
    cost      analytic multiply-add model
    report    statistics, LaTeX tables, figures (SVG), manifest.json
    figures   the figures again (--outdir, --ext svg,pgf)

Seeds 1..N are the evaluation seeds; nothing in the filter is fitted to them.
Uncertainty: paired bootstrap over seeds.

Run: python3 tools/carma_shared_modes_study.py {smoke,run,robust,traces,dataset,gates,timing,size,cost,report,figures,all} [--jobs 4]
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

sys.path.insert(0, str(Path(__file__).resolve().parent))
import carma3_analysis  # noqa: E402  (analytic cost model of the CARMA3 core)

ROOT = Path(__file__).resolve().parents[1]
TDIR = ROOT / "tests" / "kalman_carma_shared_modes"
SIM = TDIR / "kalman_carma_shared_modes-sim"        # dataset replay
BENCH = TDIR / "kalman_carma_shared_modes-bench"    # paired Monte Carlo
OUT = ROOT / "reports" / "results" / "kalman_carma_shared_modes"
RAW = OUT / "raw_runs.csv"
ROB = OUT / "robustness_runs.csv"
SMOKE = OUT / "smoke_runs.csv"
TIM = OUT / "timing_runs.csv"
DATASET = OUT / "dataset_summary.csv"
SIZE = OUT / "esp32s3_size.csv"
COST = OUT / "cost_model.csv"
BARS = OUT / "csm_regression_bars.json"
TRACE_DIR = OUT / "traces"
FIG_DIR = OUT / "figures"
DOC_DIR = ROOT / "doc" / "kalman_carma_shared_modes"

J29 = "jonswap_fp0.29_hs0.8_g3.3"
# scenario -> (group, description)
SCENARIOS = {
    J29: ("fixed", "directionally spread JONSWAP (cos^8)"),
    f"{J29}+sp10000": ("coupling", "long-crested JONSWAP (strongest coupling)"),
    f"{J29}+sp1": ("coupling", "broadly spread JONSWAP (cos^2, weak coupling)"),
    f"{J29}+rel90": ("coupling", "beam sea (phase and RAO change the coupling)"),
    f"{J29}+rel180": ("coupling", "head sea"),
    f"{J29}+hscale0.5": ("coupling", "heave RAO x0.5 (heave weaker than the surge it couples to)"),
    f"{J29}+hscale1.5": ("coupling", "heave RAO x1.5"),
    "regular_f0.25_h1": ("fixed", "regular monochromatic wave"),
    "pm_fp0.25_hs1": ("fixed", "Pierson-Moskowitz"),
    "swell": ("fixed", "long-period swell"),
    "jonswap_fp0.45_hs0.4_g3.3": ("fixed", "short-period wind waves"),
    "mixed_swell_wind": ("fixed", "mixed wind sea + swell"),
    "bimodal": ("fixed", "two systems, different frequencies"),
    "crossing": ("fixed", "two crossing systems (90 deg)"),
    "small_amplitude": ("fixed", "very low sea state (5 cm)"),
    "changing_fp": ("gradual", "peak frequency 0.22 -> 0.36 Hz over 600 s"),
    "abrupt_change": ("rapid", "wind sea -> long higher sea in 10 s (period, direction, height)"),
    "dir_change": ("rapid", "direction 30 -> 120 deg over 60 s"),
    f"{J29}+on400w60+off900w60": ("rapid", "wave onset and cessation"),
    f"{J29}+turns": ("disturbance", "heading changes (90 deg zigzag)"),
    f"{J29}+maneuver": ("disturbance", "vessel manoeuvring"),
    f"{J29}+cnoise": ("disturbance", "coloured accelerometer noise (4 mg, 2 s)"),
    f"{J29}+nvar": ("disturbance", "accelerometer noise variance cycling 1x-2x"),
    f"{J29}+jit": ("disturbance", "+-25 % sample-interval jitter"),
    f"{J29}+magdist+scal": ("disturbance", "magnetic disturbances, calibration residuals"),
    "bias_drift": ("disturbance", "accelerometer bias ramp and oscillation"),
}
FILTERS = ["A41@C", "CSM"]
BASE = "A41@C"         # the CARMA3 baseline every change is judged against
LABEL = {"A41@C": "A41 (21 states, CARMA3 baseline)", "CSM": "CSM (27 states, shared modes)"}
COLORS = {"A41@C": "#555555", "CSM": "#1f77b4", "CSM:noconf": "#9edae5", "truth": "#000000"}

ROB_SEAS = [J29, "mixed_swell_wind", f"{J29}+sp10000"]
ROB_WRONG = ["flip", "rot90", "ph45", "amp2"]
ROB_FILTERS = ",".join(["A41@C", "CSM/ext-true"] + [f"CSM{c}/ext-{w}" for w in ROB_WRONG for c in ("", ":noconf")])
ROB_TURNS = f"{J29}+turns"
ROB_TURN_FILTERS = "A41@C,CSM,CSM/fz400/keel"

TRACE_SEAS = ["changing_fp", "abrupt_change", "dir_change", f"{J29}+on400w60+off900w60"]
TRACE_FILTERS = "A41@C,CSM"

DATASET_VARIANTS = ["CSM", "CSM_OFF"]


def require(exe=None):
    for x in ([exe] if exe else [SIM, BENCH]):
        if not x.exists():
            sys.exit(f"missing {x.relative_to(ROOT)}: make -C tests/kalman_carma_shared_modes build")


def read_rows(path):
    if not path.exists():
        return []
    with path.open() as f:
        return list(csv.DictReader(f))


def cmd(sc, seed, filters, score_start=300.0):
    return ["--level", "full", "--scenario", sc, "--seed", str(seed), "--filters", filters,
            "--score-start", f"{score_start:g}"]


def run_rows(cmds, out_path, jobs, extra=("--extended", "--coupling")):
    require()
    extra = list(extra)
    hdr = subprocess.run([str(BENCH), "--header"] + extra, capture_output=True, text=True, check=True).stdout.strip()
    new_file = not out_path.exists()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    failures = 0
    with out_path.open("a") as out:
        if new_file:
            out.write(hdr + "\n")
        with ThreadPoolExecutor(max_workers=jobs) as ex:
            futs = {ex.submit(subprocess.run, [str(BENCH)] + extra + c, capture_output=True, text=True): c for c in cmds}
            for k, fut in enumerate(as_completed(futs), 1):
                c, p = futs[fut], fut.result()
                if p.returncode != 0:
                    failures += 1
                    print(f"ERROR {' '.join(c)}: rc={p.returncode} {p.stderr[-400:]}", flush=True)
                    continue
                out.write(p.stdout)
                out.flush()
                if k % 20 == 0 or k == len(cmds):
                    print(f"  {k}/{len(cmds)}", flush=True)
    return failures


def pending(path, work):
    have = {(r["scenario"], r["seed"]) for r in read_rows(path)}
    return [w for w in work if (w[0], str(w[1])) not in have]


def paired(df, sc, cand, col, base=BASE, n_boot=4000):
    """Relative change 1 - cand / base of the seed means (+ better), with a
    paired bootstrap 95 % interval over seeds."""
    a = df[(df.scenario == sc) & (df["filter"] == base)].set_index("seed")[col]
    b = df[(df.scenario == sc) & (df["filter"] == cand)].set_index("seed")[col]
    idx = a.index.intersection(b.index)
    a, b = a.loc[idx].values.astype(float), b.loc[idx].values.astype(float)
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    if len(a) < 3:
        return None
    rng = np.random.default_rng(20261106)
    k = rng.integers(0, len(a), (n_boot, len(a)))
    bs = 1.0 - b[k].mean(axis=1) / a[k].mean(axis=1)
    lo, hi = np.percentile(bs, [2.5, 97.5])
    return {"base": float(a.mean()), "cand": float(b.mean()), "change": float(1.0 - b.mean() / a.mean()),
            "ci": [float(lo), float(hi)], "n": int(len(a))}


def git_rev():
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    except OSError:
        return "unknown"


def tex(t):
    return t.replace("_", "\\_").replace("%", "\\%").replace("&", "\\&").replace("->", "$\\to$")


# ---------------------------------------------------------------------------
# stages
# ---------------------------------------------------------------------------

def stage_smoke(args):
    if SMOKE.exists():
        SMOKE.unlink()
    run_rows([cmd(f"{J29}+T600", 1, "A41@C,CSM,CSM:noconf,CSM:off", 200.0)], SMOKE, 1)
    rows = {r["filter"]: r for r in read_rows(SMOKE)}

    def rel(a, b, k):
        x, y = float(rows[a][k]), float(rows[b][k])
        return abs(x - y) / max(abs(y), 1e-12)

    # Coupling off, the filter IS the CARMA3 baseline in law, through adaptation
    # (float rounding only: carma_shared_modes_realization-test [3] checks it in double).
    for k in ("rms_x", "rms_y", "rms_z", "nees_pos3", "nees_pv6"):
        assert rel("CSM:off", "A41@C", k) < 1e-4, ("CSM:off must reproduce A41@C", k, rel("CSM:off", "A41@C", k))
    for f in ("CSM", "CSM:noconf"):
        assert rows[f]["failed"] == "0", f
    print("smoke ok:", {f: round(1000 * float(r["rms_h"]), 2) for f, r in rows.items()})


def stage_run(args):
    work = [(sc, s) for sc in SCENARIOS for s in range(1, args.seeds + 1)]
    todo = pending(RAW, work)
    print(f"run: {len(todo)} of {len(work)} rows to run", flush=True)
    run_rows([cmd(sc, s, ",".join(FILTERS)) for sc, s in todo], RAW, args.jobs)


def stage_robust(args):
    work = [(sc, s, ROB_FILTERS) for sc in ROB_SEAS for s in range(1, args.robust_seeds + 1)]
    work += [(ROB_TURNS, s, ROB_TURN_FILTERS) for s in range(1, args.robust_seeds + 1)]
    have = {(r["scenario"], r["seed"]) for r in read_rows(ROB)}
    todo = [w for w in work if (w[0], str(w[1])) not in have]
    print(f"robust: {len(todo)} of {len(work)} rows to run", flush=True)
    run_rows([cmd(sc, s, f) for sc, s, f in todo], ROB, args.jobs)


def stage_traces(args):
    TRACE_DIR.mkdir(parents=True, exist_ok=True)
    cmds = []
    for sc in TRACE_SEAS:
        tag = re.sub(r"[^A-Za-z0-9_.+-]", "_", sc)
        for s in range(1, args.trace_seeds + 1):
            p = TRACE_DIR / f"{tag}_s{s}_ts.csv"
            if p.exists():
                continue
            cmds.append(cmd(sc, s, TRACE_FILTERS) + ["--timeseries", str(p), "--decimate", "20"])
    if cmds:
        require()
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = [ex.submit(subprocess.run, [str(BENCH), "--extended", "--coupling"] + c, capture_output=True, text=True)
                for c in cmds]
        for f in futs:
            p = f.result()
            if p.returncode != 0:
                print("ERROR trace:", p.stderr[-300:])
    # Reduce to per-scenario ensemble RMS envelopes (the raw series are large
    # and not committed).
    import pandas as pd
    for sc in TRACE_SEAS:
        tag = re.sub(r"[^A-Za-z0-9_.+-]", "_", sc)
        files = sorted(TRACE_DIR.glob(f"{tag}_s*_ts.csv"))
        if not files:
            continue
        acc, n, t = {}, 0, None
        for p in files:
            d = pd.read_csv(p)
            if t is None:
                t = d.t.values
            m = min(len(t), len(d))
            for fl in TRACE_FILTERS.split(","):
                for ax in "xyz":
                    e = (d[f"{ax}_{fl}"].values[:m] - d[f"{ax}_true"].values[:m]) ** 2
                    acc[(fl, ax)] = acc.get((fl, ax), 0) + e
            for ax in "xyz":
                acc[("truth", ax)] = acc.get(("truth", ax), 0) + d[f"{ax}_true"].values[:m] ** 2
            t = t[:m]
            n += 1
        # Ensemble mean square per sample, averaged over 1 s blocks (the
        # figures smooth over 10 s): the committed envelope stays small.
        blk = max(1, int(round(1.0 / float(np.median(np.diff(t))))))
        nb = len(t) // blk
        out = {"t": t[:nb * blk].reshape(nb, blk).mean(axis=1)}
        for (fl, ax), v in acc.items():
            out[f"{fl}|{ax}"] = np.sqrt(v[:nb * blk].reshape(nb, blk).mean(axis=1) / n)
        pd.DataFrame(out).to_csv(TRACE_DIR / f"{tag}_envelope.csv", index=False, float_format="%.4g")
        print(f"traces {sc}: {n} seeds", flush=True)


# ---------------------------------------------------------------------------
# regression gates of the dataset replay (the quantities print_summary_and_fail_if_needed scores)
# ---------------------------------------------------------------------------

GATES = OUT / "dataset_gates.csv"
TFG_SIM = ROOT / "tests" / "kalman_tfg_carma" / "tfg_carma3-sim"
# label -> (simulator, environment)
GATE_RUNS = {
    "CSM": (SIM, {"W3D_CSM_VARIANT": "CSM"}),
    "T41": (TFG_SIM, {}),                                      # TFG-CARMA3, its gated default (Legacy excitation)
    "T41@C": (TFG_SIM, {"W3D_CARMA3_EXCITATION": "C"}),       # ... Evidence excitation (like for like; diagnostic)
}
GATE_KEYS = ["z_pct_hs", "pct_3d", "roll_deg", "pitch_deg", "yaw_deg", "acc_z_bias_pct", "acc_bias_3d_pct",
             "gyro_bias_3d_pct"]


def parse_gates(text):
    """Per-record gate quantities from a replay's summary blocks."""
    rows, cur = [], None
    num = r"([-+0-9.eE]+|nan|inf)"
    for line in text.splitlines():
        m = re.match(r"=== Last 900 s RMS summary for \S*w3d_(\S+?)_fusion_", line)
        if m:
            cur = {"record": m.group(1)}
            continue
        if cur is None:
            continue
        if line.startswith("XYZ RMS (%Hs)"):
            cur["z_pct_hs"] = float(re.search(r"Z=" + num + "%", line).group(1))
        elif line.startswith("3D RMS (m)"):
            cur["pct_3d"] = float(re.search(r"3D % of max \|disp_ref\|_3D = " + num + "%", line).group(1))
        elif line.startswith("Angles RMS (deg)"):
            for k, n in (("roll_deg", "Roll"), ("pitch_deg", "Pitch"), ("yaw_deg", "Yaw")):
                cur[k] = float(re.search(n + "=" + num, line).group(1))
        elif line.startswith("Bias error RMS (% of max TRUE bias) (acc)"):
            cur["acc_z_bias_pct"] = float(re.search(r" Z=" + num + "%", line).group(1))
            cur["acc_bias_3d_pct"] = float(re.search(r"\|3D\|=" + num + "%", line).group(1))
        elif line.startswith("Bias error RMS (% of max TRUE bias) (gyro)"):
            cur["gyro_bias_3d_pct"] = float(re.search(r"\|3D\|=" + num + "%", line).group(1))
        elif line.startswith("QUALITY_GATE:"):
            m = re.search(r"PASS=(\d) REASON=(\S+)", line)
            cur["pass"], cur["reason"] = int(m.group(1)), m.group(2)
            rows.append(cur)
            cur = None
    return rows


def fit_bars(rows, margin=0.005):
    """The repository's convention: each bar is the worst default-seed value
    plus about 0.5 %, rounded up at the fourth significant digit."""
    def up(x):
        e = 10 ** (np.floor(np.log10(abs(x))) - 3)
        return float(np.ceil(x * (1 + margin) / e) * e)
    js = [r for r in rows if r["record"].startswith("jonswap")]
    pm = [r for r in rows if r["record"].startswith("pmstokes")]
    w = {}
    for name, sub, key in (("err_limit_percent_z_jonswap", js, "z_pct_hs"), ("err_limit_percent_z_pmstokes", pm, "z_pct_hs"),
                           ("err_limit_yaw_deg", rows, "yaw_deg"), ("err_limit_roll_deg", rows, "roll_deg"),
                           ("err_limit_pitch_deg", rows, "pitch_deg"), ("err_limit_percent_3d_jonswap", js, "pct_3d"),
                           ("err_limit_percent_3d_pmstokes", pm, "pct_3d"), ("acc_z_bias_percent", rows, "acc_z_bias_pct"),
                           ("bias_3d_percent", rows, "acc_bias_3d_pct"), ("gyro_bias_3d_percent", rows, "gyro_bias_3d_pct")):
        worst = max(sub, key=lambda r: r[key])
        w[name] = {"worst": worst[key], "record": worst["record"], "bar": up(worst[key])}
    return w


def stage_gates(args):
    """Default-seed replays of the eight scored records: CSM and TFG-CARMA3
    (T41, Legacy and Evidence excitation); the gate quantities and CSM's bars."""
    rows = []
    for label, (exe, extra) in GATE_RUNS.items():
        if not exe.exists():
            sys.exit(f"missing {exe.relative_to(ROOT)}: make -C {exe.parent.relative_to(ROOT)} build")
        if args.logs:
            text = (Path(args.logs) / f"replay_{label.replace('@', '')}.log").read_text()
        else:
            env = {k: v for k, v in os.environ.items() if not k.startswith("W3D_")}
            env.update(extra, W3D_COLLECT_ALL_GATES="1", W3D_WRITE_TIMESERIES="0")
            text = subprocess.run([str(exe)], cwd=exe.parent, capture_output=True, text=True, env=env).stdout
        got = parse_gates(text)
        if len(got) != 8:
            sys.exit(f"gates: {label}: {len(got)} scored records parsed (expected 8)")
        for r in got:
            r["filter"] = label
        rows += got
        print(f"gates {label}: {sum(r['pass'] for r in got)}/8 records within the simulator's bars", flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    with GATES.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["filter", "record"] + GATE_KEYS + ["pass", "reason"])
        w.writeheader()
        w.writerows(rows)
    bars = fit_bars([r for r in rows if r["filter"] == "CSM"])
    BARS.write_text(json.dumps(bars, indent=1))
    for k, v in bars.items():
        print(f"  {k:32s} worst {v['worst']:.6g} ({v['record']}) -> bar {v['bar']:.6g}")


def stage_dataset(args):
    require(SIM)
    rows = []
    for variant in DATASET_VARIANTS:
        env = {k: v for k, v in os.environ.items() if not k.startswith("W3D_")}
        env.update(W3D_CSM_VARIANT=variant, W3D_WRITE_TIMESERIES="0")
        p = subprocess.run([str(SIM)], cwd=SIM.parent, capture_output=True, text=True, env=env)
        for line in p.stdout.splitlines():
            if not line.startswith("VALIDATION_METRICS family="):
                continue
            kv = dict(tok.split("=", 1) for tok in line.split()[1:] if "=" in tok)
            row = {"variant": variant, "record": kv.get("input", "").replace("wave_data_", "").replace(".csv", "")}
            for k in ("disp_x_rms_m", "disp_y_rms_m", "disp_z_rms_m", "disp_3d_rms_m", "disp_z_pct_hs",
                      "disp_3d_pct_refmax", "roll_rms_deg", "pitch_rms_deg", "yaw_rms_deg", "accel_3d_rms_mps2",
                      "accel_bias_3d_rms_mps2", "gyro_bias_3d_rms_radps"):
                try:
                    row[k] = float(kv[k])
                except (KeyError, ValueError):
                    row[k] = float("nan")
            rows.append(row)
        print(f"dataset {variant}: rc={p.returncode}", flush=True)
    if rows:
        OUT.mkdir(parents=True, exist_ok=True)
        with DATASET.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)


CORE_TIM = OUT / "timing_core.csv"
TIMING_EXE = TDIR / "kalman_carma_shared_modes-timing"


def stage_timing(args):
    """Host timing, one job on an otherwise idle machine: the complete wrapper
    (bench runtime column, three seeds) and the bare QMEKF cores with and
    without x86 flush-to-zero (kalman_carma_shared_modes-timing)."""
    if TIM.exists():
        TIM.unlink()
    run_rows([cmd(J29, s, f) for s in range(1, 4) for f in FILTERS], TIM, 1)
    if not TIMING_EXE.exists():
        subprocess.run(["make", "-C", str(TDIR), TIMING_EXE.name], check=True)
    rows = []
    for extra in ([], ["--ftz"], ["--no-commits"], ["--no-commits", "--ftz"]):
        p = subprocess.run([str(TIMING_EXE)] + extra, capture_output=True, text=True, check=True)
        rows += list(csv.DictReader(p.stdout.splitlines()))
    with CORE_TIM.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"timing: {TIM.relative_to(ROOT)}, {CORE_TIM.relative_to(ROOT)}")


def stage_size(args):
    tc = args.toolchain or os.environ.get("XTENSA_TOOLCHAIN_BIN")
    if not tc:
        sys.exit("size: --toolchain <Espressif xtensa bin dir> (esp-14.2.0_20241119)")
    probe = TDIR / "arduino" / "csm_size_probe" / "build.sh"
    p = subprocess.run(["bash", str(probe), tc], capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit(f"size probe failed: {p.stderr[-800:]}")
    OUT.mkdir(parents=True, exist_ok=True)
    SIZE.write_text(p.stdout)
    print(p.stdout)


# ---------------------------------------------------------------------------
# cost model (analytic multiply-adds; carma3_analysis.cost_model, extended)
# ---------------------------------------------------------------------------

def cost_rows():
    """Per-sample multiply-adds of predict + accelerometer update (the CARMA3
    core's structured implementation, carma3_analysis.cost_model), the jitter
    series of the motion discretization, and the per-commit work amortized
    over the 0.1 s model cadence at 200 Hz (20 samples)."""
    out = []
    for name, nx, blocks, kind in (("A41 (21)", 21, [4, 4, 4], "two4"), ("CSM (27)", 27, [6, 6, 6], "pair6")):
        c = carma3_analysis.cost_model(nx, blocks, 9, "carma3")
        nl = sum(blocks)
        q_add = (nl * (nl + 1)) // 2 * (2 if kind == "pair6" else 1)
        if kind == "two4":            # horizontal and vertical model, 4 states each
            series = sum(3 * m * m + 2 * 3 * m * (m + 1) // 2 for m in (4, 4))
            disc = 2 * ((3 + 7) * 64 + (2 * 3 + 7) * 64 * 2)
            recoord = 0
        else:                          # one 6-state model, symmetric + antisymmetric Q series
            m, kf = 6, 5
            series = kf * m * m + 2 * 2 * kf * m * (m + 1) // 2
            disc = (m - 1 + 7) * m ** 3 + (2 * (m - 1) + 7) * m ** 3 * 3
            # The core applies every axis's 6 x 6 transform densely to the rows and
            # the columns (Kalman3D_Wave_CARMA3::recoordinate_motion_), the
            # diagonal colouring rescale included.
            recoord = 3 * 2 * m * m * nx
        per_commit = disc + recoord
        out.append({"filter": name, "states": nx, "P_floats": nx * nx, "P_bytes": 4 * nx * nx,
                    "predict_madds": c["predict_madds"] + q_add, "acc_update_madds": c["acc_update_madds"],
                    "per_sample_madds": c["predict_madds"] + q_add + c["acc_update_madds"],
                    "jitter_series_madds": series, "per_commit_madds": per_commit,
                    "commit_amortized_madds": per_commit // 20,
                    "total_madds_per_sample": c["predict_madds"] + q_add + c["acc_update_madds"] + per_commit // 20})
    return out


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------

GROUP_ORDER = ["fixed", "coupling", "gradual", "rapid", "disturbance"]


def _frame(path):
    import pandas as pd
    if not path.exists():
        return pd.DataFrame()
    d = pd.read_csv(path)
    d = d.copy()
    d["rms3"] = np.sqrt(d.rms_x ** 2 + d.rms_y ** 2 + d.rms_z ** 2)
    return d


def compute_stats():
    """Every number the article quotes, from the committed CSVs."""
    df = _frame(RAW)
    st = {"n_rows": {}, "scenarios": {}, "groups": {}, "pooled": {}, "robust": {}, "dataset": {}, "cost": cost_rows()}
    if not df.empty:
        st["n_rows"]["run"] = int(len(df))
        for sc, (grp, desc) in SCENARIOS.items():
            sub = df[df.scenario == sc]
            if sub.empty:
                continue
            r = {"group": grp, "desc": desc, "seeds": int(sub.seed.nunique()), "filters": {}}
            for f in FILTERS:
                g = sub[sub["filter"] == f]
                if g.empty:
                    continue
                e = {"failed": int(g.failed.sum())}
                for c in ("rms_x", "rms_y", "rms_z", "rms_h", "rms3", "rmsv_z", "rmsv_h", "peak_z", "nees_pos3", "nees_pv6",
                          "nis", "cov95_x", "cov95_y", "cov95_z", "roll_rms_deg", "pitch_rms_deg", "yaw_rms_deg",
                          "bg_rms_dph", "ba_rms_mg", "cond_max", "runtime_ns"):
                    if c in g:
                        e[c] = float(g[c].astype(float).mean())
                for c in ("nees_pos3", "nees_pv6"):
                    e[c + "_median"] = float(g[c].astype(float).median())
                if f != BASE:
                    for c in ("rms_h", "rms_z", "rms3", "nees_pos3"):
                        pr = paired(df, sc, f, c)
                        if pr:
                            e["vsA41_" + c] = pr
                r["filters"][f] = e
            st["scenarios"][sc] = r
        for grp in GROUP_ORDER:
            scs = [sc for sc, (g, _) in SCENARIOS.items() if g == grp and sc in st["scenarios"]]
            vals = {c: [st["scenarios"][sc]["filters"]["CSM"]["vsA41_" + c]["change"] for sc in scs
                        if "vsA41_" + c in st["scenarios"][sc]["filters"].get("CSM", {})] for c in ("rms_h", "rms_z", "rms3")}
            gr = {c: (float(np.mean(v)) if v else float("nan")) for c, v in vals.items()}
            gr.update({c + "_min": (float(np.min(v)) if v else float("nan")) for c, v in vals.items()})
            gr.update({c + "_max": (float(np.max(v)) if v else float("nan")) for c, v in vals.items()})
            st["groups"][grp] = {"scenarios": scs, "CSM": gr}
        for f in FILTERS:
            g = df[df["filter"] == f]
            st["pooled"][f] = {c: float(g[c].astype(float).median()) for c in ("nees_pos3", "nees_pv6", "nis")}
            st["pooled"][f].update({c + "_mean": float(g[c].astype(float).mean()) for c in ("cov95_x", "cov95_y", "cov95_z")})
            st["pooled"][f]["failed"] = int(g.failed.sum())
            st["pooled"][f]["cond_max_median"] = float(g.cond_max.astype(float).median())
            st["pooled"][f]["runtime_ns_median"] = float(g.runtime_ns.astype(float).median())
        # Seas where CSM is better / worse than A41 (paired change of the seed means).
        for c in ("rms_z", "rms_h", "rms3"):
            ch = {sc: r["filters"]["CSM"]["vsA41_" + c]["change"] for sc, r in st["scenarios"].items()
                  if "vsA41_" + c in r["filters"].get("CSM", {})}
            st["pooled"]["CSM"]["better_" + c] = sum(1 for v in ch.values() if v > 0)
            st["pooled"]["CSM"]["worst_" + c] = min(ch.items(), key=lambda kv: kv[1]) if ch else None
            st["pooled"]["CSM"]["best_" + c] = max(ch.items(), key=lambda kv: kv[1]) if ch else None
    rob = _frame(ROB)
    if not rob.empty:
        st["n_rows"]["robust"] = int(len(rob))
        for sc in sorted(rob.scenario.unique()):
            sub = rob[rob.scenario == sc]
            st["robust"][sc] = {}
            for f in sorted(sub["filter"].unique()):
                g = sub[sub["filter"] == f]
                st["robust"][sc][f] = {c: float(g[c].astype(float).mean()) for c in ("rms_h", "rms_z", "rms3", "nees_pos3")}
                st["robust"][sc][f]["nees_pos3_median"] = float(g.nees_pos3.astype(float).median())
                st["robust"][sc][f]["n"] = int(len(g))
    if DATASET.exists():
        import pandas as pd
        d = pd.read_csv(DATASET)
        for rec in sorted(d.record.unique()):
            st["dataset"][rec] = {r.variant: {"z_pct_hs": float(r.disp_z_pct_hs), "pct_3d": float(r.disp_3d_pct_refmax),
                                              "yaw_deg": float(r.yaw_rms_deg), "ba": float(r.accel_bias_3d_rms_mps2),
                                              "bg": float(r.gyro_bias_3d_rms_radps)}
                                  for r in d[d.record == rec].itertuples()}
    if GATES.exists():
        g = read_rows(GATES)
        st["gates"] = {}
        for f in GATE_RUNS:
            rows = [r for r in g if r["filter"] == f]
            st["gates"][f] = {"records": {r["record"]: {k: float(r[k]) for k in GATE_KEYS} for r in rows},
                              "passed": sum(int(r["pass"]) for r in rows)}
        if BARS.exists():
            st["csm_bars"] = json.loads(BARS.read_text())
    if SIZE.exists():
        st["esp32s3"] = {r["variant"]: r for r in read_rows(SIZE)}
    if TIM.exists():
        t = _frame(TIM)
        st["host_timing"] = {f: {"mean_ns": float(g.runtime_ns.mean()), "p99_ns": float(g.runtime_p99_ns.mean())}
                             for f, g in t.groupby("filter")}
    if CORE_TIM.exists():
        import pandas as pd
        t = pd.read_csv(CORE_TIM)
        st["core_timing"] = {f"{r.core}|commits{r.commits}|ftz{r.ftz}": {
            "total_us": float(g.total_us.median()), "time_update_us": float(g.time_update_us.median()),
            "acc_mag_us": float(g.acc_mag_us.median()), "states": int(r.states)}
            for r, g in ((g.iloc[0], g) for _, g in t.groupby(["core", "commits", "ftz"]))}
    an = OUT / "analysis.json"
    if an.exists():
        st["analysis"] = json.loads(an.read_text())
    return st


def _pct(r, key="change"):
    return f"{100 * r[key]:+.2f}" if r else "--"


def _ci(r):
    return f"{100 * r['change']:+.2f} [{100 * r['ci'][0]:+.2f}, {100 * r['ci'][1]:+.2f}]" if r else "--"


def short(sc):
    return (sc.replace("jonswap_fp0.29_hs0.8_g3.3", "J29").replace("jonswap_fp0.45_hs0.4_g3.3", "wind sea")
            .replace("_", " ").replace("+", " +"))


def write_tex(st):
    """LaTeX parts (doc/kalman_carma_shared_modes/csm-*.tex-part)."""
    parts = {}
    # per-scenario table: RMS of A41 and CSM, paired changes vs A41
    lines = [r"\begin{tabular}{@{}l r r r r r@{}}", r"\toprule",
             r"sea & A41 3-D (mm) & CSM 3-D (mm) & \multicolumn{3}{c}{CSM change vs A41 (\%, + better)} \\",
             r"\cmidrule(l){4-6} & & & heave & horizontal & 3-D \\", r"\midrule"]
    last = None
    for sc, r in st["scenarios"].items():
        if r["group"] != last:
            if last is not None:
                lines.append(r"\addlinespace")
            last = r["group"]
        fl = r["filters"]
        if "CSM" not in fl or BASE not in fl:
            continue
        row = [tex(short(sc)), f"{1000 * fl[BASE]['rms3']:.1f}", f"{1000 * fl['CSM']['rms3']:.1f}"]
        for c in ("rms_z", "rms_h", "rms3"):
            row.append(_ci(fl["CSM"].get("vsA41_" + c)))
        lines.append(" & ".join(row) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    parts["csm-scenarios.tex-part"] = "\n".join(lines) + "\n"
    # groups
    lines = [r"\begin{tabular}{@{}l r r r@{}}", r"\toprule",
             r"group (seas) & heave & horizontal & 3-D \\", r"\midrule"]
    for grp in GROUP_ORDER:
        g = st["groups"].get(grp)
        if not g or not g["scenarios"]:
            continue
        c = g["CSM"]
        lines.append(" & ".join([f"{grp} ({len(g['scenarios'])})"] +
                                [f"{100 * c[k]:+.2f} ({100 * c[k + '_min']:+.1f} .. {100 * c[k + '_max']:+.1f})"
                                 for k in ("rms_z", "rms_h", "rms3")]) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    parts["csm-groups.tex-part"] = "\n".join(lines) + "\n"
    # consistency
    lines = [r"\begin{tabular}{@{}l r r r r r@{}}", r"\toprule",
             r"filter & NEES$_p$ & NEES$_{pv}$ & NIS & 95\,\% cover $z$ & failures \\", r"\midrule"]
    for f in FILTERS:
        p = st["pooled"].get(f)
        if p:
            lines.append(" & ".join([tex(f), f"{p['nees_pos3']:.3f}", f"{p['nees_pv6']:.3f}", f"{p['nis']:.3f}",
                                     f"{100 * p['cov95_z_mean']:.1f}", str(p["failed"])]) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    parts["csm-consistency.tex-part"] = "\n".join(lines) + "\n"
    # robustness
    if st["robust"]:
        lines = [r"\begin{tabular}{@{}l l r r r r@{}}", r"\toprule",
                 r"sea & prior & CSM 3-D (mm) & no conf.\ 3-D & CSM NEES$_p$ & no conf.\ NEES$_p$ \\", r"\midrule"]
        for sc, fl in st["robust"].items():
            if "turns" in sc:
                continue
            a41 = fl.get(BASE)
            if a41:
                lines.append(" & ".join([tex(short(sc)), "A41 (no coupling)", f"{1000 * a41['rms3']:.1f}", "",
                                         f"{a41['nees_pos3_median']:.2f}", ""]) + r" \\")
            for w in ["true"] + ROB_WRONG:
                a, b = fl.get(f"CSM/ext-{w}"), fl.get(f"CSM:noconf/ext-{w}")
                if a:
                    lines.append(" & ".join(["", tex(w), f"{1000 * a['rms3']:.1f}",
                                             f"{1000 * b['rms3']:.1f}" if b else "--",
                                             f"{a['nees_pos3_median']:.2f}",
                                             f"{b['nees_pos3_median']:.2f}" if b else "--"]) + r" \\")
            lines.append(r"\addlinespace")
        lines += [r"\bottomrule", r"\end{tabular}"]
        parts["csm-robust.tex-part"] = "\n".join(lines) + "\n"
    # dataset
    if st["dataset"]:
        vs = DATASET_VARIANTS
        lines = [r"\begin{tabular}{@{}l " + "r " * len(vs) + r"@{}}", r"\toprule",
                 "record & " + " & ".join(tex(v) for v in vs) + r" \\", r"\midrule"]
        for rec, fl in st["dataset"].items():
            lines.append(tex(rec.replace("_L", " L").split(" L")[0]) + " & "
                         + " & ".join(f"{fl[v]['z_pct_hs']:.2f} / {fl[v]['pct_3d']:.2f}" if v in fl else "--" for v in vs) + r" \\")
        lines += [r"\bottomrule", r"\end{tabular}"]
        parts["csm-dataset.tex-part"] = "\n".join(lines) + "\n"
    # dataset regression gates: worst of the eight scored records
    if st.get("gates"):
        gt = st["gates"]
        cats = [("heave \\% $H_s$, JONSWAP", "z_pct_hs", "jonswap", "err_limit_percent_z_jonswap"),
                ("heave \\% $H_s$, PM-Stokes", "z_pct_hs", "pmstokes", "err_limit_percent_z_pmstokes"),
                ("3-D \\%, JONSWAP", "pct_3d", "jonswap", "err_limit_percent_3d_jonswap"),
                ("3-D \\%, PM-Stokes", "pct_3d", "pmstokes", "err_limit_percent_3d_pmstokes"),
                ("roll (deg)", "roll_deg", "", "err_limit_roll_deg"), ("pitch (deg)", "pitch_deg", "", "err_limit_pitch_deg"),
                ("yaw (deg)", "yaw_deg", "", "err_limit_yaw_deg"), ("acc.\\ $z$ bias \\%", "acc_z_bias_pct", "", "acc_z_bias_percent"),
                ("acc.\\ 3-D bias \\%", "acc_bias_3d_pct", "", "bias_3d_percent"),
                ("gyro 3-D bias \\%", "gyro_bias_3d_pct", "", "gyro_bias_3d_percent")]
        lines = [r"\begin{tabular}{@{}l r r r r r@{}}", r"\toprule",
                 r"gate (worst of 8) & CSM & bar & TFG-CARMA3 & TFG, Evidence & CSM/TFG \\", r"\midrule"]

        def worst(f, k, fam):
            return max(v[k] for rec, v in gt[f]["records"].items() if rec.startswith(fam))
        for name, k, fam, bar in cats:
            wc, wt = worst("CSM", k, fam), worst("T41", k, fam)
            lines.append(" & ".join([name, f"{wc:.4g}", f"{st.get('csm_bars', {}).get(bar, {}).get('bar', float('nan')):.4g}",
                                     f"{wt:.4g}", f"{worst('T41@C', k, fam):.4g}",
                                     f"{100 * (wc / wt - 1):+.1f}\\,\\%"]) + r" \\")
        lines += [r"\bottomrule", r"\end{tabular}"]
        parts["csm-gates.tex-part"] = "\n".join(lines) + "\n"
    # cost
    lines = [r"\begin{tabular}{@{}l r r r r r@{}}", r"\toprule",
             r"filter & $n_x$ & $P$ (B) & madds/sample & madds/commit & total/sample \\", r"\midrule"]
    for r in st["cost"]:
        lines.append(" & ".join([tex(r["filter"]), str(r["states"]), str(r["P_bytes"]), f"{r['per_sample_madds']:,}",
                                 f"{r['per_commit_madds']:,}", f"{r['total_madds_per_sample']:,}"]).replace(",", r"\,")
                     + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    parts["csm-cost.tex-part"] = "\n".join(lines) + "\n"
    if st.get("esp32s3"):
        lines = [r"\begin{tabular}{@{}l r r r r r@{}}", r"\toprule",
                 r"variant & filter object (B) & code (B) & rodata (B) & largest frame (B) & update chain (B) \\",
                 r"\midrule"]
        for v, r in st["esp32s3"].items():
            lines.append(" & ".join([tex(v), r["sizeof_qmekf"], f"{int(r['code_bytes']):,}".replace(",", r"\,"),
                                     f"{int(r['rodata_bytes']):,}".replace(",", r"\,"), r["max_stack_frame"],
                                     r["update_chain_stack"]]) + r" \\")
        lines += [r"\bottomrule", r"\end{tabular}"]
        parts["csm-esp32s3.tex-part"] = "\n".join(lines) + "\n"
    DOC_DIR.mkdir(parents=True, exist_ok=True)
    for name, body in parts.items():
        (DOC_DIR / name).write_text("% generated by tools/carma_shared_modes_study.py report; do not edit\n" + body)
    return parts


def figures(st, outdir, exts):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd
    if "pgf" in exts:
        plt.rcParams.update({"pgf.texsystem": "lualatex", "pgf.rcfonts": False, "font.family": "serif"})
    plt.rcParams.update({"font.size": 8, "axes.titlesize": 8, "legend.fontsize": 7})
    outdir.mkdir(parents=True, exist_ok=True)

    def save(fig, name):
        for e in exts:
            fig.savefig(outdir / f"{name}.{e}", bbox_inches="tight")
        plt.close(fig)

    # 1. paired changes vs A41, per scenario
    scs = [sc for sc in st["scenarios"] if "CSM" in st["scenarios"][sc]["filters"]]
    if scs:
        fig, axes = plt.subplots(1, 3, figsize=(7.2, 5.6), sharey=True)
        y = np.arange(len(scs))
        for ax, col, title in ((axes[0], "rms_z", "heave"), (axes[1], "rms_h", "horizontal"), (axes[2], "rms3", "3-D")):
            ch, lo, hi = [], [], []
            for sc in scs:
                r = st["scenarios"][sc]["filters"]["CSM"].get("vsA41_" + col)
                ch.append(100 * r["change"] if r else np.nan)
                lo.append(100 * (r["change"] - r["ci"][0]) if r else 0)
                hi.append(100 * (r["ci"][1] - r["change"]) if r else 0)
            ax.errorbar(ch, y, xerr=[lo, hi], fmt="o", ms=3, lw=0.8, color=COLORS["CSM"])
            ax.axvline(0, color="k", lw=0.6)
            ax.set_xlabel(f"{title} RMS: CSM vs A41\n(%, + better)")
            ax.grid(alpha=0.3, axis="x")
        axes[0].set_yticks(y)
        axes[0].set_yticklabels([short(s) for s in scs])
        axes[0].invert_yaxis()
        fig.tight_layout()
        save(fig, "csm_paired_changes")
    # 2. absolute 3-D RMS (log bar chart)
    if scs:
        fig, ax = plt.subplots(figsize=(7.2, 2.8))
        x = np.arange(len(scs))
        w = 0.38
        for k, f in enumerate(FILTERS):
            vals = [1000 * st["scenarios"][sc]["filters"].get(f, {}).get("rms3", np.nan) for sc in scs]
            ax.bar(x + (k - 0.5) * w, vals, w, color=COLORS[f], label=LABEL[f])
        ax.set_yscale("log")
        ax.set_ylabel("3-D position RMS (mm)")
        ax.set_xticks(x)
        ax.set_xticklabels([short(s) for s in scs], rotation=60, ha="right", fontsize=6)
        h, lab = ax.get_legend_handles_labels()
        fig.legend(h, lab, ncol=2, frameon=False, fontsize=6, loc="upper center", bbox_to_anchor=(0.5, 1.04))
        ax.grid(alpha=0.3, axis="y")
        fig.tight_layout(rect=(0, 0, 1, 0.95))
        save(fig, "csm_rms3")
    # 3. transient envelopes (from 150 s: the proxy start-up precedes Live)
    changes = {"changing_fp": [(300.0, 900.0)], "abrupt_change": [(595.0, 605.0)], "dir_change": [(570.0, 630.0)],
               f"{J29}+on400w60+off900w60": [(370.0, 430.0), (870.0, 930.0)]}
    envs = [(sc, TRACE_DIR / f"{re.sub(r'[^A-Za-z0-9_.+-]', '_', sc)}_envelope.csv") for sc in TRACE_SEAS]
    envs = [(sc, p) for sc, p in envs if p.exists()]
    if envs:
        fig, axes = plt.subplots(len(envs), 2, figsize=(7.2, 1.9 * len(envs)), squeeze=False)
        handles = []
        for i, (sc, p) in enumerate(envs):
            d = pd.read_csv(p)
            d = d[d.t.values >= 150.0]
            t = d.t.values
            k = max(1, int(round(10.0 / float(np.median(np.diff(t))))))   # 10 s moving RMS of the 1 s envelope

            def smooth(v, k=k):
                return np.sqrt(np.convolve(np.nan_to_num(v) ** 2, np.ones(k) / k, mode="same"))
            for f in TRACE_FILTERS.split(","):
                ln, = axes[i, 0].plot(t, 1000 * smooth(d[f"{f}|z"].values), color=COLORS[f], lw=0.8, label=LABEL[f])
                if i == 0:
                    handles.append(ln)
                hz = np.sqrt(d[f"{f}|x"].values ** 2 + d[f"{f}|y"].values ** 2)
                axes[i, 1].plot(t, 1000 * smooth(hz), color=COLORS[f], lw=0.8)
            axes[i, 0].set_ylabel(f"{short(sc)}\nheave (mm)", fontsize=7)
            axes[i, 1].set_ylabel("horizontal (mm)", fontsize=7)
            for a in axes[i]:
                for t0, t1 in changes.get(sc, []):
                    a.axvspan(t0, t1, color="#bbbbbb", alpha=0.3, lw=0)
                a.set_yscale("log")
                a.set_ylim(bottom=1e-1)
                a.grid(alpha=0.3)
                a.set_xlim(t[0], t[-1])
        axes[-1, 0].set_xlabel("time (s)")
        axes[-1, 1].set_xlabel("time (s)")
        fig.legend(handles, [h.get_label() for h in handles], loc="upper center", ncol=2, frameon=False, fontsize=6,
                   bbox_to_anchor=(0.5, 1.03))
        fig.tight_layout(rect=(0, 0, 1, 0.96))
        save(fig, "csm_transients")
    # 4. consistency
    df = _frame(RAW)
    if not df.empty:
        fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.4))
        data = [df[df["filter"] == f].nees_pos3.astype(float).clip(upper=50).values for f in FILTERS]
        axes[0].boxplot(data, showfliers=False)
        axes[0].axhline(3, color="k", lw=0.6, ls=":")
        axes[0].set_xticks(range(1, len(FILTERS) + 1))
        axes[0].set_xticklabels(FILTERS, fontsize=6)
        axes[0].set_ylabel(r"position NEES (3 dof)")
        a = df[df["filter"] == BASE].set_index(["scenario", "seed"]).nees_pos3
        b = df[df["filter"] == "CSM"].set_index(["scenario", "seed"]).nees_pos3
        idx = a.index.intersection(b.index)
        axes[1].scatter(a.loc[idx], b.loc[idx], s=3, color=COLORS["CSM"], alpha=0.6)
        lim = [0.3, 30]
        axes[1].plot(lim, lim, color="k", lw=0.6)
        axes[1].set_xscale("log")
        axes[1].set_yscale("log")
        axes[1].set_xlim(lim)
        axes[1].set_ylim(lim)
        axes[1].set_xlabel("A41 NEES$_p$")
        axes[1].set_ylabel("CSM NEES$_p$")
        fig.tight_layout()
        save(fig, "csm_consistency")
    # 5. robustness
    seas = [s for s in st["robust"] if "turns" not in s]
    if seas:
        fig, axes = plt.subplots(1, len(seas), figsize=(7.2, 2.4), squeeze=False)
        for i, sc in enumerate(seas):
            fl = st["robust"][sc]
            labels = ["true"] + ROB_WRONG
            x = np.arange(len(labels))
            for k, (conf, col) in enumerate((("", COLORS["CSM"]), (":noconf", COLORS["CSM:noconf"]))):
                vals = [1000 * fl.get(f"CSM{conf}/ext-{w}", {}).get("rms3", np.nan) for w in labels]
                axes[0, i].bar(x + (k - 0.5) * 0.38, vals, 0.38, color=col, label=f"CSM{conf}")
            a41 = 1000 * fl.get(BASE, {}).get("rms3", np.nan)
            axes[0, i].axhline(a41, color=COLORS[BASE], lw=0.8, ls="--", label="A41")
            axes[0, i].set_xticks(x)
            axes[0, i].set_xticklabels(labels, fontsize=6)
            vals = [1000 * v.get("rms3", np.nan) for k, v in fl.items() if "ext-" in k] + [a41]
            axes[0, i].set_ylim(0.85 * np.nanmin(vals), 1.05 * np.nanmax(vals))
            axes[0, i].set_title(short(sc), fontsize=7)
            axes[0, i].grid(alpha=0.3, axis="y")
        axes[0, 0].set_ylabel("3-D RMS (mm)")
        h, lab = axes[0, 0].get_legend_handles_labels()
        fig.legend(h, lab, fontsize=6, frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.06))
        fig.tight_layout(rect=(0, 0, 1, 0.93))
        save(fig, "csm_robustness")
    # 6. cost
    cost = st["cost"]
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.2))
    names = [c["filter"].split(" ")[0] for c in cost]
    cols = [COLORS["A41@C"], COLORS["CSM"]]
    axes[0].bar(names, [c["total_madds_per_sample"] for c in cost], color=cols)
    axes[0].set_ylabel("multiply-adds / sample")
    axes[1].bar(names, [c["P_bytes"] for c in cost], color=cols)
    axes[1].set_ylabel("covariance P (bytes)")
    if st.get("esp32s3"):
        e = st["esp32s3"]
        vv = [v for v in ("A41", "CSM") if v in e]
        axes[2].bar(vv, [int(e[v]["sizeof_qmekf"]) for v in vv], color=cols[:len(vv)])
        axes[2].set_ylabel("ESP32-S3 filter object (B)")
    for a in axes:
        a.tick_params(axis="x", labelsize=6)
        a.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    save(fig, "csm_cost")


def stage_report(args):
    st = compute_stats()
    (OUT / "summary.json").write_text(json.dumps(st, indent=1, default=float))
    write_tex(st)
    man = {"git_revision": git_rev(), "scenarios": {k: v[1] for k, v in SCENARIOS.items()},
           "filters": FILTERS, "robust_filters": ROB_FILTERS, "turn_filters": ROB_TURN_FILTERS,
           "trace_filters": TRACE_FILTERS, "dataset_variants": DATASET_VARIANTS,
           "seeds": {"run": args.seeds, "robust": args.robust_seeds, "traces": args.trace_seeds},
           "n_rows": st["n_rows"], "duration_s": "per scenario (default 1200 s), scored from 300 s",
           "result_files": sorted(p.name for p in OUT.glob("*.csv"))}
    (OUT / "manifest.json").write_text(json.dumps(man, indent=1))
    figures(st, FIG_DIR, ["svg"])
    print(f"report: {OUT.relative_to(ROOT)} (summary.json, manifest.json, figures/), tex parts in "
          f"{DOC_DIR.relative_to(ROOT)}")


def stage_figures(args):
    st = json.loads((OUT / "summary.json").read_text()) if (OUT / "summary.json").exists() else compute_stats()
    figures(st, Path(args.outdir) if args.outdir else FIG_DIR, args.ext.split(","))


def stage_cost(args):
    rows = cost_rows()
    OUT.mkdir(parents=True, exist_ok=True)
    with COST.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(json.dumps(rows, indent=1))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("stage", choices=["smoke", "run", "robust", "traces", "dataset", "gates", "timing", "size",
                                      "cost", "report", "figures", "all"])
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--seeds", type=int, default=30)
    ap.add_argument("--robust-seeds", type=int, default=20)
    ap.add_argument("--trace-seeds", type=int, default=10)
    ap.add_argument("--toolchain", default=None, help="size stage: Espressif xtensa toolchain bin directory")
    ap.add_argument("--outdir", default=None, help="figures stage: output directory")
    ap.add_argument("--logs", default=None, help="gates stage: parse replay_<label>.log files from this directory")
    ap.add_argument("--ext", default="svg", help="figures stage: comma-separated formats (svg, pgf, png, pdf)")
    args = ap.parse_args(argv)
    OUT.mkdir(parents=True, exist_ok=True)
    stages = {"smoke": stage_smoke, "run": stage_run, "robust": stage_robust, "traces": stage_traces,
              "dataset": stage_dataset, "gates": stage_gates, "timing": stage_timing, "size": stage_size,
              "cost": stage_cost, "report": stage_report, "figures": stage_figures}
    if args.stage == "all":
        for s in ("cost", "smoke", "run", "robust", "traces", "dataset", "gates", "timing", "report"):
            stages[s](args)
    else:
        stages[args.stage](args)


if __name__ == "__main__":
    main()
