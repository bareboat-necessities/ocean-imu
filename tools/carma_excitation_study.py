#!/usr/bin/env python3
"""Wave-excitation study of the CARMA research filters (docs/carma-wave-excitation.md).

Paired runs of tests/kalman_carma/carma3_bench: every filter of a row sees the
same truth, sensor stream and seed.  Filters are the CARMA wrappers in their
excitation modes (A Legacy, B0 floor removed, B debiased, C evidence and its
ablations Cnf / Cng), on the quaternion MEKF (A41...) and the two-frame Lie
group (T41...), with the shipping OU-III as reference.

    tune       tuning seeds 101-110: evidence-threshold and noise-model
               variants of C on the tuning conditions
    run        validation seeds 1-30 (10 or 20 on secondary grids): experiments
               E1-E8 (--exp to select); resumable
    traces     seed-1 time series kept for the figures
    report     summary.md, paired_stats.json and the figures

Run: python3 tools/carma_excitation_study.py {tune,run,traces,report} [--jobs 4]
"""
from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BENCH = ROOT / "tests" / "kalman_carma" / "carma3_bench"
OUT = ROOT / "reports" / "results" / "carma_excitation"
RAW = OUT / "raw_runs.csv"
TUNE_RAW = OUT / "tuning_runs.csv"
TRANS = OUT / "transitions.csv"
PSD = OUT / "psd_calm.csv"
TRACE_DIR = OUT / "traces"

CORE = ["OU3", "A41@A", "A41@B0", "A41@B", "A41@C", "T41@A", "T41@C"]
ABL = ["A41@Cnf", "A41@Cng", "T41@B"]
FULL = CORE + ABL
NOISE = ["A41@B/nv0.5", "A41@B/nv2", "A41@C/nv0.5", "A41@C/nv2", "T41@C/nv0.5", "T41@C/nv2"]
ORDERS = ["A30@A", "A30@C", "A52@A", "A52@C", "T52@A", "T52@C"]

J = "jonswap_fp0.29_hs{h}_g3.3"
HS = ["0.005", "0.01", "0.02", "0.05", "0.1", "0.2"]
SWELL_F = ["0.125", "0.1", "0.08", "0.0625", "0.05", "0.04"]   # 8 .. 25 s
SWELL_H = ["0.05", "0.2", "0.5"]

# (experiment, scenario, n_seeds, filters, needs time series, extended)
def experiments():
    g = []
    for sc in ["stationary+T2400", "calm+T2400", "stationary+scal+T2400", "stationary+calbias+T2400"]:
        g.append(("E1", sc, 30, FULL, "psd", True))
    for h in HS:
        g.append(("E2", J.format(h=h), 30, FULL, None, False))
        g.append(("E2", f"jonswap_fp0.15_hs{h}_g3.3", 10, CORE, None, False))
        g.append(("E2", f"jonswap_fp0.45_hs{h}_g3.3", 10, CORE, None, False))
        g.append(("E2", f"pm_fp0.25_hs{h}", 10, CORE, None, False))
    for f in SWELL_F:
        for h in SWELL_H:
            g.append(("E3", f"jonswap_fp{f}_hs{h}_g5", 10, CORE, None, False))
    for w in ("10", "120"):
        g.append(("E4", f"{J.format(h='0.05')}+on400w{w}", 30, FULL, "trans", False))
        g.append(("E4", f"{J.format(h='0.8')}+on400w{w}", 30, FULL, "trans", False))
        g.append(("E4", f"{J.format(h='0.8')}+off600w{w}", 30, FULL, "trans", False))
        g.append(("E4", f"jonswap_fp0.29_hs1_g3.3+lev0w0.001x0.1+lev600w{w}x1", 30, FULL, "trans", False))
        g.append(("E4", f"jonswap_fp0.29_hs1_g3.3+lev600w{w}x0.1", 30, FULL, "trans", False))
        rw = "60" if w == "120" else w
        g.append(("E4", f"{J.format(h='0.5')}+T1500+on300w{rw}+off600w{rw}+on900w{rw}+off1200w{rw}", 30, FULL,
                  "trans", False))
    for base in ["stationary", J.format(h="0.05"), J.format(h="0.8")]:
        for m in ["", "+cnoise", "+nvar", "+engine1500", "+engine2400", "+temp"]:
            g.append(("E5", base + m, 10, CORE + NOISE, None, False))
    for m in ["", "+cnoise"]:
        g.append(("E5", "jonswap_fp0.08_hs0.2_g5" + m, 10, CORE + NOISE, None, False))
    for sc in ["calm", "calm+U2.5", "calm+maneuver", "calm+drift", "calm+anchor", "calm+U2.5+turns",
               f"{J.format(h='0.8')}+hscale0.2", f"{J.format(h='0.2')}+hscale0.2"]:
        g.append(("E6", sc, 20, CORE, None, True))
    for sc in ["mixed_swell_wind", "bimodal", "broadband", "abrupt_change", "changing_fp", J.format(h="0.8"),
               "pm_fp0.25_hs1"]:
        g.append(("E7", sc, 20, CORE, None, False))
    for sc in ["stationary", J.format(h="0.8"), J.format(h="0.02"), "jonswap_fp0.08_hs0.2_g5", "mixed_swell_wind"]:
        g.append(("E8", sc, 10, ORDERS, None, False))
    return g


TUNE_FILTERS = ["A41@A", "A41@B", "A41@C", "A41@C/u0.25", "A41@C/u1", "A41@C/z1.5_3", "A41@C/z3_6",
                "A41@C/zf4_8", "A41@C/zf16_32", "A41@C/k1.75_2.25", "A41@C/k2.5_3",
                "A41@C/nv0.5", "A41@C/nv0.5/u1", "A41@C/nv2"]
TUNE_SCEN = ["stationary", "calm+maneuver", J.format(h="0.01"), J.format(h="0.02"), J.format(h="0.05"),
             J.format(h="0.2"), J.format(h="0.8"), "jonswap_fp0.08_hs0.2_g5", "jonswap_fp0.05_hs0.05_g5",
             "mixed_swell_wind", "bimodal", f"{J.format(h='0.05')}+on400w10", f"{J.format(h='0.8')}+off600w10"]


def require(path):
    if not path.exists():
        sys.exit(f"missing {path.relative_to(ROOT)}: make -C tests/kalman_carma carma3_bench")


def header(extended):
    a = [str(BENCH), "--header", "--excitation"] + (["--extended"] if extended else [])
    return subprocess.run(a, capture_output=True, text=True, check=True).stdout.strip()


def done_keys(path):
    if not path.exists():
        return set()
    with path.open() as f:
        return {(r["exp"], r["scenario"], r["seed"]) for r in csv.DictReader(f)}


# ---------------------------------------------------------------------------
# Time-series derived metrics
# ---------------------------------------------------------------------------

def envelope_steps(scenario):
    """Transitions of the envelope modifiers: (t_c, width, level_from, level_to)."""
    level = 1.0
    steps = []
    first = True
    for m in scenario.split("+")[1:]:
        if m.startswith("on") or m.startswith("off") or m.startswith("lev"):
            on, off = m.startswith("on"), m.startswith("off")
            body = m[2:] if on else m[3:]
            tc, rest = body.split("w", 1)
            if on or off:
                w, to = rest, (1.0 if on else 0.0)
            else:
                w, x = rest.split("x", 1)
                to = float(x)
            if first and on:
                level = 0.0
            first = False
            steps.append((float(tc), float(w), level, to))
            level = to
    return steps


def transition_metrics(ts, filters, scenario, seed, exp):
    """Detection delay, excitation rise/fall, heave transient of every step."""
    t = ts["t"]
    rows = []
    for k, (tc, w, lf, lt) in enumerate(envelope_steps(scenario)):
        if w < 0.01:
            continue   # the initial level, not a transition
        t_s = tc - 0.5 * w
        nxt = [s for s in envelope_steps(scenario) if s[0] - 0.5 * s[1] > t_s + 1e-6]
        t_end = (nxt[0][0] - 0.5 * nxt[0][1]) if nxt else t[-1]
        kind = "onset" if lt > lf else "decay"
        for f in filters:
            if f"pw_{f}" not in ts:
                continue
            pw, sig, z, zt = ts[f"pw_{f}"], ts[f"sig_{f}"], ts[f"z_{f}"], ts["z_true"]
            pre = (t >= t_s - 120) & (t < t_s)
            post = (t >= t_end - 120) & (t < t_end)
            win = (t >= t_s) & (t < t_end)
            s_pre = np.nanmedian(sig[pre]) if pre.any() else np.nan
            s_post = np.nanmedian(sig[post]) if post.any() else np.nan
            delay = rise = np.nan
            tw = t[win]
            if kind == "onset":
                hit = np.nonzero(pw[win] >= 0.5)[0]
                delay = tw[hit[0]] - t_s if hit.size else np.inf
                thr = s_pre + 0.9 * (s_post - s_pre)
                hit = np.nonzero(sig[win] >= thr)[0]
                rise = tw[hit[0]] - t_s if hit.size else np.inf
            else:
                if lt == 0.0:
                    hit = np.nonzero(pw[win] < 0.5)[0]
                    delay = tw[hit[0]] - t_s if hit.size else np.inf
                thr = s_post + 0.1 * (s_pre - s_post)
                hit = np.nonzero(sig[win] <= thr)[0]
                rise = tw[hit[0]] - t_s if hit.size else np.inf
            e = z - zt
            def rms(a, b):
                m = (t >= a) & (t < b)
                return float(np.sqrt(np.nanmean(e[m] ** 2))) if m.any() else np.nan
            def rms_true(a, b):
                m = (t >= a) & (t < b)
                return float(np.sqrt(np.nanmean(zt[m] ** 2))) if m.any() else np.nan
            def peak(a, b):
                m = (t >= a) & (t < b)
                return float(np.nanmax(np.abs(e[m]))) if m.any() else np.nan
            ex = ts[f"x_{f}"] - ts["x_true"]
            mh = (t >= t_s) & (t < min(t_s + 120, t_end))
            rows.append(dict(exp=exp, scenario=scenario, seed=seed, filter=f, step=k, kind=kind, t_start=t_s,
                             width=w, level_from=lf, level_to=lt, detect_s=delay, sigma_settle_s=rise,
                             sig_pre=s_pre, sig_post=s_post,
                             rms_z_0_60=rms(t_s, min(t_s + 60, t_end)),
                             rms_z_60_180=rms(t_s + 60, min(t_s + 180, t_end)),
                             peak_z_0_120=peak(t_s, min(t_s + 120, t_end)),
                             rms_x_0_120=float(np.sqrt(np.nanmean(ex[mh] ** 2))) if mh.any() else np.nan,
                             true_z_0_60=rms_true(t_s, min(t_s + 60, t_end)),
                             true_z_60_180=rms_true(t_s + 60, min(t_s + 180, t_end)),
                             sz_end=float(ts[f"sz_{f}"][win][-1]) if win.any() else np.nan))
    return rows


def welch(x, fs, nper):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size < nper:
        return None, None
    step = nper // 2
    win = np.hanning(nper)
    scale = 1.0 / (fs * (win ** 2).sum())
    acc, n = None, 0
    for i in range(0, x.size - nper + 1, step):
        seg = x[i:i + nper] - x[i:i + nper].mean()
        p = np.abs(np.fft.rfft(seg * win)) ** 2 * scale
        p[1:-1] *= 2.0
        acc = p if acc is None else acc + p
        n += 1
    return np.fft.rfftfreq(nper, 1.0 / fs), acc / n


def read_ts(path):
    with open(path) as f:
        hdr = f.readline().strip().split(",")
    data = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)
    return {h: data[:, i] for i, h in enumerate(hdr)}


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

def run_one(exp, sc, seed, filters, tsmode, extended):
    args = [str(BENCH), "--level", "full", "--scenario", sc, "--seed", str(seed), "--filters", ",".join(filters),
            "--excitation"] + (["--extended"] if extended else [])
    tmp = None
    if tsmode:
        fd, tmp = tempfile.mkstemp(suffix=".csv", dir=os.environ.get("TMPDIR"))
        os.close(fd)
        args += ["--timeseries", tmp, "--decimate", "20" if tsmode == "trans" else "20"]
    p = subprocess.run(args, capture_output=True, text=True)
    extra_trans, extra_psd = [], []
    if p.returncode == 0 and tmp:
        ts = read_ts(tmp)
        if tsmode == "trans":
            extra_trans = transition_metrics(ts, filters, sc, seed, exp)
        elif tsmode == "psd":
            fs = 10.0
            m = ts["t"] >= 300.0
            for f in filters:
                for ax in "xyz":
                    key = f"{ax}_{f}"
                    if key not in ts:
                        continue
                    est = ts[key][m]
                    fr, pp = welch(est, fs, 8192)
                    if fr is not None:
                        extra_psd.append((exp, sc, seed, f, ax, fr, pp))
    if tmp:
        os.unlink(tmp)
    return p, extra_trans, extra_psd


def run_grid(work, raw, jobs, extended_any=True):
    require(BENCH)
    OUT.mkdir(parents=True, exist_ok=True)
    hdr = "exp," + header(True)
    cols = hdr.split(",")
    std_cols = ("exp," + header(False)).split(",")
    new = not raw.exists()
    tnew = not TRANS.exists()
    pnew = not PSD.exists()
    done = done_keys(raw)
    work = [w for w in work if (w[0], w[1], str(w[2])) not in done]
    print(f"{len(work)} runs to do", flush=True)
    with raw.open("a") as out, TRANS.open("a") as tout, PSD.open("a") as pout:
        if new:
            out.write(hdr + "\n")
        tw = None
        with ThreadPoolExecutor(max_workers=jobs) as ex:
            futs = {ex.submit(run_one, *w): w for w in work}
            n = 0
            for fut in as_completed(futs):
                w = futs[fut]
                p, trans, psd = fut.result()
                n += 1
                if p.returncode != 0:
                    print(f"ERROR {w[0]} {w[1]} seed {w[2]}: {p.stderr[-300:]}", flush=True)
                    continue
                for line in p.stdout.strip().splitlines():
                    vals = line.split(",")
                    if not w[5]:
                        # Standard + excitation columns: pad the extended block.
                        d = dict(zip(std_cols[1:], vals))
                        vals = [d.get(c, "") for c in cols[1:]]
                    out.write(w[0] + "," + ",".join(vals) + "\n")
                out.flush()
                if trans:
                    if tw is None:
                        tw = csv.DictWriter(tout, fieldnames=list(trans[0].keys()))
                        if tnew:
                            tw.writeheader()
                    for r in trans:
                        tw.writerow(r)
                    tout.flush()
                for (e, sc, seed, f, ax, fr, pp) in psd:
                    if pnew:
                        pout.write("exp,scenario,seed,filter,axis,df," + ",".join(f"{x:.6g}" for x in fr) + "\n")
                        pnew = False
                    pout.write(f"{e},{sc},{seed},{f},{ax},{fr[1]:.6g}," + ",".join(f"{x:.4e}" for x in pp) + "\n")
                pout.flush()
                if n % 20 == 0 or n == len(work):
                    print(f"  {n}/{len(work)}", flush=True)


def cmd_tune(a):
    work = [("T", sc, s, TUNE_FILTERS, None, False) for sc in TUNE_SCEN for s in range(101, 111)]
    run_grid(work, TUNE_RAW, a.jobs)


def cmd_run(a):
    sel = set(a.exp.split(",")) if a.exp else None
    work = []
    for (e, sc, n, filt, ts, ext) in experiments():
        if sel and e not in sel:
            continue
        for s in range(1, n + 1):
            work.append((e, sc, s, filt, ts, ext))
    run_grid(work, RAW, a.jobs)


TRACES = [
    ("calm", "stationary+T2400", ["OU3", "A41@A", "A41@B", "A41@C", "T41@A", "T41@C"]),
    ("onset", f"{J.format(h='0.05')}+on400w10", FULL),
    ("repeated", f"{J.format(h='0.5')}+T1500+on300w10+off600w10+on900w10+off1200w10", FULL),
    ("decay", f"{J.format(h='0.8')}+off600w120", FULL),
    ("swell", "jonswap_fp0.05_hs0.2_g5", CORE),
    ("small", J.format(h="0.02"), FULL),
    ("maneuver", "calm+maneuver", CORE),
]


def cmd_traces(a):
    require(BENCH)
    TRACE_DIR.mkdir(parents=True, exist_ok=True)
    for name, sc, filt in TRACES:
        out = TRACE_DIR / f"{name}.csv"
        fd, tmp = tempfile.mkstemp(suffix=".csv")
        os.close(fd)
        subprocess.run([str(BENCH), "--level", "full", "--scenario", sc, "--seed", "1", "--filters", ",".join(filt),
                        "--excitation", "--timeseries", tmp, "--decimate", "40"], capture_output=True, check=True)
        ts = read_ts(tmp)
        os.unlink(tmp)
        keep = ["t", "x_true", "z_true"] + [k for k in ts if any(k.startswith(p + "_") for p in
                                                                   ("x", "z", "sig", "pw", "zw", "cf", "var", "var0",
                                                                    "fop", "sz"))
                                             and k not in ("x_true", "z_true", "x_total")]
        with out.open("w") as f:
            f.write(",".join(keep) + "\n")
            for i in range(len(ts["t"])):
                f.write(",".join(f"{ts[k][i]:.6g}" for k in keep) + "\n")
        print(f"wrote {out.relative_to(ROOT)}")


def cmd_report(a):
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import carma_excitation_report
    carma_excitation_report.main(OUT)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["tune", "run", "traces", "report"])
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--exp", default="")
    a = ap.parse_args()
    {"tune": cmd_tune, "run": cmd_run, "traces": cmd_traces, "report": cmd_report}[a.cmd](a)


if __name__ == "__main__":
    main()
