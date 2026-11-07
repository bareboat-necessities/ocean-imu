#!/usr/bin/env python3
"""Report of the coupled three-axis CARMA benchmark (EXPERIMENTAL research).

Reads reports/results/carma3_coupled/{raw_runs,ablation_runs,timing_runs}.csv, spectra/*.csv, traces/*.csv, esp32s3_size.csv and
analysis.json (tools/carma3_coupled_benchmark.py, tools/carma3_coupled_analysis.py)
and writes summary.md, paired_stats.json and the figures.

Paired statistics, as tools/carma3_report.py: change = 1 - mean RMS
(candidate) / mean RMS (baseline A41) over paired seeds, 95 % bootstrap
interval (10 000 resamples of seeds, fixed generator seed), win rate (share of
seeds on which the candidate is better) and the Wilcoxon signed-rank p.
Positive change = better than the baseline.

Run: python3 tools/carma3_coupled_report.py
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "reports" / "results" / "carma3_coupled"
sys.path.insert(0, str(Path(__file__).resolve().parent))
import carma3_coupled_benchmark as bench  # noqa: E402

BASE = "A41"
LABEL = {"A41": "Baseline CARMA3 (A41)", "COR:est": "Correlated, estimated", "COR:orc": "Correlated, oracle",
         "COR:est/xy": "Correlated x-y only, estimated",
         "MOD:est/nz": "Complex-OU modal, directional (numerator), estimated",
         "MOD:est": "Complex-OU modal, real coupling, estimated", "MOD:orc/nz": "Complex-OU modal, directional, oracle",
         "MOD:orc/nz/flip": "Complex-OU modal, directional, oracle reversed (180 deg)"}
COLORS = {"A41": "#555555", "COR:est": "#1f77b4", "COR:orc": "#9ecae1", "MOD:est/nz": "#ff7f0e"}
WAVE_EXCLUDE = {"calm", "stationary"}


def load(path):
    return pd.read_csv(path) if path.exists() else pd.DataFrame()


def rms3(df):
    return np.sqrt(df["rms_x"] ** 2 + df["rms_y"] ** 2 + df["rms_z"] ** 2)


def paired(df, sc, cand, col, base=BASE, n_boot=10000):
    a = df[(df.scenario == sc) & (df["filter"] == base)].set_index("seed")[col]
    b = df[(df.scenario == sc) & (df["filter"] == cand)].set_index("seed")[col]
    idx = a.index.intersection(b.index)
    a, b = a.loc[idx].values.astype(float), b.loc[idx].values.astype(float)
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    if len(a) < 3:
        return None
    rng = np.random.default_rng(20261105)
    imp = 1.0 - b.mean() / a.mean()
    bs = []
    for _ in range(n_boot):
        k = rng.integers(0, len(a), len(a))
        bs.append(1.0 - b[k].mean() / a[k].mean())
    lo, hi = np.percentile(bs, [2.5, 97.5])
    try:
        from scipy.stats import wilcoxon
        p = float(wilcoxon(a, b).pvalue) if np.any(a != b) else 1.0
    except Exception:
        p = float("nan")
    return {"n": int(len(a)), "base": float(a.mean()), "cand": float(b.mean()), "change": float(imp),
            "ci": [float(lo), float(hi)], "win": float(np.mean(b < a)), "wilcoxon_p": p}


def fmt_change(r):
    if r is None:
        return "--"
    s = f"{100 * r['change']:+.1f} % [{100 * r['ci'][0]:+.1f}, {100 * r['ci'][1]:+.1f}]"
    if 0.0 < r["win"] < 1.0:
        s += f", win {100 * r['win']:.0f} %"
    return s


# ---------------------------------------------------------------------------
def final_tables(md, stats):
    df = load(bench.RAW)
    if df.empty:
        return df
    filters = [f for f in bench.FINAL_FILTERS.split(",")]
    stats["final"] = {}
    md.append("## Final comparison (complete INS, seeds 1-50)\n")
    md.append("RMS error of the wave displacement in % of incident H_s (mm for calm water and the still IMU), "
              "mean over seeds; in brackets the paired change against the baseline with its 95 % bootstrap "
              "interval (positive = better); better on every seed unless a win rate is given.\n")
    for axis, col in (("x (north)", "x"), ("y (east)", "y"), ("z (heave)", "z")):
        md.append(f"### Axis {axis}\n")
        md.append("| Scenario | n | " + " | ".join(filters) + " |")
        md.append("|---" * (2 + len(filters)) + "|")
        for sc, purpose in bench.SCENARIOS.items():
            sub = df[df.scenario == sc]
            if sub.empty:
                continue
            calm = sc in WAVE_EXCLUDE
            metric = f"rms_{col}" if calm else f"pct_{col}"
            cells = []
            for f in filters:
                v = sub[sub["filter"] == f][metric].mean() * (1000.0 if calm else 1.0)
                if f == BASE:
                    cells.append(f"{v:.2f}")
                else:
                    r = paired(df, sc, f, metric)
                    stats["final"].setdefault(sc, {}).setdefault(f, {})[metric] = r
                    cells.append(f"{v:.2f} ({fmt_change(r)})")
            n = sub[sub["filter"] == BASE].shape[0]
            md.append(f"| {purpose} | {n} | " + " | ".join(cells) + " |")
        md.append("")
    # Consistency and other metrics.
    md.append("### Consistency, velocity, acceleration, attitude (means over seeds)\n")
    md.append("NEES of [p, v] (6 dof; 6 if consistent, larger = overconfident), accelerometer NIS (3 dof; 3 if "
              "consistent), 95 % coverage of the position intervals (x, z), whitened-innovation autocorrelation at "
              "1 s (z), horizontal velocity RMS (mm/s), vertical acceleration RMS (mm/s^2), yaw RMS (deg).\n")
    cols = [("nees_pv6", "NEES6", "{:.2f}"), ("nis", "NIS", "{:.3f}"), ("cov95_x", "cov x", "{:.3f}"),
            ("cov95_z", "cov z", "{:.3f}"), ("wi_ac1_z", "innov ac z", "{:.2f}"), ("rmsv_h", "v_h", "{:.1f}"),
            ("rmsa_z", "a_z", "{:.1f}"), ("yaw_rms_deg", "yaw", "{:.2f}")]
    md.append("| Scenario | filter | " + " | ".join(c[1] for c in cols) + " |")
    md.append("|---" * (2 + len(cols)) + "|")
    for sc in [bench.J29, f"{bench.J29}+sp10000", "mixed_swell_wind", "bimodal", "crossing", "dir_change",
               f"{bench.J29}+turns", "calm"]:
        sub = df[df.scenario == sc]
        if sub.empty:
            continue
        for f in filters:
            s = sub[sub["filter"] == f]
            vals = []
            for c, _, fm in cols:
                v = s[c].mean()
                if c in ("rmsv_h", "rmsa_z"):
                    v *= 1000.0
                vals.append(fm.format(v))
            md.append(f"| {bench.SCENARIOS[sc]} | {f} | " + " | ".join(vals) + " |")
    md.append("")
    return df


def pooled(df, filters, metric_fn, scen):
    """Pooled change over wave scenarios: mean over scenarios of mean over seeds;
    bootstrap resamples seeds jointly across scenarios."""
    rng = np.random.default_rng(7)
    seeds = sorted(df.seed.unique())
    out = {}
    tab = {}
    for f in filters:
        m = []
        for sc in scen:
            s = df[(df.scenario == sc) & (df["filter"] == f)].set_index("seed")
            m.append(metric_fn(s).reindex(seeds).values)
        tab[f] = np.array(m)        # scenarios x seeds
    for f in filters:
        a, b = tab[BASE], tab[f]
        ok = np.isfinite(a) & np.isfinite(b)
        val = np.nanmean(np.where(ok, b, np.nan))
        if f == BASE:
            out[f] = {"value": float(val)}
            continue
        def ratio(k, ok=ok, a=a, b=b):
            return np.nanmean(np.where(ok[:, k], b[:, k], np.nan)) / np.nanmean(np.where(ok[:, k], a[:, k], np.nan))
        r0 = ratio(np.arange(len(seeds)))
        bs = [ratio(rng.integers(0, len(seeds), len(seeds))) for _ in range(4000)]
        out[f] = {"value": float(val), "change": float(1 - r0),
                  "ci": [float(1 - np.percentile(bs, 97.5)), float(1 - np.percentile(bs, 2.5))]}
    return out


def summary_table(md, stats, df, spec_acc):
    if df.empty:
        return
    filters = bench.FINAL_FILTERS.split(",")
    scen = [s for s in bench.SCENARIOS if s not in WAVE_EXCLUDE and not df[df.scenario == s].empty]
    px = pooled(df, filters, lambda s: s["pct_x"], scen)
    py = pooled(df, filters, lambda s: s["pct_y"], scen)
    pz = pooled(df, filters, lambda s: s["pct_z"], scen)
    ne = pooled(df, filters, lambda s: s["nees_pv6"], scen)
    ni = pooled(df, filters, lambda s: s["nis"], scen)
    tim = load(bench.TIMING_RAW)
    size = load(OUT / "esp32s3_size.csv")
    states = {f: int(df[df["filter"] == f]["states"].iloc[0]) for f in filters}
    md.append("## Summary (section 12 of the study brief)\n")
    md.append(f"Pooled over the {len(scen)} wave scenarios (every final scenario except calm water and the still "
              "IMU): mean RMS in % of incident H_s, change against the baseline with a 95 % bootstrap interval "
              "(seeds resampled jointly across scenarios).  NEES/NIS: pooled means (6 / 3 if consistent).  "
              "Cross-spectral accuracy: energy-weighted relative error of the estimated displacement "
              "cross-spectral matrix against the truth's (spectra stage, 10 seeds x 5 seas; lower is better).  "
              "Runtime: host mean / 99th percentile per update (one job, x86-64), ESP32-S3 object code and "
              "filter sizeof (Xtensa GCC 14.2.0 -Os).  Completed experiments only.\n")
    md.append("| Model | States | X RMS | Y RMS | Z RMS | NEES / NIS | Cross-spectral accuracy | Runtime host (us) | ESP32-S3 code / filter (B) |")
    md.append("|---|---|---|---|---|---|---|---|---|")
    stats["summary"] = {}
    sz_name = {"A41": "A41", "COR:est": "COR", "COR:orc": "COR", "MOD:est/nz": "MOD"}
    for f in filters:
        def cell(d, f=f):
            return f"{d[f]['value']:.2f}" + ("" if f == BASE else f" ({100 * d[f]['change']:+.1f} % [{100 * d[f]['ci'][0]:+.1f}, {100 * d[f]['ci'][1]:+.1f}])")
        rt = "--"
        if not tim.empty and f in set(tim["filter"]):
            t = tim[(tim["filter"] == f) & (tim.scenario == bench.J29)]
            rt = f"{t['runtime_ns'].mean() / 1000:.1f} / {t['runtime_p99_ns'].mean() / 1000:.1f}"
        sz = "--"
        if not size.empty:
            r = size[size.variant == sz_name.get(f, "")]
            if not r.empty:
                sz = f"{int(r.code_bytes.iloc[0]):,} / {int(r.sizeof_filter.iloc[0]):,}"
        csa = spec_acc.get(f)
        md.append(f"| {LABEL.get(f, f)} | {states[f]} | {cell(px)} | {cell(py)} | {cell(pz)} | "
                  f"{ne[f]['value']:.2f} / {ni[f]['value']:.3f} | {('%.3f' % csa) if csa is not None else '--'} | {rt} | {sz} |")
        stats["summary"][f] = {"states": states[f], "x": px[f], "y": py[f], "z": pz[f], "nees_pv6": ne[f], "nis": ni[f],
                               "cross_spectral_error": csa}
    md.append("")


def ablation_tables(md, stats):
    df = load(bench.ABL_RAW)
    if df.empty:
        return df
    filters = bench.ABL_FILTERS.split(",")
    stats["ablation"] = {}
    md.append("## Ablations (seeds 1-30)\n")
    md.append("3-D RMS (x, y, z) change against the baseline (positive = better) with 95 % bootstrap intervals; "
              "NEES_pv6 in brackets.  Oracle rows use the true wave acceleration (labelled); flip is the "
              "deliberately reversed oracle coupling.\n")
    df = df.copy()
    df["rms3"] = rms3(df)
    md.append("| filter | " + " | ".join(bench.SCENARIOS.get(s, s) for s in bench.ABL_SCENARIOS) + " |")
    md.append("|---" * (1 + len(bench.ABL_SCENARIOS)) + "|")
    for f in filters:
        cells = []
        for sc in bench.ABL_SCENARIOS:
            sub = df[(df.scenario == sc) & (df["filter"] == f)]
            if f == BASE:
                cells.append(f"{1000 * sub.rms3.mean():.1f} mm ({sub.nees_pv6.mean():.1f})")
                continue
            r = paired(df, sc, f, "rms3")
            stats["ablation"].setdefault(sc, {})[f] = r
            cells.append(f"{fmt_change(r)} ({sub.nees_pv6.mean():.1f})")
        md.append(f"| {LABEL.get(f, f)} | " + " | ".join(cells) + " |")
    md.append("")
    return df


# ---------------------------------------------------------------------------
def spectra_load(sc):
    p = bench.SPEC_DIR / f"{sc.replace('+', '_')}.csv"
    if not p.exists():
        return None
    df = pd.read_csv(p)
    f = df["f_hz"].values
    out = {}
    for c in df.columns[1:]:
        name, key, part = c.split("|")
        out.setdefault((name, key), np.zeros(len(f), complex))
        if part == "re":
            out[(name, key)] += df[c].values
        else:
            out[(name, key)] += 1j * df[c].values
    return f, out


def cross_spectral_accuracy():
    """Energy-weighted relative error of the estimated displacement cross-spectral
    matrix (upper triangle incl. diagonal) against the truth, 0.05-1 Hz,
    averaged over the spectra seas."""
    acc = {}
    for sc in bench.SPEC_SCENARIOS:
        r = spectra_load(sc)
        if r is None:
            continue
        f, S = r
        band = (f >= 0.05) & (f <= 1.0)
        keys = ["00", "01", "02", "11", "12", "22"]
        den = sum(np.trapezoid(np.abs(S[("truth", k)][band]) ** 2, f[band]) for k in keys)
        for flt in bench.SPEC_FILTERS:
            num = sum(np.trapezoid(np.abs(S[(flt, k)][band] - S[("truth", k)][band]) ** 2, f[band]) for k in keys)
            acc.setdefault(flt, []).append(math.sqrt(num / den))
    return {k: float(np.mean(v)) for k, v in acc.items()}


def plots(stats, df, abl, outdir=OUT, ext="svg"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.3, "svg.hashsalt": "carma3-coupled"})
    figs = []

    # (3)-(6) spectra
    for sc in bench.SPEC_SCENARIOS:
        r = spectra_load(sc)
        if r is None:
            continue
        f, S = r
        m = (f > 0.03) & (f < 1.2)
        fig, axs = plt.subplots(3, 4, figsize=(13, 8.2))
        for i, ax in enumerate(axs[0, :3]):
            k = f"{i}{i}"
            ax.loglog(f[m], np.real(S[("truth", k)][m]), color="k", lw=2, label="truth")
            for flt in bench.SPEC_FILTERS:
                ax.loglog(f[m], np.real(S[(flt, k)][m]), color=COLORS.get(flt), lw=1, label=flt)
            ax.set_title(f"{'xyz'[i]} displacement PSD (m^2/Hz)")
        ax = axs[0, 3]
        for i in range(3):
            for flt in ("A41", "MOD:est/nz"):
                ax.loglog(f[m], np.real(S[(flt, f"err{i}")][m]), color=COLORS.get(flt), ls=["-", "--", ":"][i], lw=1,
                          label=f"{flt} {'xyz'[i]}")
        ax.set_title("error PSD (m^2/Hz)")
        ax.legend(fontsize=6)
        for c, (i, j, lab) in enumerate(((0, 2, "x-z"), (1, 2, "y-z"), (0, 1, "x-y"))):
            k = f"{i}{j}"
            axs[1, c].semilogx(f[m], np.abs(S[("truth", k)][m]), color="k", lw=2)
            axs[2, c].semilogx(f[m], np.angle(S[("truth", k)][m], deg=True), color="k", lw=2)
            for flt in bench.SPEC_FILTERS:
                axs[1, c].semilogx(f[m], np.abs(S[(flt, k)][m]), color=COLORS.get(flt), lw=1)
                axs[2, c].semilogx(f[m], np.angle(S[(flt, k)][m], deg=True), color=COLORS.get(flt), lw=1)
            axs[1, c].set_yscale("log")
            axs[1, c].set_title(f"|S_{lab}| (m^2/Hz)")
            axs[2, c].set_title(f"{lab} cross-spectral phase (deg)")
            axs[2, c].set_ylim(-180, 180)
        ax = axs[1, 3]
        for name in ["truth"] + bench.SPEC_FILTERS:
            k = "02"
            cxy = np.abs(S[(name, k)]) ** 2 / (np.real(S[(name, "00")]) * np.real(S[(name, "22")]))
            ax.semilogx(f[m], cxy[m], color="k" if name == "truth" else COLORS.get(name), lw=2 if name == "truth" else 1, label=name)
        ax.set_title("x-z coherence")
        ax.set_ylim(0, 1.05)
        ax.legend(fontsize=6)
        ax = axs[2, 3]
        for name in ["truth"] + bench.SPEC_FILTERS:
            cxy = np.abs(S[(name, "12")]) ** 2 / (np.real(S[(name, "11")]) * np.real(S[(name, "22")]))
            ax.semilogx(f[m], cxy[m], color="k" if name == "truth" else COLORS.get(name), lw=2 if name == "truth" else 1)
        ax.set_title("y-z coherence")
        ax.set_ylim(0, 1.05)
        axs[0, 0].legend(fontsize=6)
        for ax in axs[-1]:
            ax.set_xlabel("frequency (Hz)")
        fig.suptitle(f"Truth vs estimated displacement spectra, {bench.SCENARIOS.get(sc, sc)} (10 seeds, Welch)")
        fig.tight_layout()
        name = f"spectra_{sc.replace('+', '_')}.{ext}"
        fig.savefig(outdir / name, dpi=150)
        plt.close(fig)
        figs.append(name)

    # (7) RMS change by sea state and model
    if not df.empty:
        filters = [f for f in bench.FINAL_FILTERS.split(",") if f != BASE]
        scen = [s for s in bench.SCENARIOS if s not in WAVE_EXCLUDE]
        fig, axs = plt.subplots(1, 3, figsize=(14, 6.2), sharey=True)
        for a, col in enumerate("xyz"):
            ax = axs[a]
            y = np.arange(len(scen))
            w = 0.8 / len(filters)
            for k, f in enumerate(filters):
                v, lo, hi = [], [], []
                for sc in scen:
                    r = stats["final"].get(sc, {}).get(f, {}).get(f"pct_{col}")
                    v.append(100 * r["change"] if r else np.nan)
                    lo.append(100 * (r["change"] - r["ci"][0]) if r else 0)
                    hi.append(100 * (r["ci"][1] - r["change"]) if r else 0)
                ax.barh(y + (k - len(filters) / 2 + 0.5) * w, v, w, xerr=[lo, hi], color=COLORS.get(f), label=f,
                        error_kw={"lw": 0.6})
            ax.axvline(0, color="k", lw=0.8)
            ax.set_yticks(y)
            ax.set_yticklabels([bench.SCENARIOS[s] for s in scen], fontsize=7)
            ax.set_title(f"{col}: change of RMS vs baseline (%), + = better")
            ax.invert_yaxis()
        axs[0].legend(fontsize=7, loc="lower left")
        fig.tight_layout()
        fig.savefig(outdir / f"rms_change_by_sea.{ext}", dpi=150)
        plt.close(fig)
        figs.append(f"rms_change_by_sea.{ext}")

        # (8) NEES / NIS
        fig, axs = plt.subplots(1, 2, figsize=(13, 5.5), sharey=True)
        allf = bench.FINAL_FILTERS.split(",")
        for a, (col, ideal, lab) in enumerate((("nees_pv6", 6.0, "NEES of [p, v] (6 dof)"), ("nis", 3.0, "accelerometer NIS (3 dof)"))):
            ax = axs[a]
            y = np.arange(len(bench.SCENARIOS))
            w = 0.8 / len(allf)
            for k, f in enumerate(allf):
                v = [df[(df.scenario == sc) & (df["filter"] == f)][col].mean() for sc in bench.SCENARIOS]
                ax.barh(y + (k - len(allf) / 2 + 0.5) * w, v, w, color=COLORS.get(f), label=f)
            ax.axvline(ideal, color="k", ls="--", lw=0.8)
            ax.set_xscale("log")
            ax.set_yticks(y)
            ax.set_yticklabels(list(bench.SCENARIOS.values()), fontsize=7)
            ax.set_title(lab + " (dashed: consistent)")
            ax.invert_yaxis()
        axs[0].legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(outdir / f"consistency_nees_nis.{ext}", dpi=150)
        plt.close(fig)
        figs.append(f"consistency_nees_nis.{ext}")

    # (11) sensitivity to incorrect coupling (ablations)
    if abl is not None and not abl.empty:
        filters = [f for f in bench.ABL_FILTERS.split(",") if f != BASE]
        fig, ax = plt.subplots(figsize=(12, 6))
        y = np.arange(len(filters))
        w = 0.8 / len(bench.ABL_SCENARIOS)
        for k, sc in enumerate(bench.ABL_SCENARIOS):
            v = []
            for f in filters:
                r = stats["ablation"].get(sc, {}).get(f)
                v.append(100 * r["change"] if r else np.nan)
            ax.barh(y + (k - len(bench.ABL_SCENARIOS) / 2 + 0.5) * w, v, w, label=bench.SCENARIOS.get(sc, sc))
        ax.axvline(0, color="k", lw=0.8)
        ax.set_yticks(y)
        ax.set_yticklabels([LABEL.get(f, f) for f in filters], fontsize=8)
        ax.invert_yaxis()
        ax.set_xlabel("3-D RMS change vs baseline (%), + = better")
        ax.legend(fontsize=7)
        ax.set_title("Co-spectrum only, real modal coupling, oracles and a reversed oracle (30 seeds)")
        fig.tight_layout()
        fig.savefig(outdir / f"ablation_sensitivity.{ext}", dpi=150)
        plt.close(fig)
        figs.append(f"ablation_sensitivity.{ext}")

    figs += trace_plots(plt, outdir, ext)

    # (12) execution time and memory
    tim = load(bench.TIMING_RAW)
    size = load(OUT / "esp32s3_size.csv")
    if not tim.empty or not size.empty:
        fig, axs = plt.subplots(1, 3, figsize=(13, 3.8))
        if not tim.empty:
            labs = bench.TIMING_FILTERS
            for k, sc in enumerate((bench.J29, f"{bench.J29}+jit1", f"{bench.J29}+jit")):
                t = tim[tim.scenario == sc]
                v = [t[t["filter"] == f]["runtime_ns"].mean() / 1000 for f in labs]
                p = [t[t["filter"] == f]["runtime_p99_ns"].mean() / 1000 for f in labs]
                x = np.arange(len(labs)) + (k - 1) * 0.27
                axs[0].bar(x, v, 0.27, label={bench.J29: "fixed dt", f"{bench.J29}+jit1": "+-1 % (1 us)", f"{bench.J29}+jit": "+-25 %"}[sc])
                axs[1].bar(x, p, 0.27)
            for ax, t in zip(axs[:2], ("host mean update (us)", "host 99th percentile update (us)")):
                ax.set_xticks(np.arange(len(labs)))
                ax.set_xticklabels(labs)
                ax.set_title(t)
                ax.set_yscale("log")
            axs[0].legend(fontsize=7)
        if not size.empty:
            s = size[size.variant.isin(["OU3", "A41", "COR", "MOD"])]
            x = np.arange(len(s))
            axs[2].bar(x - 0.2, s.code_bytes / 1000, 0.4, label="code (kB)")
            axs[2].bar(x + 0.2, s.sizeof_filter / 1000, 0.4, label="filter object (kB)")
            ax2 = axs[2].twinx()
            ax2.plot(x, s.max_stack_frame / 1000, "kD", label="largest stack frame (kB)")
            ax2.set_ylabel("largest stack frame (kB)")
            axs[2].set_xticks(x)
            axs[2].set_xticklabels(s.variant)
            axs[2].set_title("ESP32-S3 (Xtensa GCC 14.2, -Os)")
            axs[2].legend(fontsize=7, loc="upper left")
        fig.tight_layout()
        fig.savefig(outdir / f"cost_time_memory.{ext}", dpi=150)
        plt.close(fig)
        figs.append(f"cost_time_memory.{ext}")
    return figs


def trace_plots(plt, outdir=OUT, ext="svg"):
    figs = []
    def tload(sc):
        p = bench.TRACE_DIR / f"{sc.replace('+', '_')}_seed1.csv"
        return pd.read_csv(p) if p.exists() else None
    names = ["A41", "COR:est", "MOD:est/nz"]
    t = tload(bench.J29)
    if t is not None:
        # (1) true vs estimated X/Y/Z
        w = (t.t >= 700) & (t.t <= 760)
        fig, axs = plt.subplots(3, 1, figsize=(11, 7), sharex=True)
        for a, c in enumerate("xyz"):
            axs[a].plot(t.t[w], t[f"{c}_true"][w], color="k", lw=2, label="truth")
            for n in names:
                axs[a].plot(t.t[w], t[f"{c}_{n}"][w], color=COLORS.get(n), lw=1, label=n)
            axs[a].set_ylabel(f"{c} (m)")
        axs[0].legend(ncol=4, fontsize=7)
        axs[-1].set_xlabel("time (s)")
        fig.suptitle("Directionally spread JONSWAP 0.29 Hz, 0.8 m, seed 1: truth vs estimates (NED, z down)")
        fig.tight_layout()
        fig.savefig(outdir / f"timeseries_xyz.{ext}", dpi=150)
        plt.close(fig)
        figs.append(f"timeseries_xyz.{ext}")
        # (2) 3-D trajectory
        fig = plt.figure(figsize=(10, 4.6))
        w = (t.t >= 700) & (t.t <= 730)
        for k, (n, title) in enumerate((("A41", "baseline"), ("MOD:est/nz", "modal directional coupling, estimated"))):
            ax = fig.add_subplot(1, 2, k + 1, projection="3d")
            ax.plot(t.x_true[w], t.y_true[w], -t.z_true[w], color="k", lw=1.8, label="truth")
            ax.plot(t[f"x_{n}"][w], t[f"y_{n}"][w], -t[f"z_{n}"][w], color=COLORS[n], lw=1, label=n)
            ax.set_xlabel("north (m)"); ax.set_ylabel("east (m)"); ax.set_zlabel("up (m)")
            ax.set_title(title)
            ax.legend(fontsize=7)
        fig.suptitle("Three-dimensional wave displacement, 30 s, JONSWAP seed 1")
        fig.tight_layout()
        fig.savefig(outdir / f"trajectory_3d.{ext}", dpi=150)
        plt.close(fig)
        figs.append(f"trajectory_3d.{ext}")

    def rolling(t, n, c, win=30.0):
        dt = np.median(np.diff(t.t.values))
        e2 = (t[f"{c}_{n}"] - t[f"{c}_true"]) ** 2
        k = max(1, int(round(win / dt)))
        return np.sqrt(e2.rolling(k, min_periods=k // 2).mean())

    # (9) heading change and wave-direction change; (10) onset / cessation
    for sc, title, fname in ((f"{bench.J29}+turns", "90 deg zigzag heading changes every 120 s", f"heading_change.{ext}"),
                             ("dir_change", "wave direction 30 -> 120 deg over 570-630 s", f"direction_change.{ext}"),
                             (f"{bench.J29}+on400w60+off900w60", "wave onset at 400 s, cessation at 900 s", f"onset_cessation.{ext}")):
        t = tload(sc)
        if t is None:
            continue
        fig, axs = plt.subplots(3, 1, figsize=(11, 7.5), sharex=True)
        for n in names:
            h = np.sqrt(rolling(t, n, "x") ** 2 + rolling(t, n, "y") ** 2)
            axs[0].plot(t.t, 1000 * h, color=COLORS.get(n), lw=1, label=n)
            axs[1].plot(t.t, 1000 * rolling(t, n, "z"), color=COLORS.get(n), lw=1)
            if f"phxz_{n}" in t.columns and n != BASE:
                axs[2].plot(t.t, t[f"phxz_{n}"], color=COLORS.get(n), lw=1, label=f"{n} applied x-z phase")
        if "turns" in sc and "dyaw_A41" in t.columns:
            ax2 = axs[2].twinx()
            ax2.plot(t.t, t["dyaw_A41"], color="#555555", lw=0.6, ls=":")
            ax2.set_ylabel("yaw error A41 (deg)")
        axs[0].set_ylabel("horizontal 30 s RMS (mm)")
        axs[1].set_ylabel("heave 30 s RMS (mm)")
        axs[2].set_ylabel("x-z phase of Sigma (deg)")
        axs[0].legend(ncol=5, fontsize=7)
        axs[2].legend(fontsize=6)
        axs[-1].set_xlabel("time (s)")
        fig.suptitle(f"{title} (seed 1)")
        fig.tight_layout()
        fig.savefig(outdir / fname, dpi=150)
        plt.close(fig)
        figs.append(fname)
    return figs


def analysis_section(md):
    p = OUT / "analysis.json"
    if not p.exists():
        return
    a = json.loads(p.read_text())
    c = a["complexification"]
    md.append("## Analytic results (tools/carma3_coupled_analysis.py)\n")
    md.append(f"* Complexified C41J: real-part PSD equals C41J to {c['max_rel_error_realpart_psd_vs_C41J']:.1e} "
              f"(relative), [p, v, a] stationary covariance to {c['realpart_stationary_cov_pva_max_rel_dev_vs_C41J']:.1e}; "
              "mirror leakage |H(-iw)|^2/|H(iw)|^2 at kappa = 1: " +
              ", ".join(f"{k} Hz {v:.3f}" for k, v in c["mirror_leakage_Hminus_over_Hplus"].items()) + ".")
    s = a["system_checks"]
    md.append(f"* Kronecker identity (Experiment A): F {s['kronecker_F_err']:.1e}, Q {s['kronecker_Q_rel_err']:.1e}; "
              f"Hermitian Kronecker (Experiment B): Q {s['hermitian_kronecker_Q_rel_err']:.1e}; SHM with a real "
              f"excitation reproduces COR's cross-spectral matrix to {s['shm_real_sigma_vs_cor_spectrum_rel_err']:.1e}.")
    md.append("\nSteady-state (DARE, white accelerometer noise only, oracle parameters) horizontal / heave RMS (mm):\n")
    names = ["IND", "COR", "MODN", "MODN_wrong180", "SHM", "SHM_diag", "MCA", "MCA_diag", "SHR", "SHM_wrong180", "SHM_wrong90",
             "K2", "K2_diag"]
    md.append("| sea | " + " | ".join(names) + " |")
    md.append("|---" * (1 + len(names)) + "|")
    for sea, r in a["steady_state"].items():
        m = r["models"]["hnoise_x1"]
        cells = []
        for n in names:
            if n in m:
                cells.append(f"{1000 * math.hypot(m[n]['rms_m'][0], m[n]['rms_m'][1]):.1f} / {1000 * m[n]['rms_m'][2]:.1f}")
            else:
                cells.append("--")
        md.append(f"| {sea} | " + " | ".join(cells) + " |")
    md.append("\nRepresentable share of the truth's cross-spectral energy (real-only = COR, constant phase = SHM) "
              "and phase at the energy peak:\n")
    md.append("| sea | x-z real / const-phase | y-z real / const-phase | x-z phase (deg) |")
    md.append("|---|---|---|---|")
    for sea, r in a["steady_state"].items():
        rp = r["representable"]
        md.append(f"| {sea} | {rp['xz']['real_part_energy_fraction']:.2f} / {rp['xz']['constant_phase_fit_energy_fraction']:.2f} | "
                  f"{rp['yz']['real_part_energy_fraction']:.2f} / {rp['yz']['constant_phase_fit_energy_fraction']:.2f} | "
                  f"{rp['xz']['phase_deg_at_peak']:.0f} |")
    if "modal" in a:
        mo = a["modal"]
        eq = mo["equivalence_vs_kappa"]
        worst = {k: max(r[k] for r in eq) for k in eq[0] if k != "kappa"}
        ind = mo["independent_driving"]
        md.append("\nComplex-OU modal form of C41J (kappa 0.2-2.5): worst transfer-function error "
                  f"{worst['tf_rel_err']:.1e}, stationary covariance {worst['stationary_rel_err']:.1e}, F {worst['F_err']:.1e}, "
                  f"Q {worst['Q_rel_err']:.1e}, autocovariance {worst['autocov_rel_err']:.1e}; cond(T) <= {worst['cond_T']:.0f}. "
                  f"Residue sums: sum Re r = {mo['sum_re_r']:.1e}, sum Re(r p) = {mo['sum_re_rp']:.1e}, sum Im r = {mo['sum_im_r']:.4f}. "
                  f"Independently driven modes: Q_pp slope {ind['independent']['Qpp_log2_slope']:.2f} (exact: "
                  f"{ind['common']['Qpp_log2_slope']:.2f}), displacement variance x{ind['independent']['var_p_ratio_vs_common']:.3f}, "
                  f"acceleration PSD high-frequency log2 slope {ind['independent']['a_psd_hf_log2_slope']:+.2f} per octave "
                  f"(exact: {ind['common']['a_psd_hf_log2_slope']:+.2f}).")
    idn = a["identifiability"]
    md.append(f"\nIdentifiability (Whittle Fisher information, {idn['T_s']:.0f} s of three-axis acceleration, JONSWAP):\n")
    md.append("| family | parameters | numerical rank | smallest / largest FIM eigenvalue |")
    md.append("|---|---|---|---|")
    for k in ("COR", "SHM", "MCA", "MCARMA41_unrestricted"):
        v = idn[k]
        md.append(f"| {k} | {v['n_params']} | {v['numerical_rank']} | {v['fim_eig_min']:.2e} / {v['fim_eig_max']:.2e} |")
    md.append("")


# ---------------------------------------------------------------------------
# LaTeX table parts for doc/kalman_carma/kalman_carma_coupled.tex
# ---------------------------------------------------------------------------

TEX_DIR = ROOT / "doc" / "kalman_carma"
SHORT = {"A41": "Baseline", "COR:est": "A est.", "COR:orc": "A oracle", "MOD:est/nz": "Modal dir. est."}


def tex_escape(t):
    return (t.replace("\\", "/").replace("_", "\\_").replace("%", "\\%").replace("&", "\\&")
            .replace("->", "$\\to$").replace("^2", "$^2$").replace("^8", "$^8$"))


def fit_tables(text):
    """Wrap every tabular in an adjustbox no wider than its float."""
    out, wide = [], False
    for line in text.split("\n"):
        if line.startswith("\\begin{table*}"):
            wide = True
        elif line.startswith("\\begin{table}"):
            wide = False
        if line.startswith("\\begin{tabular}"):
            out.append("\\begin{adjustbox}{max width=" + ("\\textwidth}" if wide else "\\linewidth}"))
        out.append(line)
        if line.startswith("\\end{tabular}"):
            out.append("\\end{adjustbox}")
    return "\n".join(out)


def write_tex_parts(stats, df, abl):
    TEX_DIR.mkdir(parents=True, exist_ok=True)
    a = json.loads((OUT / "analysis.json").read_text()) if (OUT / "analysis.json").exists() else None
    if a:
        L = ["\\begin{table}[t]", "\\caption{Share of the truth's $x$--$z$ / $y$--$z$ cross-spectral energy a real "
             "(A) or constant-phase (B) coefficient can represent, and the $x$--$z$ phase at the energy peak.}",
             "\\label{tab:repr}", "\\centering\\small", "\\begin{tabular}{lccc}", "\\toprule",
             "Sea & $x$--$z$ A / B & $y$--$z$ A / B & phase \\\\", "\\midrule"]
        for sea, r in a["steady_state"].items():
            rp = r["representable"]
            L.append(f"{tex_escape(sea)} & {rp['xz']['real_part_energy_fraction']:.2f} / {rp['xz']['constant_phase_fit_energy_fraction']:.2f} & "
                     f"{rp['yz']['real_part_energy_fraction']:.2f} / {rp['yz']['constant_phase_fit_energy_fraction']:.2f} & "
                     f"\\SI{{{rp['xz']['phase_deg_at_peak']:.0f}}}{{\\degree}} \\\\")
        L += ["\\bottomrule", "\\end{tabular}", "\\end{table}"]
        (TEX_DIR / "carma-coupled-repr-generated.tex-part").write_text(fit_tables("\n".join(L)) + "\n")
        names = [("IND", "Base"), ("COR", "A"), ("MODN", "Modal dir."), ("MODN_wrong180", "M. rev."), ("SHM", "B"),
                 ("SHM_diag", "B diag"), ("MCA", "C"), ("MCA_diag", "C diag"), ("SHR", "B red."), ("SHM_wrong180", "B rev."),
                 ("K2", "2 modes")]
        L = ["\\begin{table*}[t]", "\\caption{Steady-state RMS error, horizontal / heave (mm), white accelerometer "
             "noise only, oracle parameters, against the true cross-spectral matrix. \\emph{diag}: the same marginals "
             "without cross terms; \\emph{red.}: state-reduced shared mode; \\emph{rev.}: coupling of the reversed wave "
             "direction; \\emph{Modal dir.}: complex-OU modal form with per-axis numerator zeros (the kinematically admissible "
             "phase coupling); \\emph{2 modes}: oracle two-mode model (two-system seas).}", "\\label{tab:steady}",
             "\\centering\\small", "\\begin{tabular}{l" + "c" * len(names) + "}", "\\toprule",
             "Sea & " + " & ".join(n for _, n in names) + " \\\\", "\\midrule"]
        for sea, r in a["steady_state"].items():
            m = r["models"]["hnoise_x1"]
            cells = []
            for k, _ in names:
                cells.append(f"{1000 * math.hypot(m[k]['rms_m'][0], m[k]['rms_m'][1]):.1f}/{1000 * m[k]['rms_m'][2]:.1f}" if k in m else "--")
            L.append(f"{tex_escape(sea)} & " + " & ".join(cells) + " \\\\")
        L += ["\\bottomrule", "\\end{tabular}", "\\end{table*}"]
        (TEX_DIR / "carma-coupled-steady-generated.tex-part").write_text(fit_tables("\n".join(L)) + "\n")
        if "modal" in a:
            mo = a["modal"]
            ind = mo["independent_driving"]
            L = ["\\begin{table}[t]", "\\caption{Complex-OU modal form of C41J against the baseline realization: "
                 "relative errors of the transfer function, stationary covariance, discrete $\\mat Q$ (under "
                 "$\\mat T$) and displacement autocovariance, and the condition number of $\\mat T$. Last rows: "
                 "independently driven modes.}", "\\label{tab:modal}", "\\centering\\small",
                 "\\begin{tabular}{lccccc}", "\\toprule",
                 "$\\kappa$ & $H$ & $\\mat P_\\infty$ & $\\mat Q$ & autocov. & cond $\\mat T$ \\\\", "\\midrule"]
            rows = mo["equivalence_vs_kappa"]
            for r in rows[::3] + ([rows[-1]] if (len(rows) - 1) % 3 else []):
                L.append(f"{r['kappa']:.2f} & \\num{{{r['tf_rel_err']:.1e}}} & \\num{{{r['stationary_rel_err']:.1e}}} & "
                         f"\\num{{{r['Q_rel_err']:.1e}}} & \\num{{{r['autocov_rel_err']:.1e}}} & {r['cond_T']:.0f} \\\\")
            L += ["\\midrule", "\\multicolumn{6}{l}{$Q_{pp}$ log-slope in $h$: common noise "
                  f"{ind['common']['Qpp_log2_slope']:.2f}, independent {ind['independent']['Qpp_log2_slope']:.2f}}} \\\\",
                  "\\multicolumn{6}{l}{displacement variance, independent / common: "
                  f"{ind['independent']['var_p_ratio_vs_common']:.3f}}} \\\\",
                  "\\multicolumn{6}{l}{acceleration PSD log$_2$ slope per octave at 5--10\\,Hz: "
                  f"{ind['common']['a_psd_hf_log2_slope']:+.2f} / {ind['independent']['a_psd_hf_log2_slope']:+.2f}}} \\\\",
                  "\\bottomrule", "\\end{tabular}", "\\end{table}"]
            (TEX_DIR / "carma-coupled-modal-generated.tex-part").write_text(fit_tables("\n".join(L)) + "\n")
    if df is None or df.empty:
        return
    filters = bench.FINAL_FILTERS.split(",")
    L = []
    summ = stats.get("summary", {})
    if summ:
        L += ["\\begin{table*}[t]", "\\caption{Summary over the wave scenarios (50 seeds each): mean RMS in \\% of "
              "incident $H_s$ with the change against the baseline and its 95\\,\\% bootstrap interval (positive = "
              "better), pooled NEES of $[\\vct p,\\vct v]$ (6 = consistent) and NIS (3), cross-spectral error "
              "(lower = better), host update time (mean / 99th percentile) and ESP32-S3 object code / filter object.}",
              "\\label{tab:summary}", "\\centering\\small", "\\begin{tabular}{lccccccc}", "\\toprule",
              "Model & States & $x$ RMS & $y$ RMS & $z$ RMS & NEES/NIS & CS err. & host \\si{\\micro s}; ESP32 B \\\\", "\\midrule"]
        tim = load(bench.TIMING_RAW)
        size = load(OUT / "esp32s3_size.csv")
        sz_name = {"A41": "A41", "COR:est": "COR", "COR:orc": "COR", "MOD:est/nz": "MOD"}
        for f in filters:
            r = summ[f]
            def cell(d, f=f):
                if f == BASE:
                    return f"{d['value']:.2f}"
                return f"{d['value']:.2f} ({100 * d['change']:+.1f} [{100 * d['ci'][0]:+.1f}, {100 * d['ci'][1]:+.1f}])"
            rt = "--"
            if not tim.empty and f in set(tim["filter"]):
                t = tim[(tim["filter"] == f) & (tim.scenario == bench.J29)]
                rt = f"{t['runtime_ns'].mean() / 1000:.1f}/{t['runtime_p99_ns'].mean() / 1000:.0f}"
            if not size.empty:
                q = size[size.variant == sz_name.get(f, "")]
                if not q.empty:
                    rt += f"; {int(q.code_bytes.iloc[0]) / 1000:.0f}k/{int(q.sizeof_filter.iloc[0]) / 1000:.1f}k"
            csa = r.get("cross_spectral_error")
            L.append(f"{tex_escape(LABEL.get(f, f))} & {r['states']} & {cell(r['x'])} & {cell(r['y'])} & {cell(r['z'])} & "
                     f"{r['nees_pv6']['value']:.2f}/{r['nis']['value']:.3f} & {('%.3f' % csa) if csa is not None else '--'} & {rt} \\\\")
        L += ["\\bottomrule", "\\end{tabular}", "\\end{table*}", ""]
    # per-sea 3-D change table
    L += ["\\begin{table*}[t]", "\\caption{Per-sea RMS (baseline, \\% of $H_s$; mm for calm water and the still IMU) "
          "and paired change of the other models, $x$ / $y$ / $z$, \\% (positive = better; 50 seeds).}",
          "\\label{tab:persea}", "\\centering\\scriptsize", "\\begin{tabular}{l" + "c" * len(filters) + "}", "\\toprule",
          "Sea & " + " & ".join(SHORT.get(f, f) for f in filters) + " \\\\", "\\midrule"]
    for sc, purpose in bench.SCENARIOS.items():
        sub = df[df.scenario == sc]
        if sub.empty:
            continue
        calm = sc in WAVE_EXCLUDE
        cells = []
        for f in filters:
            if f == BASE:
                v = [sub[sub["filter"] == f][(f"rms_{c}" if calm else f"pct_{c}")].mean() * (1000 if calm else 1) for c in "xyz"]
                cells.append("/".join(f"{x:.1f}" if calm else f"{x:.2f}" for x in v))
            else:
                v = []
                for c in "xyz":
                    r = stats["final"].get(sc, {}).get(f, {}).get(f"rms_{c}" if calm else f"pct_{c}")
                    v.append(f"{100 * r['change']:+.0f}" if r else "--")
                cells.append("/".join(v))
        L.append(f"{tex_escape(purpose)} & " + " & ".join(cells) + " \\\\")
    L += ["\\bottomrule", "\\end{tabular}", "\\end{table*}", ""]
    # ablations
    if abl is not None and not abl.empty:
        af = [f for f in bench.ABL_FILTERS.split(",") if f != BASE]
        L += ["\\begin{table*}[t]", "\\caption{Ablations (30 seeds): paired 3-D RMS change against the baseline, \\% "
              "(positive = better), and NEES$_6$ in brackets.}", "\\label{tab:ablation}", "\\centering\\scriptsize",
              "\\begin{tabular}{l" + "c" * len(bench.ABL_SCENARIOS) + "}", "\\toprule",
              "Variant & " + " & ".join(tex_escape(bench.SCENARIOS.get(s, s)) for s in bench.ABL_SCENARIOS) + " \\\\", "\\midrule"]
        for f in af:
            cells = []
            for sc in bench.ABL_SCENARIOS:
                r = stats["ablation"].get(sc, {}).get(f)
                ne = abl[(abl.scenario == sc) & (abl["filter"] == f)].nees_pv6.mean()
                cells.append(f"{100 * r['change']:+.1f} ({ne:.1f})" if r else "--")
            L.append(f"{tex_escape(LABEL.get(f, f))} & " + " & ".join(cells) + " \\\\")
        L += ["\\bottomrule", "\\end{tabular}", "\\end{table*}", ""]
    (TEX_DIR / "carma-coupled-results-generated.tex-part").write_text(fit_tables("\n".join(L)) + "\n")


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--figures-only", action="store_true", help="write only the figures (for plots/kalman_carma)")
    ap.add_argument("--outdir", default=str(OUT))
    ap.add_argument("--ext", default="svg", choices=["svg", "png", "pdf"])
    args = ap.parse_args(argv)
    if args.figures_only:
        stats = {}
        quiet = []
        df = final_tables(quiet, stats)
        abl = ablation_tables(quiet, stats)
        figs = plots(stats, df, abl, Path(args.outdir), args.ext)
        print(f"wrote {len(figs)} figures to {args.outdir}")
        return 0
    stats = {}
    md = ["# Coupled three-axis CARMA priors: generated results\n",
          "Generated by tools/carma3_coupled_report.py; do not edit by hand.  Simulation only (no hardware).\n"]
    spec_acc = cross_spectral_accuracy()
    stats["cross_spectral_accuracy"] = spec_acc
    df = final_tables(md, stats)
    summary_table(md, stats, df, spec_acc)
    abl = ablation_tables(md, stats)
    analysis_section(md)
    figs = plots(stats, df, abl)
    md.append("## Figures\n")
    for fn in figs:
        md.append(f"![{fn}]({fn})")
    write_tex_parts(stats, df, abl)
    (OUT / "summary.md").write_text("\n".join(md) + "\n")
    (OUT / "paired_stats.json").write_text(json.dumps(stats, indent=1))
    print(f"wrote {OUT / 'summary.md'} ({len(figs)} figures)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
