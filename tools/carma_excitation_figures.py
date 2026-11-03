"""Figures of the CARMA wave-excitation study (tools/carma_excitation_report.py).

Colours follow one entity through every figure (the validated categorical order:
evidence C blue, debiased B orange, TFG evidence aqua, floor-removed B0 yellow,
C without the evidence weight magenta, C without the frequency hold green, TFG
legacy violet, legacy A red; OU-III grey).  Legacy filters are dashed.  The
tables of summary.md carry every plotted number.
"""
from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path

import numpy as np

from carma_excitation_report import LABEL, aligned, bootstrap_ci, fnum, load_trace, rms_h, style

C = {"A41@C": "#2a78d6", "A41@B": "#eb6834", "T41@C": "#1baf7a", "A41@B0": "#eda100", "A41@Cng": "#e87ba4",
     "A41@Cnf": "#008300", "T41@A": "#6250d6", "A41@A": "#e34948", "OU3": "#8a8984", "T41@B": "#e87ba4"}
LS = {"A41@A": "--", "T41@A": "--", "OU3": ":"}
ORDER = ["A41@C", "A41@B", "T41@C", "A41@B0", "A41@Cng", "A41@Cnf", "T41@A", "A41@A", "OU3"]
INK, INK2 = "#0b0b0b", "#52514e"


def hs_of(sc):
    m = re.search(r"_hs([0-9.]+)", sc)
    return float(m.group(1)) if m else 0.0


def fp_of(sc):
    m = re.search(r"_fp([0-9.]+)", sc)
    return float(m.group(1)) if m else np.nan


def save(fig, svgd, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    svgd.save_svg(fig, path, bbox_inches="tight")


def line(ax, x, y, f, **kw):
    ax.plot(x, y, color=C.get(f, INK2), ls=LS.get(f, "-"), lw=kw.pop("lw", 1.6), label=LABEL.get(f, f), **kw)


def fig_calm(r, plt, svgd, out):
    tr = out / "traces" / "calm.csv"
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.6), gridspec_kw={"width_ratios": [1.6, 1]})
    if tr.exists():
        t = load_trace(tr)
        for f in ["A41@A", "T41@A", "OU3", "A41@B", "A41@C", "T41@C"]:
            k = f"z_{f}"
            if k in t:
                line(axs[0], t["t"], 1e3 * (t[k] - t["z_true"]), f, lw=1.1)
        axs[0].set_xlabel("time (s)")
        axs[0].set_ylabel("heave estimate, true heave 0 (mm)")
        axs[0].set_title("Still IMU, seed 1", loc="left", fontsize=10, color=INK)
        axs[0].set_xlim(0, t["t"][-1])
        axs[0].set_ylim(-250, 250)
        axs[0].legend(fontsize=7, ncol=3, frameon=False, loc="upper right")
    style(axs[0])
    d = r.g.get(("E1", "stationary+T2400"))
    if d:
        fl = [f for f in ["A41@A", "T41@A", "OU3", "A41@B0", "A41@B", "A41@C", "T41@C"] if f in d]
        seeds, vz = aligned(d, fl, lambda x: fnum(x, "rms_z"))
        _, vh = aligned(d, fl, rms_h)
        y = np.arange(len(fl))
        for i, f in enumerate(fl):
            for off, v, a in ((-0.18, vz[f], 1.0), (0.18, vh[f], 0.45)):
                m = 1e3 * np.mean(v)
                lo, hi = (1e3 * np.array(bootstrap_ci(v)))
                axs[1].barh(i + off, max(m, 1e-3), height=0.32, color=C.get(f, INK2), alpha=a, edgecolor="none")
                axs[1].errorbar(max(m, 1e-3), i + off, xerr=[[max(m - lo, 0)], [max(hi - m, 0)]], color=INK, lw=0.8,
                                capsize=2)
                axs[1].text(max(hi, m, 1e-3) * 1.3, i + off, f"{m:.2f}" if m < 10 else f"{m:.0f}", va="center",
                            fontsize=7, color=INK2)
        axs[1].set_yticks(y)
        axs[1].set_yticklabels([LABEL.get(f, f) for f in fl], fontsize=8)
        axs[1].set_xscale("log")
        axs[1].set_xlim(1e-3, 2e3)
        axs[1].set_xlabel("RMS (mm), log; dark heave, light horizontal")
        axs[1].set_title(f"2400 s, t >= 300 s, {len(seeds)} seeds, 95 % CI", loc="left", fontsize=10, color=INK)
    style(axs[1])
    fig.suptitle("Phantom displacement in still water", x=0.01, ha="left", fontsize=12, color=INK)
    save(fig, svgd, out / "figures" / "fig01_calm_phantom.svg")
    plt.close(fig)


def fig_psd(r, plt, svgd, out):
    path = out / "psd_calm.csv"
    if not path.exists():
        return
    with path.open() as f:
        rd = csv.reader(f)
        hdr = next(rd)
        freqs = np.array([float(x) for x in hdr[6:]])
        acc = defaultdict(list)
        for row in rd:
            if row[1] != "stationary+T2400":
                continue
            acc[(row[3], row[4])].append(np.array([float(x) for x in row[6:]]))
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.8), sharey=True)
    for ax, axis, name in ((axs[0], "z", "heave"), (axs[1], "x", "surge")):
        for f in ["A41@A", "T41@A", "OU3", "A41@B", "A41@C", "T41@C"]:
            v = acc.get((f, axis))
            if not v:
                continue
            p = np.mean(v, axis=0)
            m = freqs > 0
            line(ax, freqs[m], np.maximum(p[m], 1e-14), f, lw=1.3)
        ax.axvspan(freqs[1], 0.05, color="#eeeeea", zorder=0)
        ax.text(0.0012, 2e-13, "below 0.05 Hz", fontsize=8, color=INK2)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(1e-3, 5)
        ax.set_ylim(1e-13, 1e2)
        ax.set_xlabel("frequency (Hz)")
        ax.set_title(f"{name}", loc="left", fontsize=10, color=INK)
        style(ax)
    axs[0].set_ylabel("PSD of the estimate (m^2/Hz)")
    axs[0].legend(fontsize=7, frameon=False, loc="upper right")
    fig.suptitle("Still-water displacement PSD, mean over seeds (true displacement is zero)", x=0.01, ha="left",
                 fontsize=12, color=INK)
    save(fig, svgd, out / "figures" / "fig02_calm_psd.svg")
    plt.close(fig)


def fig_excitation(r, plt, svgd, out):
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    pts = defaultdict(lambda: ([], []))
    for (e, sc), d in r.g.items():
        if e not in ("E1", "E2", "E3"):
            continue
        for f in ["A41@A", "A41@B", "A41@C"]:
            if f not in d:
                continue
            for x in d[f].values():
                a = fnum(x, "true_az_rms")
                s = fnum(x, "exc_sig_med")
                if e == "E1":
                    a = 3e-4   # still water, drawn at the left edge
                if np.isfinite(a) and np.isfinite(s):
                    pts[f][0].append(a)
                    pts[f][1].append(s)
    for f in ["A41@A", "A41@B", "A41@C"]:
        if pts[f][0]:
            ax.scatter(pts[f][0], np.maximum(pts[f][1], 1e-6), s=10, color=C[f], alpha=0.6, edgecolors="none",
                       label=LABEL[f], marker="o" if f != "A41@A" else "s")
    xx = np.logspace(-3.6, 0.5, 20)
    ax.plot(xx, xx, color=INK2, lw=0.8, ls=":")
    ax.text(1.2, 1.5, "1:1", fontsize=8, color=INK2)
    ax.axhline(0.05, color=C["A41@A"], lw=0.6, ls="--")
    ax.text(4e-4, 0.06, "legacy floor 0.05", fontsize=8, color=INK2)
    ax.axhline(1e-5, color=C["A41@C"], lw=0.6, ls="--")
    ax.text(4e-4, 1.3e-5, "numerical floor 1e-5", fontsize=8, color=INK2)
    ax.axvline(3e-4, color="#d6d6d1", lw=6, zorder=0)
    ax.text(3.3e-4, 2e-6, "still water", fontsize=8, color=INK2, rotation=90)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(2e-4, 4)
    ax.set_ylim(1e-6, 4)
    ax.set_xlabel("true vertical wave acceleration RMS (m/s^2)")
    ax.set_ylabel("CARMA driving sigma, median (m/s^2)")
    ax.legend(fontsize=8, frameon=False, loc="lower right")
    ax.set_title("CARMA excitation against the true wave energy (E1-E3, every run)", loc="left", fontsize=10,
                 color=INK)
    style(ax)
    save(fig, svgd, out / "figures" / "fig03_excitation_vs_energy.svg")
    plt.close(fig)


def envelope(t, sc):
    import carma_excitation_study as S
    e = np.ones_like(t)
    lvl = 1.0
    steps = S.envelope_steps(sc)
    if steps and steps[0][2] == 0.0:
        lvl = 0.0
    e[:] = lvl
    for (tc, w, lf, lt) in steps:
        x = np.clip((t - (tc - 0.5 * w)) / w, 0, 1)
        sm = x ** 3 * (10 - 15 * x + 6 * x * x)
        m = t >= tc - 0.5 * w
        e[m] = lf + (lt - lf) * sm[m]
    return e


def fig_detection(r, plt, svgd, out):
    tr = out / "traces" / "repeated.csv"
    if not tr.exists():
        return
    import carma_excitation_study as S
    sc = [x for x in S.TRACES if x[0] == "repeated"][0][1]
    t = load_trace(tr)
    fig, axs = plt.subplots(3, 1, figsize=(10, 6.6), sharex=True)
    env = envelope(t["t"], sc)
    for ax in axs:
        ax.fill_between(t["t"], 0, 1, where=env > 0.01, color="#eeeeea", transform=ax.get_xaxis_transform(), lw=0)
    axs[0].plot(t["t"], t["pw_A41@C"], color=C["A41@C"], lw=1.4, label="evidence weight p")
    axs[0].plot(t["t"], t["cf_A41@C"], color=C["A41@Cnf"], lw=1.0, label="period confidence c_f")
    axs[0].set_ylabel("weight")
    axs[0].set_ylim(-0.05, 1.1)
    axs[0].legend(fontsize=8, frameon=False, loc="center right")
    axs[1].plot(t["t"], np.clip(t["zw_A41@C"], 1e-2, None), color=C["A41@C"], lw=1.2)
    axs[1].axhline(2, color=INK2, lw=0.7, ls="--")
    axs[1].axhline(4, color=INK2, lw=0.7, ls="--")
    axs[1].text(5, 2.3, "z0 = 2", fontsize=7, color=INK2)
    axs[1].text(5, 4.6, "z1 = 4", fontsize=7, color=INK2)
    axs[1].set_yscale("log")
    axs[1].set_ylabel("detection statistic z")
    for f in ["A41@A", "A41@B", "A41@C"]:
        axs[2].plot(t["t"], np.maximum(t[f"fop_{f}"], 1e-3), color=C[f], ls=LS.get(f, "-"), lw=1.2, label=LABEL[f])
    axs[2].set_ylabel("CARMA frequency (Hz)")
    axs[2].set_xlabel("time (s); shaded: waves present (J29, H_s 0.5 m)")
    axs[2].legend(fontsize=8, frameon=False, loc="upper right")
    for ax in axs:
        style(ax)
    fig.suptitle("Wave-detection confidence through repeated calm / wave / calm transitions, seed 1", x=0.01,
                 ha="left", fontsize=12, color=INK)
    save(fig, svgd, out / "figures" / "fig04_detection_confidence.svg")
    plt.close(fig)


def fig_variance(r, plt, svgd, out):
    tr = out / "traces" / "repeated.csv"
    if not tr.exists():
        return
    t = load_trace(tr)
    fig, ax = plt.subplots(figsize=(10, 3.8))
    ax.plot(t["t"], t["var_A41@C"], color=INK, lw=1.2, label="measured band variance (tuner)")
    ax.plot(t["t"], t["var0_A41@C"], color=C["A41@B0"], lw=1.6, label="propagated sensor-noise variance E[V | 0]")
    ax.plot(t["t"], np.maximum(t["sig_A41@C"] ** 2, 1e-12), color=C["A41@C"], lw=1.4,
            label="evidence-mode CARMA sigma^2 (p max(V - E[V|0], 0), smoothed)")
    ax.plot(t["t"], np.maximum(t["sig_A41"] ** 2, 1e-12), color=C["A41@A"], ls="--", lw=1.2,
            label="legacy CARMA sigma^2")
    ax.set_yscale("log")
    ax.set_ylim(1e-11, 3)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("(m/s^2)^2")
    ax.legend(fontsize=7, frameon=False, loc="lower left", ncol=2)
    ax.set_title("Variance channel: measured, sensor-noise and physical wave variance, seed 1", loc="left",
                 fontsize=10, color=INK)
    style(ax)
    save(fig, svgd, out / "figures" / "fig05_variance_channel.svg")
    plt.close(fig)


def fig_hs(r, plt, svgd, out):
    fam = "jonswap_fp0.29_hs{}_g3.3"
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
    hs_list = [0.0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.8]
    for ax, key, name in ((axs[0], lambda x: fnum(x, "rms_z"), "heave"), (axs[1], rms_h, "horizontal")):
        series = defaultdict(lambda: ([], [], []))
        truth = []
        for hs in hs_list:
            if hs == 0.0:
                d = r.g.get(("E1", "stationary+T2400"))
            elif hs == 0.8:
                d = r.g.get(("E7", fam.format("0.8")))
            else:
                d = r.g.get(("E2", fam.format(f"{hs:g}")))
            if not d:
                continue
            xpos = hs if hs > 0 else 0.002
            for f in ["A41@A", "T41@A", "OU3", "A41@B0", "A41@B", "A41@C", "T41@C"]:
                if f not in d:
                    continue
                v = np.array([key(x) for x in d[f].values()])
                lo, hi = bootstrap_ci(v)
                series[f][0].append(xpos)
                series[f][1].append(1e3 * np.mean(v))
                series[f][2].append((1e3 * (np.mean(v) - lo), 1e3 * (hi - np.mean(v))))
            if name == "heave":
                tv = np.array([fnum(x, "std_z") for x in d[next(iter(d))].values()])
            else:
                tv = np.array([np.hypot(fnum(x, "std_x"), fnum(x, "std_y")) for x in d[next(iter(d))].values()])
            truth.append((xpos, 1e3 * np.mean(tv)))
        for f in ["A41@A", "T41@A", "OU3", "A41@B0", "A41@B", "A41@C", "T41@C"]:
            x, y, e = series[f]
            if x:
                e = np.array(e).T
                ax.errorbar(x, y, yerr=e, color=C[f], ls=LS.get(f, "-"), lw=1.4, marker="o", ms=4, capsize=2,
                            label=LABEL[f])
        if truth:
            tx, ty = zip(*truth)
            ax.plot(tx, np.maximum(ty, 1e-3), color=INK, lw=1.0, ls="-.", label="zero-output error (true RMS)")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xticks([0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.8])
        ax.set_xticklabels(["0", "0.005", "0.01", "0.02", "0.05", "0.1", "0.2", "0.8"])
        ax.set_xlabel("H_s (m), JONSWAP f_p 0.29 Hz; 0 = still IMU")
        ax.set_ylabel(f"{name} RMS error (mm)")
        ax.set_ylim(1e-2, 2e2)
        ax.set_title(name, loc="left", fontsize=10, color=INK)
        style(ax)
    axs[1].legend(fontsize=7, frameon=False, loc="lower right")
    fig.suptitle("RMS displacement error against significant wave height (mean, 95 % CI over seeds)", x=0.01,
                 ha="left", fontsize=12, color=INK)
    save(fig, svgd, out / "figures" / "fig06_error_vs_hs.svg")
    plt.close(fig)


def fig_transitions(r, plt, svgd, out):
    tr = out / "traces" / "repeated.csv"
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.0), gridspec_kw={"width_ratios": [1.5, 1]})
    if tr.exists():
        t = load_trace(tr)
        dt = np.median(np.diff(t["t"]))
        n = max(1, int(round(10.0 / dt)))
        k = np.ones(n) / n
        for f in ["A41@A", "A41@B", "A41@C"]:
            e2 = (t[f"z_{f}"] - t["z_true"]) ** 2
            e2 = np.where(np.isfinite(e2), e2, 0.0)
            axs[0].plot(t["t"], 1e3 * np.sqrt(np.convolve(e2, k, mode="same")), color=C[f], ls=LS.get(f, "-"),
                        lw=1.2, label=LABEL[f])
        axs[0].plot(t["t"], 1e3 * np.sqrt(np.convolve(t["z_true"] ** 2, k, mode="same")), color=INK, lw=0.8,
                    ls="-.", label="true heave (10 s RMS)")
        axs[0].set_yscale("log")
        axs[0].set_ylim(1e-2, 1e3)
        axs[0].set_xlabel("time (s)")
        axs[0].set_ylabel("heave error, 10 s RMS (mm)")
        axs[0].legend(fontsize=7, frameon=False, loc="upper left", ncol=2)
        axs[0].set_title("Repeated transitions, seed 1", loc="left", fontsize=10, color=INK)
    style(axs[0])
    rows = [t for t in r.trans if t["filter"] == "A41@C" and t["kind"] == "onset"]
    groups = defaultdict(list)
    for t in rows:
        lab = (t["scenario"].replace("jonswap_fp0.29_", "").replace("_g3.3", "").split("+T")[0].split("+")[0]
               + f" w{float(t['width']):g}" + f" {float(t['level_from']):g}->{float(t['level_to']):g}")
        groups[lab].append(fnum(t, "detect_s"))
    labs = sorted(groups)
    if labs:
        data = [np.array([x for x in groups[l] if np.isfinite(x)]) for l in labs]
        axs[1].boxplot(data, vert=False, widths=0.5, medianprops={"color": C["A41@C"]},
                       flierprops={"markersize": 3})
        axs[1].set_yticks(range(1, len(labs) + 1))
        axs[1].set_yticklabels(labs, fontsize=7)
        axs[1].set_xlabel("detection delay after the onset starts (s)")
        axs[1].set_title("Evidence mode: detection delay, all seeds", loc="left", fontsize=10, color=INK)
    style(axs[1])
    fig.suptitle("Detection and recovery during sea-state transitions", x=0.01, ha="left", fontsize=12, color=INK)
    save(fig, svgd, out / "figures" / "fig07_transitions.svg")
    plt.close(fig)


def fig_swell(r, plt, svgd, out):
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.0))
    by = defaultdict(dict)
    for (e, sc), d in r.g.items():
        if e == "E3":
            by[hs_of(sc)][1.0 / fp_of(sc)] = d
    marks = {0.05: "o", 0.2: "s", 0.5: "^"}
    for hs in sorted(by):
        Ts = sorted(by[hs])
        for f in ["A41@A", "A41@B", "A41@C"]:
            err, amp = [], []
            for T in Ts:
                d = by[hs][T]
                if f not in d:
                    err.append(np.nan)
                    amp.append(np.nan)
                    continue
                v = np.array([fnum(x, "rms_z") / fnum(x, "std_z") for x in d[f].values()])
                a = np.array([fnum(x, "exc_pz_rms") / fnum(x, "std_z") for x in d[f].values()])
                err.append(100 * np.nanmean(v))
                amp.append(np.nanmean(a))
            kw = dict(color=C[f], ls=LS.get(f, "-"), marker=marks.get(hs, "o"), ms=4, lw=1.2)
            axs[0].plot(Ts, err, label=f"{LABEL[f]}, H_s {hs:g} m", **kw)
            axs[1].plot(Ts, amp, **kw)
    axs[0].axhline(100, color=INK, lw=0.8, ls="-.")
    axs[0].text(8.2, 110, "zero output", fontsize=8, color=INK2)
    axs[0].set_yscale("log")
    axs[0].set_xlabel("peak period (s)")
    axs[0].set_ylabel("heave RMS error (% of true heave RMS)")
    axs[0].legend(fontsize=6, frameon=False, ncol=1, loc="upper left")
    axs[1].axhline(1.0, color=INK, lw=0.8, ls="-.")
    axs[1].set_xlabel("peak period (s)")
    axs[1].set_ylabel("estimated / true heave RMS")
    axs[1].set_ylim(0, 3)
    for ax in axs:
        style(ax)
    fig.suptitle("Long-period swell (JONSWAP gamma 5): error and amplitude preservation", x=0.01, ha="left",
                 fontsize=12, color=INK)
    save(fig, svgd, out / "figures" / "fig08_swell.svg")
    plt.close(fig)


def fig_noise(r, plt, svgd, out):
    conds = ["stationary", "jonswap_fp0.29_hs0.05_g3.3", "jonswap_fp0.29_hs0.8_g3.3", "jonswap_fp0.08_hs0.2_g5"]
    names = ["still IMU", "J29 H_s 0.05 m", "J29 H_s 0.8 m", "swell 12.5 s, 0.2 m"]
    fl = ["A41@A", "A41@B/nv0.5", "A41@B", "A41@B/nv2", "A41@C/nv0.5", "A41@C", "A41@C/nv2"]
    col = {"A41@A": C["A41@A"], "A41@B/nv0.5": C["A41@B"], "A41@B": C["A41@B"], "A41@B/nv2": C["A41@B"],
           "A41@C/nv0.5": C["A41@C"], "A41@C": C["A41@C"], "A41@C/nv2": C["A41@C"]}
    hatch = {"/nv0.5": "////", "/nv2": "\\\\\\\\"}
    fig, axs = plt.subplots(1, 4, figsize=(13, 3.8))
    for ax, sc, nm in zip(axs, conds, names):
        d = r.g.get(("E5", sc))
        if not d:
            continue
        f2 = [f for f in fl if f in d]
        seeds, vz = aligned(d, f2, lambda x: fnum(x, "rms_z"))
        for i, f in enumerate(f2):
            m = 1e3 * np.mean(vz[f])
            lo, hi = 1e3 * np.array(bootstrap_ci(vz[f]))
            h = next((v for k, v in hatch.items() if f.endswith(k)), None)
            ax.bar(i, max(m, 1e-3), color=col[f], hatch=h, edgecolor="#fcfcfb", lw=0.5, width=0.8)
            ax.errorbar(i, max(m, 1e-3), yerr=[[max(m - lo, 0)], [max(hi - m, 0)]], color=INK, lw=0.8, capsize=2)
        ax.set_xticks(range(len(f2)))
        ax.set_xticklabels([f.replace("A41@", "").replace("A41@A", "A") for f in f2], fontsize=7, rotation=45)
        ax.set_yscale("log")
        ax.set_title(nm, loc="left", fontsize=10, color=INK)
        style(ax)
    axs[0].set_ylabel("heave RMS error (mm)")
    fig.suptitle("Sensor-noise mis-calibration: model variance x0.5 (////) and x2 (\\\\\\\\), 10 seeds", x=0.01,
                 ha="left", fontsize=12, color=INK)
    save(fig, svgd, out / "figures" / "fig09_noise_miscalibration.svg")
    plt.close(fig)


def fig_tfg(r, plt, svgd, out):
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.2))
    for ax, key, nm in ((axs[0], lambda x: fnum(x, "rms_z"), "heave"), (axs[1], rms_h, "horizontal")):
        for (cand, base, col) in (("T41@C", "A41@C", C["T41@C"]), ("T41@A", "A41@A", C["T41@A"])):
            xs, ys = [], []
            for (e, sc), d in r.g.items():
                if cand in d and base in d:
                    for s in set(d[cand]) & set(d[base]):
                        xs.append(key(d[base][s]))
                        ys.append(key(d[cand][s]))
            if xs:
                ax.scatter(1e3 * np.array(xs), 1e3 * np.array(ys), s=7, color=col, alpha=0.5, edgecolors="none",
                           label=f"{LABEL[cand]} vs {LABEL[base]} ({len(xs)} runs)")
        lim = [1e-3, 1e3]
        ax.plot(lim, lim, color=INK2, lw=0.8, ls=":")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(*lim)
        ax.set_ylim(*lim)
        ax.set_xlabel(f"quaternion MEKF {nm} RMS (mm)")
        ax.set_ylabel(f"two-frame Lie group {nm} RMS (mm)")
        ax.set_title(nm, loc="left", fontsize=10, color=INK)
        style(ax)
    axs[0].legend(fontsize=7, frameon=False, loc="upper left")
    fig.suptitle("Q-MEKF CARMA3 against TFG-CARMA3, same excitation mode, every paired run", x=0.01, ha="left",
                 fontsize=12, color=INK)
    save(fig, svgd, out / "figures" / "fig10_qmekf_vs_tfg.svg")
    plt.close(fig)


def all_figures(r, plt, svgd, out):
    for fn in (fig_calm, fig_psd, fig_excitation, fig_detection, fig_variance, fig_hs, fig_transitions, fig_swell,
               fig_noise, fig_tfg):
        fn(r, plt, svgd, out)
