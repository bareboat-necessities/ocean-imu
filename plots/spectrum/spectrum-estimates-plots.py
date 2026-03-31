#!/usr/bin/env python3
import glob
import os
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pgf import FigureCanvasPgf

mpl.backend_bases.register_backend('pgf', FigureCanvasPgf)
plt.rcParams.update({
    "pgf.texsystem": "xelatex",
    "font.family": "serif",
    "text.usetex": True,
    "pgf.rcfonts": False,
    "pgf.preamble": "\n".join([
        r"\usepackage{fontspec}",
        r"\usepackage{unicode-math}",
        r"\usepackage{amsmath}",
        r"\setmainfont{DejaVu Serif}",
        r"\setmathfont{Latin Modern Math}",
        r"\providecommand{\mathdefault}[1]{#1}",
    ]),
})

height_groups = {
    "low": 0.27,
    "medium": 1.50,
    "high": 8.50,
}

fname_re = re.compile(
    r"^spectrum_estimate_(?P<wtype>[a-z]+)_H(?P<H>[0-9.]+)_L(?P<L>[0-9.]+)_A(?P<A>[-0-9.]+)_P(?P<P>[-0-9.]+)\.csv$"
)


def parse_filename(fname):
    m = fname_re.match(os.path.basename(fname))
    if not m:
        return None
    return {
        "wtype": m.group("wtype"),
        "H": float(m.group("H")),
        "A": float(m.group("A")),
    }


def save_all(fig, base):
    fig.savefig(f"{base}.pgf", bbox_inches="tight", backend="pgf")
    fig.savefig(f"{base}.svg", bbox_inches="tight", dpi=150)
    with mpl.rc_context({"text.usetex": False}):
        fig.savefig(f"{base}.png", bbox_inches="tight", dpi=300)


def spectrum_cols(df):
    freq_map = {}
    for c in df.columns:
        if c.startswith("S_eta_est_f") and "_Hz=" in c:
            idx = int(c.split("_f")[1].split("_Hz=")[0])
            freq_map[idx] = float(c.split("_Hz=")[1])
    idxs = sorted(freq_map.keys())
    est_cols = [f"S_eta_est_f{i}_Hz={freq_map[i]}" for i in idxs]
    ref_cols = [f"S_eta_ref_f{i}" for i in idxs]
    freqs = np.array([freq_map[i] for i in idxs])
    return freqs, est_cols, ref_cols


def make_plots(fname, group, meta):
    df = pd.read_csv(fname)
    if df.empty:
        return

    freqs, est_cols, ref_cols = spectrum_cols(df)
    last = df.iloc[-1]

    base = f"spectrum_estimates_{meta['wtype']}_{group}"
    title = fr"{meta['wtype'].capitalize()} estimator vs reference ($H_s={meta['H']:.2f}\,\mathrm{{m}}$)"

    fig, ax = plt.subplots()
    ax.plot(freqs, last[est_cols].to_numpy(dtype=float), label="Estimator", linewidth=2)
    ax.plot(freqs, last[ref_cols].to_numpy(dtype=float), label="Reference", linewidth=2, linestyle="--")
    ax.set_xlabel(r"Frequency $f$ [Hz]")
    ax.set_ylabel(r"$S_{\eta}(f)$")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend()
    save_all(fig, f"{base}_spectrum")
    plt.close(fig)

    fig2, ax2 = plt.subplots()
    ax2.plot(df["time"], df["Hs"], label=r"$H_s$", linewidth=2)
    ax2.plot(df["time"], 1.0 / np.maximum(df["Fp"], 1e-9), label=r"$T_p=1/f_p$", linewidth=2)
    ax2.set_xlabel("Time [s]")
    ax2.set_ylabel("Estimate")
    ax2.set_title(title)
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    save_all(fig2, f"{base}_timeseries")
    plt.close(fig2)


if __name__ == "__main__":
    files = glob.glob("spectrum_estimate_*.csv")
    if not files:
        print("No spectrum_estimate_*.csv files found.")
        raise SystemExit(0)

    for fname in sorted(files):
        meta = parse_filename(fname)
        if not meta:
            continue
        for group, target_H in height_groups.items():
            if abs(meta["H"] - target_H) < 1e-3:
                print(f"Processing {fname} as {group} ...")
                make_plots(fname, group, meta)
