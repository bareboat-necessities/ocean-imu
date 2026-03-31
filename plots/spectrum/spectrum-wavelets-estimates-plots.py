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
    r"^spectrum_wavelets_"
    r"(?P<wtype>[A-Za-z0-9]+)_"
    r"H(?P<H>[-0-9.]+)"
    r"(?:_L(?P<L>[-0-9.]+))?"
    r"(?:_A(?P<A>[-0-9.]+))?"
    r"(?:_P(?P<P>[-0-9.]+))?"
    r"_N(?P<noise>[-0-9.]+)"
    r"_B(?P<bias>[-0-9.]+)"
    r"\.csv$"
)


def parse_filename(fname):
    m = fname_re.match(os.path.basename(fname))
    if not m:
        return None
    return {
        "wtype": m.group("wtype"),
        "H": float(m.group("H")),
        "A": float(m.group("A")) if m.group("A") is not None else 0.0,
    }


def save_all(fig, base):
    fig.savefig(f"{base}.pgf", bbox_inches="tight", backend="pgf")
    fig.savefig(f"{base}.svg", bbox_inches="tight", dpi=150)
    with mpl.rc_context({"text.usetex": False}):
        fig.savefig(f"{base}.png", bbox_inches="tight", dpi=300)


def make_plots(fname, group, meta):
    df = pd.read_csv(fname)
    if df.empty:
        return

    base = f"spectrum_wavelets_estimates_{meta['wtype']}_{group}"
    title = fr"{meta['wtype'].capitalize()} wavelet estimator vs reference ($H_s={meta['H']:.2f}\,\mathrm{{m}}$)"

    fig, ax = plt.subplots()
    ax.plot(df["freq_hz"], df["S_eta_hz"], label="Wavelet estimator", linewidth=2)
    ax.plot(df["freq_hz"], df["S_ref_interp"], label="Reference", linewidth=2, linestyle="--")
    ax.set_xlabel(r"Frequency $f$ [Hz]")
    ax.set_ylabel(r"$S_{\eta}(f)$")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend()
    save_all(fig, f"{base}_spectrum")
    plt.close(fig)

    fig2, ax2 = plt.subplots()
    ax2.plot(df["freq_hz"], df["CumVar_est"], label=r"Estimated cumulative variance", linewidth=2)
    ax2.plot(df["freq_hz"], df["CumVar_ref"], label=r"Reference cumulative variance", linewidth=2, linestyle="--")
    ax2.set_xlabel("Frequency [Hz]")
    ax2.set_ylabel("Cumulative variance")
    ax2.set_title(title)
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    save_all(fig2, f"{base}_timeseries")
    plt.close(fig2)


if __name__ == "__main__":
    files = glob.glob("spectrum_wavelets_*.csv")
    if not files:
        print("No spectrum_wavelets_*.csv files found.")
        raise SystemExit(0)

    for fname in sorted(files):
        meta = parse_filename(fname)
        if not meta:
            continue
        for group, target_H in height_groups.items():
            if abs(meta["H"] - target_H) < 1e-3:
                print(f"Processing {fname} as {group} ...")
                make_plots(fname, group, meta)
