#!/usr/bin/env python3
"""Heave over the last 30 s: literature baselines against the OU families.

Reads the time series the simulators write with W3D_WRITE_TIMESERIES=1:
  tests/kalman_ou_iii/*_fusion_ou3.csv, tests/kalman_ou_ii/*_fusion_ou2.csv,
  tests/kalman_tfg/*_fusion_tfg.csv and tests/heave_baselines/*_heave_*.csv
(the latter for pii, godhavn, godhavn_ext, richter_zd and kuchler, per
levelling front end).  All share one noise
realisation per record, so the traces are directly comparable.

For each JONSWAP / PM-Stokes record at the low, medium and high heights it
writes heave_baselines_<wave>_<group>.svg (reference and every estimate,
plus the error) and <method>_<wave>_<group>_zkin.svg for each baseline
(heave and heave rate against the reference); --png adds PNG previews.
"""

import argparse
import re
import sys
from pathlib import Path

import matplotlib as mpl
import numpy as np
import pandas as pd

sys.path.append(str(Path(__file__).resolve().parents[1]))
from svg_determinism import configure_svg, save_svg  # noqa: E402

mpl.use("Agg")
configure_svg(mpl)
import matplotlib.pyplot as plt  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
PLOT_WINDOW_S = 30.0
PLOT_STEP = 4  # 200 Hz -> 50 Hz
GROUPS = {"low": 0.27, "medium": 1.50, "high": 8.50}
WAVES = ("jonswap", "pmstokes")

SOURCES = {
    "ou3": ("kalman_ou_iii", "_fusion_ou3.csv", "OU-III"),
    "ou2": ("kalman_ou_ii", "_fusion_ou2.csv", "OU-II"),
    "tfg": ("kalman_tfg", "_fusion_tfg.csv", "TFG"),
    "pii": ("heave_baselines", "_heave_pii.csv", "PII"),
    "godhavn": ("heave_baselines", "_heave_godhavn.csv", "Godhavn (published cutoff law)"),
    "godhavn_ext": ("heave_baselines", "_heave_godhavn_ext.csv", "Godhavn + bias term"),
    "richter_zd": ("heave_baselines", "_heave_richter_zd.csv", "Richter zero-displacement"),
    "kuchler": ("heave_baselines", "_heave_kuchler.csv", "Küchler EKF"),
}
BASELINES = ("godhavn", "godhavn_ext", "richter_zd", "kuchler", "pii")
FRONTENDS = {
    "mahony": ("", "baselines levelled by the shipped Mahony"),
    "mahony_slow": ("_mahony_slow", "baselines levelled by the slow IMU-only Mahony"),
    "truth": ("_truth", "baselines levelled with the true attitude"),
}
NAME = re.compile(r"w3d_(?P<wave>[a-z]+)_H(?P<h>[0-9.]+)_")


FRONTEND = "mahony"


def find(tests_dir: Path, method: str, wave: str, height: float):
    sub, suffix, _ = SOURCES[method]
    if method in BASELINES:
        suffix = suffix.replace(".csv", FRONTENDS[FRONTEND][0] + ".csv")
    # The glob ends in the full suffix, so *_heave_godhavn.csv never matches
    # *_heave_godhavn_ext.csv and the Mahony series never match the others.
    for path in sorted((tests_dir / sub).glob(f"w3d_{wave}_H*{suffix}")):
        m = NAME.match(path.name)
        if m and abs(float(m.group("h")) - height) < 1e-6:
            return path
    return None


def tail(path: Path) -> pd.DataFrame:
    cols = ["time", "disp_ref_z", "disp_est_z", "vel_ref_z", "vel_est_z"]
    df = pd.read_csv(path, usecols=cols)
    df = df[df["time"] >= df["time"].max() - PLOT_WINDOW_S]
    return df.iloc[::PLOT_STEP].reset_index(drop=True)


def rms(x) -> float:
    return float(np.sqrt(np.mean(np.square(x))))


WRITE_PNG = False


def save(fig, base: Path):
    fig.tight_layout()
    save_svg(fig, base.with_suffix(".svg"), bbox_inches="tight")
    if WRITE_PNG:
        fig.savefig(base.with_suffix(".png"), dpi=110, bbox_inches="tight")
    plt.close(fig)


def overview(dfs: dict, wave: str, group: str, height: float, out: Path):
    t = dfs["ou3"]["time"].to_numpy() if "ou3" in dfs else next(iter(dfs.values()))["time"].to_numpy()
    ref = next(iter(dfs.values()))
    ref_z = np.interp(t, ref["time"], ref["disp_ref_z"])
    est = {m: np.interp(t, df["time"], df["disp_est_z"]) for m, df in dfs.items()}

    def label(m):
        return f"{SOURCES[m][2]} ({rms(est[m] - ref_z) * 100:.1f} cm)"

    fig, axes = plt.subplots(4, 1, figsize=(10, 10.5), sharex=True)
    levelling = FRONTENDS[FRONTEND][1]
    fig.suptitle(f"{wave.upper()} Hs = {height:g} m: heave, last {PLOT_WINDOW_S:.0f} s; "
                 f"{levelling} (window RMS error in legend)")
    panels = (
        (axes[0], ("ou3", "ou2", "tfg", "pii"), "OU families and PII" if "pii" in est else "OU families"),
        (axes[1], ("richter_zd", "godhavn_ext", "kuchler"), "Literature baselines"),
        (axes[2], ("godhavn",), "Godhavn, published cutoff law (own scale)"),
    )
    styles = {"ou3": "-", "ou2": "-.", "tfg": (0, (3, 1, 1, 1)), "pii": ":",
              "godhavn": "--", "godhavn_ext": "--", "richter_zd": "-", "kuchler": "-."}
    for ax, methods, title in panels:
        ax.plot(t, ref_z, color="black", linewidth=1.6, label="Reference")
        for m in methods:
            if m in est:
                ax.plot(t, est[m], linewidth=1.1, linestyle=styles[m], label=label(m))
        ax.set_title(title, fontsize=9, loc="left")
        ax.set_ylabel("Heave [m]")
        ax.grid(True)
        ax.legend(loc="upper right", fontsize=7, ncol=2)
    for m in ("ou3", "ou2", "tfg", "pii", "richter_zd", "godhavn_ext", "kuchler"):
        if m in est:
            axes[3].plot(t, est[m] - ref_z, linewidth=0.9, linestyle=styles[m], label=SOURCES[m][2])
    axes[3].set_title("Error (estimate - reference); published-law Godhavn omitted", fontsize=9, loc="left")
    axes[3].set_ylabel("Error [m]")
    axes[3].set_xlabel("Time [s]")
    axes[3].grid(True)
    axes[3].legend(loc="upper right", fontsize=7, ncol=3)
    save(fig, out / f"heave_baselines{FRONTENDS[FRONTEND][0]}_{wave}_{group}")


def zkin(df: pd.DataFrame, method: str, wave: str, group: str, height: float, out: Path):
    fig, axes = plt.subplots(2, 1, figsize=(10, 5.0), sharex=True)
    fig.suptitle(f"{SOURCES[method][2]}: {wave.upper()} Hs = {height:g} m (Z-axis, last "
                 f"{PLOT_WINDOW_S:.0f} s)")
    for ax, prefix, ylabel in ((axes[0], "disp", "Disp Z [m]"), (axes[1], "vel", "Vel Z [m/s]")):
        ax.plot(df["time"], df[f"{prefix}_ref_z"], label="Ref")
        ax.plot(df["time"], df[f"{prefix}_est_z"], label="Est", linestyle="--")
        ax.set_ylabel(ylabel)
        ax.grid(True)
        ax.legend(loc="upper right", fontsize=8)
    axes[-1].set_xlabel("Time [s]")
    save(fig, out / f"{method}{FRONTENDS[FRONTEND][0]}_{wave}_{group}_zkin")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--tests-dir", type=Path, default=REPO / "tests")
    ap.add_argument("--output-dir", type=Path, default=Path("."))
    ap.add_argument("--png", action="store_true", help="also write PNG previews")
    ap.add_argument("--frontend", choices=tuple(FRONTENDS), default="mahony",
                    help="which levelling the baseline series used")
    args = ap.parse_args()
    global WRITE_PNG, FRONTEND
    WRITE_PNG = args.png
    FRONTEND = args.frontend
    args.output_dir.mkdir(parents=True, exist_ok=True)

    made = 0
    for wave in WAVES:
        for group, height in GROUPS.items():
            dfs = {}
            for m in SOURCES:
                path = find(args.tests_dir, m, wave, height)
                if path is not None:
                    dfs[m] = tail(path)
            if not any(m in dfs for m in BASELINES):
                print(f"skip {wave} {group}: no baseline time series")
                continue
            overview(dfs, wave, group, height, args.output_dir)
            for m in BASELINES:
                if m in dfs:
                    zkin(dfs[m], m, wave, group, height, args.output_dir)
            print(f"{wave} {group}: {', '.join(sorted(dfs))}")
            made += 1
    if made == 0:
        sys.exit("no time series found; run the simulators with W3D_WRITE_TIMESERIES=1")


if __name__ == "__main__":
    main()
