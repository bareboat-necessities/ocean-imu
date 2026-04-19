#!/usr/bin/env python3
from pathlib import Path
import csv
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.use("pgf")
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

SCRIPT_DIR = Path(__file__).resolve().parent
TEST_DIR = SCRIPT_DIR.parent.parent / "tests" / "wave_sim"
CSV_PATH = TEST_DIR / "vertical_accel_drift_timeseries.csv"


def short_name(name: str) -> str:
    stem = name.replace("wave_data_", "").replace(".csv", "")
    return stem.replace("_", r"\_")


def read_rows():
    with CSV_PATH.open(newline="") as handle:
        return list(csv.DictReader(handle))


def unique_files(rows, wave_type: str):
    names = sorted({r["wave_file"] for r in rows if r["wave_type"] == wave_type})
    return names


def float_col(rows, key):
    return [float(r[key]) for r in rows]


def make_plot(rows, wave_type: str, out_name: str) -> None:
    files = unique_files(rows, wave_type)
    if not files:
        print(f"no rows for {wave_type}")
        return

    n = len(files)
    fig, axes = plt.subplots(n, 1, figsize=(14, max(3.3 * n, 5.0)), sharex=False)
    if n == 1:
        axes = [axes]

    for ax, file in zip(axes, files):
        case = [r for r in rows if r["wave_type"] == wave_type and r["wave_file"] == file]
        t = float_col(case, "time_s")
        ref = float_col(case, "disp_ref_z_m")
        dint = float_col(case, "disp_int_z_m")
        ddet = float_col(case, "disp_int_detrended_z_m")

        ax.plot(t, ref, color="black", linewidth=1.2, linestyle="--", label="reference $z$")
        ax.plot(t, dint, color="#d62728", linewidth=1.0, label="double-integrated noisy $a_z$")
        ax.plot(t, ddet, color="#1f77b4", linewidth=1.0, label="double-integrated + detrender")
        ax.set_title(short_name(file), fontsize=10)
        ax.set_ylabel("z [m]")
        ax.grid(True, alpha=0.35)

    axes[-1].set_xlabel("time [s]")
    axes[0].legend(ncol=3, fontsize=8, loc="upper center")
    fig.tight_layout()

    pgf_path = SCRIPT_DIR / out_name
    svg_path = pgf_path.with_suffix(".svg")
    fig.savefig(pgf_path, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {pgf_path}")


def main() -> None:
    rows = read_rows()
    make_plot(rows, "pmstokes", "pmstokes_vertical_accel_drift.pgf")
    make_plot(rows, "jonswap", "jonswap_vertical_accel_drift.pgf")


if __name__ == "__main__":
    main()
