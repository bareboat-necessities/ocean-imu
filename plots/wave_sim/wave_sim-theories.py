#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy import special


AREA_LABELS = [
    "Beyond breaking / invalid",
    "Linear",
    "Stokes II",
    "Stokes III",
    "Stokes IV",
    "Stokes V",
    "Cnoidal",
    "Cnoidal III",
    "Cnoidal V",
    "Solitary-wave side",
    "Stream-function / numerical",
]


def configure_pgf():
    mpl.rcParams.update({
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


def longest_path(paths):
    if not paths:
        return None
    return max(paths, key=lambda p: (np.nanmax(p[:, 0]) - np.nanmin(p[:, 0]), len(p)))


def contour_level_segments(xgrid2d, ygrid2d, zgrid2d, level=0.01):
    fig, ax = plt.subplots()
    cs = ax.contour(xgrid2d, ygrid2d, zgrid2d, levels=[level])
    segs = [np.asarray(seg) for seg in cs.allsegs[0] if len(seg) > 1]
    plt.close(fig)
    return segs


def interp_log_curve(x, y, x_new):
    mask = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    x = np.asarray(x)[mask]
    y = np.asarray(y)[mask]
    if len(x) < 2:
        return np.full_like(x_new, np.nan, dtype=float)

    order = np.argsort(x)
    x = x[order]
    y = y[order]

    xu, idx = np.unique(x, return_index=True)
    yu = y[idx]
    if len(xu) < 2:
        return np.full_like(x_new, np.nan, dtype=float)

    out = np.full_like(x_new, np.nan, dtype=float)
    valid = (x_new >= xu.min()) & (x_new <= xu.max())
    out[valid] = 10 ** np.interp(np.log10(x_new[valid]), np.log10(xu), np.log10(yu))
    return out


def fenton_cnoidal_curve(m, h=10.0, hvals=None):
    """
    Return x=h/L, y=H/L curve for a chosen elliptic modulus m using the same
    Fenton (1999) series form already present in the user's MATLAB code.
    """
    if hvals is None:
        hvals = np.arange(0.01, 7.8 + 1e-12, 0.001)

    Km, Em = special.ellipk(m), special.ellipe(m)
    em = Em / Km
    Hoverh = hvals / h
    L = (
        4 * Km / np.sqrt(3 * Hoverh)
        * (
            1
            + Hoverh * (5 / 8 - 3 / 2 * em)
            + Hoverh**2 * (-21 / 128 + 1 / 16 * em + 3 / 8 * em**2)
            + Hoverh**3 * (
                20127 / 179200
                - 409 / 6400 * em
                + 7 / 64 * em**2
                + 1 / 16 * em**3
            )
            + Hoverh**4 * (
                -1575087 / 28672000
                + 1086367 / 1792000 * em
                - 2679 / 25600 * em**2
                + 13 / 128 * em**3
                + 3 / 128 * em**4
            )
        )
    ) * h
    return h / L, hvals / L


def callout(ax, text, xy, xytext):
    ax.annotate(
        text,
        xy=xy,
        xytext=xytext,
        textcoords="data",
        ha="center",
        va="center",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="0.35", alpha=0.92),
        arrowprops=dict(arrowstyle="->", color="0.25", lw=1.0, shrinkA=4, shrinkB=4),
        zorder=10,
    )


def build_figure():
    # -------------------------------------------------------
    # Same computational setup as the user's MATLAB workflow
    # -------------------------------------------------------
    kh = (10.0 ** np.arange(-3, 0.8001, 0.02))[:, None]
    kH = (10.0 ** np.arange(-5, -0.23 + 1e-12, 0.02)) * 2.0
    HOverL, hOverL = np.meshgrid(kH / (2 * np.pi), kh[:, 0] / (2 * np.pi))

    ax_range = [0.01, 1.0, 1e-5, 0.3]
    UR_lim = 26.0

    sigma = np.tanh(kh)
    B31 = (3 + 8 * sigma**2 - 9 * sigma**4) / (16 * sigma**4)
    B33 = (27 - 9 * sigma**2 + 9 * sigma**4 - 3 * sigma**6) / (64 * sigma**6)

    alpha1 = np.cosh(2 * kh)
    B51 = (
        121 * alpha1**5 + 263 * alpha1**4 + 376 * alpha1**3
        - 1999 * alpha1**2 + 2509 * alpha1 - 1108
    ) / (192 * (alpha1 - 1) ** 5)

    B53 = (
        (57 * alpha1**7 + 204 * alpha1**6 - 53 * alpha1**5
         - 782 * alpha1**4 - 741 * alpha1**3 - 52 * alpha1**2
         + 371 * alpha1 + 186) * 9
    ) / (((3 * alpha1 + 2) * (alpha1 - 1) ** 6) * 128)

    B55 = (
        (300 * alpha1**8 + 1579 * alpha1**7 + 3176 * alpha1**6
         + 2949 * alpha1**5 + 1188 * alpha1**4 + 675 * alpha1**3
         + 1326 * alpha1**2 + 827 * alpha1 + 130) * 5
    ) / (((alpha1 - 1) ** 6) * (12 * alpha1**2 + 11 * alpha1 + 2) * 384)

    # Stokes-region solve
    stokes_mask = HOverL <= (UR_lim * hOverL ** 3)
    C = np.where(stokes_mask, np.pi * HOverL, np.nan)
    A = B31 + B33
    B = B51 + B53 + B55

    ka = np.where(np.isnan(C), np.nan, C.copy())
    for _ in range(30):
        f = ka + A * ka**3 + B * ka**5 - C
        fp = 1 + 3 * A * ka**2 + 5 * B * ka**4
        step = f / fp
        ka_new = ka - step
        ka_new = np.where(ka_new > 0, ka_new, ka / 2)
        ka = np.where(np.isnan(ka), np.nan, ka_new)

    eta1 = np.ones_like(ka)
    eta2 = ka / 4 * (3 - sigma**2) / sigma**3
    eta3 = ka**2 * (B31 + B33)
    sigma1 = 24 * (3 * np.cosh(2 * kh) + 2) * (np.cosh(2 * kh) - 1) ** 4 * np.sinh(2 * kh)
    eta4 = ka**3 / sigma1 * (
        (60 * alpha1**6 + 232 * alpha1**5 - 118 * alpha1**4 - 989 * alpha1**3
         - 607 * alpha1**2 + 352 * alpha1 + 260)
        + (24 * alpha1**6 + 116 * alpha1**5 + 214 * alpha1**4 + 188 * alpha1**3
           + 133 * alpha1**2 + 101 * alpha1 + 34)
    )
    eta5 = ka**4 * (B51 + B53 + B55)

    eta_ratio2 = eta2 / eta1
    eta_ratio3 = eta3 / (eta1 + eta2)
    eta_ratio4 = eta4 / (eta1 + eta2 + eta3)
    eta_ratio5 = eta5 / (eta1 + eta2 + eta3 + eta4)

    # Order-boundary contours for legend/overlay only
    A2_paths = contour_level_segments(hOverL, np.where(stokes_mask, HOverL, np.nan), eta_ratio2, 0.01)
    A3_paths = contour_level_segments(hOverL, np.where(stokes_mask, HOverL, np.nan), eta_ratio3, 0.01)
    A4_paths = contour_level_segments(hOverL, np.where(stokes_mask, HOverL, np.nan), eta_ratio4, 0.01)
    A5_paths = contour_level_segments(hOverL, np.where(stokes_mask, HOverL, np.nan), eta_ratio5, 0.01)

    # Breaking line
    lambda_over_h = 2 * np.pi / kh[:, 0]
    HoverL_break = (
        kh[:, 0]
        * (0.141063 * lambda_over_h
           + 0.0095721 * lambda_over_h**2
           + 0.0077829 * lambda_over_h**3)
        / (1
           + 0.0788340 * lambda_over_h
           + 0.0317567 * lambda_over_h**2
           + 0.0093407 * lambda_over_h**3)
        / (2 * np.pi)
    )
    y_break_2d = HoverL_break[:, None]

    # Cnoidal family boundaries
    x1d = hOverL[:, 0]
    x_m05_raw, y_m05_raw = fenton_cnoidal_curve(0.5)
    x_m96_raw, y_m96_raw = fenton_cnoidal_curve(0.96)
    x_msol_raw, y_msol_raw = fenton_cnoidal_curve(1 - 4e-8)

    y_m05 = interp_log_curve(x_m05_raw, y_m05_raw, x1d)[:, None]
    y_m96 = interp_log_curve(x_m96_raw, y_m96_raw, x1d)[:, None]
    y_msol = interp_log_curve(x_msol_raw, y_msol_raw, x1d)[:, None]

    # Optional numerical/stream-function limit line from the MATLAB script
    y_num = (0.4 * x1d)
    y_num[(x1d > 0.20004) | (x1d < 0.024)] = np.nan

    # -------------------------------------------------------
    # Disjoint classification: one category per grid cell
    # -------------------------------------------------------
    # Codes:
    # 0 beyond breaking / invalid
    # 1 linear
    # 2 stokes ii
    # 3 stokes iii
    # 4 stokes iv
    # 5 stokes v
    # 6 cnoidal
    # 7 cnoidal iii
    # 8 cnoidal v
    # 9 solitary-wave side
    # 10 stream-function / numerical fallback
    category = np.full(HOverL.shape, 10, dtype=int)

    # Above breaking is invalid
    category[HOverL > y_break_2d] = 0

    # Stokes side (UR <= 26) below breaking
    stokes_valid = stokes_mask & (HOverL <= y_break_2d)
    category[stokes_valid] = 5  # default deepest Stokes order below breaking
    category[stokes_valid & (eta_ratio5 < 0.01)] = 4
    category[stokes_valid & (eta_ratio4 < 0.01)] = 3
    category[stokes_valid & (eta_ratio3 < 0.01)] = 2
    category[stokes_valid & (eta_ratio2 < 0.01)] = 1

    # Shallow-water / non-Stokes side
    nonstokes = (~stokes_mask) & (HOverL <= y_break_2d)

    # Start with stream-function/numerical fallback for all non-Stokes points,
    # then overwrite where cnoidal/solitary boundaries classify them.
    category[nonstokes] = 10
    category[nonstokes & (HOverL < y_m05)] = 6
    category[nonstokes & (HOverL >= y_m05) & (HOverL < y_m96)] = 7
    category[nonstokes & (HOverL >= y_m96) & (HOverL < y_msol)] = 8
    category[nonstokes & (HOverL >= y_msol)] = 9

    # -------------------------------------------------------
    # Plot
    # -------------------------------------------------------
    fig, ax = plt.subplots(figsize=(12, 8.5))

    cmap_base = plt.get_cmap("tab20")
    colors = [cmap_base(i) for i in [14, 1, 3, 5, 7, 9, 11, 13, 15, 17, 19]]
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(np.arange(-0.5, len(AREA_LABELS) + 0.5, 1), cmap.N)

    pm = ax.pcolormesh(hOverL, HOverL, category, cmap=cmap, norm=norm, shading="nearest", alpha=0.30)

    # Boundary overlays
    for paths in [A2_paths, A3_paths, A4_paths, A5_paths]:
        p = longest_path(paths)
        if p is not None:
            ax.plot(p[:, 0], p[:, 1], "-", linewidth=1.8, color="k")

    line_break, = ax.plot(x1d, HoverL_break, "-.", linewidth=2.2, label="Breaking limit")
    line_ur26, = ax.plot(x1d, 26 * x1d**3, "--", linewidth=2.0, label=r"$U_r = 26$")
    line_m05, = ax.plot(x_m05_raw, y_m05_raw, ":", linewidth=2.2, label=r"$m = 0.5$")
    line_m96, = ax.plot(x_m96_raw, y_m96_raw, ":", linewidth=2.2, label=r"$m = 0.96$")
    line_msol, = ax.plot(x_msol_raw, y_msol_raw, "-", linewidth=2.2, label=r"$m = 1 - 4\times 10^{-8}$")
    line_num, = ax.plot(x1d, y_num, "-.", linewidth=2.0, label="Stream-function / numerical limit")
    line_h005 = ax.axvline(0.05, linestyle="--", linewidth=1.2, color="0.4")
    line_h05 = ax.axvline(0.5, linestyle="--", linewidth=1.2, color="0.4")

    ax.text(1.25e-2, 6.4e-5, r"$U_r = 26$", rotation=45)
    ax.text(2.1e-2, 1.8e-4, r"$m = 0.5$", rotation=45)
    ax.text(1.8e-2, 1.2e-2, r"$m = 1-4\times 10^{-8}$", rotation=56)

    # Callouts for disjoint areas
    callout(ax, "Linear", xy=(0.22, 2.1e-3), xytext=(0.34, 1.0e-3))
    callout(ax, "Stokes II", xy=(0.20, 1.35e-2), xytext=(0.33, 8.0e-3))
    callout(ax, "Stokes III", xy=(0.24, 3.8e-2), xytext=(0.39, 2.8e-2))
    callout(ax, "Stokes IV", xy=(0.30, 6.6e-2), xytext=(0.48, 5.6e-2))
    callout(ax, "Stokes V", xy=(0.14, 7.5e-2), xytext=(0.12, 1.20e-1))

    callout(ax, "Cnoidal", xy=(0.040, 3.5e-4), xytext=(0.026, 7.0e-4))
    callout(ax, "Cnoidal III", xy=(0.060, 4.0e-3), xytext=(0.026, 3.0e-3))
    callout(ax, "Cnoidal V", xy=(0.095, 2.4e-2), xytext=(0.030, 2.1e-2))
    callout(ax, "Solitary-wave\nside", xy=(0.017, 8.0e-3), xytext=(0.013, 1.9e-2))
    callout(ax, "Stream-function /\nnumerical", xy=(0.11, 3.8e-2), xytext=(0.23, 1.05e-1))
    callout(ax, "Beyond breaking /\ninvalid", xy=(0.11, 1.7e-1), xytext=(0.24, 2.1e-1))

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(ax_range[0], ax_range[1])
    ax.set_ylim(ax_range[2], ax_range[3])
    ax.set_xlabel(r"$h/L$")
    ax.set_ylabel("H/L")
    ax.set_title("Wave-theory applicability chart with disjoint regions")
    ax.grid(True, which="major", alpha=0.25)

    area_handles = [
        Patch(facecolor=colors[i], edgecolor="none", alpha=0.30, label=AREA_LABELS[i])
        for i in range(len(AREA_LABELS))
    ]
    leg_area = ax.legend(
        handles=area_handles,
        title="Applicability areas",
        loc="lower right",
        bbox_to_anchor=(0.995, 0.01),
        borderaxespad=0.4,
        fontsize=9,
        title_fontsize=10,
        framealpha=0.95,
    )
    ax.add_artist(leg_area)

    line_legend_handles = [
        Line2D([0], [0], color="k", linewidth=1.8, label=r"A2/A1 = 1%"),
        Line2D([0], [0], color="k", linewidth=1.8, label=r"A3/(A1+A2) = 1%"),
        Line2D([0], [0], color="k", linewidth=1.8, label=r"A4/(A1+A2+A3) = 1%"),
        Line2D([0], [0], color="k", linewidth=1.8, label=r"A5/(A1+A2+A3+A4) = 1%"),
        Line2D([0], [0], color=line_break.get_color(), linestyle="-.", linewidth=2.2, label="Breaking limit"),
        Line2D([0], [0], color=line_ur26.get_color(), linestyle="--", linewidth=2.0, label=r"$U_r = 26$"),
        Line2D([0], [0], color=line_m05.get_color(), linestyle=":", linewidth=2.2, label=r"$m = 0.5$"),
        Line2D([0], [0], color=line_m96.get_color(), linestyle=":", linewidth=2.2, label=r"$m = 0.96$"),
        Line2D([0], [0], color=line_msol.get_color(), linestyle="-", linewidth=2.2, label=r"$m = 1 - 4\times10^{-8}$"),
        Line2D([0], [0], color=line_num.get_color(), linestyle="-.", linewidth=2.0, label="Stream-function / numerical limit"),
        Line2D([0], [0], color="0.4", linestyle="--", linewidth=1.2, label=r"$h/L = 0.05,\ 0.5$"),
    ]
    leg_lines = ax.legend(
        handles=line_legend_handles,
        title="Contours and boundaries",
        loc="lower right",
        bbox_to_anchor=(0.70, 0.01),
        borderaxespad=0.4,
        fontsize=9,
        title_fontsize=10,
        framealpha=0.95,
    )
    ax.add_artist(leg_lines)

    fig.tight_layout()
    return fig


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate a cleaned wave-theory applicability chart with disjoint regions."
    )
    parser.add_argument(
        "--formats",
        nargs="+",
        default=["png"],
        choices=["png", "svg", "pgf"],
        help="Output formats to save. Example: --formats png svg pgf",
    )
    parser.add_argument(
        "--output-dir",
        default=".",
        help="Directory where output files will be written.",
    )
    parser.add_argument(
        "--basename",
        default="wave_theory_applicability_clean",
        help="Base filename without extension.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if "pgf" in args.formats:
        configure_pgf()

    fig = build_figure()

    saved = []
    for fmt in args.formats:
        out = output_dir / f"{args.basename}.{fmt}"
        fig.savefig(out, bbox_inches="tight", dpi=220 if fmt == "png" else None)
        saved.append(out)

    plt.close(fig)

    print("Saved files:")
    for path in saved:
        print(path)


if __name__ == "__main__":
    main()
