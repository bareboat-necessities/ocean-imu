#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy import special, ndimage


AREA_LABELS = [
    "Beyond breaking / invalid",
    "Linear (Airy)",
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

AREA_COLORS = [
    "#d0d0d0",
    "#cfe1f2",
    "#efdcc4",
    "#d8ead2",
    "#f4cccc",
    "#d9d2e9",
    "#fff2cc",
    "#ead1dc",
    "#d9d9d9",
    "#d9ead3",
    "#cfe2f3",
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


def interior_point_visible(xgrid2d, ygrid2d, category2d, code, xlim, ylim):
    vis = (xgrid2d >= xlim[0]) & (xgrid2d <= xlim[1]) & (ygrid2d >= ylim[0]) & (ygrid2d <= ylim[1])
    mask = (category2d == code) & vis
    if not np.any(mask):
        return None
    dist = ndimage.distance_transform_edt(mask)
    idx = np.unravel_index(np.argmax(dist), dist.shape)
    return float(xgrid2d[idx]), float(ygrid2d[idx])


def fenton_cnoidal_curve(m, h=10.0, hvals=None):
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


def interp_curve_to_x(x1d, xraw, yraw):
    mask = np.isfinite(xraw) & np.isfinite(yraw) & (xraw > 0) & (yraw > 0)
    x = np.asarray(xraw)[mask]
    y = np.asarray(yraw)[mask]
    if len(x) < 2:
        return np.full_like(x1d, np.nan, dtype=float)

    order = np.argsort(x)
    x = x[order]
    y = y[order]

    out = np.full_like(x1d, np.nan, dtype=float)
    valid = (x1d >= x.min()) & (x1d <= x.max())
    out[valid] = 10 ** np.interp(np.log10(x1d[valid]), np.log10(x), np.log10(y))
    return out


def first_visible_point(x, y, frac=0.3):
    mask = np.isfinite(x) & np.isfinite(y)
    if not np.any(mask):
        return None
    xv = np.asarray(x)[mask]
    yv = np.asarray(y)[mask]
    idx = max(0, min(len(xv) - 1, int(frac * (len(xv) - 1))))
    return float(xv[idx]), float(yv[idx])


def callout(ax, text, xy, xytext):
    if xy is None:
        return
    ax.annotate(
        text,
        xy=xy,
        xytext=xytext,
        textcoords="data",
        ha="center",
        va="center",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="0.35", alpha=0.96),
        arrowprops=dict(arrowstyle="->", color="0.25", lw=1.0, shrinkA=4, shrinkB=4),
        zorder=12,
    )


def endpoint_label(ax, x, y, text, dx, dy):
    ax.plot([x], [y], marker='o', markersize=3.0, color='0.2', zorder=8)
    ax.annotate(
        text, xy=(x, y), xytext=(dx, dy), textcoords='offset points',
        fontsize=8, ha='left', va='center',
        bbox=dict(boxstyle='round,pad=0.16', fc='white', ec='0.45', alpha=0.92),
        arrowprops=dict(arrowstyle='->', color='0.3', lw=0.85, shrinkA=2, shrinkB=3),
        zorder=12
    )


def build_figure():
    kh = (10.0 ** np.arange(-3, 0.8001, 0.02))[:, None]
    kH = (10.0 ** np.arange(-5, -0.23 + 1e-12, 0.02)) * 2.0
    HOverL, hOverL = np.meshgrid(kH / (2 * np.pi), kh[:, 0] / (2 * np.pi))
    x1d = hOverL[:, 0]

    xlim = (0.01, 1.0)
    ylim = (1e-5, 0.3)
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

    lambda_over_h = 2 * np.pi / kh[:, 0]
    y_break = (
        kh[:, 0]
        * (
            0.141063 * lambda_over_h
            + 0.0095721 * lambda_over_h**2
            + 0.0077829 * lambda_over_h**3
        )
        / (
            1
            + 0.0788340 * lambda_over_h
            + 0.0317567 * lambda_over_h**2
            + 0.0093407 * lambda_over_h**3
        )
        / (2 * np.pi)
    )
    y_break_2d = y_break[:, None]

    x_m05_raw, y_m05_raw = fenton_cnoidal_curve(0.5)
    x_m96_raw, y_m96_raw = fenton_cnoidal_curve(0.96)
    x_msol_raw, y_msol_raw = fenton_cnoidal_curve(1 - 4e-8)

    y_m05 = interp_curve_to_x(x1d, x_m05_raw, y_m05_raw)[:, None]
    y_m96 = interp_curve_to_x(x1d, x_m96_raw, y_m96_raw)[:, None]
    y_msol = interp_curve_to_x(x1d, x_msol_raw, y_msol_raw)[:, None]

    category = np.full(HOverL.shape, 10, dtype=int)
    category[HOverL > y_break_2d] = 0

    stokes_valid = stokes_mask & (HOverL <= y_break_2d)
    category[stokes_valid] = 5
    category[stokes_valid & (eta_ratio5 < 0.01)] = 4
    category[stokes_valid & (eta_ratio4 < 0.01)] = 3
    category[stokes_valid & (eta_ratio3 < 0.01)] = 2
    category[stokes_valid & (eta_ratio2 < 0.01)] = 1

    nonstokes = (~stokes_mask) & (HOverL <= y_break_2d)
    category[nonstokes & np.isfinite(y_m05) & (HOverL < y_m05)] = 6
    category[nonstokes & np.isfinite(y_m05) & np.isfinite(y_m96) & (HOverL >= y_m05) & (HOverL < y_m96)] = 7
    category[nonstokes & np.isfinite(y_m96) & np.isfinite(y_msol) & (HOverL >= y_m96) & (HOverL < y_msol)] = 8
    category[nonstokes & np.isfinite(y_msol) & (HOverL >= y_msol) & (HOverL <= y_break_2d)] = 9

    A2_paths = contour_level_segments(hOverL, np.where(stokes_mask, HOverL, np.nan), eta_ratio2, 0.01)
    A3_paths = contour_level_segments(hOverL, np.where(stokes_mask, HOverL, np.nan), eta_ratio3, 0.01)
    A4_paths = contour_level_segments(hOverL, np.where(stokes_mask, HOverL, np.nan), eta_ratio4, 0.01)
    A5_paths = contour_level_segments(hOverL, np.where(stokes_mask, HOverL, np.nan), eta_ratio5, 0.01)

    fig, ax = plt.subplots(figsize=(12, 8.3))
    cmap = ListedColormap(AREA_COLORS)
    norm = BoundaryNorm(np.arange(-0.5, len(AREA_LABELS) + 0.5, 1), cmap.N)

    ax.pcolormesh(hOverL, HOverL, category, cmap=cmap, norm=norm, shading="nearest", zorder=1)
    ax.contour(hOverL, HOverL, category, levels=np.arange(0.5, len(AREA_LABELS) - 0.5, 1.0),
               colors="white", linewidths=0.7, zorder=2)

    for paths in [A2_paths, A3_paths, A4_paths, A5_paths]:
        p = longest_path(paths)
        if p is not None:
            ax.plot(p[:, 0], p[:, 1], "-", linewidth=1.8, color="k", zorder=4)

    line_break, = ax.plot(x1d, y_break, "-.", linewidth=2.2, zorder=5, label="Breaking waves line")
    line_ur26, = ax.plot(x1d, 26 * x1d**3, "--", linewidth=2.0, zorder=5, label=r"$U_r = 26$")

    ur26_raw_m05 = 26 * x_m05_raw**3
    ur26_raw_m96 = 26 * x_m96_raw**3
    ur26_raw_msol = 26 * x_msol_raw**3
    y_break_on_m05 = np.interp(x_m05_raw, x1d, y_break, left=np.nan, right=np.nan)
    y_break_on_m96 = np.interp(x_m96_raw, x1d, y_break, left=np.nan, right=np.nan)
    y_break_on_msol = np.interp(x_msol_raw, x1d, y_break, left=np.nan, right=np.nan)

    m05_mask = np.isfinite(x_m05_raw) & np.isfinite(y_m05_raw) & (y_m05_raw >= ur26_raw_m05) & (y_m05_raw <= y_break_on_m05)
    m96_mask = np.isfinite(x_m96_raw) & np.isfinite(y_m96_raw) & (y_m96_raw >= ur26_raw_m96) & (y_m96_raw <= y_break_on_m96)
    msol_mask = np.isfinite(x_msol_raw) & np.isfinite(y_msol_raw) & (y_msol_raw >= ur26_raw_msol) & (y_msol_raw <= y_break_on_msol)

    x_m05_plot = np.where(m05_mask, x_m05_raw, np.nan)
    y_m05_plot = np.where(m05_mask, y_m05_raw, np.nan)
    x_m96_plot = np.where(m96_mask, x_m96_raw, np.nan)
    y_m96_plot = np.where(m96_mask, y_m96_raw, np.nan)
    x_msol_plot = np.where(msol_mask, x_msol_raw, np.nan)
    y_msol_plot = np.where(msol_mask, y_msol_raw, np.nan)

    line_m05, = ax.plot(x_m05_plot, y_m05_plot, ":", linewidth=2.2, zorder=5, label=r"Cnoidal line $(m=0.5)$")
    line_m96, = ax.plot(x_m96_plot, y_m96_plot, ":", linewidth=2.2, zorder=5, label=r"$m = 0.96$")
    line_msol, = ax.plot(x_msol_plot, y_msol_plot, "-", linewidth=2.2, zorder=5, label=r"$m = 1 - 4\times10^{-8}$")

    # depth regime labels only
    y_top = 0.275
    ax.text(np.sqrt(0.01 * 0.05), y_top, "Shallow", fontsize=10, ha="center",
            bbox=dict(boxstyle='round,pad=0.16', fc='white', ec='0.5', alpha=0.92))
    ax.text(np.sqrt(0.05 * 0.5), y_top, "Intermediate", fontsize=10, ha="center",
            bbox=dict(boxstyle='round,pad=0.16', fc='white', ec='0.5', alpha=0.92))
    ax.text(np.sqrt(0.5 * 1.0), y_top, "Deep", fontsize=10, ha="center",
            bbox=dict(boxstyle='round,pad=0.16', fc='white', ec='0.5', alpha=0.92))
    for xpos, txt in [(0.05, r"$h/L=0.05$"), (0.5, r"$h/L=0.5$")]:
        ax.plot([xpos, xpos], [0.235, 0.262], color="0.45", lw=1.0, zorder=6)
        ax.text(xpos, 0.228, txt, ha="center", va="top", fontsize=8,
                bbox=dict(boxstyle='round,pad=0.12', fc='white', ec='0.5', alpha=0.9))

    # line labels
    ax.annotate("Breaking waves line",
                xy=(0.19, float(np.interp(0.19, x1d, y_break))),
                xytext=(0.34, 0.175),
                fontsize=9,
                bbox=dict(boxstyle='round,pad=0.16', fc='white', ec='0.4', alpha=0.94),
                arrowprops=dict(arrowstyle='->', color='0.3', lw=0.9),
                zorder=12)
    p_cnoidal_line = first_visible_point(x_m05_plot, y_m05_plot, frac=0.35)
    if p_cnoidal_line is not None:
        ax.annotate("Cnoidal line",
                    xy=p_cnoidal_line,
                    xytext=(0.022, 2.0e-4),
                    fontsize=9,
                    bbox=dict(boxstyle='round,pad=0.16', fc='white', ec='0.4', alpha=0.94),
                    arrowprops=dict(arrowstyle='->', color='0.3', lw=0.9),
                    zorder=12)

    valid96 = np.isfinite(x_m96_plot) & np.isfinite(y_m96_plot)
    validsol = np.isfinite(x_msol_plot) & np.isfinite(y_msol_plot)
    if np.any(valid96):
        endpoint_label(ax, x_m96_plot[valid96][-1], y_m96_plot[valid96][-1], "Cnoidal III ends", 12, 8)
    if np.any(validsol):
        endpoint_label(ax, x_msol_plot[validsol][-1], y_msol_plot[validsol][-1], "Cnoidal V ends", 14, -8)

    # callouts anchored to visible interiors with cleaner layout
    callout(ax, "Linear\n(Airy wave theory)", interior_point_visible(hOverL, HOverL, category, 1, xlim, ylim), (0.52, 1.3e-3))
    callout(ax, "Stokes II", interior_point_visible(hOverL, HOverL, category, 2, xlim, ylim), (0.38, 8.5e-3))
    callout(ax, "Stokes III", interior_point_visible(hOverL, HOverL, category, 3, xlim, ylim), (0.44, 2.9e-2))
    callout(ax, "Stokes IV", interior_point_visible(hOverL, HOverL, category, 4, xlim, ylim), (0.56, 5.8e-2))
    callout(ax, "Stokes V", interior_point_visible(hOverL, HOverL, category, 5, xlim, ylim), (0.11, 1.12e-1))
    # no separate "Cnoidal" area callout because that region is not visible in this plotting window
    callout(ax, "Cnoidal III", interior_point_visible(hOverL, HOverL, category, 7, xlim, ylim), (0.022, 2.7e-3))
    callout(ax, "Cnoidal V", interior_point_visible(hOverL, HOverL, category, 8, xlim, ylim), (0.030, 1.8e-2))
    callout(ax, "Solitary-wave\nside", interior_point_visible(hOverL, HOverL, category, 9, xlim, ylim), (0.014, 6.0e-2))
    callout(ax, "Stream-function /\nnumerical", interior_point_visible(hOverL, HOverL, category, 10, xlim, ylim), (0.13, 4.5e-2))
    callout(ax, "Beyond breaking /\ninvalid", interior_point_visible(hOverL, HOverL, category, 0, xlim, ylim), (0.08, 0.18))

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel(r"$h/L$")
    ax.set_ylabel("H/L")
    ax.set_title("Wave-theory applicability chart")
    ax.grid(True, which="major", alpha=0.25)

    # legends: only visible region classes in this plotting window
    vis = (hOverL >= xlim[0]) & (hOverL <= xlim[1]) & (HOverL >= ylim[0]) & (HOverL <= ylim[1])
    present_codes = [i for i in range(len(AREA_LABELS)) if np.any((category == i) & vis)]
    area_handles = [
        Patch(facecolor=AREA_COLORS[i], edgecolor="none", label=AREA_LABELS[i])
        for i in present_codes
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
        Line2D([0], [0], color=line_break.get_color(), linestyle="-.", linewidth=2.2, label="Breaking waves line"),
        Line2D([0], [0], color=line_ur26.get_color(), linestyle="--", linewidth=2.0, label=r"$U_r = 26$"),
        Line2D([0], [0], color=line_m05.get_color(), linestyle=":", linewidth=2.2, label=r"Cnoidal line $(m=0.5)$"),
        Line2D([0], [0], color=line_m96.get_color(), linestyle=":", linewidth=2.2, label=r"$m = 0.96$"),
        Line2D([0], [0], color=line_msol.get_color(), linestyle="-", linewidth=2.2, label=r"$m = 1 - 4\times10^{-8}$"),
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

    return fig


def parse_args():
    parser = argparse.ArgumentParser(description="Generate the tidied wave-theory applicability chart.")
    parser.add_argument("--formats", nargs="+", default=["png"], choices=["png", "svg", "pgf"])
    parser.add_argument("--output-dir", default=".")
    parser.add_argument("--basename", default="wave_theory_applicability_tidy")
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
