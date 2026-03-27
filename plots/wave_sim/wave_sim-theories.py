#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import special
from matplotlib.patches import Patch
from matplotlib.lines import Line2D


def configure_pgf():
    """
    Enable PGF/XeLaTeX export styling.
    This is only activated if --formats includes pgf.
    """
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


def callout(ax, text, xy, xytext):
    ax.annotate(
        text,
        xy=xy,
        xytext=xytext,
        textcoords='data',
        ha='center',
        va='center',
        fontsize=9,
        bbox=dict(boxstyle='round,pad=0.22', fc='white', ec='0.35', alpha=0.92),
        arrowprops=dict(arrowstyle='->', color='0.25', lw=1.0, shrinkA=4, shrinkB=4),
        zorder=8
    )


def build_figure():
    # -------------------------------------------------------
    # Core calculations translated from the MATLAB workflow
    # -------------------------------------------------------
    kh = (10.0 ** np.arange(-3, 0.8001, 0.02))[:, None]
    kH = (10.0 ** np.arange(-5, -0.23 + 1e-12, 0.02)) * 2.0
    HOverL, hOverL = np.meshgrid(kH / (2 * np.pi), kh[:, 0] / (2 * np.pi))

    ax_range = [0.01, 1.0, 1e-5, 0.3]
    UR_lim = 26.0

    HOverL_masked = HOverL.copy()
    HOverL_masked[HOverL > (UR_lim * hOverL ** 3)] = np.nan

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

    A = B31 + B33
    B = B51 + B53 + B55
    C = np.pi * HOverL_masked

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

    A2_paths = contour_level_segments(hOverL, HOverL_masked, eta_ratio2, 0.01)
    A3_paths = contour_level_segments(hOverL, HOverL_masked, eta_ratio3, 0.01)
    A4_paths = contour_level_segments(hOverL, HOverL_masked, eta_ratio4, 0.01)
    A5_paths = contour_level_segments(hOverL, HOverL_masked, eta_ratio5, 0.01)

    p2 = longest_path(A2_paths)
    p3 = longest_path(A3_paths)
    p4 = longest_path(A4_paths)
    p5 = longest_path(A5_paths)

    x = np.logspace(np.log10(ax_range[0]), np.log10(ax_range[1]), 900)
    y_min = np.full_like(x, ax_range[2], dtype=float)

    y_a2 = interp_log_curve(p2[:, 0], p2[:, 1], x) if p2 is not None else np.full_like(x, np.nan)
    y_a3 = interp_log_curve(p3[:, 0], p3[:, 1], x) if p3 is not None else np.full_like(x, np.nan)
    y_a4 = interp_log_curve(p4[:, 0], p4[:, 1], x) if p4 is not None else np.full_like(x, np.nan)
    y_a5 = interp_log_curve(p5[:, 0], p5[:, 1], x) if p5 is not None else np.full_like(x, np.nan)

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
    y_break = interp_log_curve(hOverL[:, 0], HoverL_break, x)

    y_ur1 = x**3
    y_ur10 = 10 * x**3
    y_ur26 = 26 * x**3
    y_ur1 = np.where(y_ur1 <= y_break, y_ur1, np.nan)
    y_ur10 = np.where(y_ur10 <= y_break, y_ur10, np.nan)
    y_ur26 = np.where(y_ur26 <= y_break, y_ur26, np.nan)

    h = 10.0
    H = np.arange(0.01, 7.8 + 1e-12, 0.001)

    def fenton_cnoidal_curve(m):
        Km, Em = special.ellipk(m), special.ellipe(m)
        em = Em / Km
        Hoverh = H / h
        L = (
            4 * Km / np.sqrt(3 * Hoverh)
            * (1
               + Hoverh * (5/8 - 3/2 * em)
               + Hoverh**2 * (-21/128 + 1/16 * em + 3/8 * em**2)
               + Hoverh**3 * (20127/179200 - 409/6400 * em + 7/64 * em**2 + 1/16 * em**3)
               + Hoverh**4 * (-1575087/28672000 + 1086367/1792000 * em
                              - 2679/25600 * em**2 + 13/128 * em**3 + 3/128 * em**4))
        ) * h
        return h / L, H / L

    x_m96_raw, y_m96_raw = fenton_cnoidal_curve(0.96)
    x_msol_raw, y_msol_raw = fenton_cnoidal_curve(1 - 4e-8)

    y_m96 = interp_log_curve(x_m96_raw, y_m96_raw, x)
    y_msol = interp_log_curve(x_msol_raw, y_msol_raw, x)

    y_num = 0.4 * x
    y_num[(x > 0.20004) | (x < 0.024)] = np.nan

    fig, ax = plt.subplots(figsize=(12, 8.5))

    fills = []
    fills.append(ax.fill_between(x, y_min, y_a2, where=np.isfinite(y_a2), alpha=0.22))
    fills.append(ax.fill_between(x, y_a2, y_a3, where=np.isfinite(y_a2) & np.isfinite(y_a3), alpha=0.22))
    fills.append(ax.fill_between(x, y_a3, y_a4, where=np.isfinite(y_a3) & np.isfinite(y_a4), alpha=0.22))
    fills.append(ax.fill_between(x, y_a4, y_a5, where=np.isfinite(y_a4) & np.isfinite(y_a5), alpha=0.22))
    fills.append(ax.fill_between(x, y_a5, y_break, where=np.isfinite(y_a5) & np.isfinite(y_break), alpha=0.22))

    fills.append(ax.fill_between(x, y_ur1, y_ur10, where=np.isfinite(y_ur1) & np.isfinite(y_ur10), alpha=0.22))
    fills.append(ax.fill_between(x, y_ur10, y_m96, where=np.isfinite(y_ur10) & np.isfinite(y_m96), alpha=0.22))
    fills.append(ax.fill_between(x, y_m96, y_msol, where=np.isfinite(y_m96) & np.isfinite(y_msol), alpha=0.22))
    fills.append(ax.fill_between(x, y_msol, y_break, where=np.isfinite(y_msol) & np.isfinite(y_break), alpha=0.22))
    fills.append(ax.fill_between(x, y_num, y_break, where=np.isfinite(y_num) & np.isfinite(y_break), alpha=0.18))

    if p2 is not None:
        ax.plot(p2[:, 0], p2[:, 1], '-', linewidth=1.8, color='k')
    if p3 is not None:
        ax.plot(p3[:, 0], p3[:, 1], '-', linewidth=1.8, color='k')
    if p4 is not None:
        ax.plot(p4[:, 0], p4[:, 1], '-', linewidth=1.8, color='k')
    if p5 is not None:
        ax.plot(p5[:, 0], p5[:, 1], '-', linewidth=1.8, color='k')

    ax.plot(hOverL[:, 0], HoverL_break, '-.', linewidth=2.2)
    ax.plot(x, y_ur26, '--', linewidth=2.0)
    ax.plot(x, y_ur10, '--', linewidth=2.0)
    ax.plot(x, y_ur1, '--', linewidth=2.0)
    ax.plot(x_m96_raw, y_m96_raw, ':', linewidth=2.2)
    ax.plot(x_msol_raw, y_msol_raw, '-', linewidth=2.2)
    ax.plot(x, y_num, '-.', linewidth=2.0)
    ax.axvline(0.05, linestyle='--', linewidth=1.2, color='0.4')
    ax.axvline(0.5, linestyle='--', linewidth=1.2, color='0.4')

    ax.text(1.22e-2, 6.3e-5, r'$U_r = 26$', rotation=45)
    ax.text(2.05e-2, 6.0e-5, r'$U_r = 10$', rotation=45)
    ax.text(2.70e-2, 1.5e-5, r'$U_r = 1$', rotation=45)

    callout(ax, 'Linear', xy=(0.22, 2.1e-3), xytext=(0.34, 1.0e-3))
    callout(ax, 'Stokes II', xy=(0.20, 1.35e-2), xytext=(0.33, 8.0e-3))
    callout(ax, 'Stokes III', xy=(0.24, 3.8e-2), xytext=(0.39, 2.8e-2))
    callout(ax, 'Stokes IV', xy=(0.30, 6.6e-2), xytext=(0.48, 5.6e-2))
    callout(ax, 'Stokes V', xy=(0.14, 7.5e-2), xytext=(0.12, 1.20e-1))

    callout(ax, 'Cnoidal\n(low-order)', xy=(0.030, 8.5e-5), xytext=(0.015, 1.7e-4))
    callout(ax, 'Cnoidal III', xy=(0.060, 3.0e-3), xytext=(0.026, 2.8e-3))
    callout(ax, 'Cnoidal V', xy=(0.090, 3.1e-2), xytext=(0.030, 2.0e-2))
    callout(ax, 'Solitary-wave\nside', xy=(0.017, 8.0e-3), xytext=(0.013, 1.9e-2))
    callout(ax, 'Numerical /\nstream-function', xy=(0.12, 4.5e-2), xytext=(0.22, 1.05e-1))

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(ax_range[0], ax_range[1])
    ax.set_ylim(ax_range[2], ax_range[3])
    ax.set_xlabel(r'$h/L$')
    ax.set_ylabel('H/L')
    ax.set_title('Wave-theory applicability chart with colored regions and callouts')
    ax.grid(True, which='major', alpha=0.25)

    theory_labels = [
        'Linear', 'Stokes II', 'Stokes III', 'Stokes IV', 'Stokes V',
        'Cnoidal (low-order)', 'Cnoidal III', 'Cnoidal V',
        'Solitary-wave side', 'Numerical / stream-function'
    ]
    theory_handles = [
        Patch(facecolor=fills[i].get_facecolor()[0], edgecolor='none',
              alpha=0.22 if i < 9 else 0.18, label=theory_labels[i])
        for i in range(len(theory_labels))
    ]
    leg1 = ax.legend(handles=theory_handles, title='Colored applicability areas',
                     loc='lower right', bbox_to_anchor=(0.995, 0.01),
                     borderaxespad=0.4, fontsize=9, title_fontsize=10, framealpha=0.95)
    ax.add_artist(leg1)

    line_legend_handles = [
        Line2D([0], [0], color='k', linewidth=1.8, label=r'A2/A1 = 1%'),
        Line2D([0], [0], color='k', linewidth=1.8, label=r'A3/(A1+A2) = 1%'),
        Line2D([0], [0], color='k', linewidth=1.8, label=r'A4/(A1+A2+A3) = 1%'),
        Line2D([0], [0], color='k', linewidth=1.8, label=r'A5/(A1+A2+A3+A4) = 1%'),
        Line2D([0], [0], color='C0', linestyle='-.', linewidth=2.2, label='Breaking limit'),
        Line2D([0], [0], color='C1', linestyle='--', linewidth=2.0, label=r'$U_r = 26$'),
        Line2D([0], [0], color='C2', linestyle='--', linewidth=2.0, label=r'$U_r = 10$'),
        Line2D([0], [0], color='C3', linestyle='--', linewidth=2.0, label=r'$U_r = 1$'),
        Line2D([0], [0], color='C4', linestyle=':', linewidth=2.2, label=r'$m = 0.96$'),
        Line2D([0], [0], color='C5', linestyle='-', linewidth=2.2, label=r'$m = 1 - 4\times10^{-8}$'),
        Line2D([0], [0], color='C6', linestyle='-.', linewidth=2.0, label='Numerical / stream-function limit'),
        Line2D([0], [0], color='0.4', linestyle='--', linewidth=1.2, label=r'$h/L = 0.05,\ 0.5$'),
    ]
    leg2 = ax.legend(handles=line_legend_handles, title='Contours and boundaries',
                     loc='lower right', bbox_to_anchor=(0.73, 0.01),
                     borderaxespad=0.4, fontsize=9, title_fontsize=10, framealpha=0.95)
    ax.add_artist(leg2)

    fig.tight_layout()
    return fig


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate the wave-theory applicability chart with colored regions and callouts."
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
        default="wave_theory_applicability_callouts",
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
