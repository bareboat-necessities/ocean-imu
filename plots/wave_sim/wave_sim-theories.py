import numpy as np
import matplotlib.pyplot as plt
from scipy import special


def reproduce_wave_theory_chart():
    # ----------------------------
    # Direct Python translation of the MATLAB workflow
    # ----------------------------
    kh = (10.0 ** np.arange(-3, 0.8001, 0.02))[:, None]   # column vector
    kH = (10.0 ** np.arange(-5, -0.23 + 1e-12, 0.02)) * 2.0
    HOverL, hOverL = np.meshgrid(kH / (2 * np.pi), kh[:, 0] / (2 * np.pi))

    ax_range = [0.01, 1.0, 1e-5, 0.3]
    UR_lim = 26.0  # Stokes wave applicable in the range with Ursell <= 26

    # Mask Stokes region beyond UR = 26
    HOverL_masked = HOverL.copy()
    HOverL_masked[HOverL > (UR_lim * hOverL ** 3)] = np.nan

    # ----------------------------
    # Coefficients
    # ----------------------------
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

    # ----------------------------
    # Solve for ka from:
    # 2*pi*H/L - ka*(2 + 2*(B31+B33)*ka^2 + 2*(B51+B53+B55)*ka^4) = 0
    # which simplifies to:
    # pi*H/L = ka + A*ka^3 + B*ka^5
    # ----------------------------
    A = B31 + B33
    B = B51 + B53 + B55
    C = np.pi * HOverL_masked

    ka = np.where(np.isnan(C), np.nan, C.copy())  # small-wave initial guess

    # Newton iteration
    for _ in range(30):
        f = ka + A * ka**3 + B * ka**5 - C
        fp = 1 + 3 * A * ka**2 + 5 * B * ka**4
        step = f / fp
        ka_new = ka - step
        ka_new = np.where(ka_new > 0, ka_new, ka / 2)
        ka = np.where(np.isnan(ka), np.nan, ka_new)

    # ----------------------------
    # Free-surface order contributions evaluated at crest (theta = 0)
    # ----------------------------
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

    # ----------------------------
    # Helper functions
    # ----------------------------
    def contour_level_segments(Z, level=0.01):
        fig, ax = plt.subplots()
        cs = ax.contour(hOverL, HOverL_masked, Z, levels=[level])
        segs = [np.asarray(seg) for seg in cs.allsegs[0] if len(seg) > 1]
        plt.close(fig)
        return segs

    def plot_contourf(Z, title):
        fig, ax = plt.subplots(figsize=(7, 5.5))
        cs = ax.contourf(hOverL, HOverL_masked, Z, levels=[0.01, 0.05, 0.1, 0.2])
        fig.colorbar(cs, ax=ax)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(ax_range[0], ax_range[1])
        ax.set_ylim(ax_range[2], ax_range[3])
        ax.set_xlabel("h/L")
        ax.set_ylabel("H/L")
        ax.set_title(title)
        fig.tight_layout()
        return fig, ax

    # ----------------------------
    # Individual contour figures
    # ----------------------------
    plot_contourf(eta_ratio2, "A2 / A1 at crest")
    plot_contourf(eta_ratio3, "A3 / (A1 + A2) at crest")
    plot_contourf(eta_ratio4, "A4 / (A1 + A2 + A3) at crest")
    plot_contourf(eta_ratio5, "A5 / (A1 + A2 + A3 + A4) at crest")

    # Extract 1% contour lines
    A2_paths = contour_level_segments(eta_ratio2, 0.01)
    A3_paths = contour_level_segments(eta_ratio3, 0.01)
    A4_paths = contour_level_segments(eta_ratio4, 0.01)
    A5_paths = contour_level_segments(eta_ratio5, 0.01)

    # ----------------------------
    # Combined plot
    # ----------------------------
    fig, ax = plt.subplots(figsize=(8, 6.2))

    for seg in A2_paths:
        ax.plot(seg[:, 0], seg[:, 1], 'k')
    for seg in A3_paths:
        ax.plot(seg[:, 0], seg[:, 1], 'k')
    for seg in A4_paths:
        ax.plot(seg[:, 0], seg[:, 1], 'k')
    for seg in A5_paths:
        ax.plot(seg[:, 0], seg[:, 1], 'k')

    # Breaking lines by Fenton
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
    ax.plot(hOverL[:, 0], HoverL_break, 'b-.', linewidth=2)

    # Ursell number lines
    Ur10 = 10 * hOverL[:, 0] ** 3
    Ur10[Ur10 > HoverL_break] = np.nan
    ax.plot(hOverL[:, 0], Ur10, 'b', linewidth=2)

    Ur1 = 1 * hOverL[:, 0] ** 3
    Ur1[Ur1 > HoverL_break] = np.nan
    ax.plot(hOverL[:, 0], Ur1, 'b--', linewidth=2)

    Ur26 = 26 * hOverL[:, 0] ** 3
    Ur26[Ur26 > HoverL_break] = np.nan
    ax.plot(hOverL[:, 0], Ur26, 'b--', linewidth=2)

    # ----------------------------
    # Plot m curves based on Fenton's cnoidal wave solutions
    # ----------------------------
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

    x96, y96 = fenton_cnoidal_curve(0.96)
    ax.plot(x96, y96, 'r:', linewidth=2)

    xsol, ysol = fenton_cnoidal_curve(1 - 4e-8)
    ax.plot(xsol, ysol, 'r', linewidth=2)

    # Numerical solutions region as analytical solution becomes problematic
    sol_lim = hOverL[:, 0] * 0.4
    sol_lim[hOverL[:, 0] > 0.20004] = np.nan
    sol_lim[hOverL[:, 0] < 0.024] = np.nan
    ax.plot(hOverL[:, 0], sol_lim, 'k-.', linewidth=2)

    # Vertical reference lines
    ax.axvline(0.5, color='k', linestyle='--', linewidth=1)
    ax.axvline(0.05, color='k', linestyle='--', linewidth=1)

    # Annotation
    xt = [1.2e-2, 2.1e-2, 2.7e-2]
    yt = [6e-5, 6e-5, 1.4e-5]
    labels = [r'$U_r = 26$', r'$U_r = 10,\ m \approx 0.5$', r'$U_r = 1$']
    for x, y, txt in zip(xt, yt, labels):
        ax.text(x, y, txt, rotation=45, color='b')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(ax_range[0], ax_range[1])
    ax.set_ylim(ax_range[2], ax_range[3])
    ax.set_xlabel(r'$h/L$')
    ax.set_ylabel('H/L')
    ax.set_title('Combined applicability chart reproduced from MATLAB workflow')
    ax.tick_params(length=8)
    ax.set_axisbelow(False)

    plt.show()


if __name__ == "__main__":
    reproduce_wave_theory_chart()
  
