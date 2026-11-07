#!/usr/bin/env python3
"""Mathematical analysis of the shared-mode motion block's minimal realization.

Companion of docs/kalman-carma-shared-modes.md and
doc/kalman_carma_shared_modes/kalman_carma_shared_modes-math.tex.  Double
precision, independent of the C++ code (the C++ identities are tested in
tests/kalman_carma_shared_modes/carma_shared_modes_realization-test.cpp).

  1. observability of the complexified C41J mode per axis, realified to eight
     states: rank 6, the unobservable subspace lies in the imaginary half and
     equals {0} x span(Re w, Im w), w = [B(r), r B(r), r^2 B(r), 1], D(r) = 0
  2. its motion with kappa (principal angles between U(kappa) and U(1))
  3. the six-state cascade C41J * V, V = [N_r, N_i] / D co-inner: transfer
     functions against the complexified mode, lossless identities, kappa
     covariance
  4. minimal stochastic order of the three-axis law: Hankel singular values
     of the positive-real part, n = 4 rank(Sigma) + 2 rank(Im Sigma)
  5. conditioning (stationary correlation matrices) of the realizations

Writes reports/results/kalman_carma_shared_modes/analysis.json and
figures/analysis_*.{svg,pgf}.

Run: python3 tools/carma_shared_modes_analysis.py [--ext svg,pgf]
"""
from __future__ import annotations

import sys as _export_sys
from pathlib import Path as _ExportPath
_export_sys.path.insert(0, str(_ExportPath(__file__).resolve().parents[1] / "plots"))
from vector_export import save_figure


import argparse
import json
from pathlib import Path

import numpy as np
import scipy.linalg as sl

import importlib.util as _ilu  # noqa: E402
from pathlib import Path as _Path  # noqa: E402

# Uniform decimation of long lines in every saved figure, SVG and PGF alike (plots/figure_decimation.py).
_fd_spec = _ilu.spec_from_file_location("figure_decimation", _Path(__file__).resolve().parents[1] / "plots" / "figure_decimation.py")
_fd = _ilu.module_from_spec(_fd_spec)
_fd_spec.loader.exec_module(_fd)
_fd.install()

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "reports" / "results" / "kalman_carma_shared_modes"
FIG = OUT / "figures"
TWO_PI = 2.0 * np.pi
# CARMA41_JONSWAP_SHAPE (src/kalman_carma/CarmaVerticalModel.h).
SHAPE = dict(f1=0.306945, z1=0.104191, f2=0.815318, z2=0.478188, fz=0.0)


# ---------------------------------------------------------------------------
# models
# ---------------------------------------------------------------------------

def base_A(k):
    """C41J in [p, v, a, l] (carma_dynamics<4>)."""
    w1, z1 = TWO_PI * k * SHAPE["f1"], SHAPE["z1"]
    w2, z2 = TWO_PI * k * SHAPE["f2"], SHAPE["z2"]
    wz0 = TWO_PI * k * SHAPE["fz"]
    d = wz0 - 2 * z2 * w2
    A = np.zeros((4, 4))
    A[0, 1] = A[1, 2] = 1
    A[2] = [d * w1 * w1, 2 * d * z1 * w1 - w1 * w1, d - 2 * z1 * w1, -w2 * w2 - d * wz0]
    A[3] = [w1 * w1, 2 * z1 * w1, 1, -wz0]
    return A


def polys(k):
    """P+ = s^2 + b1 s + b0 (upper-half-plane poles of P), D = s^2 + d1 s + d0."""
    w1, z1 = TWO_PI * k * SHAPE["f1"], SHAPE["z1"]
    w2, z2 = TWO_PI * k * SHAPE["f2"], SHAPE["z2"]
    p1 = complex(-z1 * w1, w1 * np.sqrt(1 - z1 * z1))
    p2 = complex(-z2 * w2, w2 * np.sqrt(1 - z2 * z2))
    b1, b0 = -(p1 + p2), p1 * p2
    pp = lambda w: (1j * w) ** 2 + b1 * 1j * w + b0  # noqa: E731
    q = lambda w: abs(pp(w)) ** 2 + abs(pp(-w)) ** 2  # noqa: E731
    gamma = 0.5 * q(0.0)
    wr = max(w1, w2)
    beta = (0.5 * q(wr) - wr ** 4 - gamma) / wr ** 2
    d0 = np.sqrt(gamma)
    d1 = np.sqrt(beta + 2 * d0)
    return b1, b0, d1, d0


def complexified_A8(k):
    """The complexified mode zeta = B / (P+ D) xi, realified: [Re(p, v, a, l), Im(p, v, a, l)]."""
    b1, b0, d1, d0 = polys(k)
    wz0 = TWO_PI * k * SHAPE["fz"]
    c3, c2, c1, c0 = b1 + d1, b0 + b1 * d1 + d0, b0 * d1 + b1 * d0, b0 * d0
    Ac = np.zeros((4, 4), complex)
    Ac[0, 1] = Ac[1, 2] = Ac[2, 3] = 1
    Ac[3] = [-c0, -c1, -c2, -c3]
    cv = np.array([wz0, 1, 0, 0], complex)
    T = np.vstack([cv, cv @ Ac, cv @ Ac @ Ac, [1, 0, 0, 0]])
    Az = T @ Ac @ np.linalg.inv(T)
    return np.block([[Az.real, -Az.imag], [Az.imag, Az.real]])


def unit_g2(k):
    W = np.zeros((4, 4))
    W[2, 2] = 1
    return 1.0 / sl.solve_continuous_lyapunov(base_A(k), -W)[2, 2]


def coinner(k):
    """Input-normal realization of V = [N_r, N_i] / D (C++: coinner_input_normal)."""
    b1, b0, d1, d0 = polys(k)
    A0 = np.array([[-d1, 1.0], [-d0, 0.0]])
    B0 = np.array([[b1.real - d1, b1.imag], [b0.real - d0, b0.imag]])
    X = sl.solve_continuous_lyapunov(A0, -B0 @ B0.T)
    L = np.linalg.cholesky(X)
    Li = np.linalg.inv(L)
    A, B, C = Li @ A0 @ L, Li @ B0, np.array([[1.0, 0.0]]) @ L
    Q, R = np.linalg.qr(B)
    O = np.diag(np.sign(np.diag(R))) @ Q.T
    return O @ A @ O.T, O @ B, C @ O.T


def cascade6(k):
    Aw, Bw, Cw = coinner(k)
    A = np.zeros((6, 6))
    A[:4, :4] = base_A(k)
    A[2, 4:] = Cw[0]
    A[4:, 4:] = Aw
    G = np.zeros((6, 2))
    G[2, 0] = 1
    G[4:] = Bw
    return A, G


def obsv(A, C):
    return np.vstack([C @ np.linalg.matrix_power(A, i) for i in range(A.shape[0])])


# ---------------------------------------------------------------------------
# 1-5: structure
# ---------------------------------------------------------------------------

def unobservable(k):
    A8 = complexified_A8(k)
    ca = np.zeros((1, 8))
    ca[0, 2] = 1
    # scale-aware rank: normalize the observability rows
    O = obsv(A8, ca)
    O = O / np.linalg.norm(O, axis=1, keepdims=True)
    s = np.linalg.svd(O, compute_uv=False)
    U = sl.null_space(O, rcond=1e-10)
    return U, s


def analytic_U(k):
    """{0} x span(Re w, Im w), w = [B(r), r B(r), r^2 B(r), 1] at the roots of D."""
    _, _, d1, d0 = polys(k)
    wz0 = TWO_PI * k * SHAPE["fz"]
    r = np.roots([1.0, d1, d0])
    vecs = []
    for rr in r:
        B = rr + wz0
        w = np.array([B, rr * B, rr * rr * B, 1.0])
        vecs += [np.real(w), np.imag(w)]
    S = sl.orth(np.array(vecs).T)
    return np.vstack([np.zeros((4, S.shape[1])), S])


def tf(A, G, c, w):
    return c @ np.linalg.solve(1j * w * np.eye(A.shape[0]) - A, G)


def noise_K(S):
    """Per-axis-pair 2 x 2 noise blocks [[Re, -Im], [Im, Re]] (the C++ convention)."""
    K = np.zeros((6, 6))
    for a in range(3):
        for b in range(3):
            K[2 * a:2 * a + 2, 2 * b:2 * b + 2] = [[S[a, b].real, -S[a, b].imag], [S[a, b].imag, S[a, b].real]]
    return K


def hankel_order(A1, G1, c, S, tol=1e-7):
    A = np.kron(np.eye(3), A1)
    G = np.kron(np.eye(3), G1)
    P = sl.solve_continuous_lyapunov(A, -(G @ noise_K(S) @ G.T))
    P = 0.5 * (P + P.T)
    C = np.kron(np.eye(3), c[None, :])
    B = P @ C.T                                          # positive-real part Z(s) = C (sI - A)^-1 P C^T
    Wc = sl.solve_continuous_lyapunov(A, -B @ B.T)
    Wo = sl.solve_continuous_lyapunov(A.T, -C.T @ C)
    h = np.sort(np.sqrt(np.abs(np.linalg.eigvals(Wc @ Wo))))[::-1]
    return int((h > h[0] * tol).sum()), h / h[0]


def corr_condition(P):
    d = np.sqrt(np.diag(P))
    e = np.linalg.eigvalsh(P / np.outer(d, d))
    return float(e.max() / e.min())


# ---------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--ext", default="svg,pgf")
    ap.add_argument("--figures-only", action="store_true", help="draw the figures from the committed analysis.json")
    ap.add_argument("--outdir", default=None, help="figure directory (default: the results' figures/)")
    args = ap.parse_args(argv)
    outdir = Path(args.outdir) if args.outdir else FIG
    if args.figures_only:
        figures(json.loads((OUT / "analysis.json").read_text()), args.ext.split(","), outdir)
        return
    res = {}
    # 1. observability
    U1, s = unobservable(1.0)
    Ua = analytic_U(1.0)
    res["observability_singular_values_normalized"] = s.tolist()
    res["unobservable_dim"] = int(U1.shape[1])
    res["unobservable_real_half_norm"] = float(np.abs(U1[:4]).max())
    res["unobservable_vs_analytic_angle_deg"] = float(np.degrees(sl.subspace_angles(U1, Ua)).max())
    b1, b0, d1, d0 = polys(1.0)
    res["D_roots_kappa1"] = [complex(r).real for r in np.roots([1, d1, d0])]
    res["A8_eigenvalues_kappa1"] = [[float(e.real), float(e.imag)] for e in np.sort_complex(np.linalg.eigvals(complexified_A8(1.0)))]
    # 2. motion with kappa
    kap = np.exp(np.linspace(np.log(0.5), np.log(2.0), 61))
    ang = np.array([np.degrees(sl.subspace_angles(U1, analytic_U(k))) for k in kap])
    res["unobservable_angle_vs_kappa"] = {"kappa": kap.tolist(), "angle1_deg": ang[:, 0].tolist(),
                                          "angle2_deg": ang[:, 1].tolist()}
    res["unobservable_angle_per_percent_kappa_deg"] = float(np.degrees(sl.subspace_angles(U1, analytic_U(1.01))).max())
    # 3. cascade
    A8 = complexified_A8(1.0)
    G8 = np.zeros((8, 2))
    G8[2, 0] = G8[6, 1] = 1
    A6, G6 = cascade6(1.0)
    e = 0.0
    for idx in range(3):
        c8 = np.eye(8)[idx]
        c6 = np.eye(6)[idx]
        for w in np.linspace(0.05, 30, 300):
            h8 = tf(A8, G8, c8, w)
            e = max(e, float(np.abs(h8 - tf(A6, G6, c6, w)).max() / np.abs(h8).max()))
    res["cascade_tf_rel_error_vs_complexified"] = e
    Aw, Bw, Cw = coinner(1.0)
    res["coinner_lossless_residual"] = float(np.abs(Aw + Aw.T + Bw @ Bw.T).max() + np.abs(Cw + Bw[:, 0]).max())
    res["coinner_kappa_covariance_residual"] = float(max(np.abs(coinner(k)[0] - k * Aw).max() for k in (0.5, 2.0)))
    # 4. stochastic order
    rng = np.random.default_rng(1)
    X = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
    Sh = X @ X.conj().T
    Sw = np.array([[1, 0.3, 0.6j], [0.3, 0.5, 0.2j], [-0.6j, -0.2j, 1.0]])
    v = np.array([0.8, 0.5, 1j])
    cases = {"Hermitian, full rank": Sh, "wave-like (quadrature x-z, y-z)": Sw, "real (in-phase coupling)": Sh.real.astype(complex),
             "diagonal (coupling off)": np.diag([1.0, 2.0, 0.5]).astype(complex),
             "rank 1 (one long-crested wave)": np.outer(v, v.conj())}
    cp6 = np.eye(6)[0]
    orders = {}
    for name, S in cases.items():
        n, h = hankel_order(A6, G6, cp6, S)
        orders[name] = {"rank_Sigma": int(np.linalg.matrix_rank(S, 1e-9)), "rank_Im_Sigma": int(np.linalg.matrix_rank(S.imag, 1e-9)),
                        "predicted": int(4 * np.linalg.matrix_rank(S, 1e-9) + 2 * np.linalg.matrix_rank(S.imag, 1e-9)),
                        "hankel_order": n, "hankel_sv": h[:18].tolist()}
    res["stochastic_minimal_order"] = orders
    # 5. conditioning
    g2 = unit_g2(1.0)
    W8 = np.zeros((8, 8))
    W8[2, 2] = W8[6, 6] = g2
    P8 = sl.solve_continuous_lyapunov(A8, -W8)
    PC = sl.solve_continuous_lyapunov(A6, -g2 * G6 @ G6.T)
    W4 = np.zeros((4, 4))
    W4[2, 2] = g2
    P4 = sl.solve_continuous_lyapunov(base_A(1.0), -W4)
    res["stationary_correlation_condition"] = {"complexified, realified (8)": corr_condition(P8),
                                               "cascade (6)": corr_condition(PC), "A41 (4)": corr_condition(P4)}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "analysis.json").write_text(json.dumps(res, indent=1))
    print(json.dumps({k: v for k, v in res.items() if k != "unobservable_angle_vs_kappa"}, indent=1)[:4000])
    figures(res, args.ext.split(","), outdir)


def figures(res, exts, outdir=FIG):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if "pgf" in exts:
        plt.rcParams.update({"pgf.texsystem": "lualatex", "pgf.rcfonts": False, "font.family": "serif"})
    outdir.mkdir(parents=True, exist_ok=True)

    def save(fig, name):
        for e in exts:
            save_figure(fig, outdir / f"{name}.{e}", bbox_inches="tight")
        plt.close(fig)

    d = res["unobservable_angle_vs_kappa"]
    fig, ax = plt.subplots(figsize=(3.6, 2.7))
    ax.semilogx(d["kappa"], d["angle1_deg"], color="#d62728", label="largest principal angle")
    ax.semilogx(d["kappa"], d["angle2_deg"], color="#1f77b4", label="smallest")
    ax.set_xlabel(r"frequency scale $\kappa$")
    ax.set_ylabel(r"angle $U(\kappa)$ vs $U(1)$ (deg)")
    ax.legend(fontsize=7, frameon=False)
    ax.grid(alpha=0.3)
    from matplotlib.ticker import FixedLocator, NullLocator, ScalarFormatter
    ax.xaxis.set_major_locator(FixedLocator([0.5, 0.75, 1.0, 1.5, 2.0]))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_major_formatter(ScalarFormatter())
    fig.tight_layout()
    save(fig, "analysis_unobservable_subspace")

    o = res["stochastic_minimal_order"]
    fig, ax = plt.subplots(figsize=(7.2, 2.7))
    for name, v in o.items():
        h = np.array(v["hankel_sv"])
        ax.semilogy(np.arange(1, len(h) + 1), np.maximum(h, 1e-17), marker="o", ms=3,
                    label=f"{name}: {v['hankel_order']} (4 rk + 2 rk Im = {v['predicted']})")
    ax.axhline(1e-7, color="k", lw=0.6, ls=":")
    ax.set_xlabel("index")
    ax.set_ylabel("normalized Hankel s.v.")
    ax.set_xticks(range(1, 19))
    ax.legend(fontsize=6.5, frameon=False, loc="lower left")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    save(fig, "analysis_stochastic_order")


if __name__ == "__main__":
    main()
