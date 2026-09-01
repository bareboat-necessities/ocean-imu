#!/usr/bin/env python3
"""Outward full-matrix H prefix propagation for OU-III P5.

This is the matrix-valued numerical stage that follows the scalar joint-cell
screen.  It carries one source-complete normal-Live H cell through the shipping
prediction/S/accelerometer/magnetometer order with an actual 18x18 interval
covariance matrix.  At each accepted correction it rebuilds the implemented
H, R, innovation covariance, K, residual and K*r cell, applies the Joseph
posterior enclosure, then applies the immediate left-error reset congruence.
The physical attitude correction passed to the signed Cayley primitive is

    d = -E_theta K r,

not an independently chosen correction-norm bound.

The source parameter tuple is kept *joint*: tau, sigma_aw and R_S live in one
invariant cell for the complete one-second word.  This is deliberately broader
than the P3 fine partition but it is source complete under arbitrary tuner
movement inside the validated applied-parameter invariant; no lower endpoint
from one source cell is paired with an upper endpoint from another.

Two exact nonlinear reductions are retained:

* magnetometer radial residual has zero K action and only H_theta d_eff is
  propagated;
* accelerometer finite-angle eta is inserted as the exact effective a_w input.

The full covariance propagation is entrywise interval arithmetic.  Whenever a
3x3 fixed-pivot inverse is certified, that inverse is used directly.  If an
interval innovation box is too wide for fixed pivots, the code does *not* use a
floating-point inverse: it uses only the rigorous SPD fact S>=R to enclose
S^-1 entrywise.  That path is intentionally less sharp and is reported because
it often identifies the next subdivision target.

Accepted/rejected and due/not-due branches are hulled only after the accepted
shipping map (Joseph + immediate reset) has been evaluated.  Consequently the
matrix cell remains source complete without selecting favorable rejections.
The producer is fail-closed: a signed Cayley denominator crossing or a deployed
correction outside the currently validated correction-Cayley range is reported
as a numerical obstruction, never replaced by 1-|a||c|/4.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
from typing import Sequence

from ou3_interval import (
    Interval,
    hull,
    matrix_add,
    matrix_identity,
    matrix_mul,
    matrix_transpose,
)
from ou3_interval_linear_algebra import (
    IntervalPivotError,
    matrix_inverse_gauss_jordan,
    matrix_symmetric_hull,
)
import ou3_full_process_ucc as PROCESS
import ou3_p4_group_algebra as GROUP
import ou3_p5_effective_vector_input as VEFF
import ou3_p5_go_live_covariance_stage as GOLIVE
import ou3_p5_heading_handoff_contract as HEADING
import ou3_p5_signed_cayley_cell as SIGNED
import ou3_source_reachable_matrix_p3 as P3CELL
import ou3_validated_transcendentals as VT
import ou3_vector_uco_certificate as VECTOR

REPO = Path(__file__).resolve().parents[1]
MEKF = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
N = 18
TH = range(0, 3)
BG = range(3, 6)
V = range(6, 9)
P = range(9, 12)
SS = range(12, 15)
AW = range(15, 18)


def I(x: float) -> Interval:
    return Interval.point(float(x))


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def _zero(rows: int, cols: int) -> list[list[Interval]]:
    return [[I(0.0) for _ in range(cols)] for _ in range(rows)]


def _box(a: float) -> Interval:
    a = abs(float(a))
    return Interval(down(-a), up(a))


def _vec_box(a: float) -> list[Interval]:
    return [_box(a), _box(a), _box(a)]


def _mat_hull(A, B):
    if len(A) != len(B) or (A and len(A[0]) != len(B[0])):
        raise ValueError("matrix hull shape mismatch")
    return [[hull(A[i][j], B[i][j]) for j in range(len(A[0]))] for i in range(len(A))]


def _vec_hull(a: Sequence[Interval], b: Sequence[Interval]) -> list[Interval]:
    if len(a) != len(b):
        raise ValueError("vector hull shape mismatch")
    return [hull(x, y) for x, y in zip(a, b)]


def _mat_vec(A, x: Sequence[Interval]) -> list[Interval]:
    if not A or len(A[0]) != len(x):
        raise ValueError("matrix/vector shape mismatch")
    out = []
    for row in A:
        y = I(0.0)
        for a, b in zip(row, x):
            y = y + a*b
        out.append(y)
    return out


def _vec_add(a: Sequence[Interval], b: Sequence[Interval]) -> list[Interval]:
    if len(a) != len(b):
        raise ValueError("vector add shape mismatch")
    return [x+y for x, y in zip(a, b)]


def _vec_sub(a: Sequence[Interval], b: Sequence[Interval]) -> list[Interval]:
    if len(a) != len(b):
        raise ValueError("vector sub shape mismatch")
    return [x-y for x, y in zip(a, b)]


def _vec_neg(a: Sequence[Interval]) -> list[Interval]:
    return [-x for x in a]


def _norm_upper(v: Sequence[Interval]) -> float:
    s = 0.0
    for x in v:
        z = x.abs_upper()
        s = up(s + up(z*z))
    return up(math.sqrt(s))


def _intersect(a: Interval, b: Interval) -> Interval:
    lo = max(a.lo, b.lo)
    hi = min(a.hi, b.hi)
    if lo > hi:
        raise RuntimeError(f"rigorous interval enclosures became disjoint: {a} vs {b}")
    return Interval(lo, hi)


def _psd_tighten(Pm):
    """Use covariance PSD and Cauchy-Schwarz to tighten an entrywise box."""
    Pm = matrix_symmetric_hull(Pm)
    dhi = []
    for i in range(N):
        if not math.isfinite(Pm[i][i].hi) or Pm[i][i].hi < 0.0:
            raise RuntimeError("covariance diagonal upper is invalid")
        dhi.append(up(max(0.0, Pm[i][i].hi)))
    out = [[Pm[i][j] for j in range(N)] for i in range(N)]
    for i in range(N):
        out[i][i] = _intersect(out[i][i], Interval(0.0, dhi[i]))
        for j in range(i+1, N):
            b = up(math.sqrt(up(dhi[i]*dhi[j])))
            z = _intersect(hull(out[i][j], out[j][i]), Interval(-b, b))
            out[i][j] = z
            out[j][i] = z
    return out


def _posterior_loewner_entry_box(Pm):
    """Entrywise enclosure of any 0<=P+<=P using only prior diagonal uppers."""
    dhi = [up(max(0.0, Pm[i][i].hi)) for i in range(N)]
    out = _zero(N, N)
    for i in range(N):
        out[i][i] = Interval(0.0, dhi[i])
        for j in range(i+1, N):
            b = up(math.sqrt(up(dhi[i]*dhi[j])))
            out[i][j] = Interval(-b, b)
            out[j][i] = out[i][j]
    return out


def _reset_matrix(dx_theta: Sequence[Interval]):
    if len(dx_theta) != 3:
        raise ValueError("reset correction must be length three")
    x, y, z = dx_theta
    G = matrix_identity(N)
    h = I(0.5)
    G[0][1] = -h*z; G[0][2] = h*y
    G[1][0] = h*z;  G[1][2] = -h*x
    G[2][0] = -h*y; G[2][1] = h*x
    return G


def _reset_covariance(Pm, dx_theta):
    G = _reset_matrix(dx_theta)
    return _psd_tighten(matrix_mul(matrix_mul(G, Pm), matrix_transpose(G)))


def _innovation(Pm, H, R):
    PHt = matrix_mul(Pm, matrix_transpose(H))
    S = matrix_symmetric_hull(matrix_add(matrix_mul(H, PHt), R))
    return PHt, S


def _spd_inverse_enclosure(S, R):
    """Certified inverse: fixed pivots first, otherwise use only S>=R>0."""
    try:
        return matrix_inverse_gauss_jordan(S), "FIXED_PIVOT_INTERVAL_GAUSS_JORDAN"
    except (IntervalPivotError, ZeroDivisionError):
        rmin = min(R[i][i].lo for i in range(3))
        if not rmin > 0.0:
            raise RuntimeError("measurement covariance lost positive diagonal floor")
        b = up(1.0/rmin)
        out = _zero(3, 3)
        for i in range(3):
            out[i][i] = Interval(0.0, b)
            for j in range(i+1, 3):
                out[i][j] = Interval(-b, b)
                out[j][i] = out[i][j]
        return out, "SPD_S_GE_R_SPECTRAL_ENTRY_ENCLOSURE"


def _shipping_joseph(Pm, K, S, PHt):
    """Natural interval extension of shipping P-KCP-(KCP)' + KSK'."""
    out = [[Pm[i][j] for j in range(N)] for i in range(N)]
    for i in range(N):
        for j in range(i, N):
            kcp_ij = I(0.0)
            kcp_ji = I(0.0)
            for l in range(3):
                kcp_ij = kcp_ij + K[i][l]*PHt[j][l]
                kcp_ji = kcp_ji + K[j][l]*PHt[i][l]
            ksk = I(0.0)
            for a in range(3):
                for b in range(3):
                    ksk = ksk + K[i][a]*S[a][b]*K[j][b]
            z = Pm[i][j] - kcp_ij - kcp_ji + ksk
            out[i][j] = z
            out[j][i] = z
    # Both this natural extension and the exact Kalman Loewner inequality contain
    # the real posterior.  Their intersection is therefore a valid tightening.
    loose = _posterior_loewner_entry_box(Pm)
    for i in range(N):
        for j in range(N):
            out[i][j] = _intersect(out[i][j], loose[i][j])
    return _psd_tighten(out)


def _measurement_cell(Pm, H, R, r):
    PHt, S = _innovation(Pm, H, R)
    Sinv, inverse_backend = _spd_inverse_enclosure(S, R)
    K = matrix_mul(PHt, Sinv)
    dx = _mat_vec(K, r)
    Pj = _shipping_joseph(Pm, K, S, PHt)
    Pr = _reset_covariance(Pj, dx[0:3])
    return {
        "P_accepted": Pr,
        "K": K,
        "S": S,
        "r": list(r),
        "dx": dx,
        "inverse_backend": inverse_backend,
    }


def _measurement_branch_hull(Pm, e, c, H, R, r, *, allow_rejected: bool):
    cell = _measurement_cell(Pm, H, R, r)
    dx = cell["dx"]
    # E=R_true R_hat^T and R_hat+=Q(dx) R_hat imply E+=E Q(dx)^-1.
    # Invert a left product to reuse either validated SIGNED/QCOMP backend:
    # (Q(dx) E^-1)^-1 = E Q(dx)^-1.
    signed = dict(SIGNED.compose_cell(_vec_neg(c), dx[0:3]))
    signed["c_plus"] = _vec_neg(signed["c_plus"])
    if "a" in signed:
        signed["a"] = _vec_neg(signed["a"])
    signed["error_rotation_convention"] = "R_true R_hat^T"
    signed["physical_error_correction_side"] = "right"
    c_acc = signed["c_plus"]
    e_acc = list(e)
    for i in range(3, N):
        e_acc[i] = e[i] - dx[i]
    if allow_rejected:
        Pout = _psd_tighten(_mat_hull(Pm, cell["P_accepted"]))
        eout = _vec_hull(e, e_acc)
        cout = _vec_hull(c, c_acc)
    else:
        Pout, eout, cout = cell["P_accepted"], e_acc, c_acc
    return Pout, eout, cout, cell, signed


def _predict_estimator_mean(xhat, F):
    """Source mean at a prediction; local attitude coordinates stay reset."""
    out = [I(0.0) for _ in range(N)]
    for i in range(3, N):
        for j in range(3, N):
            out[i] = out[i] + F[i][j] * xhat[j]
    return out


def _update_estimator_mean(xhat, cell, *, allow_rejected: bool, bias_limit=None):
    """Carry the estimator mean separately from its physical truth error."""
    accepted = list(xhat)
    for i in range(3, N):
        accepted[i] = xhat[i] + cell["dx"][i]
    for i in range(3):
        accepted[i] = I(0.0)
    if bias_limit is not None and N == 21:
        # Euclidean projection scales each component towards zero. Retain
        # both the unchanged interior and the projected exterior branches.
        for i in range(18, 21):
            accepted[i] = _intersect(hull(I(0.0), accepted[i]), _box(bias_limit))
    return _vec_hull(xhat, accepted) if allow_rejected else accepted


def _predicted_force_upper(xhat, domain):
    """Bound the nominal force in H_a from the nominal a_w, not true a_w."""
    g = float(domain["startup"]["gravity_mps2"])
    return up(g + _norm_upper([xhat[i] for i in AW]))


def _source_cell() -> dict:
    s = P3CELL.source_schedule()
    return {
        "dt_s": float(s["dt_s"]),
        "tau_s": Interval.outward_bounds(*map(float, s["tau_applied_invariant_s"])),
        "sigma_aw_mps2": Interval.outward_bounds(*map(float, s["sigma_aw_applied_safety"])),
        "R_S_filter_std": Interval.outward_bounds(*map(float, s["R_S_applied_invariant"])),
        "pseudo_period_s": Interval.outward_bounds(
            max(float(s["pseudo_min_s"]), float(s["pseudo_ratio"])*float(s["tau_applied_invariant_s"][0])),
            min(float(s["pseudo_max_s"]), float(s["pseudo_ratio"])*float(s["tau_applied_invariant_s"][1])),
        ),
    }


def _source_pb0() -> float:
    text = MEKF.read_text(encoding="utf-8")
    m = re.search(r"T\s+Pq0\s*=\s*T\([^)]*\),\s*T\s+Pb0\s*=\s*T\(([0-9.eE+-]+)\)", text)
    if not m:
        raise RuntimeError("cannot extract H-mode gyro-bias covariance seed")
    return float(m.group(1))


def _rotation_step_box(rate_rad_s: float, h: float):
    theta = up(rate_rad_s*h)
    s = VT.sin_point(theta)
    co = VT.cos_point(theta)
    sin_hi = max(abs(s.lo), abs(s.hi))
    one_minus_c = up(1.0-co.lo)
    off = up(sin_hi + 0.5*one_minus_c)
    R = _zero(3, 3)
    for i in range(3):
        R[i][i] = Interval(co.lo, 1.0)
        for j in range(3):
            if i != j:
                R[i][j] = Interval(-off, off)
    B = _zero(3, 3)
    for i in range(3):
        B[i][i] = Interval(down(h*co.lo), up(h))
        for j in range(3):
            if i != j:
                B[i][j] = Interval(down(-h*off), up(h*off))
    return R, B


def _transition_and_Q(src: dict, domain: dict):
    h = float(src["dt_s"])
    tau = src["tau_s"]
    sigma = src["sigma_aw_mps2"]
    x = Interval.outward_bounds(h/tau.hi, h/tau.lo)
    alpha = VT.exp_interval(-x)
    em1 = VT.expm1_interval(-x)
    phi_va = -tau*em1
    phi_pa = tau.square()*(x+em1)
    phi_Sa = (tau.square()*tau)*(I(0.5)*x.square()-x-em1)

    F = matrix_identity(N)
    rate = float(domain["normal_live"]["body_rate_norm_upper_deg_s"])*math.pi/180.0
    Rstep, Bstep = _rotation_step_box(rate, h)
    for i in range(3):
        for j in range(3):
            F[i][j] = Rstep[i][j]
            F[i][3+j] = Bstep[i][j]
    for ax in range(3):
        iv, ip, iS, ia = 6+ax, 9+ax, 12+ax, 15+ax
        F[iv][iv] = I(1.0); F[iv][ia] = phi_va
        F[ip][iv] = I(h); F[ip][ip] = I(1.0); F[ip][ia] = phi_pa
        F[iS][iv] = I(0.5*h*h); F[iS][ip] = I(h); F[iS][iS] = I(1.0); F[iS][ia] = phi_Sa
        F[ia][ia] = Interval(alpha.lo, min(1.0, alpha.hi))

    Q = _zero(N, N)
    proc = PROCESS.build()["source_constants"]
    qg = float(proc["gyro_noise_density_rad_sqrt_s"])**2
    qb = float(proc["gyro_bias_rw_variance_density"])
    bb = up(qb*h*h*h/3.0)
    cross = up(qb*h*h/2.0)
    for i in range(3):
        Q[i][i] = Interval(down(qg*h), up(qg*h+bb))
        Q[3+i][3+i] = Interval.outward_bounds(qb*h, qb*h)
        for j in range(3):
            if i != j:
                Q[i][j] = Interval(-bb, bb)
            Q[i][3+j] = Interval(-cross, cross)
            Q[3+j][i] = Q[i][3+j]

    qscaled = P3CELL.step_scaled_q(x)
    D = [sigma*I(h), sigma*I(h*h), sigma*I(h*h*h), sigma]
    qaxis = [[D[i]*qscaled[i][j]*D[j] for j in range(4)] for i in range(4)]
    groups = (6, 9, 12, 15)
    for ax in range(3):
        ids = [g+ax for g in groups]
        for i in range(4):
            for j in range(4):
                Q[ids[i]][ids[j]] = qaxis[i][j]
    return F, _psd_tighten(Q), Rstep


def _initial_covariance(src: dict, domain_path: Path):
    go = GOLIVE.build(domain_path)["goLive_H_covariance_seed"]
    heading = HEADING.build(domain_path)
    # Timeout-gauged covariance contains the normal-gauged seed in Loewner
    # order, so one full H covariance family covers both gauged startup nodes.
    tilt = float(go["attitude_covariance_seed"]["tilt_variance"])
    yaw = float(go["attitude_covariance_seed"]["gauged_yaw_variance"])
    # The timeout heading contract bounds attitude more broadly, but the filter
    # covariance seed itself is the same configured gauged yaw sigma.
    _ = heading["gauged_timeout_subbranch"]
    Pm = _zero(N, N)
    for i in range(3):
        # Body orientation of the world-down yaw axis is source varying.  The
        # entire anisotropic attitude seed is enclosed entrywise by eigenvalues.
        Pm[i][i] = Interval(down(tilt), up(yaw))
        for j in range(i+1, 3):
            b = up(0.5*(yaw-tilt))
            Pm[i][j] = Interval(-b, b); Pm[j][i] = Pm[i][j]
    pb = _source_pb0()
    for i in BG:
        Pm[i][i] = Interval.outward_bounds(pb, pb)
    for i in V:
        Pm[i][i] = Interval.outward_bounds(go["P_vv_variance_per_axis"], go["P_vv_variance_per_axis"])
    for i in P:
        Pm[i][i] = Interval.outward_bounds(go["P_pp_variance_per_axis"], go["P_pp_variance_per_axis"])
    for i in SS:
        Pm[i][i] = Interval.outward_bounds(go["P_SS_variance_per_axis"], go["P_SS_variance_per_axis"])
    awv = src["sigma_aw_mps2"].square()
    for i in AW:
        Pm[i][i] = awv
    return _psd_tighten(Pm)


def _initial_error(domain: dict):
    b = domain["startup"]["physical_handoff_coordinate_bounds"]
    e = [I(0.0) for _ in range(N)]
    for idxs, key in (
        (BG, "gyro_bias_error_norm_upper_rad_s"),
        (V, "velocity_error_norm_upper_mps"),
        (P, "position_error_norm_upper_m"),
        (SS, "integral_displacement_error_norm_upper_m_s"),
        (AW, "latent_acceleration_error_norm_upper_mps2"),
    ):
        a = float(b[key])
        for i in idxs:
            e[i] = _box(a)
    return e


def _predict_error(e, F):
    out = list(e)
    # Gyro-bias error is a deterministic bounded input in P5; stochastic RW is
    # composed later.  Translation follows the exact H-mode linear chain.
    lin = list(V)+list(P)+list(SS)+list(AW)
    for i in lin:
        y = I(0.0)
        for j in lin:
            y = y + F[i][j]*e[j]
        out[i] = y
    return out


def _predict_c(c, Rstep, domain: dict, h: float):
    transported = _mat_vec(Rstep, c)
    bg = float(domain["startup"]["physical_handoff_coordinate_bounds"]["gyro_bias_error_norm_upper_rad_s"])
    wdist = float(domain["startup"]["effective_deterministic_gyro_transport_disturbance_upper_rad_s"])
    th = up(h*(bg+wdist))
    half = up(0.5*th)
    sn = VT.sin_point(half)
    co = VT.cos_point(half)
    if not co.lo > 0.0:
        raise RuntimeError("one-step prediction input reaches Cayley antipode")
    ca = up(2.0*sn.hi/co.lo)
    return GROUP.cayley_compose_left(_vec_box(ca), transported)


def _H_S():
    H = _zero(3, N)
    for i in range(3):
        H[i][12+i] = I(1.0)
    return H


def _H_acc(domain: dict, force_norm_upper=None):
    fhi = (float(domain["normal_live"]["specific_force_norm_upper_mps2"])
           if force_norm_upper is None else float(force_norm_upper))
    H = _zero(3, N)
    # -[f]_x: structural diagonal zeros are retained.
    b = _box(fhi)
    H[0][1] = b; H[0][2] = b
    H[1][0] = b; H[1][2] = b
    H[2][0] = b; H[2][1] = b
    # J_aw=R_wb.  Entrywise [-1,1] encloses every source orientation while the
    # exact orthogonality is separately consumed by the effective-input lemma.
    for i in range(3):
        for j in range(3):
            H[i][15+j] = Interval(-1.0, 1.0)
    return H


def _H_mag(domain: dict):
    mhi = float(domain["normal_live"]["magnetic_vector_norm_upper_uT"])
    H = _zero(3, N)
    b = _box(mhi)
    H[0][1] = b; H[0][2] = b
    H[1][0] = b; H[1][2] = b
    H[2][0] = b; H[2][1] = b
    return H


def _R_diag(std: float):
    r = Interval.outward_bounds(std*std, std*std)
    R = _zero(3, 3)
    for i in range(3): R[i][i] = r
    return R


def _R_S(src: dict):
    r = src["R_S_filter_std"].square()
    R = _zero(3, 3)
    factors = src.get("R_S_axis_std_factors")
    if factors is None:
        factors = P3CELL.source_rs_axis_std_factors()
    for i in range(3):
        R[i][i] = r * I(factors[i]).square()
    return R


def _acc_residual(e, c, domain: dict, q_hi: float, force_norm_upper=None):
    H = _H_acc(domain, force_norm_upper)
    fhi = (float(domain["normal_live"]["specific_force_norm_upper_mps2"])
           if force_norm_upper is None else float(force_norm_upper))
    aw_hi = max(e[i].abs_upper() for i in AW)
    eta = up(
        VEFF.accel_attitude_eta_per_vector_norm_upper(q_hi)*fhi
        + VEFF.accel_latent_cross_gain_upper(q_hi)*aw_hi
    )
    z = [I(0.0) for _ in range(N)]
    for i in range(3): z[i] = c[i]
    for i in AW: z[i] = e[i] + _box(eta)
    r = _mat_vec(H, z)
    ba = float(domain["startup"]["physical_handoff_coordinate_bounds"]["accelerometer_bias_error_norm_upper_mps2"])
    r = [x + _box(ba) for x in r]
    return H, r, eta


def _mag_effective(c, q_lo: float, q_hi: float):
    # Exact d_m=A c_perp-B alpha(c_perp x vhat), evaluated over every unit
    # direction using the box vhat_i in [-1,1], then intersected with the exact
    # nonexpansive norm consequence from the same identity.
    vhat = [Interval(-1.0, 1.0) for _ in range(3)]
    alpha = GROUP.dot(c, vhat)
    cperp = [c[i]-alpha*vhat[i] for i in range(3)]
    A = Interval.outward_bounds(
        4.0/(4.0+q_hi*q_hi),
        4.0/(4.0+q_lo*q_lo),
    )
    B = Interval.outward_bounds(
        2.0/(4.0+q_hi*q_hi),
        2.0/(4.0+q_lo*q_lo),
    )
    cp_cross_v = GROUP.cross(cperp, vhat)
    d = [A*cperp[i]-B*alpha*cp_cross_v[i] for i in range(3)]
    gain_hi = VEFF.mag_effective_coordinate_gain_upper(q_lo)
    cap = up(gain_hi*q_hi)
    return [_intersect(x, Interval(-cap, cap)) for x in d]


def _mag_residual(c, domain: dict, q_lo: float, q_hi: float):
    H = _H_mag(domain)
    d_eff = _mag_effective(c, q_lo, q_hi)
    r = _mat_vec(H, d_eff)
    return H, r, d_eff


def _serialize_vec(v):
    return [x.as_list() for x in v]


def _matrix_summary(Pm):
    diag = [Pm[i][i].as_list() for i in range(N)]
    cross = 0.0
    for i in range(N):
        for j in range(N):
            if i != j: cross = max(cross, Pm[i][j].abs_upper())
    return {
        "dimension": N,
        "diagonal_intervals": diag,
        "max_offdiagonal_abs_upper": cross,
        "theta_S_cross_abs_upper": max(Pm[i][j].abs_upper() for i in TH for j in SS),
        "theta_aw_cross_abs_upper": max(Pm[i][j].abs_upper() for i in TH for j in AW),
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain_path = Path(domain_path).resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("full H prefix domain must not be trajectory fitted")
    if domain.get("configured_runtime", {}).get("imu_lever_arm_enabled") is not False:
        raise RuntimeError("schema-1 full H prefix requires configured lever arm disabled")

    go = GOLIVE.build(domain_path)
    heading = HEADING.build(domain_path)
    veff = VEFF.build(domain_path)
    vector = VECTOR.build()
    process = PROCESS.build()
    failures = [f"goLive: {x}" for x in GOLIVE.validate(go)]
    failures += [f"heading: {x}" for x in HEADING.validate(heading)]
    failures += [f"effective-vector: {x}" for x in VEFF.validate(veff)]
    failures += [f"vector: {x}" for x in VECTOR.validate(vector)]
    failures += [f"process: {x}" for x in PROCESS.validate(process)]

    src = _source_cell()
    F, Q, Rstep = _transition_and_Q(src, domain)
    Pm = _initial_covariance(src, domain_path)
    e = _initial_error(domain)
    xhat = [I(0.0) for _ in range(N)]  # first goLive follows an untouched zero constructor mean

    q0 = float(heading["gauged_timeout_subbranch"]["full_attitude_cayley_norm_upper"])
    c = _vec_box(q0)
    q_chart = 8.0
    h = float(src["dt_s"])
    Tword = float(domain["normal_live"]["vector_pe_recurrence_window_s"])
    samples = int(math.ceil(Tword/h)) + 2
    vc = vector["configured_measurement_bounds"]
    Racc = _R_diag(float(vc["acc_measurement_std_mps2"]))
    Rmag = _R_diag(float(vc["mag_measurement_std_uT"]))
    RS = _R_S(src)

    inverse_counts = {"FIXED_PIVOT_INTERVAL_GAUSS_JORDAN": 0, "SPD_S_GE_R_SPECTRAL_ENTRY_ENCLOSURE": 0}
    first_failure = None
    max_q = _norm_upper(c)
    first_S_done = False
    first_S_sample = None
    extrema = {
        "max_cayley_norm_upper": max_q,
        "max_S_correction_norm_upper": 0.0,
        "max_acc_correction_norm_upper": 0.0,
        "max_mag_correction_norm_upper": 0.0,
        "max_acc_effective_aw_eta_norm_upper": 0.0,
    }
    last_cells = {}

    for k in range(samples):
        try:
            Pm = _psd_tighten(matrix_add(matrix_mul(matrix_mul(F, Pm), matrix_transpose(F)), Q))
            e = _predict_error(e, F)
            xhat = _predict_estimator_mean(xhat, F)
            c = _predict_c(c, Rstep, domain, h)
            qnow = _norm_upper(c)
            max_q = max(max_q, qnow)
            if qnow >= q_chart:
                raise RuntimeError(f"prediction prefix leaves q<={q_chart} chart: q={qnow}")

            # Source-complete pseudo phase: both due and not-due remain possible
            # until the branch is taken.  Evaluate the due shipping map first,
            # then hull it with the identity branch.
            HS = _H_S()
            rS = [-xhat[12+i] for i in range(3)]
            Pm, e, c, Scell, Ssigned = _measurement_branch_hull(Pm, e, c, HS, RS, rS, allow_rejected=True)
            inverse_counts[Scell["inverse_backend"]] += 1
            xhat = _update_estimator_mean(xhat, Scell, allow_rejected=True)
            ds = _norm_upper(_vec_neg(Scell["dx"][0:3]))
            extrema["max_S_correction_norm_upper"] = max(extrema["max_S_correction_norm_upper"], ds)
            last_cells["S"] = {
                "sample": k,
                "inverse_backend": Scell["inverse_backend"],
                "r": _serialize_vec(Scell["r"]),
                "d": _serialize_vec(_vec_neg(Scell["dx"][0:3])),
                "signed_denominator": Ssigned["denominator"].as_list(),
            }
            if not first_S_done:
                first_S_done = True
                first_S_sample = k

            qnow = _norm_upper(c)
            if qnow >= q_chart:
                raise RuntimeError(f"S prefix leaves q<={q_chart} chart: q={qnow}")

            Hacc, racc, eta = _acc_residual(
                e, c, domain, min(q_chart, qnow), _predicted_force_upper(xhat, domain))
            Pm, e, c, Acell, Asigned = _measurement_branch_hull(Pm, e, c, Hacc, Racc, racc, allow_rejected=True)
            inverse_counts[Acell["inverse_backend"]] += 1
            xhat = _update_estimator_mean(xhat, Acell, allow_rejected=True)
            da = _norm_upper(_vec_neg(Acell["dx"][0:3]))
            extrema["max_acc_correction_norm_upper"] = max(extrema["max_acc_correction_norm_upper"], da)
            extrema["max_acc_effective_aw_eta_norm_upper"] = max(extrema["max_acc_effective_aw_eta_norm_upper"], eta)
            last_cells["accelerometer"] = {
                "sample": k,
                "inverse_backend": Acell["inverse_backend"],
                "r": _serialize_vec(Acell["r"]),
                "d": _serialize_vec(_vec_neg(Acell["dx"][0:3])),
                "signed_denominator": Asigned["denominator"].as_list(),
                "effective_aw_eta_norm_upper": eta,
            }
            qnow = _norm_upper(c)
            if qnow >= q_chart:
                raise RuntimeError(f"accelerometer prefix leaves q<={q_chart} chart: q={qnow}")

            # The 25 Hz packet phase is asynchronous.  Allowing a packet branch
            # every IMU sample is a source-complete overapproximation; rejected
            # and not-due are both contained by the identity branch hull.
            qlo = 0.0
            Hmag, rmag, deff = _mag_residual(c, domain, qlo, min(q_chart, qnow))
            Pm, e, c, Mcell, Msigned = _measurement_branch_hull(Pm, e, c, Hmag, Rmag, rmag, allow_rejected=True)
            inverse_counts[Mcell["inverse_backend"]] += 1
            xhat = _update_estimator_mean(xhat, Mcell, allow_rejected=True)
            dm = _norm_upper(_vec_neg(Mcell["dx"][0:3]))
            extrema["max_mag_correction_norm_upper"] = max(extrema["max_mag_correction_norm_upper"], dm)
            last_cells["magnetometer"] = {
                "sample": k,
                "inverse_backend": Mcell["inverse_backend"],
                "r": _serialize_vec(Mcell["r"]),
                "d_eff": _serialize_vec(deff),
                "d": _serialize_vec(_vec_neg(Mcell["dx"][0:3])),
                "signed_denominator": Msigned["denominator"].as_list(),
                "radial_K_action_exact_zero": True,
            }
            qnow = _norm_upper(c)
            max_q = max(max_q, qnow)
            extrema["max_cayley_norm_upper"] = max_q
            if qnow >= q_chart:
                raise RuntimeError(f"magnetometer prefix leaves q<={q_chart} chart: q={qnow}")

            # Pending a_w covariance floor can be applied on the next prediction.
            # For the configured uncorrelated/isotropic a_w source target, the
            # eigenvalue map is p_i -> max(p_i,sigma^2); cross blocks are unchanged.
            if (k+1) % max(1, int(round(0.1/h))) == 0:
                u = max(src["sigma_aw_mps2"].square().hi, max(Pm[i][i].hi for i in AW))
                for i in AW:
                    Pm[i][i] = Interval(0.0, up(u))
                    for j in AW:
                        if i != j:
                            b = up(u)
                            Pm[i][j] = _intersect(Pm[i][j], Interval(-b, b))
                Pm = _psd_tighten(Pm)
        except Exception as exc:  # fail closed with exact operation/sample witness
            first_failure = {
                "sample": k,
                "operation": "prediction/S/accelerometer/magnetometer prefix",
                "reason": f"{type(exc).__name__}: {exc}",
                "cayley_norm_upper_before_failure": _norm_upper(c),
                "covariance": _matrix_summary(Pm),
            }
            break

    closed = first_failure is None and first_S_done and max_q < q_chart
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_FULL_18X18_H_PREFIX_INTERVAL_CELL_PROPAGATION",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "mode": "H",
        "dimension": N,
        "source_cell": {
            "tau_s": src["tau_s"].as_list(),
            "sigma_aw_mps2": src["sigma_aw_mps2"].as_list(),
            "R_S_filter_std": src["R_S_filter_std"].as_list(),
            "pseudo_period_s": src["pseudo_period_s"].as_list(),
            "joint_invariant_cell_over_complete_word": True,
        },
        "shipping_order_retained": ["prediction", "S_or_identity", "accelerometer_or_identity", "magnetometer_or_identity", "immediate_reset_after_each_accepted_correction"],
        "full_18x18_covariance_propagated": True,
        "H_R_S_K_r_d_eff_recomputed_in_same_prefix_cell": True,
        "shipping_Joseph_update_used": True,
        "immediate_left_error_reset_congruence_used": True,
        "physical_attitude_correction_is_minus_Etheta_Kr": True,
        "signed_cayley_primitive_consumes_actual_interval_d": True,
        "signed_a_dot_c_replaced_by_independent_abs_product": False,
        "magnetometer_radial_K_action_exact_zero": True,
        "standalone_vector_eta_penalty_used": False,
        "accelerometer_effective_aw_input_used": True,
        "source_complete_rejection_identity_hulls": True,
        "word_horizon_s": Tword,
        "word_samples_upper": samples,
        "first_S_branch_evaluated": first_S_done,
        "first_S_sample_index": first_S_sample,
        "inverse_backend_counts": inverse_counts,
        "q_chart_upper": q_chart,
        "max_reached_cayley_norm_upper": max_q,
        "smaller_source_reachable_chart_upper": up(max_q) if closed else None,
        "numerical_extrema": extrema,
        "last_prefix_cells": last_cells,
        "final_covariance": _matrix_summary(Pm),
        "first_failure": first_failure,
        "P5_FULL_H_PREFIX_MATRIX_CERTIFICATE": "PASS" if closed and not failures else "NOT_ESTABLISHED",
        "complete_q_le_8_prefix_family_closed": bool(closed and not failures),
        "next_obligation": (
            "promote the complete gauged H P5 numerical word and set N_H_words from the certified overlap"
            if closed and not failures else
            "subdivide the first reported full-matrix prefix obstruction (source/vector direction/branch phase) while retaining the same P,H,R,S,K,r,d_eff and signed-reset calculus"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    for k in (
        "source_generated_not_trajectory_fit", "full_18x18_covariance_propagated",
        "H_R_S_K_r_d_eff_recomputed_in_same_prefix_cell", "shipping_Joseph_update_used",
        "immediate_left_error_reset_congruence_used", "physical_attitude_correction_is_minus_Etheta_Kr",
        "signed_cayley_primitive_consumes_actual_interval_d", "magnetometer_radial_K_action_exact_zero",
        "accelerometer_effective_aw_input_used", "source_complete_rejection_identity_hulls",
    ):
        if d.get(k) is not True:
            failures.append(f"{k} is not true")
    for k in ("source_replay_used", "filter_changed", "signed_a_dot_c_replaced_by_independent_abs_product", "standalone_vector_eta_penalty_used"):
        if d.get(k) is not False:
            failures.append(f"{k} is not false")
    if int(d.get("dimension", 0)) != 18:
        failures.append("full H covariance dimension is not 18")
    if not isinstance(d.get("inverse_backend_counts"), dict):
        failures.append("inverse backend accounting missing")
    # Numerical non-closure is not a validation error.  The tool must remain
    # executable and fail closed with a concrete first obstruction.
    if d.get("complete_q_le_8_prefix_family_closed") is False and d.get("first_failure") is None:
        failures.append("non-closed prefix family lacks a first failure witness")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build(args.domain.resolve())
    vf = validate(out)
    out["validation_pass"] = not vf
    out["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": out["P5_FULL_H_PREFIX_MATRIX_CERTIFICATE"],
        "q8_closed": out["complete_q_le_8_prefix_family_closed"],
        "max_q": out["max_reached_cayley_norm_upper"],
        "smaller_chart": out["smaller_source_reachable_chart_upper"],
        "inverse_backends": out["inverse_backend_counts"],
        "first_failure": out["first_failure"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
