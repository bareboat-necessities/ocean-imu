#!/usr/bin/env python3
"""Validated shipping prediction primitives for the canonical SEA3 full word.

This module has no source domain and cannot promote P3.  It accepts only
already-SEA3-derived per-sample quantities and encloses the shipping prediction
matrices used by ``Kalman3D_Wave_OU_III``:

* exact attitude/gyro-bias transition with constant body rate over one sample;
* the structured attitude/gyro-bias process covariance;
* the exact implemented [v,p,S,a_w] integrated-OU transition and covariance;
* active accelerometer-bias Gauss--Markov transition/covariance.

The caller remains responsible for proving that every supplied sample belongs
to the same phase-continuous compact SEA3 word.  No tuner rectangle, arbitrary
bounded-input box, replay, or independent tau/sigma source is created here.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Sequence

from ou3_interval import (
    Interval,
    IntervalMatrix,
    hull,
    matrix_add,
    matrix_identity,
    matrix_mul,
    matrix_transpose,
)
import ou3_full_process_ucc as PROCESS
import ou3_validated_transcendentals as VT

REPO = Path(__file__).resolve().parents[2]
CORE = REPO / "src" / "kalman_ou_common" / "KalmanOUCoreMath.h"
MEKF = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"
SCHEMA = 1
QUALIFICATION = "OU3_SEA3_SHIPPING_PREDICTION_MATRIX_PRIMITIVES"
BRANCH_X = 1.0e-2
SERIES_ORDER = 7


def I(x: float) -> Interval:
    return Interval.outward_bounds(float(x), float(x))


def _zero(rows: int, cols: int) -> IntervalMatrix:
    return [[I(0.0) for _ in range(cols)] for _ in range(rows)]


def _shape(A: Sequence[Sequence[Interval]]) -> tuple[int, int]:
    r = len(A)
    c = len(A[0]) if r else 0
    if any(len(row) != c for row in A):
        raise ValueError("ragged interval matrix")
    return r, c


def _scale(A: Sequence[Sequence[Interval]], c: Interval) -> IntervalMatrix:
    return [[c * A[i][j] for j in range(len(A[i]))] for i in range(len(A))]


def _spow(x: Interval, n: int) -> Interval:
    if n < 0:
        raise ValueError("nonnegative power required")
    y = I(1.0)
    for _ in range(n):
        y = y * x
    return y


def _matrix_hull(A: IntervalMatrix, B: IntervalMatrix) -> IntervalMatrix:
    if _shape(A) != _shape(B):
        raise ValueError("matrix hull shape mismatch")
    return [[hull(A[i][j], B[i][j]) for j in range(len(A[i]))] for i in range(len(A))]


def _sym(A: IntervalMatrix) -> IntervalMatrix:
    At = matrix_transpose(A)
    half = I(0.5)
    return [[half * (A[i][j] + At[i][j]) for j in range(len(A[i]))] for i in range(len(A))]


def _skew(v: Sequence[Interval]) -> IntervalMatrix:
    if len(v) != 3:
        raise ValueError("three-vector required")
    z = I(0.0)
    return [
        [z, -v[2], v[1]],
        [v[2], z, -v[0]],
        [-v[1], v[0], z],
    ]


def _series_matrix_functions(
    omega_body: Sequence[Interval],
    t: Interval,
    order: int = SERIES_ORDER,
) -> tuple[IntervalMatrix, IntervalMatrix, IntervalMatrix]:
    """Enclose R=exp(-[w]x t), B=int R, IB=int B by a matrix Taylor series.

    The remainder bound uses ||[w]x||_2=||w||_2 and is applied entrywise.  On
    the declared Normal-Live rate cap and 5 ms sample this is extremely small,
    while retaining the component intervals supplied by the SEA3 state.
    """
    if len(omega_body) != 3 or t.lo < 0.0 or order < 3:
        raise ValueError("invalid attitude series inputs")
    W = _skew(omega_body)
    power = matrix_identity(3)

    R = _zero(3, 3)
    B = _zero(3, 3)
    IB = _zero(3, 3)

    for n in range(order + 1):
        sign = -1.0 if n & 1 else 1.0
        cn_r = I(sign / math.factorial(n)) * _spow(t, n)
        cn_b = I(sign / math.factorial(n + 1)) * _spow(t, n + 1)
        cn_i = I(sign / math.factorial(n + 2)) * _spow(t, n + 2)
        R = matrix_add(R, _scale(power, cn_r))
        B = matrix_add(B, _scale(power, cn_b))
        IB = matrix_add(IB, _scale(power, cn_i))
        power = matrix_mul(power, W)

    wnorm = math.sqrt(sum(float(x.abs_upper()) ** 2 for x in omega_body))
    z = wnorm * t.hi
    if not math.isfinite(z) or z >= 0.1:
        raise ValueError("attitude series requires ||omega||*dt < 0.1")

    def tail(first_power: int, factorial_index: int, extra_t_power: int, ratio_den: int) -> float:
        first = (
            (wnorm ** first_power)
            * (t.hi ** extra_t_power)
            / math.factorial(factorial_index)
        )
        ratio = z / float(ratio_den)
        return math.nextafter(first / (1.0 - ratio), math.inf)

    n = order + 1
    r_rem = tail(n, n, n, n + 1)
    b_rem = tail(n, n + 1, n + 1, n + 2)
    ib_rem = tail(n, n + 2, n + 2, n + 3)

    def widen(A: IntervalMatrix, r: float) -> IntervalMatrix:
        e = Interval.outward_bounds(-r, r)
        return [[A[i][j] + e for j in range(3)] for i in range(3)]

    return widen(R, r_rem), widen(B, b_rem), widen(IB, ib_rem)


def attitude_gyro_bias_F_Q(
    omega_body: Sequence[Interval],
    h: Interval,
    gyro_variance_density_xyz: Sequence[Interval],
    gyro_bias_variance_density: Interval,
) -> tuple[IntervalMatrix, IntervalMatrix]:
    """Shipping exact/structured 6x6 F_AA,Q_AA for the default gyro-bias mode."""
    if h.lo <= 0.0 or len(gyro_variance_density_xyz) != 3:
        raise ValueError("positive sample and three gyro densities required")
    if gyro_bias_variance_density.lo <= 0.0:
        raise ValueError("gyro-bias density must be positive")

    R1, B1, IB1 = _series_matrix_functions(omega_body, h)
    hm = h * I(0.5)
    Rm, Bm, _ = _series_matrix_functions(omega_body, hm)
    R0 = matrix_identity(3)
    B0 = _zero(3, 3)

    F = matrix_identity(6)
    for i in range(3):
        for j in range(3):
            F[i][j] = R1[i][j]
            F[i][3 + j] = B1[i][j]

    Qg = _zero(3, 3)
    for i in range(3):
        if gyro_variance_density_xyz[i].lo <= 0.0:
            raise ValueError("gyro density variance must be positive")
        Qg[i][i] = gyro_variance_density_xyz[i]
    Qbg = _zero(3, 3)
    for i in range(3):
        Qbg[i][i] = gyro_bias_variance_density

    isotropic = all(
        gyro_variance_density_xyz[i].lo == gyro_variance_density_xyz[0].lo
        and gyro_variance_density_xyz[i].hi == gyro_variance_density_xyz[0].hi
        for i in range(1, 3)
    )
    if isotropic:
        I_R = _zero(3, 3)
        for i in range(3):
            I_R[i][i] = gyro_variance_density_xyz[0] * h
    else:
        def rq(R: IntervalMatrix) -> IntervalMatrix:
            return matrix_mul(matrix_mul(R, Qg), matrix_transpose(R))
        I_R = _scale(
            matrix_add(matrix_add(rq(R0), _scale(rq(Rm), I(4.0))), rq(R1)),
            h / I(6.0),
        )

    def bq(Bm_: IntervalMatrix) -> IntervalMatrix:
        return matrix_mul(matrix_mul(Bm_, Qbg), matrix_transpose(Bm_))
    I_BB = _scale(
        matrix_add(matrix_add(bq(B0), _scale(bq(Bm), I(4.0))), bq(B1)),
        h / I(6.0),
    )
    Qtt = matrix_add(I_R, I_BB)
    Qtb = matrix_mul(IB1, Qbg)
    Qbb = _scale(Qbg, h)

    Q = _zero(6, 6)
    for i in range(3):
        for j in range(3):
            Q[i][j] = Qtt[i][j]
            Q[i][3 + j] = Qtb[i][j]
            Q[3 + i][j] = Qtb[j][i]
            Q[3 + i][3 + j] = Qbb[i][j]
    return F, _sym(Q)


def _transition_axis_branch(tau: Interval, h: Interval, small: bool) -> IntervalMatrix:
    if tau.lo <= 0.0 or h.lo <= 0.0:
        raise ValueError("positive tau/h required")
    x = h / tau
    alpha = VT.exp_interval(-x)
    em1 = VT.expm1_interval(-x)
    phi_va = -tau * em1
    if small:
        x2, x3 = _spow(x, 2), _spow(x, 3)
        x4, x5 = _spow(x, 4), _spow(x, 5)
        phi_pa = _spow(tau, 2) * (
            I(0.5) * x2 - I(1.0 / 6.0) * x3 + I(1.0 / 24.0) * x4
        )
        phi_Sa = _spow(tau, 3) * (
            I(1.0 / 6.0) * x3 - I(1.0 / 24.0) * x4 + I(1.0 / 120.0) * x5
        )
    else:
        phi_pa = _spow(tau, 2) * (x + em1)
        phi_Sa = _spow(tau, 3) * (I(0.5) * _spow(x, 2) - x - em1)

    Phi = _zero(4, 4)
    Phi[0][0] = I(1.0)
    Phi[0][3] = phi_va
    Phi[1][0] = h
    Phi[1][1] = I(1.0)
    Phi[1][3] = phi_pa
    Phi[2][0] = I(0.5) * _spow(h, 2)
    Phi[2][1] = h
    Phi[2][2] = I(1.0)
    Phi[2][3] = phi_Sa
    Phi[3][3] = Interval(min(alpha.lo, 1.0), min(alpha.hi, 1.0))
    return Phi


def translation_axis_transition(tau: Interval, h: Interval) -> IntervalMatrix:
    x = h / tau
    if x.hi < BRANCH_X:
        return _transition_axis_branch(tau, h, True)
    if x.lo >= BRANCH_X:
        return _transition_axis_branch(tau, h, False)
    return _matrix_hull(
        _transition_axis_branch(tau, h, True),
        _transition_axis_branch(tau, h, False),
    )


def _axis_q_small(tau: Interval, h: Interval, sigma2: Interval) -> IntervalMatrix:
    inv = I(1.0) / tau
    hp = [I(1.0)] + [_spow(h, n) for n in range(1, 10)]
    ip = [I(1.0)] + [_spow(inv, n) for n in range(1, 10)]
    s = sigma2

    qvv = s * (I(2/3)*hp[3]*ip[1]-I(1/2)*hp[4]*ip[2]+I(7/30)*hp[5]*ip[3]-I(1/12)*hp[6]*ip[4]+I(31/1260)*hp[7]*ip[5]-I(1/160)*hp[8]*ip[6]+I(127/90720)*hp[9]*ip[7])
    qvp = s * (I(1/4)*hp[4]*ip[1]-I(1/6)*hp[5]*ip[2]+I(5/72)*hp[6]*ip[3]-I(1/45)*hp[7]*ip[4]+I(17/2880)*hp[8]*ip[5]-I(41/30240)*hp[9]*ip[6])
    qva = s * (hp[2]*ip[1]-hp[3]*ip[2]+I(7/12)*hp[4]*ip[3]-I(1/4)*hp[5]*ip[4]+I(31/360)*hp[6]*ip[5]-I(1/40)*hp[7]*ip[6]+I(127/20160)*hp[8]*ip[7]-I(17/12096)*hp[9]*ip[8])
    qpp = s * (I(1/10)*hp[5]*ip[1]-I(1/18)*hp[6]*ip[2]+I(5/252)*hp[7]*ip[3]-I(1/180)*hp[8]*ip[4]+I(17/12960)*hp[9]*ip[5])
    qpa = s * (I(1/3)*hp[3]*ip[1]-I(1/3)*hp[4]*ip[2]+I(11/60)*hp[5]*ip[3]-I(13/180)*hp[6]*ip[4]+I(19/840)*hp[7]*ip[5]-I(1/168)*hp[8]*ip[6]+I(247/181440)*hp[9]*ip[7])
    qaa = s * (I(2)*hp[1]*ip[1]-I(2)*hp[2]*ip[2]+I(4/3)*hp[3]*ip[3]-I(2/3)*hp[4]*ip[4]+I(4/15)*hp[5]*ip[5]-I(4/45)*hp[6]*ip[6]+I(8/315)*hp[7]*ip[7]-I(2/315)*hp[8]*ip[8]+I(4/2835)*hp[9]*ip[9])
    qvS = s * (I(1/15)*hp[5]*ip[1]-I(1/24)*hp[6]*ip[2]+I(41/2520)*hp[7]*ip[3]-I(7/1440)*hp[8]*ip[4]+I(109/90720)*hp[9]*ip[5])
    qpS = s * (I(1/36)*hp[6]*ip[1]-I(1/72)*hp[7]*ip[2]+I(13/2880)*hp[8]*ip[3]-I(1/864)*hp[9]*ip[4])
    qSS = s * (I(1/126)*hp[7]*ip[1]-I(1/288)*hp[8]*ip[2]+I(13/12960)*hp[9]*ip[3])
    qSa = s * (I(1/12)*hp[4]*ip[1]-I(1/12)*hp[5]*ip[2]+I(2/45)*hp[6]*ip[3]-I(1/60)*hp[7]*ip[4]+I(11/2240)*hp[8]*ip[5]-I(73/60480)*hp[9]*ip[6])
    return _sym([
        [qvv, qvp, qvS, qva],
        [qvp, qpp, qpS, qpa],
        [qvS, qpS, qSS, qSa],
        [qva, qpa, qSa, qaa],
    ])


def _axis_q_exact(tau: Interval, h: Interval, sigma2: Interval) -> IntervalMatrix:
    inv = I(1.0) / tau
    x = h * inv
    a = VT.exp_interval(-x)
    a2 = a.square()
    qc = I(2.0) * sigma2 * inv
    t2, t3 = _spow(tau, 2), _spow(tau, 3)
    t4, t5 = _spow(tau, 4), _spow(tau, 5)
    t6, t7 = _spow(tau, 6), _spow(tau, 7)
    x2, x3 = _spow(x, 2), _spow(x, 3)
    x4, x5 = _spow(x, 4), _spow(x, 5)

    K00 = t3 * (-a2 + I(4)*a + I(2)*x - I(3)) / I(2)
    K01 = t4 * (a2 + I(2)*a*(x-I(1)) + x2 - I(2)*x + I(1)) / I(2)
    K03 = t2 * (a2 - I(2)*a + I(1)) / I(2)
    K11 = t5 * (-a2/I(2) - I(2)*a*x + x3/I(3) - x2 + x + I(0.5))
    K13 = t3 * (-a2 - I(2)*a*x + I(1)) / I(2)
    K33 = tau * (I(1)-a2) / I(2)
    K02 = t5 * (-I(3)*a2 + I(3)*a*(x2+I(4)) + x3 - I(3)*x2 + I(6)*x - I(9)) / I(6)
    K12 = t6 * (a2/I(2) + a*(-x2+I(2)*x-I(2))/I(2) + x4/I(8) - x3/I(2) + x2 - x + I(0.5))
    K22 = t7 * (-a2/I(2) + a*x2 + I(2)*a + x5/I(20) - x4/I(4) + I(2)*x3/I(3) - x2 + x - I(1.5))
    K23 = t4 * (a2 - a*(x2+I(2)) + I(1)) / I(2)

    return _sym([
        [qc*K00, qc*K01, qc*K02, qc*K03],
        [qc*K01, qc*K11, qc*K12, qc*K13],
        [qc*K02, qc*K12, qc*K22, qc*K23],
        [qc*K03, qc*K13, qc*K23, qc*K33],
    ])


def translation_axis_process(tau: Interval, h: Interval, sigma2: Interval) -> IntervalMatrix:
    if tau.lo <= 0.0 or h.lo <= 0.0 or sigma2.lo <= 0.0:
        raise ValueError("positive tau/h/sigma2 required")
    x = h / tau
    if x.hi < BRANCH_X:
        return _axis_q_small(tau, h, sigma2)
    if x.lo >= BRANCH_X:
        return _axis_q_exact(tau, h, sigma2)
    return _matrix_hull(
        _axis_q_small(tau, h, sigma2),
        _axis_q_exact(tau, h, sigma2),
    )


def translation_F_Q(
    tau: Interval,
    h: Interval,
    sigma_aw_std_xyz: Sequence[Interval],
) -> tuple[IntervalMatrix, IntervalMatrix]:
    if len(sigma_aw_std_xyz) != 3:
        raise ValueError("three stationary a_w standard deviations required")
    Faxis = translation_axis_transition(tau, h)
    F = _zero(12, 12)
    Q = _zero(12, 12)
    idx = (0, 3, 6, 9)
    for axis in range(3):
        s = sigma_aw_std_xyz[axis]
        if s.lo <= 0.0:
            raise ValueError("stationary a_w standard deviation must be positive")
        Qaxis = translation_axis_process(tau, h, s.square())
        for i in range(4):
            for j in range(4):
                F[idx[i] + axis][idx[j] + axis] = Faxis[i][j]
                Q[idx[i] + axis][idx[j] + axis] = Qaxis[i][j]
    return F, _sym(Q)


def active_accel_bias_F_Q(
    h: Interval,
    tau_ba: Interval,
    q_ba_density: Interval,
) -> tuple[Interval, IntervalMatrix]:
    if h.lo <= 0.0 or tau_ba.lo <= 0.0 or q_ba_density.lo <= 0.0:
        raise ValueError("positive active-bias parameters required")
    phi = VT.exp_interval(-(h / tau_ba))
    x2 = I(2.0) * h / tau_ba
    qscale = I(-0.5) * tau_ba * VT.expm1_interval(-x2)
    q = q_ba_density * qscale
    Q = _zero(3, 3)
    for i in range(3):
        Q[i][i] = q
    return phi, Q


def build() -> dict:
    core = CORE.read_text(encoding="utf-8")
    mekf = MEKF.read_text(encoding="utf-8")
    process = PROCESS.build()
    pf = PROCESS.validate(process)
    if pf:
        raise RuntimeError(f"shipping process prerequisite failed: {pf}")
    parity = {
        "integrated_OU_transition_bound_to_shipping_core": "IntegratedOUChain<T,3>" in core and "static void transition" in core,
        "integrated_OU_process_bound_to_shipping_core": "static void process_covariance" in core and "regularize_psd_if_needed<T,4>(Qd);" in core,
        "attitude_exact_transition_bound_to_shipping": "rot_and_B_from_wt_(w, Ts, Rstep, Bstep);" in mekf,
        "attitude_structured_Q_bound_to_shipping": "Q_AA.template topRightCorner<3,3>()    = Qtb;" in mekf,
        "active_ba_exact_GM_bound_to_shipping": "const T qd_scale = -T(0.5) * tau_b * std::expm1(-T(2) * Ts / tau_b);" in mekf,
        "no_source_domain_created_here": True,
    }

    h = Interval(*process["configured_runtime"]["imu_dt_outward_interval_s"])
    c = process["source_constants"]
    qg = [I(float(x) * float(x)) for x in c["gyro_noise_density_rad_sqrt_s_per_axis"]]
    qb = I(float(c["gyro_bias_rw_variance_density"]))
    zero_rate = [I(0.0), I(0.0), I(0.0)]
    Faa, Qaa = attitude_gyro_bias_F_Q(zero_rate, h, qg, qb)
    Ft, Qt = translation_F_Q(I(1.0), h, [I(0.5), I(0.5), I(0.5)])
    phi_ba, Qba = active_accel_bias_F_Q(
        h,
        I(float(c["accel_bias_tau_s"])),
        I(float(c["accel_bias_process_variance_density"])),
    )
    smoke = {
        "attitude_shape": list(_shape(Faa)),
        "attitude_Q_shape": list(_shape(Qaa)),
        "translation_shape": list(_shape(Ft)),
        "translation_Q_shape": list(_shape(Qt)),
        "active_ba_phi": phi_ba.as_list(),
        "active_ba_Q_diag": [Qba[i][i].as_list() for i in range(3)],
        "all_required_shapes_present": (
            _shape(Faa) == (6, 6)
            and _shape(Qaa) == (6, 6)
            and _shape(Ft) == (12, 12)
            and _shape(Qt) == (12, 12)
            and 0.0 < phi_ba.lo <= phi_ba.hi <= 1.0
        ),
    }
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "trajectory_replay_used": False,
        "source_domain_created_here": False,
        "arbitrary_bounded_input_source_created_here": False,
        "independent_tau_sigma_source_created_here": False,
        "consumes_only_SEA3_derived_sample_coordinates": True,
        "shipping_source_parity": parity,
        "shipping_source_parity_pass": all(parity.values()),
        "process_prerequisite": process,
        "smoke": smoke,
        "validated_matrix_primitives_ready": bool(
            all(parity.values()) and smoke["all_required_shapes_present"]
        ),
        "shipping_PSD_cleanup_noop_over_full_SEA3_family_certified": False,
        "P3_promoted": False,
        "next_obligation": (
            "prove the shipping PSD cleanup is a no-op/enclosed over the complete SEA3 family, "
            "then call these matrices from the same phase-continuous SEA3 source executor"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    for key in (
        "consumes_only_SEA3_derived_sample_coordinates",
        "shipping_source_parity_pass",
        "validated_matrix_primitives_ready",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used",
        "source_domain_created_here",
        "arbitrary_bounded_input_source_created_here",
        "independent_tau_sigma_source_created_here",
        "shipping_PSD_cleanup_noop_over_full_SEA3_family_certified",
        "P3_promoted",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if not all(d.get("shipping_source_parity", {}).values()):
        f.append("shipping source parity failed")
    if d.get("smoke", {}).get("all_required_shapes_present") is not True:
        f.append("matrix primitive smoke failed")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build()
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "qualification": d["qualification"],
        "shipping_source_parity_pass": d["shipping_source_parity_pass"],
        "smoke": d["smoke"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
