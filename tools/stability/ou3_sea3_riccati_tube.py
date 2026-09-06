#!/usr/bin/env python3
"""Quantitative SEA3 moving-Riccati tube for canonical OU-III P3.

This producer does *not* construct a source-word language.  The previous P3
route followed long histories through an 800-state tuner abstraction.  Here the
only interval cells cover the adaptive state that is active on the current
sample.  The covariance ceiling is source-uniform over a finite observation
window, so no predecessor identity is needed.

For one shipping linearized step, after all measurements that occur on that
sample,

    e+ = Phi e,
    P+ = Phi P Phi' + Omega.

A finite-memory estimator built from recurrent vector/S observations gives a
source-uniform covariance ceiling P+ <= Pbar.  In a diagonal similarity D tied
to the *current* OU scale, exact integrated-OU process arithmetic plus the
Joseph information inequality gives

    Omega >= rho D D'.

For any PSD Pbar, diagonal variance bounds u_i imply

    D^-1 Pbar D^-T <= trace(D^-1 Pbar D^-T) I
                   <= sum_i u_i / d_i^2 I.

Consequently

    Omega >= delta P+,   delta = rho / sum_i(u_i/d_i^2),

and the moving shipping-covariance metric V=e'P^-1e contracts by at least
1-delta.  Parameter changes require no history matching: P is itself the
moving metric and the current process/Joseph lower bound uses only the schedule
active on this sample.

The finite-memory ceiling intentionally uses the whole SEA3 dynamic invariant
for nuisance/process upper bounds.  Local interval cells are used only where a
current-sample lower comparison needs them.  This makes the construction valid
for time-varying adaptation without assuming a frozen tuner over the window.
"""
from __future__ import annotations

import argparse
import functools
import json
import math
import re
from pathlib import Path

from ou3_interval import (
    Interval,
    hull,
    matrix_mul,
    matrix_transpose,
    symmetric_positive_definite_ldlt,
)
from ou3_interval_linear_algebra import (
    matrix_inverse_gauss_jordan,
    matrix_symmetric_hull,
)
import ou3_full_process_ucc as PROCESS
import ou3_sea3_dynamic_source_certificate as DYNAMIC
import ou3_source_domain_contract as SOURCE
import ou3_validated_transcendentals as VT
import ou3_vector_uco_certificate as VECTOR

REPO = Path(__file__).resolve().parents[2]
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_SEA3_SOURCE_UNIFORM_MOVING_RICCATI_TUBE"
USEFUL_GATE = 1.0e-18
BRANCH_X = 1.0e-2


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def pos(x, label: str) -> float:
    y = float(x)
    if not math.isfinite(y) or y <= 0.0:
        raise RuntimeError(f"{label} must be finite positive, got {x!r}")
    return y


def I(x: float) -> Interval:
    return Interval.outward_bounds(float(x), float(x))


def ipow(x: Interval, n: int) -> Interval:
    y = Interval.point(1.0)
    for _ in range(n):
        y = y * x
    return y


def geom_edges(lo: float, hi: float, n: int) -> list[float]:
    if not (0.0 < lo < hi) or n < 1:
        raise ValueError("bad geometric grid")
    r = (hi / lo) ** (1.0 / n)
    out = [lo]
    for _ in range(n - 1):
        out.append(out[-1] * r)
    out.append(hi)
    return out


def interval_cells(edges: list[float]) -> list[Interval]:
    return [Interval.outward_bounds(edges[i], edges[i + 1]) for i in range(len(edges) - 1)]


def poly(x: Interval, terms: tuple[tuple[int, float], ...]) -> Interval:
    y = Interval.point(0.0)
    for n, c in terms:
        y = y + I(c) * ipow(x, n)
    return y


def _qbar_branch(x: Interval, small: bool) -> list[list[Interval]]:
    """Dimensionless exact integrated-OU covariance for [v,p,S,a]."""
    if small:
        qvv = poly(x, ((3,2/3),(4,-1/2),(5,7/30),(6,-1/12),(7,31/1260),(8,-1/160),(9,127/90720)))
        qvp = poly(x, ((4,1/4),(5,-1/6),(6,5/72),(7,-1/45),(8,17/2880),(9,-41/30240)))
        qva = poly(x, ((2,1),(3,-1),(4,7/12),(5,-1/4),(6,31/360),(7,-1/40),(8,127/20160),(9,-17/12096)))
        qpp = poly(x, ((5,1/10),(6,-1/18),(7,5/252),(8,-1/180),(9,17/12960)))
        qpa = poly(x, ((3,1/3),(4,-1/3),(5,11/60),(6,-13/180),(7,19/840),(8,-1/168),(9,247/181440)))
        qaa = poly(x, ((1,2),(2,-2),(3,4/3),(4,-2/3),(5,4/15),(6,-4/45),(7,8/315),(8,-2/315),(9,4/2835)))
        qvS = poly(x, ((5,1/15),(6,-1/24),(7,41/2520),(8,-7/1440),(9,109/90720)))
        qpS = poly(x, ((6,1/36),(7,-1/72),(8,13/2880),(9,-1/864)))
        qSS = poly(x, ((7,1/126),(8,-1/288),(9,13/12960)))
        qSa = poly(x, ((4,1/12),(5,-1/12),(6,2/45),(7,-1/60),(8,11/2240),(9,-73/60480)))
    else:
        a = VT.exp_interval(-x)
        a2 = a.square()
        one, two, three, four = I(1), I(2), I(3), I(4)
        x2, x3, x4, x5 = ipow(x,2), ipow(x,3), ipow(x,4), ipow(x,5)
        qvv = -a2 + four*a + two*x - three
        qvp = a2 + two*a*(x-one) + x2 - two*x + one
        qva = a2 - two*a + one
        qpp = -a2 - four*a*x + I(2/3)*x3 - two*x2 + two*x + one
        qpa = -a2 - two*a*x + one
        qaa = one - a2
        qvS = -a2 + a*(x2+four) + (x3-three*x2+I(6)*x-I(9))/three
        qpS = a2 + a*(-x2+two*x-two) + I(1/4)*x4 - x3 + two*x2 - two*x + one
        qSS = -a2 + two*a*x2 + four*a + I(1/10)*x5 - I(1/2)*x4 + I(4/3)*x3 - two*x2 + two*x - three
        qSa = a2 - a*(x2+two) + one
    return [
        [qvv,qvp,qvS,qva],
        [qvp,qpp,qpS,qpa],
        [qvS,qpS,qSS,qSa],
        [qva,qpa,qSa,qaa],
    ]


def qbar_integrated_ou(x: Interval) -> list[list[Interval]]:
    if x.lo <= 0.0 or x.hi > VT.MAX_ABS_ARGUMENT:
        raise ValueError("x=h/tau outside validated range")
    if x.hi < BRANCH_X:
        return _qbar_branch(x, True)
    if x.lo >= BRANCH_X:
        return _qbar_branch(x, False)
    left_hi = math.nextafter(BRANCH_X, -math.inf)
    families = []
    if x.lo <= left_hi:
        families.append(_qbar_branch(Interval(x.lo, left_hi), True))
    if x.hi >= BRANCH_X:
        families.append(_qbar_branch(Interval(BRANCH_X, x.hi), False))
    return [[hull(*(A[i][j] for A in families)) for j in range(4)] for i in range(4)]


def step_scaled_q(x: Interval) -> list[list[Interval]]:
    """D^-1 Q D^-T for D=diag(sigma*h,sigma*h^2,sigma*h^3,sigma)."""
    q = qbar_integrated_ou(x)
    powers = (1, 2, 3, 0)
    return [[q[i][j] / ipow(x, powers[i] + powers[j]) for j in range(4)] for i in range(4)]


def minus_rho(A, rho: float):
    return [[A[i][j] - I(rho if i == j else 0.0) for j in range(len(A))] for i in range(len(A))]


def certified_rho(A) -> float:
    ok, _ = symmetric_positive_definite_ldlt(A)
    if not ok:
        return 0.0
    hi = min(A[i][i].lo for i in range(len(A)))
    lo = 0.0
    for _ in range(44):
        mid = 0.5 * (lo + hi)
        ok, _ = symmetric_positive_definite_ldlt(minus_rho(A, mid))
        if ok:
            lo = mid
        else:
            hi = mid
    return down(lo)


def split_x_cell(x: Interval, depth: int = 0) -> list[tuple[Interval, float]]:
    rho = certified_rho(step_scaled_q(x))
    if rho > 0.0:
        return [(x, rho)]
    if depth >= 14:
        raise RuntimeError(f"cannot certify scaled OU process cell {x.as_list()}")
    mid = math.sqrt(x.lo * x.hi)
    return split_x_cell(Interval.outward_bounds(x.lo, mid), depth + 1) + split_x_cell(Interval.outward_bounds(mid, x.hi), depth + 1)


def diagonal_dominator(A) -> list[float]:
    A = matrix_symmetric_hull(A)
    out = []
    for i in range(len(A)):
        d = max(0.0, A[i][i].hi)
        for j in range(len(A)):
            if i != j:
                d = up(d + A[i][j].abs_upper())
        out.append(up(d))
    return out


@functools.lru_cache(maxsize=128)
def integrator_inverse(gap: float, spacing: float):
    g, s = pos(gap, "gap"), pos(spacing, "spacing")
    if s <= g:
        raise RuntimeError("S observation windows overlap")
    t = (
        Interval.outward_bounds(0.0, g),
        Interval.outward_bounds(s, s + g),
        Interval.outward_bounds(2 * s, 2 * s + g),
    )
    B = [[I(1), ti, I(0.5) * ti.square()] for ti in t]  # [S,p,v]
    return matrix_inverse_gauss_jordan(B)


def _axis_factors() -> list[float]:
    text = WRAPPER.read_text(encoding="utf-8")
    if "RSAdaptationLaw rs_law_ = RSAdaptationLaw::SpectralMSE;" not in text:
        raise RuntimeError("canonical SEA3 tube requires deployed SpectralMSE R_S branch")
    out = []
    for name in ("R_S_x_factor_", "R_S_y_factor_"):
        m = re.search(rf"float\s+{name}\s*=\s*([0-9.eE+-]+)f", text)
        if not m:
            raise RuntimeError(f"cannot extract deployed {name}")
        out.append(float(m.group(1)))
    return out + [1.0]


def _declared_vector_alpha6(live: dict, vector: dict) -> float:
    """Use the stronger theorem-domain PE values with the validated packet proof."""
    base = vector["operating_envelope"]
    f = pos(live["specific_force_norm_lower_mps2"], "force floor")
    m = pos(live["magnetic_vector_norm_lower_uT"], "mag floor")
    s = pos(live["vector_sine_separation_lower"], "separation floor")
    rate = pos(live["body_rate_norm_upper_deg_s"], "rate ceiling")
    if f < float(base["specific_force_norm_lower_mps2"]) or m < float(base["magnetic_vector_norm_lower_uT"]) or s < float(base["vector_sine_separation_lower"]) or rate > float(base["body_rate_norm_upper_deg_s"]):
        raise RuntimeError("declared PE does not refine vector certificate")
    vc = vector["configured_measurement_bounds"]
    ra = up(pos(vc["acc_measurement_variance_upper"], "Racc"))
    rm = up(pos(vc["mag_measurement_variance_upper"], "Rmag"))
    angular = down(s * s / up(1.0 + math.sqrt(max(0.0, 1.0 - s * s))))
    mu = down(min(f * f / ra, m * m / rm) * angular)
    dg0, dg1 = map(float, base["packet_gap_s"])
    omega = up(rate * math.pi / 180.0)
    bracket = down(1.0 - up(0.5 * omega * dg1))
    gamma = down(dg0 * bracket / pos(base["gyro_bias_time_scale_s"], "Tbg"))
    return down(mu / up(1.0 + up(2.0 / down(gamma * gamma))))


def _global_translation_upper(dynamic: dict, live: dict, axis_factors: list[float]) -> tuple[list[float], dict]:
    """Finite-memory (v,p,S,a_w) variance ceiling valid for any source motion."""
    inv = dynamic["dynamic_invariant"]
    rates = dynamic["validated_rate_and_jump_bounds"]
    h = pos(rates["dt_s"], "dt")
    tau_lo = pos(inv["tau_applied_s"][0], "tau lower")
    sigma_hi = pos(inv["sigma_aw_filter_mps2"][1], "sigma upper")
    rs_hi = pos(inv["R_S_applied"][1], "R_S upper")
    cadence_hi = pos(inv["pseudo_update_period_s"][1], "pseudo cadence upper")
    Tpe = pos(live["vector_pe_recurrence_window_s"], "PE recurrence")

    # Progress-preserving scheduler plus one configured sample is a uniform
    # firing-gap upper.  It remains valid while T_S is retargeted.
    gap = up(cadence_hi + h)
    spacing = up(max(Tpe, 2.0 * gap))
    Tobs = up(2.0 * spacing + gap)
    Tword = up(Tobs + Tpe)

    Binv = integrator_inverse(gap, spacing)
    qc_hi = up(2.0 * sigma_hi * sigma_hi / tau_lo)
    # These use only the global invariant |a_w| covariance scale and the
    # largest driving intensity, hence remain valid under time-varying tau/sigma.
    s_nuis = up(sigma_hi * sigma_hi * (Tobs ** 3 / 6.0) ** 2)
    s_proc = up(qc_hi * Tobs ** 7 / 252.0)
    rmax = up((rs_hi * max(axis_factors)) ** 2)
    rstack = up(3.0 * (rmax + s_nuis + s_proc))
    R = [[I(rstack if i == j else 0.0) for j in range(3)] for i in range(3)]
    Cspv = matrix_symmetric_hull(matrix_mul(matrix_mul(Binv, R), matrix_transpose(Binv)))
    order = (2, 1, 0)  # [S,p,v] -> [v,p,S]
    Cvps = [[Cspv[order[i]][order[j]] for j in range(3)] for i in range(3)]
    t = Interval.outward_bounds(0.0, Tword)
    F = [[I(1),I(0),I(0)],[t,I(1),I(0)],[I(0.5)*t.square(),t,I(1)]]
    Cend = matrix_symmetric_hull(matrix_mul(matrix_mul(F, Cvps), matrix_transpose(F)))
    u = diagonal_dominator(Cend)

    variances = [
        up(sigma_hi * sigma_hi * Tword * Tword + qc_hi * Tword ** 3 / 3.0),
        up(sigma_hi * sigma_hi * Tword ** 4 / 4.0 + qc_hi * Tword ** 5 / 20.0),
        up(sigma_hi * sigma_hi * Tword ** 6 / 36.0 + qc_hi * Tword ** 7 / 252.0),
        up(sigma_hi * sigma_hi),
    ]
    roots = [math.sqrt(v) for v in variances]
    total = up(sum(roots))
    noise = [up(r * total) for r in roots]
    upper = [up(u[i] + noise[i]) for i in range(3)] + [noise[3]]
    return upper, {
        "pseudo_gap_s_upper": gap,
        "observation_window_s_upper": Tobs,
        "covariance_memory_window_s_upper": Tword,
        "q_c_global_upper": qc_hi,
        "source_motion_inside_window_allowed": True,
    }


def _global_full_state_upper(dynamic: dict, live: dict, vector: dict, process: dict, alpha6: float, trans_upper: list[float], timing: dict) -> dict:
    inv = dynamic["dynamic_invariant"]
    sigma_hi = pos(inv["sigma_aw_filter_mps2"][1], "sigma upper")
    tau_lo = pos(inv["tau_applied_s"][0], "tau lower")
    vc = vector["configured_measurement_bounds"]
    ra = down(pos(vc["acc_measurement_std_mps2"], "acc std") ** 2)
    rm = down(pos(vc["mag_measurement_std_uT"], "mag std") ** 2)
    pc = process["source_constants"]
    qg = down(pos(pc["gyro_noise_density_rad_sqrt_s"], "gyro noise") ** 2)
    qb = down(pos(pc["gyro_bias_rw_variance_density"], "gyro bias process"))
    qba = down(pos(pc["accel_bias_process_variance_density"], "acc bias process"))
    tau_ba = pos(pc["accel_bias_tau_s"], "acc bias tau")
    pba = up(max(0.004 ** 2, qba * tau_ba / 2.0))
    fhi = pos(live["specific_force_norm_upper_mps2"], "force upper")
    mhi = pos(live["magnetic_vector_norm_upper_uT"], "mag upper")
    qc_hi = up(2.0 * sigma_hi * sigma_hi / tau_lo)
    pair = pos(vector["operating_envelope"]["packet_gap_s"][1], "vector packet gap")

    qab = up(3.0 * (qg * pair + qb * (pair + pair ** 3 / 3.0)))
    whitened = up(
        (fhi * fhi / ra + mhi * mhi / rm) * qab
        + (3.0 * qc_hi * pair) / ra
        + (3.0 * qba * pair) / ra
        + (6.0 * pba) / ra
    )
    u0 = up((1.0 + whitened) / alpha6)
    T = pos(timing["covariance_memory_window_s_upper"], "covariance window")
    uab_prop = up(6.0 * (qg * T + qb * (T + T ** 3 / 3.0)))
    utheta = up(2.0 * (1.0 + T * T) * u0 + uab_prop)
    ubg = up(2.0 * u0 + uab_prop)

    H = [utheta] * 3 + [ubg] * 3
    H += [trans_upper[0]] * 3 + [trans_upper[1]] * 3
    H += [trans_upper[2]] * 3 + [trans_upper[3]] * 3
    A = H + [pba] * 3
    return {"H": H, "A": A, "accel_bias_variance_upper": pba,
            "attitude_variance_upper": utheta, "gyro_bias_variance_upper": ubg}


def _local_cell(mode: str, x: Interval, rho_trans: float, sigma: Interval, rs: Interval,
                Pdiag: list[float], live: dict, vector: dict, process: dict,
                h: float, axis_factors: list[float]) -> dict:
    qtheta = pos(process["attitude_gyro_bias"]["theta_diagonal_lower"], "theta process")
    qbg = pos(process["attitude_gyro_bias"]["gyro_bias_diagonal_lower"], "bias process")
    cross = float(process["attitude_gyro_bias"]["cross_norm_upper"])
    rho_att = down(1.0 - cross / math.sqrt(qtheta * qbg))
    if rho_att <= 0.0:
        raise RuntimeError("scaled attitude/bias process comparison lost positivity")
    qba_d = pos(process["active_accelerometer_bias"]["Q_accel_bias_lambda_min_lower"], "active bias process")

    # Similarity scale is tied to the current source cell only.  Actual sigma
    # may be above sigma.lo; using the lower endpoint can only strengthen the
    # normalized process covariance lower comparison.
    sv2 = (sigma.lo * h) ** 2
    sp2 = (sigma.lo * h * h) ** 2
    sS2 = (sigma.lo * h * h * h) ** 2
    sa2 = sigma.lo * sigma.lo
    scales2 = [qtheta] * 3 + [qbg] * 3 + [sv2] * 3 + [sp2] * 3 + [sS2] * 3 + [sa2] * 3
    rho = min(rho_trans, rho_att)
    if mode == "A":
        scales2 += [qba_d] * 3
        rho = min(rho, 1.0)

    vc = vector["configured_measurement_bounds"]
    ra = down(pos(vc["acc_measurement_std_mps2"], "acc std") ** 2)
    rm = down(pos(vc["mag_measurement_std_uT"], "mag std") ** 2)
    fhi = pos(live["specific_force_norm_upper_mps2"], "force upper")
    mhi = pos(live["magnetic_vector_norm_upper_uT"], "mag upper")

    # Assimilating every possible same-sample measurement gives the *smallest*
    # optimal posterior of the process injection, hence is a valid lower bound
    # for any actual subset.  Horizontal R_S factors are std multipliers.
    rs_var_lower = down((rs.lo * min(axis_factors)) ** 2)
    betaS = up(sS2 / rs_var_lower)
    betaAcc = up((fhi * fhi * qtheta + sa2 + (qba_d if mode == "A" else 0.0)) / ra)
    betaMag = up((mhi * mhi * qtheta) / rm)
    beta = up(betaS + betaAcc + betaMag)
    rho_post = down(1.0 / up(1.0 / rho + beta))

    # PSD trace domination is intentionally used instead of the old max-diagonal
    # shortcut: lambda_max(D^-1 P D^-T) <= trace(...) is rigorous even with
    # arbitrary cross-covariances.
    scaled_trace_upper = up(sum(up(Pdiag[i] / scales2[i]) for i in range(len(Pdiag))))
    delta = down(rho_post / scaled_trace_upper)
    return {
        "mode": mode,
        "x_h_over_tau": x.as_list(),
        "sigma_aw_mps2": sigma.as_list(),
        "R_S_applied": rs.as_list(),
        "process_scaled_lambda_min_lower": rho,
        "post_measurement_scaled_Omega_lambda_min_lower": rho_post,
        "Pbar_scaled_trace_upper": scaled_trace_upper,
        "relative_Riccati_injection_margin_lower": delta,
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    live = domain["normal_live"]
    dynamic = DYNAMIC.build(path)
    df = DYNAMIC.validate(dynamic)
    vector = VECTOR.build()
    vf = VECTOR.validate(vector)
    process = PROCESS.build()
    pf = PROCESS.validate(process)
    failures = [f"dynamic: {x}" for x in df] + [f"vector: {x}" for x in vf] + [f"process: {x}" for x in pf]
    if failures:
        raise RuntimeError(f"SEA3 Riccati-tube prerequisites failed: {failures}")

    axis_factors = _axis_factors()
    alpha6 = pos(_declared_vector_alpha6(live, vector), "declared alpha6")
    trans_upper, timing = _global_translation_upper(dynamic, live, axis_factors)
    full_upper = _global_full_state_upper(dynamic, live, vector, process, alpha6, trans_upper, timing)

    inv = dynamic["dynamic_invariant"]
    h = pos(dynamic["validated_rate_and_jump_bounds"]["dt_s"], "dt")
    tau_lo, tau_hi = map(float, inv["tau_applied_s"])
    xlo, xhi = h / tau_hi, h / tau_lo
    edges = geom_edges(xlo, xhi, 24)
    if xlo < BRANCH_X < xhi:
        edges = sorted(set(edges + [BRANCH_X]))
    xcells: list[tuple[Interval, float]] = []
    for cell in interval_cells(edges):
        xcells.extend(split_x_cell(cell))

    sigma_lo, sigma_hi = map(float, inv["sigma_aw_filter_mps2"])
    rs_lo, rs_hi = map(float, inv["R_S_applied"])
    sigmas = interval_cells(geom_edges(sigma_lo, sigma_hi, 5))
    rss = interval_cells(geom_edges(rs_lo, rs_hi, 8))

    worst = {"H": None, "A": None}
    count = 0
    for x, rho_t in xcells:
        for sigma in sigmas:
            for rs in rss:
                count += 1
                for mode in ("H", "A"):
                    row = _local_cell(mode, x, rho_t, sigma, rs, full_upper[mode], live, vector, process, h, axis_factors)
                    if worst[mode] is None or row["relative_Riccati_injection_margin_lower"] < worst[mode]["relative_Riccati_injection_margin_lower"]:
                        worst[mode] = row

    modes = {}
    for mode in ("H", "A"):
        w = worst[mode]
        delta = float(w["relative_Riccati_injection_margin_lower"])
        modes[mode] = {
            "dimension": 18 if mode == "H" else 21,
            "Pbar_diagonal_variance_upper": full_upper[mode],
            "Pbar_lambda_max_trace_upper": up(sum(full_upper[mode])),
            "relative_Riccati_injection_margin_lower": delta,
            "contraction_factor_upper": up(1.0 - delta),
            "useful_margin_gate": USEFUL_GATE,
            "useful_margin_pass": delta >= USEFUL_GATE,
            "worst_current_source_cell": w,
        }

    passed = all(modes[m]["useful_margin_pass"] for m in ("H", "A"))
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_role": "QUANTITATIVE_MOVING_RICCATI_P3_TUBE",
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_shrunk": False,
        "SEA3_dynamic_source_consumed": True,
        "source_history_graph_consumed": False,
        "predecessor_path_enumeration_consumed": False,
        "P2_800_state_partition_consumed": False,
        "current_source_interval_cover_only": True,
        "time_varying_source_allowed_inside_covariance_memory_window": True,
        "metric_identity": "V=e^T P^-1 e with P the shipping Riccati covariance",
        "covariance_ceiling_argument": "finite-memory recurrent vector/S estimator, global SEA3 adaptive invariant for every nuisance/process upper",
        "process_floor_argument": "exact current-sample integrated-OU covariance plus Joseph optimal-posterior lower comparison",
        "PSD_cross_covariance_handling": "lambda_max(D^-1 Pbar D^-T) <= trace(D^-1 Pbar D^-T); no max-diagonal shortcut",
        "declared_vector_alpha6_information_lower": alpha6,
        "covariance_memory": timing,
        "cell_cover": {
            "x_cells": len(xcells),
            "sigma_cells": len(sigmas),
            "R_S_cells": len(rss),
            "joint_current_source_cells": count,
            "history_depth": 0,
        },
        "modes": modes,
        "useful_gate": USEFUL_GATE,
        "RICCATI_TUBE_PASS": passed,
        "P3_MAY_PROMOTE_FROM_THIS_TUBE": passed,
        "next_obligation": (
            "if both H/A margins pass, bind this tube into the canonical moving-Riccati P3 gate; "
            "otherwise tighten SEA3 estimator-state motion/finite-memory bounds, never reintroduce source histories"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    for key in (
        "source_generated_not_trajectory_fit", "SEA3_dynamic_source_consumed",
        "current_source_interval_cover_only", "time_varying_source_allowed_inside_covariance_memory_window",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_shrunk",
        "source_history_graph_consumed", "predecessor_path_enumeration_consumed",
        "P2_800_state_partition_consumed",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    cover = d.get("cell_cover", {})
    if int(cover.get("joint_current_source_cells", 0)) <= 0 or int(cover.get("history_depth", -1)) != 0:
        f.append("invalid current-source interval cover")
    if float(d.get("useful_gate", math.nan)) != USEFUL_GATE:
        f.append("useful gate changed")
    for mode in ("H", "A"):
        row = d.get("modes", {}).get(mode, {})
        delta = row.get("relative_Riccati_injection_margin_lower")
        if not isinstance(delta, (int, float)) or not math.isfinite(float(delta)) or float(delta) <= 0.0:
            f.append(f"{mode} Riccati margin is not finite positive")
        diag = row.get("Pbar_diagonal_variance_upper")
        if not isinstance(diag, list) or len(diag) != row.get("dimension") or any((not math.isfinite(float(x)) or float(x) <= 0.0) for x in diag):
            f.append(f"{mode} covariance ceiling is invalid")
        w = row.get("worst_current_source_cell", {})
        if not w or float(w.get("post_measurement_scaled_Omega_lambda_min_lower", 0.0)) <= 0.0:
            f.append(f"{mode} current-source process/Joseph floor missing")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "tube_pass": d["RICCATI_TUBE_PASS"],
        "cells": d["cell_cover"],
        "H_delta": d["modes"]["H"]["relative_Riccati_injection_margin_lower"],
        "A_delta": d["modes"]["A"]["relative_Riccati_injection_margin_lower"],
        "H_worst": d["modes"]["H"]["worst_current_source_cell"],
        "A_worst": d["modes"]["A"]["worst_current_source_cell"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
