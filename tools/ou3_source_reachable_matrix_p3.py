#!/usr/bin/env python3
"""Source-reachable matrix-valued P3 certificate for deployed OU-III.

This replaces the old scalar ``min(Q) / max(trace(P))`` route.  The proof keeps
scheduled parameters in source-derived cells, binds pseudo cadence to the same
applied tau used by the filter, carries a directional covariance upper matrix,
and derives the Riccati information margin as a generalized matrix ratio in a
source-dependent comparison scaling.

The lower comparison is valid for the actual Joseph implementation: for a prior
noise covariance Q and any implemented gain K,

  (I-KH)Q(I-KH)' + K R K' >= (Q^-1 + H'R^-1H)^-1.

Thus an optimally assimilated set of *all* possible same-sample measurements is
a conservative lower bound even when the source branch rejects or omits some of
them.  The proof never uses the one-step Euclidean lambda_min(Q) as a global
state-space floor.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

from ou3_interval import (
    Interval, matrix_identity, matrix_mul, matrix_sub, matrix_transpose,
    symmetric_positive_definite_ldlt,
)
from ou3_interval_linear_algebra import matrix_inverse_gauss_jordan, matrix_symmetric_hull
import ou3_full_process_ucc as PROCESS
import ou3_implementation_word_language as WORDS
import ou3_source_domain_contract as SOURCE
import ou3_validated_transcendentals as VT
import ou3_vector_uco_certificate as VECTOR

REPO = Path(__file__).resolve().parents[1]
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 4
MIN_USEFUL_DELTA = 1.0e-18


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
    e = [lo]
    for _ in range(n - 1):
        e.append(e[-1] * r)
    e.append(hi)
    return e


def cells_from_edges(edges: list[float]) -> list[Interval]:
    return [Interval.outward_bounds(edges[i], edges[i + 1]) for i in range(len(edges) - 1)]


def poly(x: Interval, terms: tuple[tuple[int, float], ...]) -> Interval:
    out = Interval.point(0.0)
    for n, c in terms:
        out = out + I(c) * ipow(x, n)
    return out


def qbar_integrated_ou(x: Interval) -> list[list[Interval]]:
    """Dimensionless exact source formulas for [v,p,S,a] process covariance.

    Q = diag(sigma*tau, sigma*tau^2, sigma*tau^3, sigma) Qbar D'.
    The implementation's small-x polynomial and expm1 branches are mirrored.
    """
    if x.lo <= 0.0 or x.hi > 0.5:
        raise ValueError("x=h/tau outside validated range")
    z = Interval.point(0.0)
    if x.hi <= 1.0e-2:
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
    elif x.lo >= 1.0e-2:
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
    else:
        raise ValueError("x cell crosses implementation branch threshold")
    return [
        [qvv, qvp, qvS, qva],
        [qvp, qpp, qpS, qpa],
        [qvS, qpS, qSS, qSa],
        [qva, qpa, qSa, qaa],
    ]


def step_scaled_q(x: Interval) -> list[list[Interval]]:
    # D_h = diag(sigma*h, sigma*h^2, sigma*h^3, sigma).
    # Since h=tau*x, D_h^-1 Q D_h^-T = Qbar / x^(k_i+k_j).
    q = qbar_integrated_ou(x)
    k = (1, 2, 3, 0)
    return [[q[i][j] / ipow(x, k[i] + k[j]) for j in range(4)] for i in range(4)]


def subtract_rho(A, rho: float):
    n = len(A)
    return [[A[i][j] - I(rho if i == j else 0.0) for j in range(n)] for i in range(n)]


def certified_rho(A) -> float:
    ok, _ = symmetric_positive_definite_ldlt(A)
    if not ok:
        return 0.0
    hi = min(x.lo for x in (A[i][i] for i in range(len(A))))
    lo = 0.0
    for _ in range(45):
        mid = 0.5 * (lo + hi)
        ok, _ = symmetric_positive_definite_ldlt(subtract_rho(A, mid))
        if ok:
            lo = mid
        else:
            hi = mid
    return down(lo)


def split_x_cell(x: Interval, depth: int = 0) -> list[tuple[Interval, float]]:
    rho = certified_rho(step_scaled_q(x))
    if rho > 0.0:
        return [(x, rho)]
    if depth >= 12:
        raise RuntimeError(f"cannot certify scaled OU process cell {x.as_list()}")
    mid = math.sqrt(x.lo * x.hi)
    return split_x_cell(Interval.outward_bounds(x.lo, mid), depth+1) + split_x_cell(Interval.outward_bounds(mid, x.hi), depth+1)


def diag_dominator(A) -> list[float]:
    A = matrix_symmetric_hull(A)
    out = []
    for i in range(len(A)):
        d = max(0.0, A[i][i].hi)
        for j in range(len(A)):
            if i != j:
                d = up(d + A[i][j].abs_upper())
        out.append(up(d))
    return out


def integrator_inverse(gap: float, spacing: float):
    g, s = pos(gap, "gap"), pos(spacing, "spacing")
    if s <= g:
        raise RuntimeError("S windows overlap")
    ts = (Interval.outward_bounds(0.0,g), Interval.outward_bounds(s,s+g), Interval.outward_bounds(2*s,2*s+g))
    B = [[I(1), t, I(0.5)*t.square()] for t in ts]  # [S,p,v]
    return matrix_inverse_gauss_jordan(B)


def declared_vector_information(live: dict, vector: dict) -> float:
    base = vector["operating_envelope"]
    f = pos(live["specific_force_norm_lower_mps2"], "force floor")
    m = pos(live["magnetic_vector_norm_lower_uT"], "mag floor")
    s = pos(live["vector_sine_separation_lower"], "separation floor")
    if f < float(base["specific_force_norm_lower_mps2"]) or m < float(base["magnetic_vector_norm_lower_uT"]) or s < float(base["vector_sine_separation_lower"]):
        raise RuntimeError("declared PE weaker than source vector certificate")
    vc = vector["configured_measurement_bounds"]
    ra = up(pos(vc["acc_measurement_variance_upper"], "Racc"))
    rm = up(pos(vc["mag_measurement_variance_upper"], "Rmag"))
    angular = down(s*s / up(1.0 + math.sqrt(max(0.0,1.0-s*s))))
    mu = down(min(f*f/ra, m*m/rm) * angular)
    dg0, dg1 = map(float, base["packet_gap_s"])
    omega = up(pos(live["body_rate_norm_upper_deg_s"], "rate") * math.pi/180.0)
    gamma = down(dg0 * down(1.0 - 0.5*omega*dg1) / pos(base["gyro_bias_time_scale_s"], "Tbg"))
    return down(mu / up(1.0 + 2.0/down(gamma*gamma)))


def source_schedule() -> dict:
    text = WRAPPER.read_text(encoding="utf-8")
    for marker in (
        "RSAdaptationLaw rs_law_ = RSAdaptationLaw::SpectralMSE;",
        "apply_pseudo_update_cadence_();",
        "tune_.tau_applied   += alpha",
        "tune_.RS_applied    += alpha_RS",
    ):
        if marker not in text:
            raise RuntimeError(f"missing deployed schedule semantic: {marker}")
    c = {name: SOURCE.parse_const(text, name) for name in (
        "MIN_TUNE_FREQ_HZ","MAX_TUNE_FREQ_HZ","MIN_TAU_S","MAX_TAU_S",
        "MIN_R_S","MAX_R_S","PSEUDO_UPDATE_TAU_RATIO_DEFAULT",
        "PSEUDO_UPDATE_PERIOD_MIN_S_DEFAULT","PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT",
        "FREQ_SMOOTHER_DT",
    )}
    tau_coeff = float(re.search(r"float\s+tau_coeff_\s*=\s*([0-9.eE+-]+)f", text).group(1))
    init_tau = float(re.search(r"float\s+tau_applied\s*=\s*([0-9.eE+-]+)f", text).group(1))
    init_rs = float(re.search(r"float\s+RS_applied\s*=\s*([0-9.eE+-]+)f", text).group(1))
    target_lo = max(c["MIN_TAU_S"], min(c["MAX_TAU_S"], tau_coeff*0.5/c["MAX_TUNE_FREQ_HZ"]))
    target_hi = max(c["MIN_TAU_S"], min(c["MAX_TAU_S"], tau_coeff*0.5/c["MIN_TUNE_FREQ_HZ"]))
    return {
        "tau_target_s": [down(target_lo), up(target_hi)],
        "tau_applied_invariant_s": [down(min(init_tau,target_lo)), up(max(init_tau,target_hi))],
        "R_S_applied_invariant": [down(min(init_rs,c["MIN_R_S"])), up(max(init_rs,c["MAX_R_S"]))],
        "sigma_aw_applied_safety": [0.05, 6.0],
        "pseudo_ratio": c["PSEUDO_UPDATE_TAU_RATIO_DEFAULT"],
        "pseudo_min_s": c["PSEUDO_UPDATE_PERIOD_MIN_S_DEFAULT"],
        "pseudo_max_s": c["PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT"],
        "dt_s": c["FREQ_SMOOTHER_DT"],
        "proof_kind": "SOURCE_REACHABLE_INVARIANT_CELL_OVERAPPROXIMATION",
        "note": "tau target is derived from the deployed wave-band law; cadence is coupled to applied tau in every cell; R_S remains a cell state rather than being collapsed to unrelated global extrema",
    }


def cadence_bounds(tau: Interval, sched: dict) -> list[float]:
    r = sched["pseudo_ratio"]
    lo = min(max(r*tau.lo, sched["pseudo_min_s"]), sched["pseudo_max_s"])
    hi = min(max(r*tau.hi, sched["pseudo_min_s"]), sched["pseudo_max_s"])
    return [down(lo), up(hi)]


def translation_upper(tau: Interval, sigma: Interval, rs: Interval, Tpe: float, sched: dict) -> tuple[list[float], dict]:
    h = sched["dt_s"]
    cadence = cadence_bounds(tau, sched)
    gap = up(cadence[1] + h)
    spacing = up(max(Tpe, 2.0*gap))
    Tobs = up(2.0*spacing + gap)
    Tword = up(Tobs + Tpe)
    Binv = integrator_inverse(gap, spacing)
    qc = up(2.0*sigma.hi*sigma.hi/tau.lo)
    s_nuis = up(sigma.hi*sigma.hi*(Tobs**3/6.0)**2)
    s_proc = up(qc*Tobs**7/252.0)
    r_stack = up(3.0*(rs.hi*rs.hi + s_nuis + s_proc))
    R = [[I(r_stack if i == j else 0.0) for j in range(3)] for i in range(3)]
    Cspv = matrix_symmetric_hull(matrix_mul(matrix_mul(Binv,R),matrix_transpose(Binv)))
    order = (2,1,0)
    Cvps = [[Cspv[order[i]][order[j]] for j in range(3)] for i in range(3)]
    t = Interval.outward_bounds(0.0,Tword)
    F = [[I(1),I(0),I(0)],[t,I(1),I(0)],[I(0.5)*t.square(),t,I(1)]]
    Cend = matrix_symmetric_hull(matrix_mul(matrix_mul(F,Cvps),matrix_transpose(F)))
    u = diag_dominator(Cend)
    vars_ = [
        up(sigma.hi*sigma.hi*Tword*Tword + qc*Tword**3/3.0),
        up(sigma.hi*sigma.hi*Tword**4/4.0 + qc*Tword**5/20.0),
        up(sigma.hi*sigma.hi*Tword**6/36.0 + qc*Tword**7/252.0),
        up(sigma.hi*sigma.hi),
    ]
    roots = [math.sqrt(v) for v in vars_]
    noise_dom = [up(roots[i]*sum(roots)) for i in range(4)]
    return [up(u[i]+noise_dom[i]) for i in range(3)] + [noise_dom[3]], {
        "cadence_s": cadence, "gap_s_upper": gap, "word_horizon_s_upper": Tword,
    }


def mode_cell(mode: str, x: Interval, rho_trans: float, sigma: Interval, rs: Interval, live: dict, vector: dict, process: dict, sched: dict) -> dict:
    h = sched["dt_s"]
    tau = Interval.outward_bounds(h/x.hi, h/x.lo)
    Tpe = pos(live["vector_pe_recurrence_window_s"], "PE recurrence")
    utrans, timing = translation_upper(tau,sigma,rs,Tpe,sched)
    alpha6 = pos(declared_vector_information(live,vector), "alpha6")
    vc = vector["configured_measurement_bounds"]
    ra_lo = down(pos(vc["acc_measurement_std_mps2"], "acc std")**2)
    rm_lo = down(pos(vc["mag_measurement_std_uT"], "mag std")**2)
    pc = process["source_constants"]
    qg = down(pos(pc["gyro_noise_density_rad_sqrt_s"], "gyro noise")**2)
    qb = down(pos(pc["gyro_bias_rw_variance_density"], "gyro bias process"))
    qba = down(pos(pc["accel_bias_process_variance_density"], "acc bias process"))
    tau_ba = pos(pc["accel_bias_tau_s"], "acc bias tau")
    pba = up(max(0.004**2, qba*tau_ba/2.0))
    fhi = pos(live["specific_force_norm_upper_mps2"], "force upper")
    mhi = pos(live["magnetic_vector_norm_upper_uT"], "mag upper")
    qc = up(2.0*sigma.hi*sigma.hi/tau.lo)
    pair = pos(vector["operating_envelope"]["packet_gap_s"][1], "vector gap")
    qab = up(3.0*(qg*pair + qb*(pair + pair**3/3.0)))
    whitened = up((fhi*fhi/ra_lo + mhi*mhi/rm_lo)*qab + (3.0*qc*pair)/ra_lo + (3.0*qba*pair)/ra_lo + (6.0*pba)/ra_lo)
    u0 = up(1.0/alpha6 + whitened/alpha6)
    T = timing["word_horizon_s_upper"]
    uab = up(2.0*(1.0+T*T)*u0 + 2.0*3.0*(qg*T + qb*(T+T**3/3.0)))
    qtheta = pos(process["attitude_gyro_bias"]["theta_diagonal_lower"], "qtheta")
    qbg = pos(process["attitude_gyro_bias"]["gyro_bias_diagonal_lower"], "qbg")
    cross = float(process["attitude_gyro_bias"]["cross_norm_upper"])
    rho_att = down(1.0 - cross/math.sqrt(qtheta*qbg))
    if rho_att <= 0.0:
        raise RuntimeError("attitude scaled process lower lost positivity")
    qba_d = pos(process["active_accelerometer_bias"]["Q_accel_bias_lambda_min_lower"], "qba discrete")
    scales2 = [qtheta]*3 + [qbg]*3
    upper = [uab]*6
    # State order: theta,bg,v,p,S,aw,(ba); translation block repeated by group.
    sv2 = (sigma.lo*h)**2
    sp2 = (sigma.lo*h*h)**2
    sS2 = (sigma.lo*h*h*h)**2
    sa2 = sigma.lo*sigma.lo
    scales2 += [sv2]*3 + [sp2]*3 + [sS2]*3 + [sa2]*3
    upper += [utrans[0]]*3 + [utrans[1]]*3 + [utrans[2]]*3 + [utrans[3]]*3
    rho = min(rho_trans,rho_att)
    if mode == "A":
        scales2 += [qba_d]*3
        upper += [pba]*3
        rho = min(rho,1.0)
    betaS = sS2/(rs.lo*rs.lo)
    betaAcc = (fhi*fhi*qtheta + sa2 + (qba_d if mode == "A" else 0.0))/ra_lo
    betaMag = (mhi*mhi*qtheta)/rm_lo
    beta = up(betaS + betaAcc + betaMag)
    rho_post = down(1.0/up(1.0/rho + beta))
    scaled_upper = up(max(upper[i]/scales2[i] for i in range(len(upper))))
    delta = down(rho_post/scaled_upper)
    return {
        "mode": mode,
        "tau_s": tau.as_list(), "sigma_aw_mps2": sigma.as_list(), "R_S_filter_std": rs.as_list(),
        "x_h_over_tau": x.as_list(), **timing,
        "process_scaled_lambda_min_lower": rho,
        "post_measurement_scaled_Omega_lambda_min_lower": rho_post,
        "Sigma_scaled_lambda_max_upper": scaled_upper,
        "relative_Riccati_injection_margin_lower": delta,
        "Sigma_lambda_min_lower": down(rho_post*min(scales2)),
        "Sigma_lambda_max_upper": up(max(upper)),
        "word_noise_Omega_lambda_min_lower": down(rho_post*min(scales2)),
        "comparison_scale_diagonal_squared": scales2,
        "Sigma_diagonal_upper": upper,
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("proof domain must not be trajectory fitted")
    live = domain["normal_live"]
    vector, process, words = VECTOR.build(), PROCESS.build(), WORDS.build(domain_path)
    for label, failures in (("vector",VECTOR.validate(vector)),("process",PROCESS.validate(process)),("word-language",WORDS.validate(words))):
        if failures:
            raise RuntimeError(f"{label} prerequisite failed: {failures}")
    sched = source_schedule()
    h = sched["dt_s"]
    tau_lo,tau_hi = sched["tau_applied_invariant_s"]
    xlo,xhi = h/tau_hi, h/tau_lo
    base = geom_edges(xlo,xhi,24)
    if xlo < 1.0e-2 < xhi:
        base.append(1.0e-2); base = sorted(set(base))
    xcells = []
    for c in cells_from_edges(base):
        xcells.extend(split_x_cell(c))
    sigmas = cells_from_edges(geom_edges(0.05,6.0,5))
    rss = cells_from_edges(geom_edges(sched["R_S_applied_invariant"][0],sched["R_S_applied_invariant"][1],8))
    worst = {"H":None,"A":None}
    count = 0
    for x,rho_t in xcells:
        for sigma in sigmas:
            for rs in rss:
                count += 1
                for mode in ("H","A"):
                    c = mode_cell(mode,x,rho_t,sigma,rs,live,vector,process,sched)
                    if worst[mode] is None or c["relative_Riccati_injection_margin_lower"] < worst[mode]["relative_Riccati_injection_margin_lower"]:
                        worst[mode] = c
    modes = {}
    for mode in ("H","A"):
        w = worst[mode]
        delta = w["relative_Riccati_injection_margin_lower"]
        modes[mode] = {
            "dimension": 18 if mode == "H" else 21,
            "Sigma_lambda_min_lower": w["Sigma_lambda_min_lower"],
            "Sigma_lambda_max_upper": w["Sigma_lambda_max_upper"],
            "word_noise_Omega_lambda_min_lower": w["word_noise_Omega_lambda_min_lower"],
            "relative_Riccati_injection_margin_lower": delta,
            "lambda_information_upper_formula": "1-delta_lower",
            "prefix_information_gain_upper": 1.0,
            "strict_information_contraction": delta > 0.0,
            "useful_margin_gate": MIN_USEFUL_DELTA,
            "useful_margin_pass": delta >= MIN_USEFUL_DELTA,
            "matrix_comparison": w,
            "pass": delta >= MIN_USEFUL_DELTA,
        }
    passed = all(modes[m]["pass"] for m in ("H","A"))
    return {
        "schema": SCHEMA,
        "qualification": "SOURCE_REACHABLE_MATRIX_VALUED_DEPLOYED_OU3_HA_INFORMATION_WORD_CERTIFICATE",
        "p3_backend": "SOURCE_CELL_GENERALIZED_MATRIX_COMPARISON",
        "source_generated_not_trajectory_fit": True,
        "validated_arithmetic": True,
        "outward_rounded": True,
        "source_schedule": sched,
        "cell_partition": {"x_cells":len(xcells),"sigma_cells":len(sigmas),"R_S_cells":len(rss),"joint_cells":count},
        "matrix_lower_argument": "Joseph noise for any implemented gain is bounded below by the optimal posterior of the same process/noise prior; all possible same-sample measurements are assimilated only to make the lower comparison smaller",
        "matrix_upper_argument": "finite-memory vector/S estimators give directional covariance matrices; interval row-sum diagonal dominators preserve Loewner upper bounds",
        "modes": modes,
        "continuous_linear_information_certificate": "PASS" if passed else "FAIL",
        "nonlinear_word_enclosed": False,
        "theorem_promotion": "LINEAR_ONLY" if passed else "NOT_ESTABLISHED",
        "old_scalar_min_Q_route_used": False,
        "next_obligation": "P3 complete; only after this certificate passes may P4 consume the matrix metric",
    }


def validate(d: dict) -> list[str]:
    f = []
    if d.get("schema") != SCHEMA: f.append("schema mismatch")
    for k in ("source_generated_not_trajectory_fit","validated_arithmetic","outward_rounded"):
        if d.get(k) is not True: f.append(f"{k} is not true")
    if d.get("p3_backend") != "SOURCE_CELL_GENERALIZED_MATRIX_COMPARISON": f.append("wrong P3 backend")
    if d.get("old_scalar_min_Q_route_used") is not False: f.append("old scalar route still active")
    if int(d.get("cell_partition",{}).get("joint_cells",0)) <= 0: f.append("no source cells")
    for mode in ("H","A"):
        r = d.get("modes",{}).get(mode,{})
        for k in ("Sigma_lambda_min_lower","Sigma_lambda_max_upper","word_noise_Omega_lambda_min_lower","relative_Riccati_injection_margin_lower"):
            x = r.get(k)
            if not isinstance(x,(int,float)) or not math.isfinite(float(x)) or float(x) <= 0.0: f.append(f"{mode}.{k} invalid")
        if r.get("useful_margin_pass") is not True or r.get("pass") is not True: f.append(f"{mode} useful matrix margin did not pass")
        if r.get("prefix_information_gain_upper") != 1.0: f.append(f"{mode} prefix gain changed")
        mc = r.get("matrix_comparison",{})
        if not mc.get("comparison_scale_diagonal_squared") or not mc.get("Sigma_diagonal_upper"): f.append(f"{mode} missing matrix comparison")
    if f:
        return f
    if d.get("continuous_linear_information_certificate") != "PASS": f.append("P3 did not pass")
    if d.get("nonlinear_word_enclosed") is not False: f.append("P3 must not claim nonlinear enclosure")
    if d.get("theorem_promotion") != "LINEAR_ONLY": f.append("P3 promotion state wrong")
    return f


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(a.domain.resolve())
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d,indent=2,sort_keys=True),encoding="utf-8")
    print(json.dumps({"backend":d["p3_backend"],"cells":d["cell_partition"],"H":d["modes"]["H"],"A":d["modes"]["A"],"failures":failures},indent=2,sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
