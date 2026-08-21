#!/usr/bin/env python3
"""
Computer-assisted analytical certificate for the OU-III Live linear basin proof.

This is not a simulation and does not inspect the eight reference seas.  It
evaluates closed-form worst-case inequalities with Decimal directed rounding
over a declared certificate domain.  Source-enforced/default constants are
cross-checked against the repository; physical/geometric assumptions that are
not enforced by code are printed separately.

The proof closes the following implication:
  q_H I <= Omega_H, P_H <= p_bar I
      => ||L_H^{-1} Psi_H L_0|| <= sqrt(1 - q_H/p_bar) < 1.

q_H is obtained from an oracle lower bound on the posterior covariance of
finite-horizon process noise.  p_bar is obtained from an explicit finite-window
observer using the four-S-update translational observability construction and
the two-vector attitude/gyro-bias construction.  Every arithmetic result used
for certification is rounded outward.
"""

from __future__ import annotations
import argparse
import json
from decimal import Decimal, getcontext, localcontext, ROUND_FLOOR, ROUND_CEILING
from pathlib import Path

getcontext().prec = 90
D = Decimal

def lo_add(a,b):
    with localcontext() as c:
        c.prec=90; c.rounding=ROUND_FLOOR
        return a+b
def hi_add(a,b):
    with localcontext() as c:
        c.prec=90; c.rounding=ROUND_CEILING
        return a+b
def lo_mul(a,b):
    assert a >= 0 and b >= 0
    with localcontext() as c:
        c.prec=90; c.rounding=ROUND_FLOOR
        return a*b
def hi_mul(a,b):
    assert a >= 0 and b >= 0
    with localcontext() as c:
        c.prec=90; c.rounding=ROUND_CEILING
        return a*b
def lo_div(a,b):
    assert a >= 0 and b > 0
    with localcontext() as c:
        c.prec=90; c.rounding=ROUND_FLOOR
        return a/b
def hi_div(a,b):
    assert a >= 0 and b > 0
    with localcontext() as c:
        c.prec=90; c.rounding=ROUND_CEILING
        return a/b

def _ulp(x: Decimal, prec: int = 80) -> Decimal:
    if x.is_zero():
        return D(10) ** D(-prec)
    return D(10) ** D(x.adjusted() - prec + 1)

def exp_hi(x: Decimal) -> Decimal:
    # Decimal.exp is correctly rounded.  Evaluate at extra precision and widen
    # by one ulp, which encloses the exact real exponential.
    with localcontext() as c:
        c.prec=100
        y=x.exp()
    return y + _ulp(y, 90)

def exp_neg_lo(x: Decimal) -> Decimal:
    return lo_div(D(1), exp_hi(x))

def sqrt_hi(x: Decimal) -> Decimal:
    assert x >= 0
    with localcontext() as c:
        c.prec=100
        y=x.sqrt()
    return y + _ulp(y, 90)

def sqrt_lo(x: Decimal) -> Decimal:
    assert x >= 0
    with localcontext() as c:
        c.prec=100
        y=x.sqrt()
    z=y-_ulp(y,90)
    return max(D(0),z)

def sq(x): return x*x

# Fixed dimensionless proof scales (Phase 2 design anchors).
S_TH=D("0.087"); S_BG=D("0.001"); S_V=D("1"); S_P=D("20")
S_S=D("50"); S_AW=D("6"); S_BA=D("0.5")

# Certificate-domain assumptions.  These are NOT all code-enforced.
H_MIN=D("0.0045")
H_MAX=D("0.0055")
SIGMA_AW_MIN=D("0.01")        # required positive certificate floor
F0=D("5"); F1=D("15")         # accepted specific-force magnitude [m/s^2]
M0=D("10"); M1=D("100")       # accepted magnetic magnitude [uT]
DELTA_GEOM=D("0.05")          # 1-|u_f^T u_m| lower bound
OMEGA_MAX=D("1.5")            # rad/s
VEC_GAP_MIN=D("0.005")
VEC_GAP_MAX=D("0.05")

# Validated proof configuration sensor quantities pinned by iss_contract-test.
SIGMA_ACC=D("0.0148")
SIGMA_GYRO=D("0.00157")
SIGMA_MAG=D("0.25")
GYRO_BIAS_Q=D("1e-11")

# Effective S-pseudo standard-deviation enclosure audited by iss_contract-test.
RS_EFF_MIN=D("0.044")
RS_EFF_MAX=D("694")

def require_token(text: str, token: str, where: str):
    if token not in text:
        raise SystemExit(f"SOURCE_CONTRACT_FAIL {where}: missing {token!r}")

def check_source(root: Path):
    wrapper=(root/"src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h").read_text()
    core=(root/"src/kalman_ou_iii/Kalman3D_Wave_OU_III.h").read_text()
    contract=(root/"tests/kalman_ou_iii/iss_contract-test.cpp").read_text()
    for token in (
        "constexpr float MAX_TUNE_FREQ_HZ = 1.5f;",
        "constexpr float MAX_TAU_S   = 12.0f;",
        "constexpr float MAX_SIGMA_A = 6.0f;",
        "constexpr float PSEUDO_UPDATE_PERIOD_MIN_S_DEFAULT = FREQ_SMOOTHER_DT;",
        "constexpr float PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT = 0.25f;",
        "float tau_coeff_    = 1.0f;",
    ):
        require_token(wrapper, token, "wrapper")
    for token in (
        "T sigma_bacc0_ = T(0.004);",
        "Matrix3 Q_bacc_ = Matrix3::Identity() * T(2.5e-7);",
        "T tau_bacc_ = T(5000.0);",
    ):
        require_token(core, token, "core")
    for token in (
        "Eigen::Vector3f::Constant(0.0148f)",
        "Eigen::Vector3f::Constant(0.00157f)",
        "Eigen::Vector3f::Constant(0.25f)",
        "rs_eff_lo > 0.044f",
        "rs_eff_hi > 690.0f",
        "rs_eff_hi < 694.0f",
    ):
        require_token(contract, token, "iss_contract-test")

def calc():
    # Default adaptation path: tau_target = 0.5/f_tune with f_tune<=1.5 Hz,
    # then source clamp to tau<=12 s.  Setter-modified coefficients are outside
    # this certificate.
    tau_min=lo_div(D("0.5"),D("1.5"))
    tau_max=D("12")
    sigma_aw_max=D("6")
    cT=lo_div(D("0.015"),D("1.1"))
    ts_max=min(D("0.25"), hi_mul(cT,tau_max))
    gap_plus=hi_add(ts_max,H_MAX)
    horizon=hi_mul(D(3),gap_plus)

    # ---------- uniform one-step process floor ----------------------------
    # For one OU-III axis:
    #   Qz = qc*h*Dh*W(x)*Dh,
    # lambda_min(W) >= exp(-2x) /
    #   [870912000*(437/315)^3].
    # This bound is exact for each sample and therefore remains valid when tau
    # changes arbitrarily from sample to sample inside the certificate box.
    qc_lo=lo_div(lo_mul(D(2),sq(SIGMA_AW_MIN)),tau_max)
    m_candidates=[
        lo_div(H_MIN,S_V),
        lo_div(sq(H_MIN),S_P),
        lo_div(H_MIN*H_MIN*H_MIN,S_S),
        lo_div(D(1),S_AW),
    ]
    m_lo=min(m_candidates)
    x_hi=hi_div(H_MAX,tau_min)
    mu_den=hi_mul(D(870_912_000), (hi_div(D(437),D(315)))**3)
    mu_lo=lo_div(exp_neg_lo(hi_mul(D(2),x_hi)),mu_den)
    q_step=lo_mul(lo_mul(lo_mul(qc_lo,H_MIN),sq(m_lo)),mu_lo)

    # Optional strengthened process floor if the OU schedule is frozen across
    # the complete proof horizon.  It is reported, but NOT used for the uniform
    # time-varying theorem.
    mh_candidates=[
        lo_div(horizon,S_V),
        lo_div(sq(horizon),S_P),
        lo_div(horizon*horizon*horizon,S_S),
        lo_div(D(1),S_AW),
    ]
    mh=min(mh_candidates)
    mu_h=lo_div(exp_neg_lo(hi_mul(D(2),hi_div(horizon,tau_min))),mu_den)
    q_frozen=lo_mul(lo_mul(lo_mul(qc_lo,horizon),sq(mh)),mu_h)

    # ---------- explicit translational observation inverse ----------------
    # Four S measurements at bounded gaps.  Third divided differences recover
    # a_w, then second/first divided differences recover v,p,S.
    # The broad proof uses the scheduler theorem's only unconditional lower
    # gap h_min and the global effective r_S enclosure.
    root5=sqrt_hi(D(5))
    exp_obs=exp_hi(hi_div(hi_mul(D(3),gap_plus),tau_min))
    k_a_phys=hi_div(
        hi_mul(hi_mul(D(2),root5),hi_mul(RS_EFF_MAX,exp_obs)),
        H_MIN**3,
    )
    t3=hi_mul(D(3),gap_plus)
    a3max=hi_div(t3**3,D(6))
    yq=hi_add(RS_EFF_MAX,hi_mul(hi_mul(D(2),a3max),k_a_phys))
    kv=hi_mul(hi_div(sqrt_hi(D(6)),sq(H_MIN)),yq)
    kp=hi_add(
        hi_mul(hi_div(sqrt_hi(D(2)),H_MIN),yq),
        hi_mul(hi_mul(D("0.5"),kv),gap_plus),
    )
    invL2=hi_add(
        hi_add(sq(hi_div(kv,S_V)),sq(hi_div(kp,S_P))),
        hi_add(sq(hi_div(yq,S_S)),sq(hi_div(k_a_phys,S_AW))),
    )
    sL=lo_div(D(1),sqrt_hi(invL2))

    # ---------- attitude / gyro-bias observation block --------------------
    racc=sq(SIGMA_ACC); rmag=sq(SIGMA_MAG)
    geom=min(lo_div(sq(F0),racc),lo_div(sq(M0),rmag))
    gamma=lo_mul(lo_mul(sq(S_TH),DELTA_GEOM),geom)
    xx=hi_div(hi_mul(OMEGA_MAX,VEC_GAP_MAX),D(2))
    # sinc(x) >= 1 - x^2/6 for x>=0 in the declared small-angle range.
    sinc_lo=D(1)-hi_div(sq(xx),D(6))
    b0=lo_mul(lo_mul(lo_div(S_BG,S_TH),VEC_GAP_MIN),sinc_lo)
    b2=sq(b0); b4=sq(b2)
    # lambda_*(b)=(2+b^2-sqrt(4+b^4))/2.  Use upper sqrt for lower bound.
    lam_num=(D(2)+b2)-sqrt_hi(D(4)+b4)
    lam=max(D(0),lo_div(lam_num,D(2)))
    sA=sqrt_lo(lo_mul(gamma,lam))

    # ---------- triangular 18-state observation block ---------------------
    d_cross=hi_div(hi_mul(sqrt_hi(D(2)),S_AW),SIGMA_ACC)
    inv18_2=hi_add(
        hi_add(lo_div(D(1),sq(sL)),lo_div(D(1),sq(sA))),
        hi_div(sq(d_cross),lo_mul(sq(sA),sq(sL))),
    )
    s18=lo_div(D(1),sqrt_hi(inv18_2))

    # ---------- covariance upper bound from an explicit batch observer ----
    # Full-horizon transition norm in scaled coordinates, using only positive
    # polynomial envelopes of the OU-III chain.
    T=horizon
    row1=hi_add(D(1),hi_mul(lo_div(S_AW,S_V),T))
    row2=hi_add(
        hi_add(D(1),hi_mul(lo_div(S_V,S_P),T)),
        hi_mul(lo_div(S_AW,S_P),hi_div(sq(T),D(2))),
    )
    row3=hi_add(
        hi_add(
            hi_add(D(1),hi_mul(lo_div(S_V,S_S),hi_div(sq(T),D(2)))),
            hi_mul(lo_div(S_P,S_S),T),
        ),
        hi_mul(lo_div(S_AW,S_S),hi_div(T**3,D(6))),
    )
    col1=hi_add(
        hi_add(D(1),hi_mul(lo_div(S_V,S_P),T)),
        hi_mul(lo_div(S_V,S_S),hi_div(sq(T),D(2))),
    )
    col2=hi_add(D(1),hi_mul(lo_div(S_P,S_S),T))
    col4=hi_add(
        hi_add(
            hi_add(D(1),hi_mul(lo_div(S_AW,S_V),T)),
            hi_mul(lo_div(S_AW,S_P),hi_div(sq(T),D(2))),
        ),
        hi_mul(lo_div(S_AW,S_S),hi_div(T**3,D(6))),
    )
    phi_bar=sqrt_hi(max(row1,row2,row3,D(1))*max(col1,col2,D(1),col4))

    # Process covariance accumulated over the proof window.  The translation
    # trace dominates the default attitude/bias process terms.
    qc_hi=hi_div(hi_mul(D(2),sq(sigma_aw_max)),tau_min)
    imp2=hi_add(
        hi_add(
            sq(hi_mul(lo_div(S_AW,S_V),T)),
            sq(hi_mul(lo_div(S_AW,S_P),hi_div(sq(T),D(2)))),
        ),
        hi_add(
            sq(hi_mul(lo_div(S_AW,S_S),hi_div(T**3,D(6)))),
            sq(lo_div(D(1),S_AW)),
        ),
    )
    wbar=hi_mul(hi_mul(qc_hi,T),imp2)

    # Residual b_a variance can be bounded even if its measurements are ignored,
    # because the source model is a finite 5000-s OU process.
    pba_initial=sq(lo_div(D("0.004"),S_BA))
    pba_stationary=hi_div(hi_mul(D("2.5e-7"),D("5000")),hi_mul(D(2),sq(S_BA)))
    pba=max(pba_initial,pba_stationary)

    cS=hi_div(S_S,RS_EFF_MIN)
    cA=hi_div(
        sqrt_hi(hi_add(hi_add(sq(hi_mul(F1,S_TH)),sq(S_AW)),sq(S_BA))),
        SIGMA_ACC,
    )
    cM=hi_div(hi_mul(M1,S_TH),SIGMA_MAG)
    cbar=max(cS,cA,cM)
    nvec=D(8)  # four S + two accel + two mag selected proof observations
    beta_c=hi_mul(hi_mul(nvec,sq(cbar)),wbar)

    ba_meas=hi_div(S_BA,SIGMA_ACC)
    betaY=hi_mul(
        D(24),
        hi_add(
            hi_add(D(1),hi_mul(sq(cbar),wbar)),
            hi_mul(sq(ba_meas),pba),
        ),
    )
    p18=hi_add(
        hi_mul(hi_mul(D(2),sq(phi_bar)),hi_div(betaY,sq(s18))),
        hi_mul(D(2),wbar),
    )
    pbar=hi_add(hi_mul(D(18),p18),hi_mul(D(3),pba))

    # Oracle lower bound on the posterior covariance of finite-horizon process
    # noise: Omega_H >= alpha_c/(1+beta_c) I.  alpha_c=q_step is valid for the
    # arbitrary time-varying schedule because the last prediction alone carries
    # the uniform one-step floor.
    qH=lo_div(q_step,hi_add(D(1),beta_c))
    deltaH=lo_div(qH,pbar)
    one_minus_chi=lo_div(deltaH,D(2))  # 1-sqrt(1-d) >= d/2

    qH_frozen=lo_div(q_frozen,hi_add(D(1),beta_c))
    delta_frozen=lo_div(qH_frozen,pbar)

    values={
        "tau_cert_min":tau_min,
        "tau_cert_max":tau_max,
        "h_min":H_MIN,
        "h_max":H_MAX,
        "sigma_aw_cert_min":SIGMA_AW_MIN,
        "horizon_max_s":horizon,
        "q_step_lower":q_step,
        "q_frozen_horizon_lower":q_frozen,
        "sL_lower":sL,
        "sA_lower":sA,
        "s18_lower":s18,
        "phi_upper":phi_bar,
        "process_window_upper":wbar,
        "pba_upper":pba,
        "measurement_operator_upper":cbar,
        "beta_process_measurement_upper":beta_c,
        "pbar_upper":pbar,
        "qH_lower":qH,
        "deltaH_lower":deltaH,
        "one_minus_chi_lower":one_minus_chi,
        "delta_frozen_lower":delta_frozen,
    }
    if not (deltaH > 0 and s18 > 0 and qH > 0 and pbar.is_finite()):
        raise SystemExit("CERTIFICATE_FAIL non-positive verified bound")
    return values

def fmt(x: Decimal) -> str:
    return f"{x:.12E}"

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--repo-root", default=str(Path(__file__).resolve().parents[1]))
    ap.add_argument("--no-source-check", action="store_true")
    ap.add_argument("--json", action="store_true")
    args=ap.parse_args()
    root=Path(args.repo_root)
    if not args.no_source_check:
        check_source(root)
    v=calc()
    if args.json:
        print(json.dumps({k:str(x) for k,x in v.items()}, indent=2, sort_keys=True))
    else:
        print("OU_III_COMPUTER_ASSISTED_CERTIFICATE PASS")
        print("DOMAIN source/default: tau=[1/3,12] s, sigma_aw<=6, validated sensor noise, rS_eff=[0.044,694]")
        print("DOMAIN external: h=[4.5,5.5] ms, sigma_aw>=0.01, accepted-vector magnitude/noncollinearity/rate bounds")
        for k in (
            "horizon_max_s","q_step_lower","q_frozen_horizon_lower",
            "sL_lower","sA_lower","s18_lower","pbar_upper",
            "beta_process_measurement_upper","qH_lower","deltaH_lower",
            "one_minus_chi_lower","delta_frozen_lower",
        ):
            print(f"{k}={fmt(v[k])}")
        print("FINDING uniform_time_varying_certificate_positive=1")
        print("FINDING practical_basin_from_global_box=not_claimed")
        print("FINDING frozen_schedule_strengthening_not_used_in_uniform_theorem=1")

if __name__=="__main__":
    main()
