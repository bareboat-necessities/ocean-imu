#!/usr/bin/env python3
"""
Computer-assisted analytical certificate for the OU-III Live basin proof.

No sea replay is used. The program evaluates closed-form worst-case inequalities
with high-precision Decimal directed rounding. Algebraic lower/upper bounds use
explicit floor/ceiling operations; transcendental results are evaluated at extra
precision and widened by one ulp in the adverse direction.
"""
from __future__ import annotations
import argparse
import json
from decimal import Decimal, getcontext, localcontext, ROUND_FLOOR, ROUND_CEILING
from pathlib import Path

getcontext().prec = 100
D = Decimal
PREC = 100


def lo_add(a,b):
    with localcontext() as c:
        c.prec=PREC; c.rounding=ROUND_FLOOR
        return a+b

def hi_add(a,b):
    with localcontext() as c:
        c.prec=PREC; c.rounding=ROUND_CEILING
        return a+b

def lo_sub(a,b):
    assert a >= b >= 0
    with localcontext() as c:
        c.prec=PREC; c.rounding=ROUND_FLOOR
        return a-b

def lo_mul(a,b):
    assert a >= 0 and b >= 0
    with localcontext() as c:
        c.prec=PREC; c.rounding=ROUND_FLOOR
        return a*b

def hi_mul(a,b):
    assert a >= 0 and b >= 0
    with localcontext() as c:
        c.prec=PREC; c.rounding=ROUND_CEILING
        return a*b

def lo_div(a,b):
    assert a >= 0 and b > 0
    with localcontext() as c:
        c.prec=PREC; c.rounding=ROUND_FLOOR
        return a/b

def hi_div(a,b):
    assert a >= 0 and b > 0
    with localcontext() as c:
        c.prec=PREC; c.rounding=ROUND_CEILING
        return a/b

def lo_sq(x): return lo_mul(x,x)
def hi_sq(x): return hi_mul(x,x)
def lo_cube(x): return lo_mul(lo_sq(x),x)
def hi_cube(x): return hi_mul(hi_sq(x),x)


def _ulp(x,prec=90):
    if x.is_zero(): return D(10) ** D(-prec)
    return D(10) ** D(x.adjusted()-prec+1)

def exp_hi(x):
    with localcontext() as c:
        c.prec=PREC+20
        y=x.exp()
    return y+_ulp(y,PREC)

def exp_neg_lo(x): return lo_div(D(1),exp_hi(x))
def sqrt_hi(x):
    assert x >= 0
    with localcontext() as c:
        c.prec=PREC+20
        y=x.sqrt()
    return y+_ulp(y,PREC)
def sqrt_lo(x):
    assert x >= 0
    with localcontext() as c:
        c.prec=PREC+20
        y=x.sqrt()
    return max(D(0),y-_ulp(y,PREC))
def ceil_decimal(x):
    n=int(x)
    return n if D(n)==x else n+1


S_TH=D("0.087"); S_BG=D("0.001"); S_V=D("1"); S_P=D("20")
S_S=D("50"); S_AW=D("6"); S_BA=D("0.5")
H_MIN=D("0.0045"); H_MAX=D("0.0055")
SIGMA_AW_MIN=D("0.05")
F0=D("5"); F1=D("15"); M0=D("10"); M1=D("100")
DELTA_GEOM=D("0.05"); OMEGA_MAX=D("1.5")
VEC_GAP_MIN=D("0.005"); VEC_GAP_MAX=D("0.05")
SIGMA_ACC=D("0.0148"); SIGMA_GYRO=D("0.00157"); SIGMA_MAG=D("0.25")
GYRO_BIAS_Q=D("1e-11"); Q_BA=D("2.5e-7"); TAU_BA=D("5000")
SIGMA_BA0=D("0.004"); RS_EFF_MIN=D("0.044"); RS_EFF_MAX=D("694")


def require_token(text,token,where):
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
        "const float sigma_floor = std::max(0.05f, band_noise_floor_sigma_());",
        "mekf_->set_aw_stationary_std(aw_std);",
    ):
        require_token(wrapper,token,"wrapper")
    for token in (
        "T sigma_bacc0_ = T(0.004);",
        "Matrix3 Q_bacc_ = Matrix3::Identity() * T(2.5e-7);",
        "T tau_bacc_ = T(5000.0);",
    ):
        require_token(core,token,"core")
    for token in (
        "Eigen::Vector3f::Constant(0.0148f)",
        "Eigen::Vector3f::Constant(0.00157f)",
        "Eigen::Vector3f::Constant(0.25f)",
        "rs_eff_lo > 0.044f",
        "rs_eff_hi > 690.0f",
        "rs_eff_hi < 694.0f",
    ):
        require_token(contract,token,"iss_contract-test")


def calc():
    tau_min=lo_div(D("0.5"),D("1.5")); tau_max=D("12")
    sigma_aw_max=D("6"); lambda_max=hi_div(D(1),tau_min)
    cT_hi=hi_div(D("0.015"),D("1.1"))
    ts_max=min(D("0.25"),hi_mul(cT_hi,tau_max))
    gap_plus=hi_add(ts_max,H_MAX)
    # Event-indexed horizon: it begins on an S update. Three bounded gaps then
    # contain the next three S updates required by the four-point inverse.
    horizon=hi_mul(D(3),gap_plus)
    n_steps=ceil_decimal(hi_div(horizon,H_MIN))

    qc_lo=lo_div(lo_mul(D(2),lo_sq(SIGMA_AW_MIN)),tau_max)
    m_step=min(
        lo_div(H_MIN,S_V),lo_div(lo_sq(H_MIN),S_P),
        lo_div(lo_cube(H_MIN),S_S),lo_div(D(1),S_AW),
    )
    x_hi=hi_div(H_MAX,tau_min)
    ratio_437=hi_div(D(437),D(315))
    mu_den=hi_mul(D(870_912_000),hi_cube(ratio_437))
    mu_step=lo_div(exp_neg_lo(hi_mul(D(2),x_hi)),mu_den)
    q_step=lo_mul(lo_mul(lo_mul(qc_lo,H_MIN),lo_sq(m_step)),mu_step)

    # Arbitrary-time-varying finite-horizon controllability of one OU-III axis.
    A1=lo_div(lo_mul(S_AW,horizon),S_V)
    A2=lo_div(lo_mul(S_AW,lo_sq(horizon)),S_P)
    A3=lo_div(lo_mul(S_AW,lo_cube(horizon)),S_S)
    invM_frob2=hi_add(
        hi_add(hi_div(D(964800),lo_sq(A1)),hi_div(D(26661600),lo_sq(A2))),
        hi_add(hi_div(D(114307200),lo_sq(A3)),D(5741)),
    )
    GROW=hi_div(D(241),D(35)); HROW=hi_div(D(19),D(20))
    energy_shape=hi_add(
        hi_div(GROW,horizon),
        hi_mul(hi_mul(hi_sq(lambda_max),horizon),HROW),
    )
    steer_energy=hi_mul(
        hi_div(hi_mul(D(2),hi_sq(S_AW)),qc_lo),
        hi_mul(energy_shape,invM_frob2),
    )
    alpha_translation=lo_div(D(1),steer_energy)

    # Complete 21-state process floor: the last sample alone gives much larger
    # floors on attitude/gyro bias and residual accelerometer bias.
    qg_lo=lo_div(lo_sq(SIGMA_GYRO),hi_sq(S_TH))
    qg_hi=hi_div(hi_sq(SIGMA_GYRO),lo_sq(S_TH))
    qbg_lo=lo_div(GYRO_BIAS_Q,hi_sq(S_BG))
    qbg_hi=hi_div(GYRO_BIAS_Q,lo_sq(S_BG))
    cb_hi=hi_div(S_BG,S_TH)
    eta_den=hi_mul(hi_mul(D(2),qbg_hi),hi_mul(hi_sq(cb_hi),hi_sq(H_MAX)))
    eta_lo=min(D(1),lo_div(lo_mul(D(3),qg_lo),eta_den))
    eta_frac_lo=lo_div(eta_lo,hi_add(D(1),eta_lo))
    alpha_attitude=lo_mul(
        H_MIN,
        min(lo_div(qg_lo,D(2)),lo_mul(qbg_lo,eta_frac_lo)),
    )

    xb=hi_div(hi_mul(D(2),H_MIN),TAU_BA)
    # 1-exp(-x) >= x-x^2/2 for x>=0; x is exact here up to the outward upper
    # representation, so subtract an upper square to preserve a lower bound.
    xb_lo=lo_div(lo_mul(D(2),H_MIN),TAU_BA)
    one_minus_exp_b=lo_sub(xb_lo,hi_div(hi_sq(xb),D(2)))
    qba_step=lo_div(
        lo_mul(lo_mul(Q_BA,lo_div(TAU_BA,D(2))),one_minus_exp_b),
        hi_sq(S_BA),
    )
    alpha_c=min(alpha_translation,alpha_attitude,qba_step)

    # Frozen-schedule comparison, not used by the theorem.
    mh=min(
        lo_div(horizon,S_V),lo_div(lo_sq(horizon),S_P),
        lo_div(lo_cube(horizon),S_S),lo_div(D(1),S_AW),
    )
    mu_h=lo_div(exp_neg_lo(hi_mul(D(2),hi_div(horizon,tau_min))),mu_den)
    q_frozen=lo_mul(lo_mul(lo_mul(qc_lo,horizon),lo_sq(mh)),mu_h)

    # Constructive four-S divided-difference inverse.
    root5=sqrt_hi(D(5)); exp_obs=exp_hi(hi_div(hi_mul(D(3),gap_plus),tau_min))
    ka_num=hi_mul(hi_mul(D(2),root5),hi_mul(RS_EFF_MAX,exp_obs))
    k_a_phys=hi_div(ka_num,lo_cube(H_MIN))
    t3=hi_mul(D(3),gap_plus); a3max=hi_div(hi_cube(t3),D(6))
    yq=hi_add(RS_EFF_MAX,hi_mul(hi_mul(D(2),a3max),k_a_phys))
    kv=hi_mul(hi_div(sqrt_hi(D(6)),lo_sq(H_MIN)),yq)
    kp=hi_add(
        hi_mul(hi_div(sqrt_hi(D(2)),H_MIN),yq),
        hi_mul(hi_mul(D("0.5"),kv),gap_plus),
    )
    invL2=hi_add(
        hi_add(hi_sq(hi_div(kv,S_V)),hi_sq(hi_div(kp,S_P))),
        hi_add(hi_sq(hi_div(yq,S_S)),hi_sq(hi_div(k_a_phys,S_AW))),
    )
    sL=lo_div(D(1),sqrt_hi(invL2))

    # Two-vector attitude/gyro-bias observation block. Instead of subtracting
    # nearly equal square roots, use det/trace for the 2x2 Gramian:
    # lambda_min >= b^2/(2+b^2).
    racc_hi=hi_sq(SIGMA_ACC); rmag_hi=hi_sq(SIGMA_MAG)
    geom=min(lo_div(lo_sq(F0),racc_hi),lo_div(lo_sq(M0),rmag_hi))
    gamma=lo_mul(lo_mul(lo_sq(S_TH),DELTA_GEOM),geom)
    xx=hi_div(hi_mul(OMEGA_MAX,VEC_GAP_MAX),D(2))
    sinc_lo=lo_sub(D(1),hi_div(hi_sq(xx),D(6)))
    b0=lo_mul(lo_mul(lo_div(S_BG,S_TH),VEC_GAP_MIN),sinc_lo)
    b2_lo=lo_sq(b0); b2_hi=hi_sq(b0)
    lam=lo_div(b2_lo,hi_add(D(2),b2_hi))
    sA=sqrt_lo(lo_mul(gamma,lam))

    d_cross=hi_div(hi_mul(sqrt_hi(D(2)),S_AW),SIGMA_ACC)
    inv18_2=hi_add(
        hi_add(hi_div(D(1),lo_sq(sL)),hi_div(D(1),lo_sq(sA))),
        hi_div(hi_sq(d_cross),lo_mul(lo_sq(sA),lo_sq(sL))),
    )
    s18=lo_div(D(1),sqrt_hi(inv18_2))

    # Explicit finite-window comparison observer.
    T=horizon
    T2_hi=hi_sq(T); T3_hi=hi_cube(T)
    row1=hi_add(D(1),hi_mul(lo_div(S_AW,S_V),T))
    row2=hi_add(
        hi_add(D(1),hi_mul(lo_div(S_V,S_P),T)),
        hi_mul(lo_div(S_AW,S_P),hi_div(T2_hi,D(2))),
    )
    row3=hi_add(
        hi_add(
            hi_add(D(1),hi_mul(lo_div(S_V,S_S),hi_div(T2_hi,D(2)))),
            hi_mul(lo_div(S_P,S_S),T),
        ),
        hi_mul(lo_div(S_AW,S_S),hi_div(T3_hi,D(6))),
    )
    col1=hi_add(
        hi_add(D(1),hi_mul(lo_div(S_V,S_P),T)),
        hi_mul(lo_div(S_V,S_S),hi_div(T2_hi,D(2))),
    )
    col2=hi_add(D(1),hi_mul(lo_div(S_P,S_S),T))
    col4=hi_add(
        hi_add(
            hi_add(D(1),hi_mul(lo_div(S_AW,S_V),T)),
            hi_mul(lo_div(S_AW,S_P),hi_div(T2_hi,D(2))),
        ),
        hi_mul(lo_div(S_AW,S_S),hi_div(T3_hi,D(6))),
    )
    phi_sq_hi=hi_mul(max(row1,row2,row3,D(1)),max(col1,col2,D(1),col4))
    phi_bar=sqrt_hi(phi_sq_hi)

    qc_hi=hi_div(hi_mul(D(2),hi_sq(sigma_aw_max)),tau_min)
    imp2=hi_add(
        hi_add(
            hi_sq(hi_mul(lo_div(S_AW,S_V),T)),
            hi_sq(hi_mul(lo_div(S_AW,S_P),hi_div(T2_hi,D(2)))),
        ),
        hi_add(
            hi_sq(hi_mul(lo_div(S_AW,S_S),hi_div(T3_hi,D(6)))),
            hi_sq(hi_div(D(1),S_AW)),
        ),
    )
    w_trans=hi_mul(hi_mul(qc_hi,T),imp2)
    w_att=hi_mul(
        D(3),
        hi_add(
            hi_mul(qg_hi,T),
            hi_mul(qbg_hi,hi_add(T,hi_mul(hi_sq(cb_hi),hi_div(T3_hi,D(3))))),
        ),
    )
    pba_initial=hi_sq(hi_div(SIGMA_BA0,S_BA))
    pba_stationary=hi_div(hi_mul(Q_BA,TAU_BA),lo_mul(D(2),lo_sq(S_BA)))
    pba=max(pba_initial,pba_stationary)
    w_ba=hi_mul(D(3),pba)
    wbar=max(w_trans,w_att,w_ba)

    cS=hi_div(S_S,RS_EFF_MIN)
    cA_sq=hi_add(
        hi_add(hi_sq(hi_mul(F1,S_TH)),hi_sq(S_AW)),hi_sq(S_BA)
    )
    cA=hi_div(sqrt_hi(cA_sq),SIGMA_ACC)
    cM=hi_div(hi_mul(M1,S_TH),SIGMA_MAG)
    cbar=max(cS,cA,cM)

    # Eight 3-vector observations = 24 scalar output components.
    ncomp=D(24)
    beta_c=hi_mul(hi_mul(ncomp,hi_sq(cbar)),wbar)
    ba_meas=hi_div(S_BA,SIGMA_ACC)
    betaY=hi_mul(
        D(24),
        hi_add(
            hi_add(D(1),hi_mul(hi_sq(cbar),wbar)),
            hi_mul(hi_sq(ba_meas),pba),
        ),
    )
    p18=hi_add(
        hi_mul(hi_mul(D(2),hi_sq(phi_bar)),hi_div(betaY,lo_sq(s18))),
        hi_mul(D(2),wbar),
    )
    pbar=hi_add(hi_mul(D(18),p18),hi_mul(D(3),pba))

    qH=lo_div(alpha_c,hi_add(D(1),beta_c))
    delta_cov=lo_div(qH,pbar)

    # Independent observable/detectable Riccati-energy closure. At the
    # homogeneous linearization the ordinary MEKF reset Jacobian is identity;
    # its first nonidentity term is quadratic in error and is charged below in
    # c_xi. Thus the correction-feedback relation is linear with no hidden
    # reset-norm factor.
    gain_R=hi_div(sqrt_hi(pbar),D(2))
    feedback_block=hi_mul(hi_mul(cbar,phi_bar),gain_R)
    Lbar=hi_mul(D(7),feedback_block)
    innovation_ratio=hi_add(D(1),hi_mul(pbar,hi_sq(cbar)))
    c_meas_energy=lo_div(
        D(1),hi_mul(innovation_ratio,hi_sq(hi_add(D(1),Lbar)))
    )
    d_b=lo_div(
        qba_step,
        hi_mul(hi_sq(pbar),hi_add(D(1),hi_div(qba_step,pbar))),
    )
    Db=hi_mul(sqrt_hi(D(2)),ba_meas)
    l11=hi_div(D(1),lo_mul(sqrt_lo(c_meas_energy),s18))
    l12=hi_div(Db,lo_mul(s18,sqrt_lo(d_b)))
    l22=hi_div(D(1),sqrt_lo(d_b))
    Lstar_frob2=hi_add(hi_add(hi_sq(l11),hi_sq(l12)),hi_sq(l22))
    g_detect=lo_div(D(1),Lstar_frob2)
    delta_detect=lo_mul(qH,g_detect)

    deltaH=max(delta_cov,delta_detect)
    # For 0<=delta<=1, 1-sqrt(1-delta) >= delta/2.
    one_minus_chi=lo_div(deltaH,D(2))

    # Nonlinear remainder and lifted Riccati radius.
    cAu=sqrt_hi(cA_sq); cMu=hi_mul(M1,S_TH); cSu=S_S
    kappa_a=hi_div(hi_mul(pbar,cAu),racc_hi)
    kappa_m=hi_div(hi_mul(pbar,cMu),rmag_hi)
    kappa_s=hi_div(hi_mul(pbar,cSu),lo_sq(RS_EFF_MIN))
    c_model=hi_add(
        hi_mul(kappa_a,hi_add(
            hi_mul(D("0.5"),hi_mul(F1,hi_sq(S_TH))),hi_mul(S_TH,S_AW))),
        hi_mul(kappa_m,hi_mul(D("0.5"),hi_mul(M1,hi_sq(S_TH)))),
    )
    a_pred=hi_mul(H_MAX,S_BG)
    a_acc=hi_mul(
        hi_mul(S_TH,kappa_a),hi_add(hi_add(hi_mul(F1,S_TH),S_AW),S_BA))
    a_mag=hi_mul(hi_mul(S_TH,kappa_m),hi_mul(M1,S_TH))
    a_s=hi_mul(hi_mul(S_TH,kappa_s),S_S)
    def group_c(a):
        # j_c<=1 on the chosen theta_c=pi/2 chart.
        return hi_mul(hi_div(D(1),S_TH),hi_mul(hi_add(S_TH,a),a))
    c_group=group_c(a_pred)
    c_group=hi_add(c_group,group_c(a_acc))
    c_group=hi_add(c_group,group_c(a_mag))
    c_group=hi_add(c_group,group_c(a_s))
    c_xi=hi_add(c_model,c_group)

    THETA_C=D("1.57079632679489661923132169163975144")
    r_xi=min(
        lo_div(THETA_C,hi_add(S_TH,a_pred)),
        lo_div(THETA_C,hi_add(S_TH,a_acc)),
        lo_div(THETA_C,hi_add(S_TH,a_mag)),
        lo_div(THETA_C,hi_add(S_TH,a_s)),
    )
    p_intermediate=lo_div(
        q_step,hi_add(D(1),hi_mul(q_step,hi_mul(D(3),hi_sq(cbar))))
    )
    cV=hi_div(hi_mul(c_xi,pbar),sqrt_lo(p_intermediate))
    C_H=hi_mul(D(n_steps),cV)
    r_metric_nl=lo_div(one_minus_chi,C_H)
    r_metric_chart=lo_div(r_xi,sqrt_hi(pbar))
    r_metric=min(r_metric_nl,r_metric_chart)

    qH_frozen=lo_div(q_frozen,hi_add(D(1),beta_c))
    delta_frozen=lo_div(qH_frozen,pbar)

    values={
        "tau_cert_min":tau_min,"tau_cert_max":tau_max,
        "h_min":H_MIN,"h_max":H_MAX,
        "sigma_aw_cert_min":SIGMA_AW_MIN,
        "horizon_max_s":horizon,"horizon_steps_upper":D(n_steps),
        "q_step_lower":q_step,
        "alpha_translation_lower":alpha_translation,
        "alpha_attitude_lower":alpha_attitude,
        "qba_step_lower":qba_step,
        "alpha_controllability_lower":alpha_c,
        "q_frozen_horizon_lower":q_frozen,
        "sL_lower":sL,"sA_lower":sA,"s18_lower":s18,
        "phi_upper":phi_bar,"process_window_upper":wbar,
        "pba_upper":pba,"measurement_operator_upper":cbar,
        "beta_process_measurement_upper":beta_c,"pbar_upper":pbar,
        "qH_lower":qH,"delta_covariance_lower":delta_cov,
        "measurement_dissipation_lower":c_meas_energy,
        "detectability_db_lower":d_b,"detectability_G_lower":g_detect,
        "delta_detectability_lower":delta_detect,
        "deltaH_lower":deltaH,"one_minus_chi_lower":one_minus_chi,
        "cxi_upper":c_xi,"rxi_lower":r_xi,"CH_upper":C_H,
        "riccati_nonlinear_radius_lower":r_metric_nl,
        "riccati_chart_radius_lower":r_metric_chart,
        "riccati_cert_radius_lower":r_metric,
        "delta_frozen_lower":delta_frozen,
    }
    required=(alpha_c,s18,qH,deltaH,one_minus_chi,r_metric)
    if not all(x>0 and x.is_finite() for x in required):
        raise SystemExit("CERTIFICATE_FAIL non-positive verified bound")
    return values


def fmt(x): return f"{x:.12E}"


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--repo-root",default=str(Path(__file__).resolve().parents[1]))
    ap.add_argument("--no-source-check",action="store_true")
    ap.add_argument("--json",action="store_true")
    args=ap.parse_args()
    if not args.no_source_check: check_source(Path(args.repo_root))
    v=calc()
    if args.json:
        print(json.dumps({k:str(x) for k,x in v.items()},indent=2,sort_keys=True)); return
    print("OU_III_COMPUTER_ASSISTED_CERTIFICATE PASS")
    print("DOMAIN source/default: tau=[1/3,12] s, sigma_aw=[0.05,6], validated sensor noise, rS_eff=[0.044,694]")
    print("DOMAIN external: h=[4.5,5.5] ms, accepted-vector magnitude/noncollinearity/rate bounds")
    for k in (
        "horizon_max_s","q_step_lower","alpha_translation_lower",
        "alpha_attitude_lower","qba_step_lower","alpha_controllability_lower",
        "q_frozen_horizon_lower","sL_lower","sA_lower","s18_lower",
        "pbar_upper","beta_process_measurement_upper","qH_lower",
        "delta_covariance_lower","delta_detectability_lower","deltaH_lower",
        "one_minus_chi_lower","cxi_upper","rxi_lower",
        "riccati_cert_radius_lower",
    ):
        print(f"{k}={fmt(v[k])}")
    print("FINDING finite_horizon_controllability_replaces_h7_primary_floor=1")
    print("FINDING full_state_process_floor_verified=1")
    print("FINDING observable_detectable_block_verified=1")
    print("FINDING lifted_nonlinear_radius_verified=1")
    print("FINDING practical_basin_from_global_box=not_useful")


if __name__=="__main__": main()
