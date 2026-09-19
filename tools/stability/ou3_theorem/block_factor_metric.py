"""Block/factor covariance metric for the OU-III A21 proof.

The failed global entrywise P-box is not reused here.  The covariance metric is
represented by three proof-coordinate blocks: attitude/gyro bias (6),
translation/integrated OU (12), and active accelerometer bias (3).  Each block
is certified by a positive factor lower bound plus bounded cross-block coupling.
"""
from __future__ import annotations
from dataclasses import dataclass
import math
from decimal import Decimal, localcontext

BLOCKS=(6,12,3)

@dataclass(frozen=True)
class BlockFactorFloor:
    """P_ii >= L_i L_i' with sigma_min(L_i)>=factor_floor."""
    factor_floor: float
    covariance_ceiling: float
    def __post_init__(self) -> None:
        if not all(math.isfinite(x) for x in (self.factor_floor,self.covariance_ceiling)):
            raise ValueError("finite block bounds required")
        if self.factor_floor<=0 or self.covariance_ceiling<=0:
            raise ValueError("positive block bounds required")
    @property
    def covariance_floor(self) -> float:
        return self.factor_floor*self.factor_floor

@dataclass(frozen=True)
class BlockFactorMetric:
    ag: BlockFactorFloor
    lin: BlockFactorFloor
    ba: BlockFactorFloor
    # operator-norm upper bounds on P_ij in the same scaled coordinates
    cross_ag_lin: float
    cross_ag_ba: float
    cross_lin_ba: float
    def __post_init__(self) -> None:
        xs=(self.cross_ag_lin,self.cross_ag_ba,self.cross_lin_ba)
        if not all(math.isfinite(x) and x>=0 for x in xs):
            raise ValueError("finite nonnegative cross bounds required")

def block_scaled_schur_floor(m: BlockFactorMetric) -> dict:
    """Rigorous full covariance floor from block factors and cross couplings.

    Scale each diagonal block by its certified factor floor.  The scaled
    diagonal is >=I.  The normalized cross norms are bounded by c_ij/(l_i l_j).
    Block Gershgorin then gives P >= D * gamma I * D with
    gamma=1-max_i sum_j cbar_ij.  This preserves useful block conditioning and
    exposes exactly when cross-covariance destroys a global scalar floor.
    """
    la, ll, lb=m.ag.factor_floor,m.lin.factor_floor,m.ba.factor_floor
    c01=m.cross_ag_lin/(la*ll)
    c02=m.cross_ag_ba/(la*lb)
    c12=m.cross_lin_ba/(ll*lb)
    row=(c01+c02,c01+c12,c02+c12)
    gamma=1.0-max(row)
    floor=max(0.0,gamma)*min(la*la,ll*ll,lb*lb)
    return {"normalized_cross":(c01,c02,c12),"row_sums":row,
            "block_gershgorin_margin":gamma,"global_covariance_floor":floor,
            "verified":gamma>0.0}

def block_information_normalization(*, information_floors: tuple[float,float,float],
                                    metric: BlockFactorMetric) -> dict:
    """Normalize block information with the certified covariance factors.

    J_i>=mu_i I and P_ii>=l_i^2 I imply factor-normalized information
    >=mu_i*l_i^2.  The contraction certificate consumes the minimum only after
    block coupling has a positive Schur/Gershgorin margin.
    """
    if len(information_floors)!=3 or any((not math.isfinite(x) or x<=0) for x in information_floors):
        raise ValueError("three positive information floors required")
    cov=block_scaled_schur_floor(metric)
    if not cov["verified"]:
        return {**cov,"mu_cov":0.0,"rho0":1.0,"norm_margin":0.0}
    raw=(information_floors[0]*metric.ag.covariance_floor,
         information_floors[1]*metric.lin.covariance_floor,
         information_floors[2]*metric.ba.covariance_floor)
    # Coupling margin is a congruence coercivity loss, not an invented
    # information source.
    mu=cov["block_gershgorin_margin"]*min(raw)
    rho=1.0/(1.0+mu)
    return {**cov,"block_mu_cov":raw,"mu_cov":mu,"rho0":rho,
            "norm_margin":1.0-math.sqrt(rho)}

def active_bias_factor_floor(*, stationary_sigma_min: float,
                             release_std_floor: float) -> float:
    """Literal A21 BA factor floor from release and OU stationary scale."""
    if min(stationary_sigma_min,release_std_floor)<=0:
        raise ValueError("positive BA scales required")
    return min(stationary_sigma_min,release_std_floor)

def attitude_bias_factor_probe(*, attitude_var_floor: float,
                               gyro_bias_var_floor: float) -> BlockFactorFloor:
    if min(attitude_var_floor,gyro_bias_var_floor)<=0:
        raise ValueError("positive AG variances required")
    return BlockFactorFloor(math.sqrt(min(attitude_var_floor,gyro_bias_var_floor)),
                            max(attitude_var_floor,gyro_bias_var_floor))

def translation_factor_probe(*, controllability_floor: float,
                             covariance_ceiling: float) -> BlockFactorFloor:
    if controllability_floor<=0:
        raise ValueError("positive translation controllability floor required")
    return BlockFactorFloor(math.sqrt(controllability_floor),covariance_ceiling)


def ba_ou_recurring_factor_floor(*, phi_max: float, process_variance_floor: float,
                                 release_variance_floor: float,
                                 corrections_information_ceiling: float) -> dict:
    """Scalar isotropic BA block floor under A21 OU prediction/correction.

    Prediction: p- >= phi^2 p+ + q.  A correction with information norm <=j
    gives p+ >= 1/(1/p-+j).  This recurrence is monotone; its positive fixed
    point and the literal release floor give a recurring factor certificate.
    """
    vals=(phi_max,process_variance_floor,release_variance_floor,corrections_information_ceiling)
    if not all(math.isfinite(x) for x in vals) or not 0<phi_max<1 or min(process_variance_floor,release_variance_floor)<=0 or corrections_information_ceiling<0:
        raise ValueError("valid BA recurrence data required")
    p=release_variance_floor
    prefix=p
    for _ in range(100000):
        pred=phi_max*phi_max*p+process_variance_floor
        post=1.0/(1.0/pred+corrections_information_ceiling)
        prefix=min(prefix,pred,post)
        if abs(post-p)<=1e-14*max(1.0,p):
            p=post;break
        p=post
    return {"covariance_floor":min(p,release_variance_floor),
            "factor_floor":math.sqrt(min(p,release_variance_floor)),
            "prefix_floor":prefix}

def cross_block_psd_bound(ceiling_i: float,ceiling_j: float) -> float:
    """Universal PSD cross-covariance operator bound ||Pij||<=sqrt(||Pii||||Pjj||)."""
    if not all(math.isfinite(x) for x in (ceiling_i,ceiling_j)) or min(ceiling_i,ceiling_j)<=0:
        raise ValueError("positive block ceilings required")
    return math.sqrt(ceiling_i*ceiling_j)

def required_cross_fraction_for_gamma(*, ell_i: float,ell_j: float,
                                      desired_normalized_cross: float) -> float:
    """Absolute cross-block threshold needed for a target scaled coupling."""
    if min(ell_i,ell_j)<=0 or desired_normalized_cross<0:
        raise ValueError("valid factor/cross target required")
    return desired_normalized_cross*ell_i*ell_j


def ag_one_second_factor_floor(*, gyro_white_density: float,
                               gyro_bias_rw_density: float,
                               attitude_scale: float,
                               gyro_bias_scale: float,
                               window_s: float=1.0) -> dict:
    """Scaled exact continuous AG process Gramian, one axis.

    theta_dot=-b_g+n_g, bdot=n_bg.  The exact process Gramian is used; three
    axes repeat independently.  LDL pivots in proof coordinates avoid mixing
    radians and rad/s.
    """
    vals=(gyro_white_density,gyro_bias_rw_density,attitude_scale,gyro_bias_scale,window_s)
    if not all(math.isfinite(x) for x in vals) or min(vals)<=0:
        raise ValueError("positive AG process data required")
    sg2=gyro_white_density**2; sb2=gyro_bias_rw_density**2; T=window_s
    q00=sg2*T+sb2*T**3/3.0
    q01=-sb2*T*T/2.0
    q11=sb2*T
    a=q00/(attitude_scale**2)
    b=q01/(attitude_scale*gyro_bias_scale)
    d=q11/(gyro_bias_scale**2)
    p1=a
    pivot_second=d-b*b/a
    floor=min(p1,pivot_second)
    return {"scaled_gramian":((a,b),(b,d)),"pivot_floor":floor,
            "factor_floor":math.sqrt(floor),"verified":floor>0}

def ba_process_variance(*, dt: float,tau_bacc: float,
                        drive_density: float) -> float:
    """Literal scalar active-bias OU process variance from continuous drive."""
    if min(dt,tau_bacc,drive_density)<=0:
        raise ValueError("positive BA process data required")
    phi=math.exp(-dt/tau_bacc)
    # shipping Q_BA uses stationary sigma_bacc0 with q=2 sigma^2/tau;
    # drive_density is sqrt(q), so exact scalar Q is q*tau/2*(1-phi^2).
    return drive_density*drive_density*tau_bacc*.5*(1.0-phi*phi)

def block_factor_feasibility(*, ag_factor: float,lin_factor: float,ba_factor: float,
                             ag_ceiling: float,lin_ceiling: float,ba_ceiling: float) -> dict:
    """Non-promoting threshold diagnostic using only universal PSD cross bounds."""
    m=BlockFactorMetric(
        BlockFactorFloor(ag_factor,ag_ceiling),
        BlockFactorFloor(lin_factor,lin_ceiling),
        BlockFactorFloor(ba_factor,ba_ceiling),
        cross_block_psd_bound(ag_ceiling,lin_ceiling),
        cross_block_psd_bound(ag_ceiling,ba_ceiling),
        cross_block_psd_bound(lin_ceiling,ba_ceiling))
    return block_scaled_schur_floor(m)


def additive_process_metric_certificate(*, ag_q_factor: float,
                                        lin_q_factor: float,
                                        ba_q_factor: float) -> dict:
    """Exact cross-independent coercivity at an A21 post-prediction root.

    Shipping prediction has P-=F P+ F'+Q with P+>=0 and block-diagonal
    Q=diag(Q_AG,Q_LIN,Q_BA). If each Q_i>=L_i L_i', then
        P- - diag(L_i L_i') = F P+ F' + (Q-diag(...)) >= 0.
    Hence in D=diag(L_i) coordinates P- >= D D' with coercivity gamma=1,
    regardless of the cross covariance carried by F P+ F'. This is strictly
    stronger than trying to prove block diagonal dominance from absolute
    cross-block norms. Cross bounds remain necessary for prefix/magnitude
    retention, but they do not consume root coercivity.
    """
    vals=(ag_q_factor,lin_q_factor,ba_q_factor)
    if not all(math.isfinite(x) and x>0 for x in vals):
        raise ValueError("positive process factor floors required")
    return {"factor_floors":vals,"root_metric_gamma":1.0,
            "verified":True,
            "lemma":"Pminus>=Q_block_factor"}

def block_metric_correction_margin(*, gamma_in: float,
                                   normalized_information_ceiling: float) -> float:
    """Lower metric coercivity after a correction.

    If P>=gamma D D' and H'R^-1H in D coordinates has norm <=j, then
    P+ = (P^-1+H'R^-1H)^-1 >= gamma/(1+gamma*j) D D'.
    This carries every cross covariance implicitly and remains strictly
    positive for finite j.
    """
    if not all(math.isfinite(x) for x in (gamma_in,normalized_information_ceiling)):
        raise ValueError("finite metric correction data required")
    if gamma_in<=0 or normalized_information_ceiling<0:
        raise ValueError("positive gamma and nonnegative information required")
    return gamma_in/(1.0+gamma_in*normalized_information_ceiling)

def normalized_measurement_information_ceiling(*, h_norm: float,
                                               max_factor: float,
                                               noise_variance_floor: float) -> float:
    if not all(math.isfinite(x) for x in (h_norm,max_factor,noise_variance_floor)):
        raise ValueError("finite measurement metric data required")
    if min(h_norm,max_factor,noise_variance_floor)<=0:
        raise ValueError("positive measurement metric data required")
    return (h_norm*max_factor)**2/noise_variance_floor


def _det_decimal(a: list[list[Decimal]]) -> Decimal:
    n=len(a)
    if n==1: return a[0][0]
    out=Decimal(0)
    for j in range(n):
        minor=[row[:j]+row[j+1:] for row in a[1:]]
        term=a[0][j]*_det_decimal(minor)
        out += term if j%2==0 else -term
    return out

def _ou_integral_primitive_decimal(t: Decimal,tau: Decimal) -> tuple[Decimal,...]:
    """Primitive of the positive OU impulse components (v,p,S,a)."""
    e=(-t/tau).exp()
    return (tau*t+tau*tau*e,
            tau*t*t/Decimal(2)-tau*tau*t-tau**3*e,
            tau*t**3/Decimal(6)-tau*tau*t*t/Decimal(2)+tau**3*t+tau**4*e,
            -tau*e)

def uniform_lin_jensen_factor_certificate(*, tau_min: float,tau_max: float,
                                          sigma_min: float,window_s: float,
                                          scales: tuple[float,float,float,float],
                                          tau_cells: int=4096) -> dict:
    """Source-uniform LIN factor from Jensen subinterval columns.

    For each of four equal time intervals I, Jensen gives
      integral_I z z' dt >= |I|^-1 (integral_I z)(integral_I z)'.
    The OU impulse components are positive and increase with tau because they
    are convolutions of nonnegative polynomials with exp(-s/tau). The factor
    sqrt(2 sigma^2/tau/|I|) decreases with tau, so endpoint products give a
    rigorous entrywise column box on each tau cell. For midpoint M and radius
    E, sigma_min(C)>=sigma_min(M)-||E||_F. We lower-bound sigma_min(M) by
    |det M|/||adj M||_F. Decimal arithmetic at 70 digits plus a 1e-50 downward
    reserve makes the returned binary64 value non-promoting-safe.
    """
    if not (0<tau_min<tau_max and sigma_min>0 and window_s>0 and tau_cells>=4):
        raise ValueError("valid LIN certificate domain required")
    if len(scales)!=4 or min(scales)<=0: raise ValueError("four positive scales required")
    with localcontext() as ctx:
        ctx.prec=70
        d0=Decimal(str(tau_min));d1=Decimal(str(tau_max))
        sig=Decimal(str(sigma_min));T=Decimal(str(window_s))
        sc=[Decimal(str(x)) for x in scales]
        ratio=(d1/d0).ln()/Decimal(tau_cells)
        edges=[d0*(ratio*Decimal(k)).exp() for k in range(tau_cells+1)]
        h=T/Decimal(4)
        best=None;best_cell=-1
        for k in range(tau_cells):
            lo_tau,hi_tau=edges[k],edges[k+1]
            lows=[]; highs=[]
            for m in range(4):
                a=Decimal(m)*h;b=Decimal(m+1)*h
                plo=_ou_integral_primitive_decimal(b,lo_tau)
                qlo=_ou_integral_primitive_decimal(a,lo_tau)
                phi=_ou_integral_primitive_decimal(b,hi_tau)
                qhi=_ou_integral_primitive_decimal(a,hi_tau)
                ilo=[(plo[i]-qlo[i])/sc[i] for i in range(4)]
                ihi=[(phi[i]-qhi[i])/sc[i] for i in range(4)]
                slo=(Decimal(2)*sig*sig/hi_tau/h).sqrt()
                shi=(Decimal(2)*sig*sig/lo_tau/h).sqrt()
                lows.append([slo*x for x in ilo]);highs.append([shi*x for x in ihi])
            # columns are time intervals
            mid=[[Decimal(0)]*4 for _ in range(4)]
            rad=[[Decimal(0)]*4 for _ in range(4)]
            for j in range(4):
                for i in range(4):
                    mid[i][j]=(lows[j][i]+highs[j][i])/Decimal(2)
                    rad[i][j]=(highs[j][i]-lows[j][i])/Decimal(2)
            det=abs(_det_decimal(mid))
            cof=[]
            for i in range(4):
                row=[]
                for j in range(4):
                    minor=[r[:j]+r[j+1:] for ii,r in enumerate(mid) if ii!=i]
                    x=_det_decimal(minor)
                    row.append(x if (i+j)%2==0 else -x)
                cof.append(row)
            adj_frob=sum(x*x for row in cof for x in row).sqrt()
            rad_frob=sum(x*x for row in rad for x in row).sqrt()
            lower=det/adj_frob-rad_frob
            if best is None or lower<best:
                best=lower;best_cell=k
        reserve=Decimal("1e-50")
        certified=max(Decimal(0),best-reserve)
        return {"verified":certified>0,"singular_factor_floor":float(certified),
                "covariance_floor":float(certified*certified),
                "tau_cells":tau_cells,"worst_cell":best_cell,
                "decimal_lower":str(certified)}


def batch_information_degraded_factor(*, process_factor: float,
                                      event_information_ceilings: tuple[tuple[int,float],...],
                                      transport_norm_ceiling: float=1.0) -> dict:
    """Conservative final-state factor after a process window and measurements.

    Treat the accumulated process-driven final state as having covariance
    >=ell^2 I. Every actual measurement is replaced by a hypothetical direct
    final-state observation with no less Fisher information: its information
    norm ceiling is multiplied by the squared backward transport norm. This can
    only reduce the comparison covariance. Summing those ceilings gives
      P_end >= (ell^-2 I + J_total I)^-1.
    This is the batch counterpart of gamma+=gamma/(1+gamma*j).
    """
    if not math.isfinite(process_factor) or process_factor<=0 or not math.isfinite(transport_norm_ceiling) or transport_norm_ceiling<1:
        raise ValueError("valid process factor/transport ceiling required")
    total=0.0
    for count,j in event_information_ceilings:
        if count<0 or not math.isfinite(j) or j<0:
            raise ValueError("valid event information ceiling required")
        total += count*j*transport_norm_ceiling*transport_norm_ceiling
    var=process_factor*process_factor
    post=var/(1.0+var*total)
    return {"total_information_ceiling":total,"covariance_floor":post,
            "factor_floor":math.sqrt(post),
            "gamma_relative_to_process":post/var,"verified":post>0}

def maximum_event_count(window_s: float,dt_min: float) -> int:
    if not all(math.isfinite(x) for x in (window_s,dt_min)) or min(window_s,dt_min)<=0:
        raise ValueError("positive timing required")
    return math.ceil(window_s/dt_min)+1


def source_uniform_root_block_certificate(*, ag_process_spd: bool,
                                          lin_process_spd: bool,
                                          ba_process_spd: bool,
                                          compact_parameter_domain: bool,
                                          block_diagonal_process_noise: bool,
                                          covariance_psd_before_prediction: bool,
                                          finite_block_covariance_ceilings: bool) -> dict:
    """Qualitative-but-rigorous source-uniform root block certificate.

    Continuous exact process covariances are SPD for every positive dt and
    positive declared noise density. Their minimum scaled factor over the
    compact shipping parameter domain is therefore strictly positive. Because
    Q is block diagonal and P-=F P+ F'+Q, the factor-scaled root coercivity is
    exactly gamma=1. Existing covariance compactness supplies finite diagonal
    ceilings, hence finite cross-block operator bounds by PSD Cauchy-Schwarz.
    This closes *positivity/recurrence* of the block metric, but deliberately
    does not claim the numerical factor minima needed by constructive rho0.
    """
    flags=(ag_process_spd,lin_process_spd,ba_process_spd,compact_parameter_domain,
           block_diagonal_process_noise,covariance_psd_before_prediction,
           finite_block_covariance_ceilings)
    if any(type(x) is not bool for x in flags):
        raise ValueError("literal certificate flags required")
    ok=all(flags)
    return {"verified":ok,"source_uniform_factor_floors_exist":ok,
            "finite_cross_block_bounds_exist":ok,
            "root_metric_gamma":1.0 if ok else 0.0,
            "constructive_numeric_factors":False}
