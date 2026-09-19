"""Block/factor covariance metric for the OU-III A21 proof.

The failed global entrywise P-box is not reused here.  The covariance metric is
represented by three proof-coordinate blocks: attitude/gyro bias (6),
translation/integrated OU (12), and active accelerometer bias (3).  Each block
is certified by a positive factor lower bound plus bounded cross-block coupling.
"""
from __future__ import annotations
from dataclasses import dataclass
import math

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
