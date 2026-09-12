"""Finite OU runtime primitive binding for ALT; not source admission.

This module removes detached transition and BA covariance operands while keeping
shipping's distinct transcendental evaluations visible. IntegratedOUChain calls
``exp(-x)`` and ``expm1(-x)`` separately; active BA prediction likewise calls
``exp(-h/tau_b)`` for the mean and ``expm1(-2h/tau_b)`` for Q_BB. At exact real
arithmetic these values satisfy em1=alpha-1 and em1_2=phi_b^2-1, but binary32
libm does not make those bit-identities automatic. The optional separate roots
below let the theorem-facing shipping layer retain that discrepancy explicitly.
Legacy component tests may omit them and recover the exact-real identities.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_prediction_graph as P
from tools.stability.ou3_alt_contraction import finite_prediction_covariance as C


@dataclass(frozen=True)
class OUDecay:
    """Shipping mean primitive after tau clamps, before source qualification."""
    h: F
    tau: F
    alpha: F
    em1: F | None = None

    def __post_init__(self):
        h, tau, alpha = P.rational(self.h), P.rational(self.tau), P.rational(self.alpha)
        if h <= 0 or tau <= 0:
            raise ValueError('positive h and tau required')
        if not 0 < alpha <= 1:
            raise ValueError('OU decay must lie in (0,1]')
        em1 = alpha-1 if self.em1 is None else P.rational(self.em1)
        if not -1 < em1 <= 0:
            raise ValueError('OU expm1 root must lie in (-1,0]')
        object.__setattr__(self, 'h', h)
        object.__setattr__(self, 'tau', tau)
        object.__setattr__(self, 'alpha', alpha)
        object.__setattr__(self, 'em1', em1)

    @property
    def x(self):
        return self.h / self.tau

    @property
    def small_coeff_branch(self):
        return abs(self.x) < F(1, 100)

    def axis_coefficients(self):
        """Literal shipping Phi graph from one h/tau and exp/expm1 roots.

        ``phi_va`` consumes make_prims().em1. The safe phi_pa/phi_Sa small-x
        branch is polynomial; the general branch independently consumes the
        same mathematical expm1(-x) value. Deployment correspondence between
        exp and expm1 remains an enclosure obligation, not an assumed bit law.
        """
        h, tau, alpha, em1, x = self.h, self.tau, self.alpha, self.em1, self.x
        phi_va = -tau * em1
        if self.small_coeff_branch:
            x2, x3 = x*x, x*x*x
            x4, x5 = x3*x, x3*x*x
            phi_pa = tau*tau * (F(1,2)*x2 - F(1,6)*x3 + F(1,24)*x4)
            phi_Sa = tau*tau*tau * (F(1,6)*x3 - F(1,24)*x4 + F(1,120)*x5)
        else:
            phi_pa = tau*tau * (x + em1)
            phi_Sa = tau*tau*tau * (F(1,2)*x*x - x - em1)
        coeff = (phi_va, phi_pa, phi_Sa, alpha, h)
        return (coeff, coeff, coeff)


@dataclass(frozen=True)
class BiasDecay:
    """Literal H18/A21 accelerometer-bias prediction branch.

    Shipping evaluates active mean and covariance transcendental terms
    separately. ``em1_2`` is expm1(-2h/tau_b). If omitted, component-level
    exact-real callers use the algebraic identity phi_b^2-1. The theorem-facing
    shipping word supplies it separately until binary32/libm error is enclosed.
    """
    active: bool
    tau_b: F
    phi_b: F
    Q_bacc: tuple
    em1_2: F | None = None

    def __post_init__(self):
        if not isinstance(self.active, bool):
            raise TypeError('literal held/active branch required')
        tau, phi = P.rational(self.tau_b), P.rational(self.phi_b)
        q = M.mat(self.Q_bacc, 3, 3)
        if q != M.transpose(q):
            raise ValueError('symmetric Q_bacc required')
        if tau < F(1,1000):
            raise ValueError('tau_b must be the post-clamp shipping value >=1e-3')
        if self.active:
            if not 0 < phi <= 1:
                raise ValueError('active phi_b must lie in (0,1]')
            em1_2 = phi*phi-1 if self.em1_2 is None else P.rational(self.em1_2)
            if not -1 < em1_2 <= 0:
                raise ValueError('active BA expm1 root must lie in (-1,0]')
        else:
            if phi != 1:
                raise ValueError('held H18 mean factor must be exactly one')
            if self.em1_2 is not None:
                raise ValueError('held H18 branch consumes no BA expm1 witness')
            em1_2 = F(0)
        object.__setattr__(self, 'tau_b', tau)
        object.__setattr__(self, 'phi_b', phi)
        object.__setattr__(self, 'Q_bacc', tuple(map(tuple, q)))
        object.__setattr__(self, 'em1_2', em1_2)

    def covariance_increment(self):
        if not self.active:
            return M.zeros(3,3)
        scale = -self.tau_b * self.em1_2 / 2
        return M.scaled(self.Q_bacc, scale)


def runtime_blocks(*, F_AA, Q_AA, ou: OUDecay, bias: BiasDecay,
                   qaxis_unit=None, sigma_aw=None, independent_qaxis=None):
    """Build the structured shipping covariance blocks without free F_LL/Q_BB."""
    coeff = ou.axis_coefficients()
    if (qaxis_unit is None) == (independent_qaxis is None):
        raise ValueError('select exactly one correlated or independent Q_LL branch')
    if qaxis_unit is not None:
        if sigma_aw is None:
            raise ValueError('correlated Q_LL requires Sigma_aw')
        qll = C.correlated_linear_process(qaxis_unit, sigma_aw)
    else:
        if sigma_aw is not None:
            raise ValueError('independent Q_LL does not consume Sigma_aw cross terms')
        qll = C.independent_linear_process(independent_qaxis)
    return C.Blocks(F_AA, Q_AA, C.linear_transition(coeff), qll,
                    bias.phi_b, bias.covariance_increment(), bias.active)


def runtime_paired_prediction(state, segment, *, gyro_body, ou: OUDecay,
                              bias: BiasDecay, F_AA, Q_AA,
                              qaxis_unit=None, sigma_aw=None,
                              independent_qaxis=None):
    """Paired finite mean/covariance prediction with shared source arguments."""
    if segment.h != ou.h:
        raise ValueError('OU runtime step detached from physical segment duration')
    if bias.active != (state.mode == 'A'):
        raise ValueError('BA runtime branch detached from H18/A21 state mode')
    blocks = runtime_blocks(F_AA=F_AA, Q_AA=Q_AA, ou=ou, bias=bias,
                            qaxis_unit=qaxis_unit, sigma_aw=sigma_aw,
                            independent_qaxis=independent_qaxis)
    phi_hat = bias.phi_b if bias.active else None
    return C.paired_prediction(state, segment, gyro_body=gyro_body,
                               axis_coefficients=ou.axis_coefficients(),
                               covariance_blocks=blocks, phi_hat=phi_hat)


def readiness():
    return {
        'mean_axis_coefficients_from_one_h_tau_exp_expm1_root': True,
        'small_safe_phi_A_branch_literal': True,
        'active_BA_mean_and_QBB_distinct_shipping_transcendentals_retained': True,
        'free_axis_coefficients_removed_at_runtime_entry': True,
        'free_phi_hat_removed_at_runtime_entry': True,
        'free_Q_BB_removed_at_runtime_entry': True,
        'exp_expm1_binary32_relation_enclosed': False,
        'alpha_exp_runtime_source_attached': False,
        'attitude_F_Q_primitives_source_attached': False,
        'Qaxis_analytic_and_regularization_source_attached': False,
        'pending_aw_covariance_inflation_attached': False,
        'periodic_S_service_attached': False,
        'finite_precision_attached': False,
        'complete_word_finite_identity': False,
        'ALT_LIVE_PASS': False,
    }
