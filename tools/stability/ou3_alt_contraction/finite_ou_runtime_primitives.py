"""Finite OU runtime primitive binding for ALT; not source admission.

This module removes three detached operands from the higher-level prediction:
per-axis mean coefficients, active ``phi_hat``, and active ``Q_BB``.  They are
all generated from the same finite step, time constant and scalar decay used by
shipping.  The real source/frontend graph must still prove that ``alpha`` is the
runtime exp(-h/tau), and must still supply attitude F/Q, Qaxis regularization,
pending a_w inflation, scheduler branches and finite precision.
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

    def __post_init__(self):
        h, tau, alpha = P.rational(self.h), P.rational(self.tau), P.rational(self.alpha)
        if h <= 0 or tau <= 0:
            raise ValueError('positive h and tau required')
        if not 0 < alpha <= 1:
            raise ValueError('OU decay must lie in (0,1]')
        object.__setattr__(self, 'h', h)
        object.__setattr__(self, 'tau', tau)
        object.__setattr__(self, 'alpha', alpha)

    @property
    def x(self):
        return self.h / self.tau

    @property
    def small_coeff_branch(self):
        return abs(self.x) < F(1, 100)

    def axis_coefficients(self):
        """Literal shipping Phi coefficients from one (h,tau,alpha) root.

        ``phi_va`` uses expm1 exactly through alpha-1.  ``phi_pa`` and
        ``phi_Sa`` use the same small-x polynomial branch as safe_phi_A_coeffs;
        otherwise they use the expm1 closed form.  All three axes share the
        same tau in the current shipping implementation.
        """
        h, tau, alpha, x = self.h, self.tau, self.alpha, self.x
        phi_va = tau * (1 - alpha)
        if self.small_coeff_branch:
            x2, x3 = x*x, x*x*x
            x4, x5 = x3*x, x3*x*x
            phi_pa = tau*tau * (F(1,2)*x2 - F(1,6)*x3 + F(1,24)*x4)
            phi_Sa = tau*tau*tau * (F(1,6)*x3 - F(1,24)*x4 + F(1,120)*x5)
        else:
            em1 = alpha - 1
            phi_pa = tau*tau * (x + em1)
            phi_Sa = tau*tau*tau * (F(1,2)*x*x - x - em1)
        coeff = (phi_va, phi_pa, phi_Sa, alpha, h)
        return (coeff, coeff, coeff)


@dataclass(frozen=True)
class BiasDecay:
    """Literal H18/A21 accelerometer-bias prediction branch.

    Active covariance uses qd_scale = tau_b/2*(1-phi_b^2), algebraically equal
    to -tau_b/2*expm1(-2h/tau_b) when phi_b=exp(-h/tau_b).  Hence mean and
    covariance cannot use independent decay factors.
    """
    active: bool
    tau_b: F
    phi_b: F
    Q_bacc: tuple

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
        elif phi != 1:
            raise ValueError('held H18 mean factor must be exactly one')
        object.__setattr__(self, 'tau_b', tau)
        object.__setattr__(self, 'phi_b', phi)
        object.__setattr__(self, 'Q_bacc', tuple(map(tuple, q)))

    def covariance_increment(self):
        if not self.active:
            return M.zeros(3,3)
        scale = self.tau_b * (1 - self.phi_b*self.phi_b) / 2
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
    """Paired finite mean/covariance prediction with shared OU/BA decay roots."""
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
        'mean_axis_coefficients_from_one_h_tau_alpha_root': True,
        'small_safe_phi_A_branch_literal': True,
        'active_BA_QBB_from_same_phi_as_mean': True,
        'free_axis_coefficients_removed_at_runtime_entry': True,
        'free_phi_hat_removed_at_runtime_entry': True,
        'free_Q_BB_removed_at_runtime_entry': True,
        'alpha_exp_runtime_source_attached': False,
        'attitude_F_Q_primitives_source_attached': False,
        'Qaxis_analytic_and_regularization_source_attached': False,
        'pending_aw_covariance_inflation_attached': False,
        'periodic_S_service_attached': False,
        'finite_precision_attached': False,
        'complete_word_finite_identity': False,
        'ALT_LIVE_PASS': False,
    }
