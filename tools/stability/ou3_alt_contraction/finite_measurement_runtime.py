"""Literal safe-LDLT control graph around the finite accepted measurement.

The exact mean/covariance measurement algebra lives in finite_core.  This layer
adds shipping's first-attempt / one-bump retry / rejection semantics.  Eigen
LDLT outcomes and the floating Frobenius norm used as ``noise_scale`` remain
explicit runtime/finite-precision witnesses; they are not inferred from rational
solvability and do not promote the theorem.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_prediction_graph as P
from tools.stability.ou3_alt_contraction import finite_core as CORE

BUMP_SCALE = F(1,10**6)


@dataclass(frozen=True)
class SafeLDLT:
    first_success: bool
    second_success: bool | None
    noise_scale: F
    machine_epsilon: F
    def __post_init__(self):
        if not isinstance(self.first_success,bool): raise TypeError('literal first LDLT outcome required')
        ns,eps=P.rational(self.noise_scale),P.rational(self.machine_epsilon)
        if ns < 0 or eps <= 0: raise ValueError('nonnegative noise scale and positive epsilon required')
        if self.first_success:
            if self.second_success is not None: raise ValueError('successful first factorization has no retry outcome')
        elif not isinstance(self.second_success,bool):
            raise TypeError('failed first factorization requires literal retry outcome')
        object.__setattr__(self,'noise_scale',ns); object.__setattr__(self,'machine_epsilon',eps)

    @property
    def bump(self):
        return max(self.machine_epsilon,BUMP_SCALE*(self.noise_scale+1))

    @property
    def accepted(self):
        return self.first_success or bool(self.second_success)

    @property
    def innovation_shift(self):
        return F(0) if self.first_success else self.bump


@dataclass(frozen=True)
class MeasurementRuntimeResult:
    state: CORE.State
    accepted: bool
    branch: SafeLDLT
    accepted_graph: CORE.Accepted | None


def measurement(state, kind, *, ldlt:SafeLDLT, **kwargs):
    """Apply exactly one shipping safe-LDLT branch.

    A double failure returns the predecessor state unchanged.  On either
    accepted branch the SAME diagonal shift is passed to finite_core, hence to
    both the inverse-free gain relation and Joseph covariance update.
    """
    if not isinstance(ldlt,SafeLDLT): raise TypeError('safe-LDLT runtime branch required')
    if not ldlt.accepted:
        return MeasurementRuntimeResult(state,False,ldlt,None)
    accepted=CORE.measurement(state,kind,innovation_shift=ldlt.innovation_shift,**kwargs)
    return MeasurementRuntimeResult(accepted.state,True,ldlt,accepted)


def readiness():
    return {
      'first_LDLT_success_branch':True,
      'single_diagonal_bump_retry_branch':True,
      'double_LDLT_failure_rejection_branch':True,
      'rejected_measurement_preserves_state_covariance':True,
      'same_retry_shift_used_by_gain_and_Joseph':True,
      'noise_scale_frobenius_runtime_source_attached':False,
      'Eigen_LDLT_outcomes_finite_precision_attached':False,
      'machine_epsilon_deployment_attached':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
