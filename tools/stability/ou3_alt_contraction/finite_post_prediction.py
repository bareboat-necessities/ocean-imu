"""Finite post-prediction prefixes for ALT; no stability/source promotion.

Shipping orders these events after the ordinary covariance prediction:
  1. optional pending a_w covariance positive-part inflation;
  2. symmetry hygiene;
  3. periodic S-service scheduler decision;
  4. if due, the S=0 measurement through shipping safe-LDLT semantics.
This module represents those finite real-arithmetic relations explicitly.
"""
from __future__ import annotations
from dataclasses import dataclass, replace
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_prediction_graph as P
from tools.stability.ou3_alt_contraction import finite_core as CORE

OFF_AW = 15


def _diag(v):
    v=P.vec(v,3); out=M.zeros(3,3)
    for i,x in enumerate(v): out[i][i]=x
    return out


@dataclass(frozen=True)
class AwFloorResult:
    state: CORE.State
    pending_after: bool
    solver_success: bool | None
    inflation_applied: bool


def aw_floor(state, *, pending, target, solver_success=None,
             eigenvectors=None, eigenvalues=None):
    """Literal positive-part a_w covariance inflation."""
    if not isinstance(pending,bool): raise TypeError('literal pending branch required')
    target=M.mat(target,3,3)
    if target != M.transpose(target): raise ValueError('symmetric floor target required')
    if not pending:
        if solver_success is not None or eigenvectors is not None or eigenvalues is not None:
            raise ValueError('non-pending branch cannot consume eigensolver operands')
        return AwFloorResult(state,False,None,False)
    if not isinstance(solver_success,bool): raise TypeError('pending branch requires literal eigensolver outcome')
    if not solver_success:
        if eigenvectors is not None or eigenvalues is not None: raise ValueError('failed eigensolver branch has no accepted eigensystem')
        return AwFloorResult(state,False,False,False)
    U=M.mat(eigenvectors,3,3); lam=P.vec(eigenvalues,3)
    if M.mm(M.transpose(U),U) != M.eye(3): raise ValueError('orthonormal self-adjoint eigensystem required')
    cov=M.mat(state.covariance,21,21)
    Paw=[row[OFF_AW:OFF_AW+3] for row in cov[OFF_AW:OFF_AW+3]]
    Paw=M.scaled(M.plus(Paw,M.transpose(Paw)),F(1,2))
    delta=M.scaled(M.plus(M.plus(target,Paw,-1),M.transpose(M.plus(target,Paw,-1))),F(1,2))
    if M.mm(delta,U) != M.mm(U,_diag(lam)): raise ValueError('eigensystem is detached from SAME floor Delta')
    clipped=[max(F(0),x) for x in lam]
    positive=M.mm(M.mm(U,_diag(clipped)),M.transpose(U)); positive=M.scaled(M.plus(positive,M.transpose(positive)),F(1,2))
    out=[row[:] for row in cov]
    for i in range(3):
        for j in range(3): out[OFF_AW+i][OFF_AW+j] += positive[i][j]
    out=M.scaled(M.plus(out,M.transpose(out)),F(1,2))
    return AwFloorResult(replace(state,covariance=tuple(map(tuple,out))),False,True,any(x>0 for x in lam))


def symmetry_hygiene(state):
    cov=M.mat(state.covariance,21,21); sym=M.scaled(M.plus(cov,M.transpose(cov)),F(1,2))
    return replace(state,covariance=tuple(map(tuple,sym)))


@dataclass(frozen=True)
class Scheduler:
    period: F
    elapsed: F
    tolerance: F = F(0)
    def __post_init__(self):
        p,e,t=map(P.rational,(self.period,self.elapsed,self.tolerance))
        if p <= 0 or e < 0 or e >= p or t < 0: raise ValueError('valid period, service credit and nonnegative tolerance required')
        object.__setattr__(self,'period',p); object.__setattr__(self,'elapsed',e); object.__setattr__(self,'tolerance',t)
    def step(self,h):
        h=P.rational(h)
        if h <= 0: raise ValueError('positive prediction duration required')
        total=self.elapsed+h
        if total+self.tolerance < self.period: return False, Scheduler(self.period,total,self.tolerance)
        if total >= self.period:
            n=total//self.period; rem=total-n*self.period
        else: rem=F(0)
        if rem < 0 or rem >= self.period: raise AssertionError('exact scheduler remainder outside [0,period)')
        return True, Scheduler(self.period,rem,self.tolerance)


@dataclass(frozen=True)
class PostPrediction:
    state: CORE.State
    floor: AwFloorResult
    scheduler: Scheduler
    S_service_due: bool


def post_prediction_prefix(state, *, h, pending_aw_floor, aw_floor_target,
                           scheduler: Scheduler, floor_solver_success=None,
                           floor_eigenvectors=None, floor_eigenvalues=None):
    floor=aw_floor(state,pending=pending_aw_floor,target=aw_floor_target,
                   solver_success=floor_solver_success,eigenvectors=floor_eigenvectors,
                   eigenvalues=floor_eigenvalues)
    clean=symmetry_hygiene(floor.state); due,next_sched=scheduler.step(h)
    return PostPrediction(clean,replace(floor,state=clean),next_sched,due)


@dataclass(frozen=True)
class ServicedPostPrediction:
    state: CORE.State
    prefix: PostPrediction
    measurement: object | None


def service_S_if_due(prefix: PostPrediction, *, R_S, ldlt=None, alpha=1,
                     radius=F(2,5)):
    """Compose scheduler decision to the literal S=0 safe-LDLT branch.

    Not-due consumes no factorization witness and is an identity suffix. A due
    event requires the declared SafeLDLT branch and uses the SAME physical S
    reference already carried by finite_core. Applied R_S is still an explicit
    runtime operand whose tuner ancestry remains open.
    """
    if not isinstance(prefix,PostPrediction): raise TypeError('post-prediction prefix required')
    from tools.stability.ou3_alt_contraction import finite_measurement_runtime as MR
    if not prefix.S_service_due:
        if ldlt is not None: raise ValueError('not-due S branch cannot consume an LDLT witness')
        return ServicedPostPrediction(prefix.state,prefix,None)
    if not isinstance(ldlt,MR.SafeLDLT): raise TypeError('due S branch requires shipping safe-LDLT witness')
    meas=MR.measurement(prefix.state,'S_zero',ldlt=ldlt,R=R_S,alpha=alpha,radius=radius)
    return ServicedPostPrediction(meas.state,prefix,meas)


def readiness():
    return {
      'pending_aw_floor_identity_branch':True,
      'pending_aw_floor_solver_failure_branch':True,
      'pending_aw_floor_positive_part_success_branch':True,
      'floor_only_changes_aw_covariance_block':True,
      'scheduler_elapsed_credit_recurrence':True,
      'scheduler_due_not_due_branches':True,
      'S_service_ldlt_branch_attached':True,
      'S_service_rejection_preserves_post_prediction_state':True,
      'floor_target_same_tuner_source_attached':False,
      'eigensolver_finite_precision_attached':False,
      'scheduler_binary_tolerance_attached':False,
      'applied_R_S_tuner_source_attached':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
