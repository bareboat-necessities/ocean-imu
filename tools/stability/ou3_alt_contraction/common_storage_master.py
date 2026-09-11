#!/usr/bin/env python3
"""Exact reduction used by the first ALT common-joint24 storage attempt.

This is theorem machinery, not a replay metric fit.  Phase 1 has already built
one physical joint24 cocycle over the COMPLETE-BRMM Normal-Live relation.  The
remaining structural issue is that H18 contains bounded, deliberately held
accelerometer-bias coordinates and the physical true-bias reference is itself a
bounded source state.  Requiring strict homogeneous contraction in those six
coordinates is therefore stronger than the target theorem.

For a word map x+ = A x + d and coercive storage M>0, ALT seeks

    x+' M x+ <= rho x'Mx + (C x)' Beta (C x) + source_supply,

where C selects only independently bounded neutral/reference coordinates.  For
fixed A,M,rho define Q=A'MA-rho M.  The coordinate Finsler/Schur lemma gives

    exists Beta>=0 : Q-C'Beta C < 0
       iff
    Z'QZ < 0,

where columns of Z span ker(C).  For the shipping joint24 state the admitted
neutral selector is exactly [e_ba(3), beta_true(3)] = coordinates 18:24.  No
motion state is moved into supply.

The lemma is valuable because it turns the first common-storage search into a
strict contraction problem on the 18-dimensional unsupplied motion subspace,
while retaining one *coercive 24-state* storage and finite bounded-bias supply.
It does not prove that a common M exists over the full source relation; that is
the next numerical/outward obligation.
"""
from __future__ import annotations

import math
import numpy as np

from tools.stability.ou3_alt_contraction import phase1_closure as PHASE1

QUALIFICATION = "OU3_ALT_COMMON_JOINT24_NEUTRAL_SUPPLY_REDUCTION_V1"
JOINT_DIM = 24
MOTION_DIM = 18
NEUTRAL_DIM = 6
NEUTRAL_INDICES = tuple(range(18, 24))


def _symmetric(a, name: str) -> np.ndarray:
    a=np.asarray(a,dtype=float)
    if a.ndim!=2 or a.shape[0]!=a.shape[1] or not np.isfinite(a).all():
        raise ValueError(name+" must be finite square")
    if not np.allclose(a,a.T,rtol=0,atol=1e-12):
        raise ValueError(name+" must be symmetric")
    return (a+a.T)/2


def coordinate_selector(n: int, supplied_indices) -> np.ndarray:
    idx=tuple(int(i) for i in supplied_indices)
    if n<=0 or len(set(idx))!=len(idx) or any(i<0 or i>=n for i in idx):
        raise ValueError("invalid supplied coordinate set")
    C=np.zeros((len(idx),n))
    for r,i in enumerate(idx): C[r,i]=1.0
    return C


def kernel_injection(n: int, supplied_indices) -> np.ndarray:
    supplied=set(int(i) for i in supplied_indices)
    keep=[i for i in range(n) if i not in supplied]
    Z=np.zeros((n,len(keep)))
    for j,i in enumerate(keep): Z[i,j]=1.0
    return Z


def joint24_neutral_selector() -> np.ndarray:
    return coordinate_selector(JOINT_DIM,NEUTRAL_INDICES)


def joint24_motion_injection() -> np.ndarray:
    return kernel_injection(JOINT_DIM,NEUTRAL_INDICES)


def unsupplied_storage_form(A, M, rho: float, supplied_indices=NEUTRAL_INDICES) -> np.ndarray:
    """Return Z'(A'MA-rho M)Z, the exact Finsler feasibility restriction."""
    A=np.asarray(A,dtype=float); M=_symmetric(M,"M")
    n=M.shape[0]
    if A.shape!=(n,n) or not np.isfinite(A).all(): raise ValueError("A shape/finite mismatch")
    if not math.isfinite(rho) or not 0.0<rho<1.0: raise ValueError("rho must be in (0,1)")
    np.linalg.cholesky(M)
    Z=kernel_injection(n,supplied_indices)
    Q=A.T@M@A-rho*M
    R=Z.T@Q@Z
    return (R+R.T)/2


def coordinate_supply_completion(A, M, rho: float, supplied_indices=NEUTRAL_INDICES, margin=1e-10):
    """Construct a finite scalar Beta=b I when the Finsler restriction is strict.

    This is a point-coefficient candidate constructor.  Rigorous source-uniform
    promotion still requires outward eigen/LDLT bounds over the universal word
    family.  The formula is the Schur complement after permuting unsupplied
    coordinates first:

      Qxx<0,
      b > lambda_max(Qcc-Qcx Qxx^{-1} Qxc).
    """
    if not math.isfinite(margin) or margin<=0: raise ValueError("positive margin required")
    A=np.asarray(A,dtype=float); M=_symmetric(M,"M"); n=M.shape[0]
    if A.shape!=(n,n): raise ValueError("A shape mismatch")
    idx=tuple(int(i) for i in supplied_indices); supplied=set(idx)
    keep=[i for i in range(n) if i not in supplied]
    Q=(A.T@M@A-rho*M); Q=(Q+Q.T)/2
    Qxx=Q[np.ix_(keep,keep)]
    eig=np.linalg.eigvalsh(Qxx)
    if eig[-1]>=-margin:
        return {"feasible":False,"projected_lambda_max":float(eig[-1]),"beta_scalar":None,"completed_lambda_max":None}
    Qxc=Q[np.ix_(keep,idx)]; Qcc=Q[np.ix_(idx,idx)]
    # Qxx is strictly negative definite here, hence nonsingular.
    schur=Qcc-Qxc.T@np.linalg.solve(Qxx,Qxc)
    lam=float(np.linalg.eigvalsh((schur+schur.T)/2)[-1]) if len(idx) else -math.inf
    beta=max(0.0,lam)+margin
    D=np.zeros((n,n))
    for i in idx: D[i,i]=beta
    completed=(Q-D+Q.T-D.T)/2
    lmax=float(np.linalg.eigvalsh(completed)[-1])
    return {"feasible":bool(lmax<0.0),"projected_lambda_max":float(eig[-1]),"beta_scalar":float(beta),"completed_lambda_max":lmax}


def bounded_additive_supply_multiplier(rho0: float, rho: float) -> float:
    """Young bound for x+=Ax+d once ||A x||_M^2 <= rho0 ||x||_M^2.

    Choosing eta=rho/rho0-1 gives
      ||Ax+d||_M^2 <= rho ||x||_M^2 + rho/(rho-rho0) ||d||_M^2.
    Thus physical-source magnitude affects the practical bound, not existence
    of the homogeneous contraction margin.  rho0=0 has multiplier 1.
    """
    if not (math.isfinite(rho0) and math.isfinite(rho)) or rho0<0 or not rho0<rho<1:
        raise ValueError("require 0 <= rho0 < rho < 1")
    return 1.0 if rho0==0 else rho/(rho-rho0)


def build() -> dict:
    phase1=PHASE1.build(); failures=PHASE1.validate(phase1)
    if failures: raise RuntimeError("Phase-1 physical master not closed: "+repr(failures))
    C=joint24_neutral_selector(); Z=joint24_motion_injection()
    return {
        "qualification":QUALIFICATION,
        "canonical_source":"COMPLETE_BRMM_NORMAL_LIVE_WORD",
        "phase1_storage_search_allowed_consumed":bool(phase1["storage_search_allowed"]),
        "joint_dimension":JOINT_DIM,"unsupplied_motion_dimension":MOTION_DIM,"bounded_neutral_dimension":NEUTRAL_DIM,
        "bounded_neutral_coordinates":"e_ba[3], beta_true[3]",
        "neutral_indices":list(NEUTRAL_INDICES),
        "selector_rank":int(np.linalg.matrix_rank(C)),"kernel_dimension":int(Z.shape[1]),
        "selector_times_kernel_zero":bool(np.array_equal(C@Z,np.zeros((NEUTRAL_DIM,MOTION_DIM)))),
        "full_joint24_storage_remains_coercive_required":True,
        "motion_coordinates_moved_to_supply":False,
        "held_bias_and_true_bias_may_enter_only_bounded_supply":True,
        "finsler_equivalence":"exists Beta>=0: Q-C.T Beta C<0 iff Z.T Q Z<0",
        "finite_scalar_supply_completion_formula_available":True,
        "physical_source_supply_magnitude_needed_for_metric_feasibility":False,
        "physical_source_supply_magnitude_needed_for_ultimate_bound":True,
        "common_M_source_uniform_search_closed":False,
        "source_uniform_outward_projected_LDLT_closed":False,
        "ALT_LIVE_PASS":False,
        "next_obligation":"search one common coercive 24x24 M and rho for which Z.T(A_w.T M A_w-rho M)Z<0 for every H18/A21 physical word in the universal relation; then outward-certify that projected family before completing finite neutral/source supplies",
    }


def validate(d: dict) -> list[str]:
    f=[]
    if d.get("qualification")!=QUALIFICATION or d.get("canonical_source")!="COMPLETE_BRMM_NORMAL_LIVE_WORD": f.append("qualification/source mismatch")
    for k in ("phase1_storage_search_allowed_consumed","selector_times_kernel_zero","full_joint24_storage_remains_coercive_required","held_bias_and_true_bias_may_enter_only_bounded_supply","finite_scalar_supply_completion_formula_available","physical_source_supply_magnitude_needed_for_ultimate_bound"):
        if d.get(k) is not True: f.append(k+" not true")
    for k in ("motion_coordinates_moved_to_supply","physical_source_supply_magnitude_needed_for_metric_feasibility","common_M_source_uniform_search_closed","source_uniform_outward_projected_LDLT_closed","ALT_LIVE_PASS"):
        if d.get(k) is not False: f.append(k+" not false")
    if d.get("joint_dimension")!=24 or d.get("unsupplied_motion_dimension")!=18 or d.get("bounded_neutral_dimension")!=6: f.append("dimension mismatch")
    if d.get("neutral_indices")!=list(NEUTRAL_INDICES) or d.get("selector_rank")!=6 or d.get("kernel_dimension")!=18: f.append("selector/kernel mismatch")
    return list(dict.fromkeys(f))
