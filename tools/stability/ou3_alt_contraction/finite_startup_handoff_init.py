"""Exact real-arithmetic core of shipping ``initialize_from_attitude`` handoff.

This module records the literal covariance/state mutation after shipping has
accepted the handoff quaternion. ALT explicitly excludes the optional wind-heel
retarget branch, so the certified deployment has B'=B and no heel transform is
present in this handoff relation.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_startup_handoff_seed as SEED
from tools.stability.ou3_alt_contraction import finite_mag_tilt_frame as TILT
from tools.stability.ou3_alt_contraction import deployment_scope as SCOPE

N=21

def R(x): return M.rational(x)
def _mat(A,n=N):
    if len(A)!=n or any(len(r)!=n for r in A): raise ValueError(f'{n}x{n} covariance required')
    return [[R(x) for x in r] for r in A]
def _vec(x,n=N): return list(M.vec(x,n))

@dataclass(frozen=True)
class UnitDownWitness:
    vector:tuple
    def __post_init__(self):
        v=tuple(M.vec(self.vector,3))
        if M.dot(v,v)!=1: raise ValueError('exact unit world-down/body witness required')
        object.__setattr__(self,'vector',v)

@dataclass(frozen=True)
class Result:
    x:tuple
    P:tuple
    q_seed:tuple
    tilt_sigma:F
    yaw_sigma:F
    gauged:bool

def down_from_unit_boat_quaternion_zero_heel(q_bw,*,scope:SCOPE.Scope):
    SCOPE.assert_certified_scope(scope)
    q=tuple(M.vec(q_bw,4))
    if M.dot(q,q)!=1: raise ValueError('exact unit boat quaternion required before zero-heel down-axis attachment')
    w,x,y,z=q
    return UnitDownWitness((2*(x*z-w*y),2*(y*z+w*x),1-2*(x*x+y*y)))

def attitude_covariance(down:UnitDownWitness,tilt_sigma,yaw_sigma):
    if not isinstance(down,UnitDownWitness): raise TypeError('unit down-axis witness required')
    st,sy=R(tilt_sigma),R(yaw_sigma)
    if st<=0 or sy<=0: raise ValueError('positive handoff attitude sigmas required')
    u=down.vector; tv=st*st; yv=sy*sy
    return [[tv*((1 if i==j else 0)-u[i]*u[j])+yv*u[i]*u[j] for j in range(3)] for i in range(3)]

def initialize(seed:SEED.Result,x_before,P_before,*,down_body:UnitDownWitness):
    if not isinstance(seed,SEED.Result): raise TypeError('startup handoff seed/control result required')
    if seed.allow_acc_bias: raise ValueError('shipping startup handoff must keep accelerometer-bias unlock closed')
    x=_vec(x_before); P=_mat(P_before); Patt=attitude_covariance(down_body,seed.tilt_sigma,seed.yaw_sigma)
    x[:3]=[F(0),F(0),F(0)]
    for i in range(3):
        for j in range(3): P[i][j]=Patt[i][j]
    for i in range(3):
        for j in range(3,6): P[i][j]=F(0); P[j][i]=F(0)
    return Result(tuple(x),tuple(tuple(r) for r in P),tuple(seed.q_seed),R(seed.tilt_sigma),R(seed.yaw_sigma),seed.gauged)

def initialize_from_seed_zero_heel(seed:SEED.Result,x_before,P_before,*,scope:SCOPE.Scope,q_norm:TILT.SqrtWitness|None=None):
    """Compose gauged or timeout/ungauged handoff through shipping normalization.

    Gauged startup already supplies an exact unit quaternion. The ungauged
    timeout path supplies the persistent proxy quaternion, so its norm is tied
    explicitly to that same predecessor. This is exact-real normalization;
    binary32/Eigen normalization correspondence remains a deployment obligation.
    """
    if not isinstance(seed,SEED.Result): raise TypeError('startup handoff seed required')
    q=tuple(M.vec(seed.q_seed,4)); q2=M.dot(q,q)
    if seed.gauged:
        if q_norm is not None: raise ValueError('gauged unit handoff consumes no normalization witness')
        if q2!=1: raise ValueError('gauged handoff seed must be exactly unit')
        qn=q
    else:
        if not isinstance(q_norm,TILT.SqrtWitness) or q_norm.radicand!=q2:
            raise TypeError('ungauged handoff requires norm witness for the same proxy quaternion')
        if q_norm.value<=F(1,10**6): raise ValueError('shipping handoff rejects tiny/degenerate proxy quaternion')
        qn=tuple(x/q_norm.value for x in q)
        if M.dot(qn,qn)!=1: raise AssertionError('ungauged normalized handoff quaternion not unit')
    normalized=SEED.Result(qn,seed.tilt_sigma,seed.yaw_sigma,seed.allow_acc_bias,seed.gauged)
    down=down_from_unit_boat_quaternion_zero_heel(qn,scope=scope)
    return initialize(normalized,x_before,P_before,down_body=down)

def initialize_from_gauged_seed_zero_heel(seed:SEED.Result,x_before,P_before,*,scope:SCOPE.Scope):
    if not isinstance(seed,SEED.Result) or not seed.gauged: raise ValueError('gauged startup seed required for exact unit handoff composition')
    return initialize_from_seed_zero_heel(seed,x_before,P_before,scope=scope)

def readiness():
    return {
      'attitude_error_zero_reset_materialized':True,
      'anisotropic_tilt_yaw_covariance_formula_materialized':True,
      'attitude_gyro_bias_cross_covariance_zero_materialized':True,
      'all_other_mean_and_covariance_entries_retained':True,
      'wind_heel_branch_excluded_by_declared_ALT_scope':True,
      'zero_heel_body_equals_body_prime':True,
      'gauged_unit_seed_to_down_axis_attached_exactly':True,
      'ungauged_proxy_seed_normalization_and_down_axis_attached_exactly':True,
      'both_gauged_and_ungauged_handoff_covariance_branches_materialized':True,
      'Eigen_normalization_and_binary32_handoff_attached':False,
      'goLive_enterLive_stage_mutation_attached':False,
      'startup_capture_closed':False,
      'ALT_STARTUP_PASS':False,
    }
