"""Exact real-arithmetic core of shipping ``initialize_from_attitude`` handoff.

This module records the literal covariance/state mutation after shipping has
accepted the handoff quaternion.  ALT explicitly excludes the optional wind-heel
retarget branch, so the certified deployment has B'=B and no heel transform is
present in this handoff relation.

* attitude error state x[0:3] is zeroed by ``set_quaternion_boat``;
* P[0:3,0:3] is replaced by the anisotropic tilt/yaw covariance about the
  world-down axis expressed directly in the physical body frame;
* attitude <-> gyro-bias cross covariance P[0:3,3:6] and transpose are zeroed;
* every other mean component and every other covariance entry is retained.

For unit down axis u,

    P_att = sigma_tilt^2 (I - u u') + sigma_yaw^2 u u'.

For a unit boat quaternion q=(w,x,y,z), B->W, zero wind heel gives

    u = R(q)^T e_z
      = (2(xz-wy), 2(yz+wx), 1-2(x^2+y^2)).

Deployment binary32 normalization remains separate and fail-closed.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_startup_handoff_seed as SEED
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
    """Exact world-down axis in body for unit q_bw under B'=B scope."""
    SCOPE.assert_certified_scope(scope)
    q=tuple(M.vec(q_bw,4))
    if M.dot(q,q)!=1:
        raise ValueError('exact unit boat quaternion required before zero-heel down-axis attachment')
    w,x,y,z=q
    u=(2*(x*z-w*y),2*(y*z+w*x),1-2*(x*x+y*y))
    return UnitDownWitness(u)


def attitude_covariance(down:UnitDownWitness,tilt_sigma,yaw_sigma):
    if not isinstance(down,UnitDownWitness): raise TypeError('unit down-axis witness required')
    st,sy=R(tilt_sigma),R(yaw_sigma)
    if st<=0 or sy<=0: raise ValueError('positive handoff attitude sigmas required')
    u=down.vector; tv=st*st; yv=sy*sy
    return [[tv*((1 if i==j else 0)-u[i]*u[j])+yv*u[i]*u[j] for j in range(3)] for i in range(3)]


def initialize(seed:SEED.Result,x_before,P_before,*,down_body:UnitDownWitness):
    if not isinstance(seed,SEED.Result): raise TypeError('startup handoff seed/control result required')
    if seed.allow_acc_bias:
        raise ValueError('shipping startup handoff must keep accelerometer-bias unlock closed')
    x=_vec(x_before); P=_mat(P_before)
    Patt=attitude_covariance(down_body,seed.tilt_sigma,seed.yaw_sigma)

    x[:3]=[F(0),F(0),F(0)]

    for i in range(3):
        for j in range(3): P[i][j]=Patt[i][j]

    for i in range(3):
        for j in range(3,6):
            P[i][j]=F(0); P[j][i]=F(0)

    return Result(tuple(x),tuple(tuple(r) for r in P),tuple(seed.q_seed),
                  R(seed.tilt_sigma),R(seed.yaw_sigma),seed.gauged)


def initialize_from_gauged_seed_zero_heel(seed:SEED.Result,x_before,P_before,*,scope:SCOPE.Scope):
    """Compose the exact gauged unit seed with the zero-heel covariance axis."""
    if not isinstance(seed,SEED.Result) or not seed.gauged:
        raise ValueError('gauged startup seed required for exact unit handoff composition')
    down=down_from_unit_boat_quaternion_zero_heel(seed.q_seed,scope=scope)
    return initialize(seed,x_before,P_before,down_body=down)


def readiness():
    return {
      'attitude_error_zero_reset_materialized':True,
      'anisotropic_tilt_yaw_covariance_formula_materialized':True,
      'attitude_gyro_bias_cross_covariance_zero_materialized':True,
      'all_other_mean_and_covariance_entries_retained':True,
      'wind_heel_branch_excluded_by_declared_ALT_scope':True,
      'zero_heel_body_equals_body_prime':True,
      'gauged_unit_seed_to_down_axis_attached_exactly':True,
      'Eigen_normalization_and_binary32_handoff_attached':False,
      'goLive_enterLive_stage_mutation_attached':False,
      'startup_capture_closed':False,
      'ALT_STARTUP_PASS':False,
    }
