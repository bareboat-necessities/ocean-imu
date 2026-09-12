"""Exact real-arithmetic core of shipping ``initialize_from_attitude`` handoff.

This module does NOT invent a generic startup reset.  It records the literal
covariance/state mutation after shipping has accepted the handoff quaternion:

* attitude error state x[0:3] is zeroed by ``set_quaternion_boat``;
* P[0:3,0:3] is replaced by the anisotropic tilt/yaw covariance about the
  world-down axis expressed in the MEKF's internal B' frame;
* attitude <-> gyro-bias cross covariance P[0:3,3:6] and transpose are zeroed;
* every other mean component and every other covariance entry is retained.

For unit down axis u,

    P_att = sigma_tilt^2 (I - u u') + sigma_yaw^2 u u'.

The source/runtime correspondence that derives ``u`` from the accepted boat
quaternion, wind-heel conversion, Eigen normalization and binary32 arithmetic is
kept separate and fail-closed.  This file proves the finite reset algebra only.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_startup_handoff_seed as SEED

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
        if M.dot(v,v)!=1: raise ValueError('exact unit world-down/body-prime witness required')
        object.__setattr__(self,'vector',v)

@dataclass(frozen=True)
class Result:
    x:tuple
    P:tuple
    q_seed:tuple
    tilt_sigma:F
    yaw_sigma:F
    gauged:bool


def attitude_covariance(down:UnitDownWitness,tilt_sigma,yaw_sigma):
    if not isinstance(down,UnitDownWitness): raise TypeError('unit down-axis witness required')
    st,sy=R(tilt_sigma),R(yaw_sigma)
    if st<=0 or sy<=0: raise ValueError('positive handoff attitude sigmas required')
    u=down.vector; tv=st*st; yv=sy*sy
    return [[tv*((1 if i==j else 0)-u[i]*u[j])+yv*u[i]*u[j] for j in range(3)] for i in range(3)]


def initialize(seed:SEED.Result,x_before,P_before,*,down_body_prime:UnitDownWitness):
    if not isinstance(seed,SEED.Result): raise TypeError('startup handoff seed/control result required')
    if seed.allow_acc_bias:
        raise ValueError('shipping startup handoff must keep accelerometer-bias unlock closed')
    x=_vec(x_before); P=_mat(P_before)
    Patt=attitude_covariance(down_body_prime,seed.tilt_sigma,seed.yaw_sigma)

    # set_quaternion_boat: the local attitude-error state is reset, not the
    # physical/navigation/bias means.
    x[:3]=[F(0),F(0),F(0)]

    # set_accel_only_attitude_covariance_(tilt_sigma,yaw_sigma).
    for i in range(3):
        for j in range(3): P[i][j]=Patt[i][j]

    # initialize_from_attitude explicitly discards only attitude/gyro-bias
    # correlation after replacing the attitude reference.
    for i in range(3):
        for j in range(3,6):
            P[i][j]=F(0); P[j][i]=F(0)

    return Result(tuple(x),tuple(tuple(r) for r in P),tuple(seed.q_seed),
                  R(seed.tilt_sigma),R(seed.yaw_sigma),seed.gauged)


def readiness():
    return {
      'attitude_error_zero_reset_materialized':True,
      'anisotropic_tilt_yaw_covariance_formula_materialized':True,
      'attitude_gyro_bias_cross_covariance_zero_materialized':True,
      'all_other_mean_and_covariance_entries_retained':True,
      'body_prime_down_axis_from_boat_quaternion_and_wind_heel_attached':False,
      'Eigen_normalization_and_binary32_handoff_attached':False,
      'goLive_enterLive_stage_mutation_attached':False,
      'startup_capture_closed':False,
      'ALT_STARTUP_PASS':False,
    }
