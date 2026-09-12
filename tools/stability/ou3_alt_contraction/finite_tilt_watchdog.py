"""Finite Live tilt-watchdog control and hard-reset relation for ALT.

SeaStateFusionFilter_OU_III evaluates the watchdog after prediction and the
accelerometer correction.  Its control state is

  * tilt_over_limit_sec,
  * tilt_reset_cooldown_sec,

with 70 deg threshold, 0.35 s hold, 3 s cooldown, and 2*dt recovery decay.
When it fires in Live it calls ``initialize_from_acc_preserve_yaw(acc_in)``.

This module proves the exact timer/control recurrence and all non-transcendental
hard-reset state/covariance effects.  The final preserve-yaw nominal quaternion
is supplied by ``PreserveYawWitness`` and checked against the resulting finite
error state, but its atan2/asin/AngleAxis binary32 ancestry is deliberately an
open obligation.  We do NOT replace it with an arbitrary attitude reset.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_core as CORE
from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_prediction_graph as P
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as SENSOR


def R(x): return P.rational(x)
THRESHOLD_DEG=F(70)
HOLD_SEC=F(35,100)
COOLDOWN_SEC=F(3)

@dataclass(frozen=True)
class State:
    over_limit:F=F(0)
    cooldown:F=F(0)
    def __post_init__(self):
        a,b=R(self.over_limit),R(self.cooldown)
        if a<0 or b<0: raise ValueError('nonnegative watchdog timers required')
        object.__setattr__(self,'over_limit',a); object.__setattr__(self,'cooldown',b)

@dataclass(frozen=True)
class StepResult:
    state:State
    fired:bool


def step(state:State,*,dt,tilt_deg):
    if not isinstance(state,State): raise TypeError('tilt watchdog state required')
    dt,tilt=R(dt),R(tilt_deg)
    if dt<=0 or tilt<0: raise ValueError('positive dt and nonnegative tilt required')
    cooldown=max(F(0),state.cooldown-dt) if state.cooldown>0 else F(0)
    if tilt>THRESHOLD_DEG:
        over=state.over_limit+dt
    else:
        over=max(F(0),state.over_limit-2*dt)
    fire=(over>=HOLD_SEC and cooldown<=0)
    if fire:
        return StepResult(State(F(0),COOLDOWN_SEC),True)
    return StepResult(State(over,cooldown),False)

@dataclass(frozen=True)
class PreserveYawWitness:
    """Named output of shipping initialize_from_acc_preserve_yaw binary32 path.

    ``q_new_hat`` is the final nominal world->internal-body quaternion used by
    the finite core. ``down_body_unit`` is the normalized world-down axis in
    that final internal body frame used to construct accel-only covariance.
    The native atan2/asin/AngleAxis correspondence that certifies these values
    remains a separate deployment obligation.
    """
    q_new_hat:tuple
    down_body_unit:tuple
    def __post_init__(self):
        q=CORE.quaternion(self.q_new_hat); u=tuple(P.vec(self.down_body_unit,3))
        if M.dot(u,u)!=1: raise ValueError('unit world-down covariance axis witness required')
        object.__setattr__(self,'q_new_hat',q); object.__setattr__(self,'down_body_unit',u)


def _attitude_covariance(u,tilt_sigma,yaw_sigma):
    u=P.vec(u,3); ts,ys=R(tilt_sigma),R(yaw_sigma)
    if ts<=0 or ys<=0: raise ValueError('positive accel-reset covariance sigmas required')
    yy=[[u[i]*u[j] for j in range(3)] for i in range(3)]
    tt=M.plus(M.eye(3),yy,-1)
    A=M.plus(M.scaled(tt,ts*ts),M.scaled(yy,ys*ys))
    return M.scaled(M.plus(A,M.transpose(A)),F(1,2))


def preserve_yaw_reset(state:CORE.State,sample:SENSOR.GuardedImuSample,witness:PreserveYawWitness,
                       *,tilt_sigma=F(35,1000),yaw_sigma=F(15708,10000)):
    """Apply shipping hard-reset effects after a certified preserve-yaw attitude.

    The guarded accelerometer is required because shipping passes the exact
    ``acc_in`` already consumed by the current sample.  Its numerical role in
    producing ``q_new_hat`` is intentionally reserved for the native
    transcendental correspondence layer; it is not ignored or replaced by a
    second sensor operand.
    """
    if not isinstance(state,CORE.State) or not isinstance(sample,SENSOR.GuardedImuSample) or not isinstance(witness,PreserveYawWitness):
        raise TypeError('finite state, same guarded sample and preserve-yaw witness required')
    if sample.physical.history_id != state.reference.history_id:
        raise ValueError('tilt reset sample detached from same physical history')
    # initialize_from_acc works on the current acc_in; require a nondegenerate
    # conditioned packet before a reset witness may be used.
    if M.dot(sample.conditioned_accel_body,sample.conditioned_accel_body) <= F(1,10**16):
        raise ValueError('degenerate accelerometer cannot drive shipping tilt reset')

    q=witness.q_new_hat
    # Recompute the finite true-minus-nominal attitude coordinate at the SAME
    # physical endpoint; hard reset changes no physical truth.
    z=list(state.z)
    z[:3]=CORE.cayley(P.quat_mul(state.reference.q_world_to_body,P.quat_conj(q)))

    cov=[list(r) for r in state.covariance]
    Patt=_attitude_covariance(witness.down_body_unit,tilt_sigma,yaw_sigma)
    for i in range(3):
        for j in range(21):
            cov[i][j]=F(0); cov[j][i]=F(0)
    for i in range(3):
        for j in range(3): cov[i][j]=Patt[i][j]
    # This is exactly initialize_from_acc's covariance effect: all non-attitude
    # marginals/cross-covariances among them survive; only attitude cross terms
    # are discarded. set_quaternion_boat then changes qref and zeros only the
    # filter's internal attitude-error bookkeeping, already represented here by
    # recomputing finite true-minus-nominal z[:3].
    return CORE.State(state.mode,tuple(z),tuple(tuple(r) for r in cov),q,state.reference)


def readiness():
    return {
      'tilt_threshold_strict_gt_70deg':True,
      'hold_035s_and_recovery_2dt_materialized':True,
      'cooldown_3s_and_fire_reset_materialized':True,
      'same_guarded_accel_required_for_reset':True,
      'hard_reset_preserves_nonattitude_nominal_states':True,
      'hard_reset_drops_all_attitude_cross_covariances':True,
      'accel_only_tilt_yaw_covariance_reseed_materialized':True,
      'preserve_yaw_atan_asin_angleaxis_binary32_attached':False,
      'watchdog_tilt_acos_binary32_attached':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
