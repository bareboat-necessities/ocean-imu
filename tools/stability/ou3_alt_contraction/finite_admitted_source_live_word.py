"""Compose the quantified admitted COMPLETE-BRMM history with the Live word.

The lower source-owning Live product checks necessary finite constraints.  This
stronger theorem-facing layer additionally carries one ``AdmittedHistory`` and
requires IMU transition k to be the kth restriction datum of that same history.
It therefore models the universal theorem in the correct direction:

    admitted history  -> finite restriction -> shipping event,

never ``finite runtime tokens -> admitted history``.

BIAS generating-history membership, sensor residual admissibility, deployment
arithmetic, and the complete 600-event composition are still open.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_admitted_brmm_restriction as ADMIT
from tools.stability.ou3_alt_contraction import finite_source_bound_live_word as WORD
from tools.stability.ou3_alt_contraction import finite_source_continuation as SOURCE
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as SENSOR


@dataclass(frozen=True)
class State:
    live_word: WORD.State
    admitted_history: ADMIT.AdmittedHistory
    def __post_init__(self):
        if not isinstance(self.live_word,WORD.State) or not isinstance(self.admitted_history,ADMIT.AdmittedHistory):
            raise TypeError('source-owning Live word and admitted physical history required')
        if self.live_word.source.root.history_id != self.admitted_history.history_id:
            raise ValueError('Live source root detached from quantified admitted history')


@dataclass(frozen=True)
class Result:
    state: State
    event: object


def from_live(live, history:ADMIT.AdmittedHistory, *, bias_family:str,
              gyro_residual_history_id:str, accel_residual_history_id:str,
              runtime:WORD.RuntimeConfig):
    if not isinstance(history,ADMIT.AdmittedHistory): raise TypeError('AdmittedHistory required')
    ref=live.live.live.mekf.reference
    root=ADMIT.source_root(history,live_origin=ref.live_origin,bias_family=bias_family)
    sensors=SOURCE.SensorDisturbanceRoot(root,gyro_residual_history_id,accel_residual_history_id)
    lower=WORD.from_live(live,root,sensors,runtime)
    return State(lower,history)


def imu_step(state:State, *, restricted:ADMIT.RestrictedSegment,
             witness:SOURCE.StepWitness, raw:SENSOR.RawImuSample, packet_id:str,
             **dynamic):
    if not isinstance(state,State): raise TypeError('admitted-history Live state required')
    if not isinstance(restricted,ADMIT.RestrictedSegment):
        raise TypeError('kth admitted-history restriction required')
    if restricted.history != state.admitted_history:
        raise ValueError('IMU restriction detached from carried admitted history')
    # This applies the universal restriction theorem and all necessary finite
    # constraints before the lower shipping word is allowed to execute.
    ADMIT.qualify_step(state.live_word.source.root,restricted,witness)
    out=WORD.imu_step(state.live_word,witness=witness,segment=restricted.segment,
                      raw=raw,packet_id=packet_id,**dynamic)
    return Result(State(out.state,state.admitted_history),out.event)


def mag_step(state:State, **kwargs):
    if not isinstance(state,State): raise TypeError('admitted-history Live state required')
    out=WORD.mag_step(state.live_word,**kwargs)
    return Result(State(out.state,state.admitted_history),out.event)


def set_hold(state:State, *, hold):
    if not isinstance(state,State): raise TypeError('admitted-history Live state required')
    out=WORD.set_hold(state.live_word,hold=hold)
    return Result(State(out.state,state.admitted_history),out.event)


def readiness():
    a=ADMIT.readiness(); lower=WORD.readiness()
    return {
      'universal_admitted_COMPLETE_BRMM_history_carried_in_Live_product':a['universal_theorem_quantifier_over_admitted_primary_history_explicit'],
      'every_IMU_event_requires_same_history_restriction_ordinal':True,
      'finite_constraints_derived_after_admitted_history_restriction':True,
      'finite_tokens_never_used_as_membership_oracle':a['arbitrary_runtime_tokens_do_not_prove_COMPLETE_BRMM_membership'],
      'magnetic_and_hold_events_preserve_admitted_history':True,
      'startup_sample_zero_equal_to_admitted_history_restriction_proved':False,
      'BIAS_generating_history_attached':False,
      'sensor_disturbance_admissibility_attached':False,
      'estimator_coefficients_same_history_closed':lower['finite_estimator_coefficients_bound_to_same_source_continuation'],
      'complete_600_step_shipping_word_composed_from_restrictions':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
      'ALT_STARTUP_PASS':False,
      'ALT_END_TO_END_PASS':False,
    }
