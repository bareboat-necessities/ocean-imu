"""Compose quantified admitted COMPLETE-BRMM and BIAS histories with Live.

The lower source-owning product checks necessary finite constraints. This
stronger theorem-facing layer carries one admitted primary wave history AND one
admitted BIAS0/1/2 history. IMU transition k must be the kth restriction datum
of both histories on the same physical segment.

``from_startup`` is the source-theorem bridge: it invokes the actual finite
startup -> interleaved Live constructor, then requires that exact fresh H18
physical reference to equal the explicit t_L restriction of the quantified
admitted COMPLETE-BRMM history and the quantified BIAS history before creating
the source-owned word. No synthetic Live root is inserted.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_admitted_brmm_restriction as ADMIT
from tools.stability.ou3_alt_contraction import finite_admitted_bias_history as ABIASS
from tools.stability.ou3_alt_contraction import finite_source_bound_live_word as WORD
from tools.stability.ou3_alt_contraction import finite_source_continuation as SOURCE
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as SENSOR
from tools.stability.ou3_alt_contraction import finite_live_interleave as INTER
from tools.stability.ou3_alt_contraction import finite_startup_live_runtime_bridge as START
from tools.stability.ou3_alt_contraction import finite_live_magnetic_word as MAGLIVE


@dataclass(frozen=True)
class State:
    live_word: WORD.State
    admitted_history: ADMIT.AdmittedHistory
    bias_history: ABIASS.AdmittedBiasHistory
    def __post_init__(self):
        if not isinstance(self.live_word,WORD.State) or not isinstance(self.admitted_history,ADMIT.AdmittedHistory):
            raise TypeError('source-owning Live word and admitted physical history required')
        if not isinstance(self.bias_history,ABIASS.AdmittedBiasHistory):
            raise TypeError('admitted BIAS history required')
        root=self.live_word.source.root
        if root.history_id != self.admitted_history.history_id:
            raise ValueError('Live source root detached from quantified admitted history')
        if root.bias_family != self.bias_history.family:
            raise ValueError('Live source root detached from quantified BIAS family')
        if self.live_word.bias_history_id != self.bias_history.history_id:
            raise ValueError('Live physical bias root detached from quantified BIAS history')


@dataclass(frozen=True)
class Result:
    state: State
    event: object


def from_live(live, history:ADMIT.AdmittedHistory, bias_history:ABIASS.AdmittedBiasHistory, *,
              gyro_residual_history_id:str, accel_residual_history_id:str,
              runtime:WORD.RuntimeConfig):
    """Component constructor; does not by itself prove startup provenance."""
    if not isinstance(history,ADMIT.AdmittedHistory): raise TypeError('AdmittedHistory required')
    if not isinstance(bias_history,ABIASS.AdmittedBiasHistory): raise TypeError('AdmittedBiasHistory required')
    ref=live.live.live.mekf.reference
    root=ADMIT.source_root(history,live_origin=ref.live_origin,bias_family=bias_history.family)
    sensors=SOURCE.SensorDisturbanceRoot(root,gyro_residual_history_id,accel_residual_history_id)
    lower=WORD.from_live(live,root,sensors,runtime)
    return State(lower,history,bias_history)


def from_startup(bridge:START.Result, magnetic:MAGLIVE.StartupState,
                 origin:ADMIT.RestrictedOrigin,
                 bias_history:ABIASS.AdmittedBiasHistory, *,
                 gyro_residual_history_id:str, accel_residual_history_id:str,
                 runtime:WORD.RuntimeConfig, proxy_q_norm, proxy_yaw_half,
                 schedule=None):
    """Bind the actual startup-produced fresh interleaved state to admitted sources."""
    if not isinstance(bridge,START.Result) or not isinstance(magnetic,MAGLIVE.StartupState):
        raise TypeError('actual startup Live bridge and magnetic startup history required')
    if not isinstance(origin,ADMIT.RestrictedOrigin) or not isinstance(bias_history,ABIASS.AdmittedBiasHistory):
        raise TypeError('admitted BRMM origin and BIAS history required')
    live=INTER.from_startup(bridge,magnetic,proxy_q_norm=proxy_q_norm,
                            proxy_yaw_half=proxy_yaw_half,schedule=schedule)
    ref=live.live.live.mekf.reference
    if origin.endpoint != ref:
        raise ValueError('startup-produced fresh H18 reference is not admitted-history t_L origin')
    if ref.bias_root != bias_history.history_id or ref.bias_family != bias_history.family:
        raise ValueError('startup-produced fresh H18 reference detached from admitted BIAS history')
    history=origin.history
    root=ADMIT.source_root(history,live_origin=ref.live_origin,bias_family=bias_history.family)
    ADMIT.qualify_origin(root,origin)
    sensors=SOURCE.SensorDisturbanceRoot(root,gyro_residual_history_id,accel_residual_history_id)
    lower=WORD.from_live(live,root,sensors,runtime)
    return State(lower,history,bias_history)


def imu_step(state:State, *, restricted:ADMIT.RestrictedSegment,
             bias_restricted:ABIASS.RestrictedBiasStep,
             witness:SOURCE.StepWitness, raw:SENSOR.RawImuSample, packet_id:str,
             **dynamic):
    if not isinstance(state,State): raise TypeError('admitted-history Live state required')
    if not isinstance(restricted,ADMIT.RestrictedSegment):
        raise TypeError('kth admitted-history restriction required')
    if not isinstance(bias_restricted,ABIASS.RestrictedBiasStep):
        raise TypeError('kth admitted-BIAS restriction required')
    if restricted.history != state.admitted_history:
        raise ValueError('IMU restriction detached from carried admitted history')
    if bias_restricted.history != state.bias_history:
        raise ValueError('BIAS restriction detached from carried admitted BIAS history')
    if restricted.ordinal != bias_restricted.ordinal or restricted.ordinal != witness.ordinal:
        raise ValueError('BRMM/BIAS/source restriction ordinals differ')
    if restricted.segment != bias_restricted.segment:
        raise ValueError('BRMM and BIAS restrictions must describe the same physical segment')
    ADMIT.qualify_step(state.live_word.source.root,restricted,witness)
    out=WORD.imu_step(state.live_word,witness=witness,segment=restricted.segment,
                      raw=raw,packet_id=packet_id,**dynamic)
    return Result(State(out.state,state.admitted_history,state.bias_history),out.event)


def mag_step(state:State, **kwargs):
    if not isinstance(state,State): raise TypeError('admitted-history Live state required')
    out=WORD.mag_step(state.live_word,**kwargs)
    return Result(State(out.state,state.admitted_history,state.bias_history),out.event)


def set_hold(state:State, *, hold):
    if not isinstance(state,State): raise TypeError('admitted-history Live state required')
    out=WORD.set_hold(state.live_word,hold=hold)
    return Result(State(out.state,state.admitted_history,state.bias_history),out.event)


def readiness():
    a=ADMIT.readiness(); b=ABIASS.readiness(); lower=WORD.readiness()
    return {
      'universal_admitted_COMPLETE_BRMM_history_carried_in_Live_product':a['universal_theorem_quantifier_over_admitted_primary_history_explicit'],
      'admitted_BIAS_history_carried_in_same_Live_product':b['BIAS0_BIAS1_BIAS2_admitted_history_quantifier_available'],
      'actual_startup_interleave_constructor_consumed_by_admitted_source_factory':True,
      'startup_fresh_H18_reference_must_equal_admitted_tL_origin':True,
      'startup_fresh_H18_reference_must_equal_admitted_BIAS_root':True,
      'synthetic_Live_root_inserted_by_source_bridge':False,
      'every_IMU_event_requires_same_BRMM_and_BIAS_restriction_ordinal':True,
      'BRMM_and_BIAS_restrictions_share_exact_physical_segment':True,
      'finite_constraints_derived_after_admitted_history_restriction':True,
      'finite_tokens_never_used_as_membership_oracle':a['arbitrary_runtime_tokens_do_not_prove_COMPLETE_BRMM_membership'],
      'magnetic_and_hold_events_preserve_admitted_histories':True,
      'startup_source_equality_bridge_available':True,
      'BIAS_generating_history_attached':True,
      'sensor_disturbance_admissibility_attached':False,
      'estimator_coefficients_same_history_closed':lower['finite_estimator_coefficients_bound_to_same_source_continuation'],
      'complete_600_step_shipping_word_composed_from_restrictions':False,
      'startup_deployment_arithmetic_closed':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
      'ALT_STARTUP_PASS':False,
      'ALT_END_TO_END_PASS':False,
    }
