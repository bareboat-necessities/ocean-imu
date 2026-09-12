"""Admitted-history magnetic edge for the ALT master word.

Magnetometer calls do not advance the 5 ms physical source. After an IMU step,
the current endpoint is already the last admitted restriction. At fresh sample
zero there is no predecessor transition, so the theorem-facing caller must
supply the explicit ``RestrictedOrigin`` of the same admitted COMPLETE-BRMM
history. The actual fresh MEKF physical reference must equal that origin exactly
before any magnetic state mutation is allowed.

Both paths retain the correlated magnetic forcing and preserve the admitted
BRMM and BIAS histories. Startup *reachability* to this fresh state and sensor
residual admissibility remain separate obligations.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_admitted_source_live_word as ADLIVE
from tools.stability.ou3_alt_contraction import finite_admitted_brmm_restriction as ABRMM
from tools.stability.ou3_alt_contraction import finite_source_bound_mag_forcing as MAG


@dataclass(frozen=True)
class Result:
    state: ADLIVE.State
    event: object
    forcing: MAG.MagForcing | None


def mag_step(state:ADLIVE.State, *, origin:ABRMM.RestrictedOrigin|None=None, **kwargs):
    if not isinstance(state,ADLIVE.State):
        raise TypeError('joint admitted BRMM/BIAS Live state required')
    if not state.live_word.source.steps:
        if not isinstance(origin,ABRMM.RestrictedOrigin):
            raise ValueError('fresh magnetic edge requires admitted-history t_L origin restriction')
        if origin.history != state.admitted_history:
            raise ValueError('fresh magnetic origin detached from carried admitted history')
        qualified=ABRMM.qualify_origin(state.live_word.source.root,origin)
        core=state.live_word.live.live.live.mekf
        if qualified.endpoint != core.reference:
            raise ValueError('fresh MEKF physical reference is not admitted-history t_L origin')
        if core.reference.bias_root != state.bias_history.history_id or core.reference.bias_family != state.bias_history.family:
            raise ValueError('fresh MEKF reference detached from admitted BIAS history')
    elif origin is not None:
        raise ValueError('post-first-IMU magnetic edge consumes no fresh origin restriction')
    before_brmm=state.admitted_history; before_bias=state.bias_history
    out=MAG.mag_step(state.live_word,**kwargs)
    nxt=ADLIVE.State(out.state,before_brmm,before_bias)
    return Result(nxt,out.word.event,out.forcing)


def readiness():
    m=MAG.readiness(); a=ABRMM.readiness()
    return {
      'post_first_IMU_magnetic_edge_uses_carried_admitted_source_endpoint':True,
      'fresh_sample_zero_magnetic_edge_requires_explicit_admitted_tL_origin':a['same_admitted_history_has_explicit_tL_origin_restriction'],
      'fresh_sample_zero_MEKF_reference_equal_to_admitted_history_origin_checked':True,
      'magnetic_event_advances_no_BRMM_or_BIAS_restriction_ordinal':True,
      'admitted_BRMM_and_BIAS_histories_preserved_across_magnetic_event':True,
      'same_event_correlated_magnetic_forcing_retained':m['magnetic_effective_residual_derived_from_same_qualified_event'],
      'startup_reachability_to_fresh_origin_state_closed':False,
      'sensor_residual_admissibility_attached':False,
      'deployment_roundoff_supply_attached':False,
      'complete_interleaved_600_step_word_composed':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
      'ALT_STARTUP_PASS':False,
      'ALT_END_TO_END_PASS':False,
    }
