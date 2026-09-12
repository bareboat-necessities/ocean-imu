"""Post-first-IMU admitted-history magnetic edge for the ALT master word.

Magnetometer calls do not advance the 5 ms physical source.  After at least one
admitted IMU restriction has established the current endpoint, this layer
executes the source-owned magnetic edge, retains its correlated forcing, and
proves that both admitted BRMM and BIAS histories survive unchanged.

The special magnetic call at fresh sample zero is intentionally rejected here:
its equality to the quantified admitted history's t_L sample is a separate
startup bridge still to be proved, and must not be inferred from matching ids.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_admitted_source_live_word as ADLIVE
from tools.stability.ou3_alt_contraction import finite_source_bound_mag_forcing as MAG


@dataclass(frozen=True)
class Result:
    state: ADLIVE.State
    event: object
    forcing: MAG.MagForcing | None


def mag_step(state:ADLIVE.State, **kwargs):
    if not isinstance(state,ADLIVE.State):
        raise TypeError('joint admitted BRMM/BIAS Live state required')
    if not state.live_word.source.steps:
        raise ValueError('sample-zero admitted-history magnetic equality is not proved here')
    before_brmm=state.admitted_history; before_bias=state.bias_history
    out=MAG.mag_step(state.live_word,**kwargs)
    nxt=ADLIVE.State(out.state,before_brmm,before_bias)
    return Result(nxt,out.word.event,out.forcing)


def readiness():
    m=MAG.readiness()
    return {
      'post_first_IMU_magnetic_edge_uses_carried_admitted_source_endpoint':True,
      'magnetic_event_advances_no_BRMM_or_BIAS_restriction_ordinal':True,
      'admitted_BRMM_and_BIAS_histories_preserved_across_magnetic_event':True,
      'same_event_correlated_magnetic_forcing_retained':m['magnetic_effective_residual_derived_from_same_qualified_event'],
      'fresh_sample_zero_magnetic_admitted_history_equality_closed':False,
      'sensor_residual_admissibility_attached':False,
      'deployment_roundoff_supply_attached':False,
      'complete_interleaved_600_step_word_composed':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
      'ALT_STARTUP_PASS':False,
      'ALT_END_TO_END_PASS':False,
    }
