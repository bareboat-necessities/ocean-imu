"""Admitted Live TuneState word with persistent machine vertical/stillness ancestry.

The lower word already carries one admitted BRMM/BIAS/ISS history, exact guarded
shipping IMU events, coherent WPE/tau/TuneState histories, and separate/FMA
machine band/statistics/sigma targets.  This layer removes two remaining free
upstream inputs from those machine frontend results:

* the adaptive-band input must be the same private-Mahony vertical successor;
* sigma still/still_time/attenuation must be the same machine tracker-LPF /
  StillnessAdapter successor.

Two LPF/stillness histories are retained because permitted contraction choices
may diverge.  Their private-Mahony predecessors and successors must remain equal:
Mahony is one deterministic source relation driven by the SAME guarded shipping
sample, not two physical histories.

The mutable tracker-LPF cutoff is still a carried local state, not yet rooted in
``RuntimeConfig``.  Likewise the machine guard arithmetic and startup-to-Live
binary32 frontend bridge remain open.  Therefore this is a composition advance,
not a source-uniform 600-edge theorem or storage certificate.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_admitted_machine_tunestate_interleaved_prefix as LOWER
from tools.stability.ou3_alt_contraction import finite_machine_vertical_stillness_source as VS
from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_source_continuation as SOURCE

QUALIFICATION='OU3_ALT_ADMITTED_MACHINE_VERTICAL_STILLNESS_INTERLEAVER_V1'


def _runtime(base:LOWER.State):
    if not isinstance(base,LOWER.State): raise TypeError('machine TuneState state required')
    return base.base.base.prefix.prefix.live.live_word.runtime


def _machine_imu_ordinal(base:LOWER.State):
    return base.machine.tau.updates-base.live_entry_machine_updates


@dataclass(frozen=True)
class State:
    base:LOWER.State
    separate_source:VS.State
    fma_source:VS.State
    entry_machine_ordinal:int
    source_steps:int=0
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.base,LOWER.State): raise TypeError('machine TuneState state required')
        if not isinstance(self.separate_source,VS.State) or not isinstance(self.fma_source,VS.State):
            raise TypeError('both machine vertical/stillness source histories required')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong admitted vertical/stillness qualification')
        if not isinstance(self.entry_machine_ordinal,int) or not isinstance(self.source_steps,int) or self.entry_machine_ordinal<0 or self.source_steps<0:
            raise ValueError('nonnegative machine source counters required')
        if _machine_imu_ordinal(self.base)-self.entry_machine_ordinal!=self.source_steps:
            raise ValueError('machine vertical/stillness count detached from TuneState IMU ordinal')
        if self.separate_source.samples!=self.fma_source.samples or self.separate_source.samples!=self.base.frontends.samples:
            raise ValueError('machine vertical/stillness sample count detached from machine frontend history')
        if self.separate_source.vertical!=self.fma_source.vertical:
            raise ValueError('separate/FMA histories cannot duplicate or diverge private-Mahony state')


@dataclass(frozen=True)
class ImuResult:
    state:State
    lower:LOWER.ImuResult
    separate_source:VS.Result
    fma_source:VS.Result
    def __post_init__(self):
        if not isinstance(self.state,State) or not isinstance(self.lower,LOWER.ImuResult):
            raise TypeError('admitted vertical/stillness result malformed')
        if not isinstance(self.separate_source,VS.Result) or not isinstance(self.fma_source,VS.Result):
            raise TypeError('both machine vertical/stillness results required')
        if self.separate_source.mahony!=self.fma_source.mahony:
            raise ValueError('compiler histories diverged in common private-Mahony source')
        if self.state.separate_source!=self.separate_source.state or self.state.fma_source!=self.fma_source.state:
            raise ValueError('persistent machine source state detached from same-event result')


@dataclass(frozen=True)
class CompleteWord:
    state:State
    lower:LOWER.CompleteWord
    def __post_init__(self):
        if not isinstance(self.state,State) or not isinstance(self.lower,LOWER.CompleteWord):
            raise TypeError('vertical/stillness complete word requires lower complete word')
        if self.lower.state!=self.state.base: raise ValueError('vertical/stillness complete word detached from lower word')
        if self.state.source_steps!=SOURCE.TRANSITIONS:
            raise ValueError('machine vertical/stillness ancestry not attached to all 600 IMU edges')


def begin(base:LOWER.State,*,separate_source:VS.State,fma_source:VS.State):
    """Begin from already-carried machine frontend states.

    This constructor deliberately does not manufacture the startup binary32
    Mahony/LPF/stillness state.  That bridge remains a separate proof obligation.
    """
    return State(base,separate_source,fma_source,_machine_imu_ordinal(base),0)


def _source_step(source,live,runtime,h,*,lpf_alpha_exp,lpf_successor,
                 still_energy_successor,still_attenuation_exp):
    return VS.step(source,live.guarded,runtime.vertical_cfg,runtime.still_cfg,
        dt=B.rn32(h),lpf_alpha_exp=lpf_alpha_exp,lpf_successor=lpf_successor,
        still_energy_successor=still_energy_successor,
        still_attenuation_exp=still_attenuation_exp)


def imu_step(state:State,*,
             separate_lpf_alpha_exp=None,separate_lpf_successor=None,
             separate_still_energy_successor,separate_still_attenuation_exp=None,
             fma_lpf_alpha_exp=None,fma_lpf_successor=None,
             fma_still_energy_successor,fma_still_attenuation_exp=None,
             **kwargs):
    if not isinstance(state,State): raise TypeError('admitted machine vertical/stillness State required')
    restricted=kwargs.get('restricted')
    if restricted is None: raise TypeError('same admitted physical restriction required')
    lower=LOWER.imu_step(state.base,**kwargs)
    live=LOWER._live_result(lower.lower)
    runtime=_runtime(state.base); h=restricted.segment.h
    sep=_source_step(state.separate_source,live,runtime,h,
        lpf_alpha_exp=separate_lpf_alpha_exp,lpf_successor=separate_lpf_successor,
        still_energy_successor=separate_still_energy_successor,
        still_attenuation_exp=separate_still_attenuation_exp)
    fma=_source_step(state.fma_source,live,runtime,h,
        lpf_alpha_exp=fma_lpf_alpha_exp,lpf_successor=fma_lpf_successor,
        still_energy_successor=fma_still_energy_successor,
        still_attenuation_exp=fma_still_attenuation_exp)
    if sep.mahony!=fma.mahony:
        raise ValueError('separate/FMA machine histories lost common private-Mahony source')
    VS.require_frontend_input(sep,lower.separate_frontend)
    VS.require_frontend_input(fma,lower.fma_frontend)
    VS.require_sigma_stillness(sep,lower.separate_sigma_join.machine)
    VS.require_sigma_stillness(fma,lower.fma_sigma_join.machine)
    nxt=State(lower.state,sep.state,fma.state,state.entry_machine_ordinal,state.source_steps+1)
    return ImuResult(nxt,lower,sep,fma)


def mag_step(state:State,**kwargs):
    if not isinstance(state,State): raise TypeError('admitted machine vertical/stillness State required')
    base,event=LOWER.mag_step(state.base,**kwargs)
    return State(base,state.separate_source,state.fma_source,state.entry_machine_ordinal,state.source_steps),event


def set_hold(state:State,*,hold):
    if not isinstance(state,State): raise TypeError('admitted machine vertical/stillness State required')
    base,event=LOWER.set_hold(state.base,hold=hold)
    return State(base,state.separate_source,state.fma_source,state.entry_machine_ordinal,state.source_steps),event


def complete(state:State): return CompleteWord(state,LOWER.complete(state.base))


def readiness():
    low=LOWER.readiness(); src=VS.readiness()
    return {
      'admitted_machine_TuneState_word_consumed':low['Live_600_step_machine_TuneState_product_attached'],
      'persistent_separate_and_FMA_machine_vertical_stillness_histories_attached':True,
      'private_Mahony_source_is_common_across_compiler_histories':True,
      'same_guarded_shipping_sample_drives_machine_private_Mahony':True,
      'machine_band_input_bound_to_same_private_Mahony_successor':True,
      'machine_sigma_stillness_bound_to_same_tracker_LPF_successor':True,
      'machine_tracker_LPF_binary32_recurrence_materialized':src['FreqInputLPF_binary32_recurrence_materialized'],
      'machine_stillness_binary32_projection_materialized':src['same_stored_LPF_output_drives_binary32_stillness_projection'],
      'complete_word_requires_machine_source_on_all_600_IMU_edges':True,
      'tracker_LPF_cutoff_runtime_ancestry_closed':False,
      'machine_guard_binary32_history_attached':False,
      'startup_machine_vertical_stillness_history_attached':False,
      'target_libm_and_compiler_profile_correspondence_closed':False,
      'source_uniform_machine_supply_bounds_closed':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,'ALT_STARTUP_PASS':False,'ALT_END_TO_END_PASS':False,
    }
