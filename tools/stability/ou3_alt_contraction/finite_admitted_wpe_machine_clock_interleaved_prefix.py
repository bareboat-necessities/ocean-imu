"""Attach full WPE moment/log machine histories to the strongest admitted Live word.

The lower clock-qualified joined product owns the literal IMU/magnetic/hold
events. This wrapper executes no second physical/filter event. On every IMU edge
it consumes the already-produced common machine Mahony vertical output and
advances the persistent dual-compiler WPE moment/log product with the same
binary32 API dt. The mode-specific log witnesses are injected into the lower
TuneState word from this single WPE product, so they cannot be spliced.

A startup constructor accepts only the WPE state preserved by the joined
startup/goLive product. ``complete`` requires exactly 600 additional WPE IMU
steps. Target libm/compiler qualification and source-uniform witness bounds
remain fail-closed.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_admitted_machine_clock_qualified_interleaved_prefix as LOWER
from tools.stability.ou3_alt_contraction import finite_admitted_machine_joined_frontend_racc_interleaved_prefix as JOIN
from tools.stability.ou3_alt_contraction import finite_startup_wpe_machine_history as STARTWPE
from tools.stability.ou3_alt_contraction import finite_startup_joined_machine_history as STARTLOW
from tools.stability.ou3_alt_contraction import finite_wpe_machine_binary32 as WPE
from tools.stability.ou3_alt_contraction import finite_source_continuation as SOURCE
from tools.stability.ou3_alt_contraction import finite_admitted_tau_interleaved_prefix as TAUJOIN

QUALIFICATION='OU3_ALT_ADMITTED_WPE_MACHINE_CLOCK_INTERLEAVER_V1'


def _logs(base:LOWER.State):
    if not isinstance(base,LOWER.State): raise TypeError('clock-qualified lower state required')
    mt=JOIN._mtune_state(base.base.base)
    return mt.base.wpe


@dataclass(frozen=True)
class State:
    base:LOWER.State
    wpe:WPE.State
    entry_wpe_samples:int
    wpe_steps:int=0
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.base,LOWER.State) or not isinstance(self.wpe,WPE.State):
            raise TypeError('clock-qualified state and full WPE machine state required')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong admitted WPE-machine qualification')
        if not isinstance(self.entry_wpe_samples,int) or not isinstance(self.wpe_steps,int) or self.entry_wpe_samples<0 or self.wpe_steps<0:
            raise ValueError('nonnegative WPE counters required')
        mt=JOIN._mtune_state(self.base.base.base)
        WPE.require_shadow_usable(self.wpe,TAUJOIN._entry_wpe(mt.base.base))
        if self.wpe.logs!=_logs(self.base):
            raise ValueError('full WPE machine logs detached from lower admitted TuneState word')
        if self.wpe.separate.samples-self.entry_wpe_samples!=self.wpe_steps or self.wpe.fma.samples-self.entry_wpe_samples!=self.wpe_steps:
            raise ValueError('full WPE machine count detached from Live WPE edge count')
        if self.base.clock_steps!=self.wpe_steps:
            raise ValueError('full WPE machine count detached from clock-qualified IMU count')


def begin(base:LOWER.State,wpe:WPE.State):
    if not isinstance(wpe,WPE.State): raise TypeError('startup-produced full WPE machine state required')
    return State(base,wpe,wpe.separate.samples,0)


def begin_from_startup(go:STARTWPE.GoLive,magnetic,origin,bias_history,**source_witnesses):
    if not isinstance(go,STARTWPE.GoLive): raise TypeError('full-WPE startup goLive result required')
    joined=STARTLOW.admitted_live(go.lower,magnetic,origin,bias_history,**source_witnesses)
    return begin(LOWER.begin(joined),go.wpe)


@dataclass(frozen=True)
class ImuResult:
    state:State
    lower:LOWER.ImuResult
    separate_wpe:object
    fma_wpe:object
    def __post_init__(self):
        if self.lower.state!=self.state.base: raise ValueError('WPE Live successor detached from lower event')


def imu_step(state:State,*,separate_wpe:WPE.ModeWitnesses,fma_wpe:WPE.ModeWitnesses,**kwargs):
    if not isinstance(state,State): raise TypeError('admitted WPE-machine Live State required')
    if 'separate_log_witness' in kwargs or 'fma_log_witness' in kwargs:
        raise TypeError('Live log witnesses are owned by the full WPE machine product')
    lower=LOWER.imu_step(state.base,separate_log_witness=separate_wpe.log,
                         fma_log_witness=fma_wpe.log,**kwargs)
    joined=lower.lower
    if joined.separate_source.band_input!=joined.fma_source.band_input:
        raise ValueError('Live compiler histories lost common Mahony WPE input')
    h=kwargs.get('machine_dt')
    if h is None: raise TypeError('same Live binary32 machine_dt required for WPE attachment')
    nxt,sr,fr,logs=WPE.step(state.wpe,dt=h,vertical_accel=joined.separate_source.band_input,
                            separate=separate_wpe,fma=fma_wpe)
    if nxt.logs!=_logs(lower.state):
        raise ValueError('machine WPE moment branch/log successor differs from lower admitted Live history')
    return ImuResult(State(lower.state,nxt,state.entry_wpe_samples,state.wpe_steps+1),lower,sr,fr)


def mag_step(state:State,**kwargs):
    if not isinstance(state,State): raise TypeError('admitted WPE-machine Live State required')
    base,event=LOWER.mag_step(state.base,**kwargs)
    return State(base,state.wpe,state.entry_wpe_samples,state.wpe_steps),event


def set_hold(state:State,*,hold):
    if not isinstance(state,State): raise TypeError('admitted WPE-machine Live State required')
    base,event=LOWER.set_hold(state.base,hold=hold)
    return State(base,state.wpe,state.entry_wpe_samples,state.wpe_steps),event


@dataclass(frozen=True)
class CompleteWord:
    state:State
    lower:LOWER.CompleteWord
    def __post_init__(self):
        if self.lower.state!=self.state.base: raise ValueError('WPE complete word detached from lower complete word')
        if self.state.wpe_steps!=SOURCE.TRANSITIONS:
            raise ValueError('full WPE machine history not attached to all 600 IMU edges')


def complete(state:State): return CompleteWord(state,LOWER.complete(state.base))


def readiness():
    sw=STARTWPE.readiness(); w=WPE.readiness()
    return {
      'strong_clock_qualified_joined_word_consumed_without_replay':True,
      'startup_produced_full_WPE_state_constructor_available':sw['goLive_preserves_full_WPE_moment_log_machine_history_by_identity'],
      'same_live_machine_Mahony_vertical_drives_full_WPE_product':w['same_machine_vertical_sample_drives_both_compiler_WPE_histories'],
      'same_mode_WPE_log_witness_injected_into_lower_TuneState_event':True,
      'MAG_and_HOLD_preserve_full_WPE_machine_history_by_identity':True,
      'complete_word_requires_full_WPE_machine_step_on_all_600_IMU_edges':True,
      'Live_600_step_WPE_machine_history_attached':True,
      'machine_vs_exact_WPE_period_branch_robustness_closed':False,
      'target_WPE_libm_and_compiler_profile_correspondence_closed':False,
      'source_uniform_WPE_machine_supply_bounds_closed':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_STARTUP_PASS':False,'ALT_LIVE_PASS':False,'ALT_END_TO_END_PASS':False,
    }
