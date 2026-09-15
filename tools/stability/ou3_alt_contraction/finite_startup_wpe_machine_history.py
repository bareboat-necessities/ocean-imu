"""Attach full binary32 WPE moment/log histories to the joined startup machine.

The lower startup product owns the one physical/frontend event. This wrapper
advances no second guard/Mahony/tuner event: it consumes the already-produced
machine Mahony vertical output and binds the dual binary32 WPE moment/log product
to the same sample. The log witnesses supplied to this product are also the
ones consumed by the lower tuner/WPE-log path.

The current lower product still expects exact-shadow and both machine log-update
branches to agree. This wrapper therefore fails closed if the newly materialized
machine moment branches would update on a different sample. Proving that branch
robustness source-uniformly is a remaining deployment obligation.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_startup_joined_machine_history as LOWER
from tools.stability.ou3_alt_contraction import finite_wpe_machine_binary32 as WPE

QUALIFICATION='OU3_ALT_STARTUP_WPE_MACHINE_HISTORY_V1'


@dataclass(frozen=True)
class State:
    base:LOWER.State
    wpe:WPE.State
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.base,LOWER.State) or not isinstance(self.wpe,WPE.State):
            raise TypeError('joined startup and full WPE machine states required')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong startup WPE-machine qualification')
        WPE.require_shadow_usable(self.wpe,self.base.base.lower.frontend.tuner.wpe)
        if self.wpe.logs!=self.base.base.lower.wpe:
            raise ValueError('full WPE machine log state detached from lower startup WPE log ledger')
        if self.wpe.separate.samples!=self.base.separate_source.samples or self.wpe.fma.samples!=self.base.fma_source.samples:
            raise ValueError('full WPE machine sample count detached from startup Mahony source history')


def initial(runtime,deployment_cfg):
    base=LOWER.initial(runtime,deployment_cfg)
    return State(base,WPE.initial(runtime.wpe_cfg))


@dataclass(frozen=True)
class StepResult:
    state:State
    lower:LOWER.StepResult
    separate_wpe:object
    fma_wpe:object
    def __post_init__(self):
        if self.state.base!=self.lower.state: raise ValueError('startup WPE-machine successor detached from lower event')


def step(state:State,raw,*,separate_wpe:WPE.ModeWitnesses,fma_wpe:WPE.ModeWitnesses,**kwargs):
    if not isinstance(state,State): raise TypeError('startup WPE-machine State required')
    if 'separate_log_witness' in kwargs or 'fma_log_witness' in kwargs:
        raise TypeError('startup log witnesses are owned by the full WPE machine product')
    lower=LOWER.step(state.base,raw,separate_log_witness=separate_wpe.log,
                     fma_log_witness=fma_wpe.log,**kwargs)
    if lower.separate_source.band_input!=lower.fma_source.band_input:
        raise ValueError('startup compiler histories lost common Mahony vertical WPE input')
    h=kwargs.get('machine_dt')
    if h is None: raise TypeError('same startup binary32 machine_dt required for WPE attachment')
    nxt,sr,fr,logs=WPE.step(state.wpe,dt=h,vertical_accel=lower.separate_source.band_input,
                            separate=separate_wpe,fma=fma_wpe)
    if nxt.logs!=lower.state.base.lower.wpe:
        raise ValueError('machine moment branch/log successor differs from lower startup WPE history')
    return StepResult(State(lower.state,nxt),lower,sr,fr)


def boundary(state:State,**kwargs):
    if not isinstance(state,State): raise TypeError('startup WPE-machine State required')
    base,event=LOWER.boundary(state.base,**kwargs)
    return State(base,state.wpe),event


@dataclass(frozen=True)
class GoLive:
    startup:State
    lower:LOWER.GoLive
    wpe:WPE.State
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.startup,State) or not isinstance(self.lower,LOWER.GoLive):
            raise TypeError('startup WPE machine and lower goLive required')
        if self.wpe!=self.startup.wpe: raise ValueError('goLive restarted full WPE machine history')
        if self.lower.startup!=self.startup.base: raise ValueError('goLive lower history detached')
        if self.lower.lower.wpe!=self.wpe.logs:
            raise ValueError('goLive WPE log state detached from full machine moment history')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong startup WPE-machine goLive qualification')


def go_live(state:State,entry,fresh,**kwargs):
    if not isinstance(state,State): raise TypeError('startup WPE-machine State required')
    return GoLive(state,LOWER.go_live(state.base,entry,fresh,**kwargs),state.wpe)


def readiness():
    w=WPE.readiness(); low=LOWER.readiness()
    return {
      'joined_guard_Mahony_startup_product_reused_without_second_event':True,
      'full_WPE_machine_histories_rooted_with_literal_startup_reset':True,
      'same_machine_Mahony_vertical_output_drives_full_WPE_moment_product':True,
      'same_mode_log_witness_consumed_by_moment_log_product_and_lower_tuner_path':True,
      'full_WPE_machine_log_state_equals_existing_lower_WPE_ledger_after_each_attached_step':True,
      'goLive_preserves_full_WPE_moment_log_machine_history_by_identity':True,
      'startup_frontend_vertical_ancestry_attached': bool(
          w['same_machine_vertical_sample_drives_both_compiler_WPE_histories'] and
          low['Cold_and_postCold_machine_sources_joined_to_one_executed_frontend_event']),
      'machine_vs_exact_WPE_period_branch_robustness_closed':False,
      'target_WPE_libm_and_compiler_profile_correspondence_closed':False,
      'source_uniform_startup_WPE_supply_bounds_closed':False,
      'Live_600_step_WPE_machine_history_attached':False,
      'storage_search_allowed':False,
      'ALT_STARTUP_PASS':False,'ALT_LIVE_PASS':False,
    }
