"""Attach full WPE moment/log machine histories to the strongest admitted Live word.

The lower clock-qualified joined product owns the single physical source history
and exact-shadow events. This wrapper carries separate and FMA CORE/control
successors through their corresponding conditional event algorithms. On every IMU edge
it consumes the already-produced common machine Mahony vertical output and
advances the persistent dual-compiler WPE moment/log product with the same
binary32 API dt. The carried machine latches/logs select frequency independently
of the exact shadow. Mode-specific log witnesses are projected into the lower TuneState
word and checked against this same moment update, so they cannot be spliced.

A startup constructor accepts only the WPE state preserved by the joined
startup/goLive product. ``complete`` requires exactly 600 additional WPE IMU
steps. Target libm/compiler qualification and source-uniform witness bounds
remain fail-closed.
"""
from __future__ import annotations
from dataclasses import dataclass, replace
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_admitted_machine_clock_qualified_interleaved_prefix as LOWER
from tools.stability.ou3_alt_contraction import finite_admitted_machine_joined_frontend_racc_interleaved_prefix as JOIN
from tools.stability.ou3_alt_contraction import finite_startup_wpe_machine_history as STARTWPE
from tools.stability.ou3_alt_contraction import finite_startup_joined_machine_history as STARTLOW
from tools.stability.ou3_alt_contraction import finite_wpe_machine_binary32 as WPE
from tools.stability.ou3_alt_contraction import finite_source_continuation as SOURCE
from tools.stability.ou3_alt_contraction import finite_machine_core_continuation as COREWORD
from tools.stability.ou3_alt_contraction import finite_wpe_uniform_bounds as SUPPLY
from tools.stability.ou3_alt_contraction import finite_live_input_contract as INPUT_DOMAIN
from tools.stability.ou3_alt_contraction import finite_frontend_uniform_bounds as FRONT_SUPPLY
from tools.stability.ou3_alt_contraction import finite_candidate_uniform_bounds as CAND_SUPPLY

QUALIFICATION='OU3_ALT_ADMITTED_WPE_MACHINE_CLOCK_INTERLEAVER_V1'


def _logs(base:LOWER.State):
    if not isinstance(base,LOWER.State): raise TypeError('clock-qualified lower state required')
    mt=JOIN._mtune_state(base.base.base)
    return mt.base.wpe


def _preword(base:LOWER.State):
    return JOIN.LOWER.LOWER._preword(base.base.base.base)


@dataclass(frozen=True)
class State:
    base:LOWER.State
    wpe:WPE.State
    entry_wpe_samples:int
    wpe_steps:int=0
    machine_core:COREWORD.State|None=None
    qualification:str=QUALIFICATION
    uniform_live_supplies:bool=False
    def __post_init__(self):
        if not isinstance(self.base,LOWER.State) or not isinstance(self.wpe,WPE.State):
            raise TypeError('clock-qualified state and full WPE machine state required')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong admitted WPE-machine qualification')
        if not isinstance(self.entry_wpe_samples,int) or not isinstance(self.wpe_steps,int) or self.entry_wpe_samples<0 or self.wpe_steps<0:
            raise ValueError('nonnegative WPE counters required')
        if self.wpe.logs!=_logs(self.base):
            raise ValueError('full WPE machine logs detached from lower admitted TuneState word')
        if self.wpe.separate.samples-self.entry_wpe_samples!=self.wpe_steps or self.wpe.fma.samples-self.entry_wpe_samples!=self.wpe_steps:
            raise ValueError('full WPE machine count detached from Live WPE edge count')
        if self.base.clock_steps!=self.wpe_steps:
            raise ValueError('full WPE machine count detached from clock-qualified IMU count')
        if not isinstance(self.machine_core,COREWORD.State):
            raise TypeError('persistent machine CORE continuation required')
        if self.machine_core.imu_steps!=self.wpe_steps:
            raise ValueError('machine CORE continuation detached from WPE IMU count')
        if self.machine_core.separate.reference!=_preword(self.base).live.live.live.mekf.reference:
            raise ValueError('machine CORE continuation detached from admitted physical endpoint')
        if type(self.uniform_live_supplies) is not bool:
            raise TypeError('literal Live supply qualification required')
        if self.uniform_live_supplies!=self.wpe.bounded_profile:
            raise ValueError('bounded startup WPE must retain physical Live supply qualification')
        if self.uniform_live_supplies:
            report=SUPPLY.physical_input_certificate(self.wpe.libm_profile)
            runtime=JOIN._runtime(self.base.base.base)
            mtune=JOIN._mtune_state(self.base.base.base)
            CAND_SUPPLY.require_config(mtune.deployment_cfg)
            CAND_SUPPLY.require_commit_config(runtime.commit_cfg)
            CAND_SUPPLY.require_state(mtune.machine)
            if self.wpe.cfg!=WPE.MOM.Config.from_shadow(runtime.wpe_cfg):
                raise ValueError('Live WPE envelope detached from carried configuration')
            if WPE.B.rn32(runtime.vertical_cfg.gravity)!=WPE.B.rn32(F(196133,20000)):
                raise ValueError('Live WPE supply detached from declared observer gravity')
            if not self.wpe.bounded_profile or self.wpe.supply_scale_exponent<report['scale_exponent']:
                raise ValueError('Live WPE envelope detached from declared MEMS input domain')
            if self.base.base.guard.weight!=0:
                raise ValueError('Live WPE physical supply requires declared dormant guard scope')
            FRONT_SUPPLY.require_config(band_cfg=runtime.band_cfg,stats_cfg=runtime.stats_cfg,
                still_cfg=runtime.still_cfg,cutoff_hz=self.base.base.separate_source.lpf.cutoff_hz,
                dt=FRONT_SUPPLY.DT,bench_noise_sigma=WPE.B.rn32(runtime.bench_noise_sigma))
            frontends=JOIN._mtune_state(self.base.base.base).frontends
            FRONT_SUPPLY.require_state(self.base.base.separate_source,frontends.separate)
            FRONT_SUPPLY.require_state(self.base.base.fma_source,frontends.fma)


def begin(base:LOWER.State,wpe:WPE.State):
    if not isinstance(wpe,WPE.State): raise TypeError('startup-produced full WPE machine state required')
    if wpe.bounded_profile:
        report=SUPPLY.physical_input_certificate(wpe.libm_profile)
        wpe=WPE.widen_supply(wpe,input_abs_upper=report['vertical_input_abs_upper'])
    return State(base,wpe,wpe.separate.samples,0,
                 COREWORD.begin_from_interleaved(_preword(base).live,
                     separate_proxy=base.base.separate_source.vertical,
                     fma_proxy=base.base.fma_source.vertical),
                 uniform_live_supplies=wpe.bounded_profile)


def begin_with_uniform_supplies(base:LOWER.State,wpe:WPE.State):
    """Require the automatically retained physical supply qualification."""
    if not isinstance(wpe,WPE.State) or not wpe.bounded_profile:
        raise TypeError('startup-produced bounded WPE history required')
    return begin(base,wpe)


def begin_from_startup(go:STARTWPE.GoLive,magnetic,origin,bias_history,*,
                       qualify_live_supplies=False,**source_witnesses):
    if not isinstance(go,STARTWPE.GoLive): raise TypeError('full-WPE startup goLive result required')
    joined=STARTLOW.admitted_live(go.lower,magnetic,origin,bias_history,**source_witnesses)
    constructor=begin_with_uniform_supplies if qualify_live_supplies else begin
    out=constructor(LOWER.begin(joined),go.wpe)
    return replace(out,machine_core=COREWORD.begin_from_goLive(_preword(out.base).live,go.lower))


@dataclass(frozen=True)
class ImuResult:
    state:State
    lower:LOWER.ImuResult
    separate_wpe:object
    fma_wpe:object
    def __post_init__(self):
        if self.lower.state!=self.state.base: raise ValueError('WPE Live successor detached from lower event')


def imu_step(state:State,*,separate_wpe:WPE.ModeWitnesses,fma_wpe:WPE.ModeWitnesses,
             separate_machine_reset=None,fma_machine_reset=None,
             separate_machine_gravity=None,fma_machine_gravity=None,**kwargs):
    if not isinstance(state,State): raise TypeError('admitted WPE-machine Live State required')
    if {'separate_log_witness','fma_log_witness','machine_wpe_entry',
        'separate_machine_predecessor','fma_machine_predecessor','machine_packet'} & set(kwargs):
        raise TypeError('Live log witnesses are owned by the full WPE machine product')
    raw=kwargs.get('raw')
    if raw is None: raise TypeError('same raw Live source packet required')
    packet=INPUT_DOMAIN.check_packet(raw)
    if state.uniform_live_supplies:
        if raw.gravity_world!=(F(0),F(0),F(196133,20000)):
            raise ValueError('Live WPE supply detached from declared physical gravity')
    lower=LOWER.imu_step(state.base,separate_log_witness=separate_wpe.log,
                         fma_log_witness=fma_wpe.log,machine_wpe_entry=state.wpe,
                         machine_packet=packet,
                         separate_machine_predecessor=state.machine_core.separate,
                         fma_machine_predecessor=state.machine_core.fma,**kwargs)
    joined=lower.lower
    if state.uniform_live_supplies:
        if (joined.guard.conditioned_acc!=joined.guard.raw_acc or joined.guard.state.weight!=0
                or joined.guard.removed_rms>state.base.base.guard_cfg.engage_lo):
            raise ValueError('Live WPE supply requires retained transparent guard event')
        mtune=JOIN.LOWER._mtune_result(joined.lower.lower)
        for mode in ('separate','fma'):
            FRONT_SUPPLY.require_event(getattr(joined,mode+'_source'),
                getattr(mtune,mode+'_frontend'),sigma_target=getattr(mtune,mode+'_sigma_join').machine)
    if joined.separate_source.band_input!=joined.fma_source.band_input:
        raise ValueError('Live compiler histories lost common Mahony WPE input')
    h=kwargs.get('machine_dt')
    if h is None: raise TypeError('same Live binary32 machine_dt required for WPE attachment')
    nxt,sr,fr,logs=WPE.step(state.wpe,dt=h,vertical_accel=joined.separate_source.band_input,
                            separate=separate_wpe,fma=fma_wpe)
    if nxt.logs!=_logs(lower.state):
        raise ValueError('machine WPE moment branch/log successor differs from lower admitted Live history')
    racc=joined.lower
    measurement=racc.lower
    carry=COREWORD.advance_imu(state.machine_core,
        separate_predecessor=measurement.separate.prediction.machine_predecessor,
        fma_predecessor=measurement.fma.prediction.machine_predecessor,
        separate_after_accel=racc.separate.accelerometer.state,
        fma_after_accel=racc.fma.accelerometer.state,
        separate_measurement=measurement.separate,fma_measurement=measurement.fma,
        separate_racc=racc.separate,fma_racc=racc.fma,
        template=_preword(lower.state).live,
        exact_guarded=JOIN.LOWER.LOWER._live_result(measurement.lower).guarded,
        machine_guard=joined.guard,dt=kwargs['restricted'].segment.h,
        machine_packet=packet,
        separate_proxy=joined.separate_source.state.vertical,
        fma_proxy=joined.fma_source.state.vertical,
        separate_reset=separate_machine_reset,fma_reset=fma_machine_reset,
        separate_gravity=separate_machine_gravity,fma_gravity=fma_machine_gravity)
    return ImuResult(State(lower.state,nxt,state.entry_wpe_samples,state.wpe_steps+1,carry,
                           uniform_live_supplies=state.uniform_live_supplies),lower,sr,fr)


def mag_step(state:State,*,separate_machine_witnesses=None,fma_machine_witnesses=None,**kwargs):
    if not isinstance(state,State): raise TypeError('admitted WPE-machine Live State required')
    base,event=LOWER.mag_step(state.base,**kwargs)
    carry,_=COREWORD.advance_mag(state.machine_core,separate=separate_machine_witnesses,
                                 fma=fma_machine_witnesses,**kwargs)
    return State(base,state.wpe,state.entry_wpe_samples,state.wpe_steps,
                 carry,
                 uniform_live_supplies=state.uniform_live_supplies),event


def set_hold(state:State,*,hold):
    if not isinstance(state,State): raise TypeError('admitted WPE-machine Live State required')
    base,event=LOWER.set_hold(state.base,hold=hold)
    carry,_=COREWORD.advance_hold(state.machine_core,hold=hold)
    return State(base,state.wpe,state.entry_wpe_samples,state.wpe_steps,
                 carry,
                 uniform_live_supplies=state.uniform_live_supplies),event


@dataclass(frozen=True)
class CompleteWord:
    state:State
    lower:LOWER.CompleteWord
    def __post_init__(self):
        if self.lower.state!=self.state.base: raise ValueError('WPE complete word detached from lower complete word')
        if self.state.wpe_steps!=SOURCE.TRANSITIONS:
            raise ValueError('full WPE machine history not attached to all 600 IMU edges')
        COREWORD.require_complete(self.state.machine_core)


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
      'physical_MEMS_to_WPE_uniform_envelope_mandatory_for_bounded_startup':True,
      'physical_Live_envelope_does_not_limit_ISS_W_or_reuse_startup_residual_cap':True,
      'conditional_Live_envelope_preserves_all_machine_WPE_history':True,
      'machine_CORE_local_successors_persist_and_next_predecessors_checked':True,
      'complete_word_requires_machine_CORE_and_control_continuation':True,
      'machine_CORE_and_control_successor_algorithms_attached':True,
      'machine_prediction_and_gravity_consume_same_rounded_API_packet':True,
      'startup_machine_goLive_aw_covariance_sync_attached':True,
      'startup_machine_nominal_attitude_and_remaining_covariance_qualified':False,
      'machine_vs_exact_WPE_period_branch_agreement_required':False,
      'independent_machine_WPE_production_and_frequency_branches_composed':True,
      'target_WPE_libm_and_compiler_profile_correspondence_closed':False,
      'source_uniform_WPE_machine_supply_bounds_closed':SUPPLY.physical_input_certificate(
          SUPPLY.ERROR_PROFILE)['source_uniform_WPE_supplies_under_declared_MEMS_prefix_closed'],
      'WPE_supply_bound_requires_no_startup_deadline':True,
      'configured_frontend_uniform_supplies_attached_to_all_bounded_Live_IMU_events':True,
      'candidate_and_commit_uniform_supplies_attached_to_bounded_Live':True,
      'WPE_supply_seed_scalar_and_dormant_guard_premises_remain_explicit':True,
      'target_exp_log_error_bounds_fit_WPE_supply_profile':SUPPLY.target_error_profile_certificate()[
          'target_exp_log_satisfy_WPE_error_profile_under_scalar_and_link_premises'],
      'pinned_WPE_exp_log_approximation_correspondence_closed':SUPPLY.target_error_profile_certificate()[
          'pinned_WPE_exp_log_approximation_correspondence_closed'],
      'pinned_WPE_sqrt_approximation_correspondence_closed':SUPPLY.target_error_profile_certificate()[
          'pinned_WPE_sqrt_approximation_correspondence_closed'],
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_STARTUP_PASS':False,'ALT_LIVE_PASS':False,'ALT_END_TO_END_PASS':False,
    }
