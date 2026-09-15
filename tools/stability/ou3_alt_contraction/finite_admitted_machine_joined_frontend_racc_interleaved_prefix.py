"""Join guard/Mahony/LPF/stillness ancestry to the SAME Racc measurement word.

Two previously independent ALT finite products meet here without executing a
second physical/filter event:

* the lower Racc/measurement product owns the one admitted IMU event and the
  machine prediction -> S service -> Racc -> accelerometer relations;
* this layer extracts that event's already-produced machine TuneState result,
  advances one persistent binary32 guard/private-Mahony/LPF/stillness history,
  and requires its outputs to equal the machine frontend and sigma-stillness
  operands that the lower event actually consumed.

Thus the adaptive-band input and sigma stillness operands can no longer be free
witnesses in the Racc/measurement word.  Separate/FMA histories may diverge only
at the tracker LPF/stillness arithmetic allowed by the lower compiler model;
the guard and private Mahony successor are common.

Configuration is also checked against the current admitted runtime relation:
the guard config is the binary32 projection of RuntimeConfig.guard_cfg and the
tracker LPF uses the represented shipping default/reset cutoff.  Mutable setter
ancestry, startup provenance, native libm/Eigen/compiler correspondence and
source-uniform supply bounds remain open.  Storage search is not authorized.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_admitted_machine_racc_supply_interleaved_prefix as LOWER
from tools.stability.ou3_alt_contraction import finite_admitted_machine_tunestate_interleaved_prefix as MTUNE
from tools.stability.ou3_alt_contraction import finite_admitted_machine_vertical_stillness_interleaved_prefix as VERTICAL
from tools.stability.ou3_alt_contraction import finite_admitted_machine_runtime_config_interleaved_prefix as CONFIG
from tools.stability.ou3_alt_contraction import finite_machine_vertical_stillness_source as VS
from tools.stability.ou3_alt_contraction import finite_machine_accel_guard_binary32 as GUARD
from tools.stability.ou3_alt_contraction import finite_source_continuation as SOURCE

QUALIFICATION='OU3_ALT_ADMITTED_MACHINE_JOINED_FRONTEND_RACC_INTERLEAVER_V1'


def _mtune_state(base:LOWER.State):
    if not isinstance(base,LOWER.State):
        raise TypeError('machine Racc state required')
    try:
        out=base.base.base.base.base.base
    except AttributeError as exc:
        raise TypeError('Racc word lost nested machine TuneState state') from exc
    if not isinstance(out,MTUNE.State):
        raise TypeError('Racc word lost nested machine TuneState state')
    return out


def _runtime(base:LOWER.State):
    return VERTICAL._runtime(_mtune_state(base))


def _validate_config_and_counts(base:LOWER.State,guard,guard_cfg,separate_source,fma_source):
    mt=_mtune_state(base); runtime=VERTICAL._runtime(mt)
    if not isinstance(guard,GUARD.State) or not isinstance(guard_cfg,GUARD.Config):
        raise TypeError('persistent common machine guard state/config required')
    if guard_cfg!=CONFIG._machine_guard_cfg(runtime):
        raise ValueError('machine guard configuration detached from carried RuntimeConfig.guard_cfg')
    if not isinstance(separate_source,VS.State) or not isinstance(fma_source,VS.State):
        raise TypeError('both persistent machine vertical/stillness histories required')
    n=mt.frontends.samples
    if guard.samples!=n or separate_source.samples!=n or fma_source.samples!=n:
        raise ValueError('joined guard/vertical/stillness sample count detached from consumed machine frontend history')
    if separate_source.vertical!=fma_source.vertical:
        raise ValueError('joined compiler histories cannot diverge private-Mahony state')
    for mode,src in (('separate',separate_source),('fma',fma_source)):
        if src.lpf.cutoff_hz!=CONFIG.DEFAULT_TRACKER_CUTOFF:
            raise ValueError(mode+' tracker LPF cutoff detached from current shipping default/reset theorem scope')
    if not CONFIG._default_tracker_cutoff_source_matches():
        raise RuntimeError('shipping tracker-LPF default/reset source shape changed')
    return mt,runtime


@dataclass(frozen=True)
class State:
    base:LOWER.State
    guard:GUARD.State
    guard_cfg:GUARD.Config
    separate_source:VS.State
    fma_source:VS.State
    entry_racc_steps:int
    source_steps:int=0
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.base,LOWER.State):
            raise TypeError('machine Racc state required')
        if self.qualification!=QUALIFICATION:
            raise ValueError('wrong joined frontend/Racc qualification')
        if not isinstance(self.entry_racc_steps,int) or not isinstance(self.source_steps,int):
            raise TypeError('integer joined-source counters required')
        if self.entry_racc_steps<0 or self.source_steps<0:
            raise ValueError('nonnegative joined-source counters required')
        if self.base.racc_steps-self.entry_racc_steps!=self.source_steps:
            raise ValueError('joined frontend ancestry count detached from Racc IMU count')
        _validate_config_and_counts(self.base,self.guard,self.guard_cfg,self.separate_source,self.fma_source)


@dataclass(frozen=True)
class ImuResult:
    state:State
    lower:LOWER.ImuResult
    guard:GUARD.Result
    separate_source:VS.Result
    fma_source:VS.Result
    def __post_init__(self):
        if not isinstance(self.state,State) or not isinstance(self.lower,LOWER.ImuResult):
            raise TypeError('joined frontend/Racc IMU result malformed')
        if not isinstance(self.guard,GUARD.Result) or not isinstance(self.separate_source,VS.Result) or not isinstance(self.fma_source,VS.Result):
            raise TypeError('joined machine guard/source results required')
        if self.lower.state!=self.state.base:
            raise ValueError('joined frontend/Racc state detached from lower Racc event')
        if self.separate_source.mahony!=self.fma_source.mahony:
            raise ValueError('joined compiler histories lost common private-Mahony source')
        if self.state.guard!=self.guard.state or self.state.separate_source!=self.separate_source.state or self.state.fma_source!=self.fma_source.state:
            raise ValueError('joined persistent machine source state detached from same-event results')


@dataclass(frozen=True)
class CompleteWord:
    state:State
    lower:LOWER.CompleteWord
    def __post_init__(self):
        if not isinstance(self.state,State) or not isinstance(self.lower,LOWER.CompleteWord):
            raise TypeError('joined frontend/Racc complete word requires lower complete word')
        if self.lower.state!=self.state.base:
            raise ValueError('joined frontend/Racc complete word detached from lower word')
        if self.state.source_steps!=SOURCE.TRANSITIONS:
            raise ValueError('joined machine frontend ancestry not attached to all 600 Racc/measurement IMU edges')


def begin(base:LOWER.State,*,guard:GUARD.State,guard_cfg:GUARD.Config,separate_source:VS.State,fma_source:VS.State):
    _validate_config_and_counts(base,guard,guard_cfg,separate_source,fma_source)
    return State(base,guard,guard_cfg,separate_source,fma_source,base.racc_steps,0)


def begin_from_goLive(base:LOWER.State,go):
    """Substitute the joined startup histories, never caller-selected Live seeds.

    This is a conditional handoff identity. It does not establish that every
    admitted startup reaches TunerReady or qualify missing deployment branches.
    """
    from tools.stability.ou3_alt_contraction import finite_startup_joined_machine_history as START
    if not isinstance(go,START.GoLive):
        raise TypeError('joined startup goLive result required')
    mt=_mtune_state(base); source=go.startup; lower=go.lower
    if _runtime(base)!=source.runtime or mt.deployment_cfg!=source.deployment_cfg:
        raise ValueError('Live runtime configuration detached from joined startup')
    if MTUNE._entry_live(mt.base)!=lower.live.state:
        raise ValueError('Live exact state detached from joined startup goLive')
    if (mt.machine,mt.base.wpe,mt.frontends)!=(lower.machine,lower.wpe,lower.frontends):
        raise ValueError('Live machine TuneState/WPE/band memory detached from joined startup')
    if (mt.separate_active,mt.fma_active)!=(lower.separate_active,lower.fma_active):
        raise ValueError('Live active parameters detached from same machine goLive commits')
    if base.racc_steps or mt.machine.tau.updates!=mt.live_entry_machine_updates:
        raise ValueError('startup attachment must precede the first Live IMU event')
    if base.separate_racc!=source.racc or base.fma_racc!=source.racc:
        raise ValueError('Live Racc state detached from bootstrap identity history')
    return begin(base,guard=source.guard,guard_cfg=source.guard_cfg,
                 separate_source=source.separate_source,fma_source=source.fma_source)


def imu_step(state:State,*,machine_dt,machine_gyro_body,machine_acc_body,
             guard_lp_alpha_exp=None,guard_lp_successors=None,guard_detect_gamma_exp=None,guard_detect_successors=None,
             guard_removed_beta_exp=None,guard_removed_ms_successor=None,guard_removed_rms_sqrt=None,guard_slew_exp=None,guard_weight_successor=None,guard_output_successor=None,
             separate_lpf_alpha_exp=None,separate_lpf_successor=None,separate_still_energy_successor,separate_still_attenuation_exp=None,
             fma_lpf_alpha_exp=None,fma_lpf_successor=None,fma_still_energy_successor,fma_still_attenuation_exp=None,
             **kwargs):
    if not isinstance(state,State):
        raise TypeError('joined frontend/Racc State required')
    restricted=kwargs.get('restricted')
    if restricted is None:
        raise TypeError('same admitted physical restriction required')

    # Execute the one admitted Racc/measurement event first. Everything below
    # is an arithmetic/source relation attached to this result, not a replay.
    lower=LOWER.imu_step(state.base,**kwargs)
    mtune=LOWER._mtune_result(lower.lower)
    live=MTUNE._live_result(mtune.lower)
    runtime=_runtime(state.base); h=restricted.segment.h

    hq=GUARD.api_scalar(h,machine_dt,'filter update dt')
    gyro=GUARD.api_vec(live.guarded.raw_gyro_body,machine_gyro_body,'filter update gyro')
    acc=GUARD.api_vec(live.guarded.raw_accel_body,machine_acc_body,'filter update accel')
    guarded=GUARD.step(state.guard,state.guard_cfg,raw_gyro=gyro,raw_acc=acc,dt=hq,
        lp_alpha_exp=guard_lp_alpha_exp,lp_successors=guard_lp_successors,
        detect_gamma_exp=guard_detect_gamma_exp,detect_successors=guard_detect_successors,
        removed_beta_exp=guard_removed_beta_exp,removed_ms_successor=guard_removed_ms_successor,
        removed_rms_sqrt=guard_removed_rms_sqrt,slew_exp=guard_slew_exp,
        weight_successor=guard_weight_successor,output_successor=guard_output_successor)

    sep=VERTICAL._source_step(state.separate_source,guarded,runtime,h,
        lpf_alpha_exp=separate_lpf_alpha_exp,lpf_successor=separate_lpf_successor,
        still_energy_successor=separate_still_energy_successor,still_attenuation_exp=separate_still_attenuation_exp)
    fma=VERTICAL._source_step(state.fma_source,guarded,runtime,h,
        lpf_alpha_exp=fma_lpf_alpha_exp,lpf_successor=fma_lpf_successor,
        still_energy_successor=fma_still_energy_successor,still_attenuation_exp=fma_still_attenuation_exp)
    if sep.mahony!=fma.mahony:
        raise ValueError('separate/FMA machine histories lost common private-Mahony source')

    # Critical join: these are the frontend/sigma objects already consumed by
    # the nested TuneState event inside the Racc/measurement word.
    VS.require_frontend_input(sep,mtune.separate_frontend)
    VS.require_frontend_input(fma,mtune.fma_frontend)
    VS.require_sigma_stillness(sep,mtune.separate_sigma_join.machine)
    VS.require_sigma_stillness(fma,mtune.fma_sigma_join.machine)

    # Strengthen the lower comparison shadow with the already-executed machine
    # guard arithmetic. No source/filter event is replayed: only the proof-side
    # Racc and held-accelerometer numerical relations are rebound to the same
    # raw packet's binary32 conditioned sample and excess RMS.
    lower=LOWER.rebind_machine_guard(state.base,lower,
        guard_excess_rms=guarded.excess_rms,conditioned_accel_body=guarded.conditioned_acc,
        separate_racc_accel_ldlt=kwargs['separate_racc_accel_ldlt'],
        fma_racc_accel_ldlt=kwargs['fma_racc_accel_ldlt'],
        separate_rao_witness=kwargs.get('separate_rao_witness'),
        separate_racc_sqrt=kwargs.get('separate_racc_sqrt'),
        fma_rao_witness=kwargs.get('fma_rao_witness'),fma_racc_sqrt=kwargs.get('fma_racc_sqrt'),
        temperature_c=kwargs['temperature_c'],restricted=restricted,
        accel_alpha=kwargs.get('accel_alpha',1),accel_radius=kwargs.get('accel_radius'))

    nxt=State(lower.state,guarded.state,state.guard_cfg,sep.state,fma.state,
              state.entry_racc_steps,state.source_steps+1)
    return ImuResult(nxt,lower,guarded,sep,fma)


def mag_step(state:State,**kwargs):
    if not isinstance(state,State):
        raise TypeError('joined frontend/Racc State required')
    base,event=LOWER.mag_step(state.base,**kwargs)
    return State(base,state.guard,state.guard_cfg,state.separate_source,state.fma_source,
                 state.entry_racc_steps,state.source_steps),event


def set_hold(state:State,*,hold):
    if not isinstance(state,State):
        raise TypeError('joined frontend/Racc State required')
    base,event=LOWER.set_hold(state.base,hold=hold)
    return State(base,state.guard,state.guard_cfg,state.separate_source,state.fma_source,
                 state.entry_racc_steps,state.source_steps),event


def complete(state:State):
    return CompleteWord(state,LOWER.complete(state.base))


def readiness():
    low=LOWER.readiness(); cfg=CONFIG.readiness()
    return {
      'machine_Racc_measurement_word_consumed':low['machine_Racc_effect_propagated_through_accelerometer'],
      'same_executed_TuneState_event_supplies_frontend_and_sigma_join':True,
      'no_second_physical_or_filter_event_executed_for_frontend_ancestry':True,
      'common_machine_guard_and_private_Mahony_history_joined_to_Racc_word':True,
      'machine_guard_excess_drives_machine_Racc_recurrence':True,
      'machine_guard_conditioned_sample_drives_machine_accelerometer_relation':True,
      'machine_band_input_bound_to_frontend_consumed_by_Racc_word':True,
      'machine_sigma_stillness_bound_to_sigma_target_consumed_by_Racc_word':True,
      'machine_guard_runtime_config_bound_in_joined_word':cfg['machine_guard_runtime_config_ancestry_closed_for_current_word'],
      'tracker_LPF_default_cutoff_bound_in_joined_word':cfg['tracker_LPF_cutoff_bound_to_default_constructor_reset_for_current_word'],
      'complete_word_requires_join_on_all_600_Racc_measurement_IMU_edges':True,
      'tracker_LPF_mutable_setter_ancestry_closed':False,
      'private_Mahony_mutable_config_setter_ancestry_closed':False,
      'startup_joined_machine_history_attached':True,
      'startup_attachment_is_conditional_not_universal_reachability':True,
      'target_libm_Eigen_and_compiler_profile_correspondence_closed':False,
      'source_uniform_machine_supply_bounds_closed':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,'ALT_STARTUP_PASS':False,'ALT_END_TO_END_PASS':False,
    }
