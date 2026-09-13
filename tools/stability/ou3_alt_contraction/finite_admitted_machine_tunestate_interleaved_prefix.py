"""Admitted BRMM/BIAS/ISS/WPE Live word with whole machine TuneState history.

This strengthens the coherent WPE-log/tau admitted product with persistent
``sigma_applied`` and ``RS_applied`` machine histories AND the separately
persisted applied ActiveParameters.  The distinction is essential: a tuner
candidate changes TuneState now, while prediction/measurement keep using the
last applied parameters until goLive or the next pending boundary.

Literal shipping order retained here is

  sample-entry candidate memory + applied parameters
    -> exact boundary inside the already-composed admitted IMU event
    -> matching whole-machine pending boundary and applied-parameter update
    -> prediction/measurements in the lower exact shadow
    -> exact tuner candidate from that same event
    -> per-mode tau/sigma/R_S machine candidate successor
    -> current-sample WPE-log successor.

The lower admitted word still owns all BRMM, BIAS and bounded ISS restrictions.
MAG/HOLD preserve WPE, TuneState and applied-machine parameters by identity.  A
completed word must advance all three candidate scalar ledgers exactly 600
Live IMU events after the one Live entry.

The lower Live dynamics still use exact-shadow applied parameters.  This layer
now exposes exact-vs-machine applied-parameter displacement at every state but
does not yet inject that displacement into prediction/measurement coefficients.
Storage therefore remains hard blocked.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_admitted_wpe_tau_interleaved_prefix as LOWER
from tools.stability.ou3_alt_contraction import finite_admitted_tau_interleaved_prefix as TAUJOIN
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as C
from tools.stability.ou3_alt_contraction import finite_tuner_deployment_config as D
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_binary32 as SIGM
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_machine_real_join as SIGJOIN
from tools.stability.ou3_alt_contraction import finite_tuner_machine_tunestate_product as PRODUCT
from tools.stability.ou3_alt_contraction import finite_tuner_machine_candidate_step as CANDSTEP
from tools.stability.ou3_alt_contraction import finite_tuner_machine_boundary as MBOUND
from tools.stability.ou3_alt_contraction import finite_tuner_common_alpha_qualification as COMMON
from tools.stability.ou3_alt_contraction import finite_startup_live_machine_tunestate_bridge as GO
from tools.stability.ou3_alt_contraction import finite_active_parameter_machine_real_join as ACTIVEJOIN
from tools.stability.ou3_alt_contraction import finite_runtime_parameters as ACTIVE
from tools.stability.ou3_alt_contraction import finite_source_continuation as SOURCE

QUALIFICATION='OU3_ALT_ADMITTED_MACHINE_TUNESTATE_INTERLEAVER_V2'


def _entry_live(base:LOWER.State):
    try: return base.base.prefix.prefix.live.live_word.live.live.live
    except AttributeError as exc: raise TypeError('admitted WPE/tau product lost Live IMU entry state') from exc


def _live_result(lower:LOWER.ImuResult):
    try: return lower.event.event.event.live
    except AttributeError as exc: raise TypeError('admitted whole-machine event lost finite Live IMU result') from exc


@dataclass(frozen=True)
class State:
    base:LOWER.State
    machine:PRODUCT.State
    deployment_cfg:D.DeploymentConfig
    separate_active:ACTIVE.ActiveParameters
    fma_active:ACTIVE.ActiveParameters
    live_entry_machine_updates:int
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.base,LOWER.State) or not isinstance(self.machine,PRODUCT.State):
            raise TypeError('admitted WPE/tau state and whole machine TuneState required')
        if not isinstance(self.deployment_cfg,D.DeploymentConfig): raise TypeError('DeploymentConfig required')
        if not isinstance(self.separate_active,ACTIVE.ActiveParameters) or not isinstance(self.fma_active,ACTIVE.ActiveParameters):
            raise TypeError('both global compiler applied ActiveParameters required')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong admitted machine TuneState qualification')
        if not isinstance(self.live_entry_machine_updates,int) or self.live_entry_machine_updates<0:
            raise ValueError('nonnegative Live-entry machine update count required')
        if self.machine.tau!=self.base.base.tau:
            raise ValueError('whole machine tau ledger detached from admitted WPE/tau product')
        if self.machine.tau.updates<self.live_entry_machine_updates:
            raise ValueError('whole machine update count moved backward')
        entry=_entry_live(self.base)
        if self.machine.pending!=entry.tuner.pending:
            raise ValueError('whole machine pending bit detached from exact Live tuner state')
        runtime=self.base.base.prefix.prefix.live.live_word.runtime
        if not COMMON._configs_match(runtime.candidate_cfg,self.deployment_cfg):
            raise ValueError('deployment tuner config detached from carried shipping runtime scalars')
        # These joins are state invariants, not a claim that the supplies are
        # already small enough for storage.  They also enforce the same Live
        # R_S activation branch on exact and machine applied states.
        ACTIVEJOIN.join(entry.active,self.separate_active,'separate')
        ACTIVEJOIN.join(entry.active,self.fma_active,'fma')

    @property
    def separate_active_join(self): return ACTIVEJOIN.join(_entry_live(self.base).active,self.separate_active,'separate')
    @property
    def fma_active_join(self): return ACTIVEJOIN.join(_entry_live(self.base).active,self.fma_active,'fma')


@dataclass(frozen=True)
class ImuResult:
    state:State
    lower:LOWER.ImuResult
    machine_boundary:MBOUND.Boundary
    machine_candidate:CANDSTEP.Result
    separate_sigma_join:SIGJOIN.Join
    fma_sigma_join:SIGJOIN.Join
    separate_active_join:ACTIVEJOIN.Join
    fma_active_join:ACTIVEJOIN.Join


@dataclass(frozen=True)
class CompleteWord:
    state:State
    lower:LOWER.CompleteWord
    def __post_init__(self):
        if not isinstance(self.state,State) or not isinstance(self.lower,LOWER.CompleteWord):
            raise TypeError('whole-machine complete word requires strong state and lower complete word')
        if self.lower.state!=self.state.base:
            raise ValueError('whole-machine complete word detached from WPE/tau complete word')
        if self.state.machine.tau.updates-self.state.live_entry_machine_updates!=SOURCE.TRANSITIONS:
            raise ValueError('whole machine TuneState did not advance exactly 600 Live IMU candidates')
        if not (self.state.machine.tau.updates==self.state.machine.sigma.updates==self.state.machine.rs.updates):
            raise ValueError('whole machine scalar ledgers lost common update count')


def begin(base:LOWER.State,machine:PRODUCT.State,deployment_cfg:D.DeploymentConfig,
          *,separate_active:ACTIVE.ActiveParameters,fma_active:ACTIVE.ActiveParameters):
    if not isinstance(machine,PRODUCT.State): raise TypeError('whole machine TuneState required')
    return State(base,machine,deployment_cfg,separate_active,fma_active,machine.tau.updates)


def begin_from_goLive(base:LOWER.State,go:GO.Result,deployment_cfg:D.DeploymentConfig):
    """Bind admitted Live state to exact whole-machine startup candidate+applied histories."""
    if not isinstance(go,GO.Result): raise TypeError('whole-machine goLive result required')
    try: embedded=base.base.prefix.prefix.live.live_word.live.live.live
    except AttributeError as exc: raise TypeError('admitted WPE/tau Live state lost goLive frontend/filter state') from exc
    if embedded!=go.live.state:
        raise ValueError('admitted Live product detached from exact whole-machine goLive state')
    if base.wpe!=go.wpe or base.base.tau!=go.machine.tau:
        raise ValueError('admitted Live WPE/tau ledgers detached from whole-machine goLive histories')
    return State(base,go.machine,deployment_cfg,go.separate_active,go.fma_active,go.machine.tau.updates)


def imu_step(state:State,*,
             separate_sigma_machine:SIGM.Target,fma_sigma_machine:SIGM.Target,
             separate_boundary_rs_sqrt_scale=None,fma_boundary_rs_sqrt_scale=None,
             separate_spectral_pow,separate_spectral_sqrt,
             fma_spectral_pow,fma_spectral_sqrt,
             separate_rs_exp_decay,fma_rs_exp_decay,
             **kwargs):
    if not isinstance(state,State): raise TypeError('admitted whole-machine State required')
    sigma_wave_sqrt=kwargs.get('sigma_wave_sqrt')
    if sigma_wave_sqrt is None: raise TypeError('same-event exact sigma-wave sqrt witness required')

    lower=LOWER.imu_step(state.base,**kwargs)
    live=_live_result(lower)
    exact_boundary=live.boundary
    runtime=state.base.base.prefix.prefix.live.live_word.runtime
    nf=(runtime.boundary_bench_noise_sigma if exact_boundary.band_noise_floor_sigma is None
        else exact_boundary.band_noise_floor_sigma)
    mb=MBOUND.imu_boundary(state.machine,runtime.commit_cfg,live=True,
        band_noise_floor_sigma=nf,
        separate_rs_sqrt_scale=separate_boundary_rs_sqrt_scale,
        fma_rs_sqrt_scale=fma_boundary_rs_sqrt_scale)
    if mb.consumed:
        sep_active=ACTIVE.ActiveParameters.from_commit(mb.separate_commit)
        fma_active=ACTIVE.ActiveParameters.from_commit(mb.fma_commit)
    else:
        sep_active=state.separate_active; fma_active=state.fma_active

    suffix=live.tuner_suffix; cand=suffix.candidate
    if not isinstance(cand,C.CandidateResult):
        raise TypeError('Live admitted event lost exact tuner candidate')
    exact_sample=C.sample_from_runtime(suffix.band,suffix.stillness,
                                       sigma_wave_sqrt=sigma_wave_sqrt)
    sj=SIGJOIN.join(exact_sample,runtime.candidate_cfg,state.deployment_cfg,separate_sigma_machine)
    fj=SIGJOIN.join(exact_sample,runtime.candidate_cfg,state.deployment_cfg,fma_sigma_machine)
    mc=CANDSTEP.step(mb.state,cand,runtime.candidate_cfg,state.deployment_cfg,
        lower.tau_step,dt=TAUJOIN.LEDGER.DT,
        separate_sigma_join=sj,fma_sigma_join=fj,
        separate_spectral_pow=separate_spectral_pow,separate_spectral_sqrt=separate_spectral_sqrt,
        fma_spectral_pow=fma_spectral_pow,fma_spectral_sqrt=fma_spectral_sqrt,
        separate_rs_exp_decay=separate_rs_exp_decay,fma_rs_exp_decay=fma_rs_exp_decay)
    nxt=State(lower.state,mc.product.state,state.deployment_cfg,sep_active,fma_active,
              state.live_entry_machine_updates)
    return ImuResult(nxt,lower,mb,mc,sj,fj,nxt.separate_active_join,nxt.fma_active_join)


def mag_step(state:State,**kwargs):
    if not isinstance(state,State): raise TypeError('admitted whole-machine State required')
    base,event=LOWER.mag_step(state.base,**kwargs)
    machine=MBOUND.mag_or_hold(state.machine)
    return State(base,machine,state.deployment_cfg,state.separate_active,state.fma_active,
                 state.live_entry_machine_updates),event


def set_hold(state:State,*,hold):
    if not isinstance(state,State): raise TypeError('admitted whole-machine State required')
    base,event=LOWER.set_hold(state.base,hold=hold)
    machine=MBOUND.mag_or_hold(state.machine)
    return State(base,machine,state.deployment_cfg,state.separate_active,state.fma_active,
                 state.live_entry_machine_updates),event


def complete(state:State): return CompleteWord(state,LOWER.complete(state.base))


def readiness():
    low=LOWER.readiness(); cand=CANDSTEP.readiness(); bound=MBOUND.readiness(); go=GO.readiness(); aj=ACTIVEJOIN.readiness()
    return {
      'admitted_BRMM_BIAS_ISS_WPE_tau_product_consumed':True,
      'whole_tau_sigma_RS_machine_TuneState_carried_in_same_Live_product':True,
      'candidate_memory_and_applied_machine_parameters_carried_separately':aj['candidate_state_not_confused_with_applied_active_state'],
      'exact_and_machine_pending_boundary_share_same_IMU_event':bound['next_boundary_common_machine_commit_attached'],
      'machine_applied_parameters_update_only_when_boundary_consumes_pending':True,
      'same_exact_tuner_suffix_anchors_sigma_and_RS_machine_candidate':cand['one_exact_frontend_candidate_anchors_both_global_compiler_histories'],
      'exact_vs_machine_applied_parameter_supplies_exposed_every_Live_state':aj['tau_stationary_Sigma_pseudo_period_and_RS_displacements_exposed'],
      'MAG_and_HOLD_preserve_WPE_TuneState_and_applied_machine_parameters_by_identity':True,
      'complete_word_requires_600_common_tau_sigma_RS_machine_updates':True,
      'whole_machine_goLive_provenance_constructor_available':go['goLive_carries_whole_machine_TuneState_product'],
      'Live_600_step_machine_TuneState_product_attached':True,
      'WPE_log_libm_correspondence_closed': bool(low['WPE_log_std_log_target_libm_correspondence_closed'] and low['WPE_log_exp_target_libm_correspondence_closed']),
      'upstream_frontend_binary32_correspondence_closed':False,
      'all_target_libm_correspondence_closed':False,
      'machine_active_parameter_displacement_injected_into_Live_coefficients':False,
      'source_uniform_machine_supply_bounds_closed':False,
      'all_event_arithmetic_witnesses_source_uniformly_qualified':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,'ALT_STARTUP_PASS':False,'ALT_END_TO_END_PASS':False,
    }
