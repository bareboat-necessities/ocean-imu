"""Guarded startup frontend with coherent WPE + whole machine TuneState history.

This strengthens ``finite_guarded_tuner_wpe_tau_deployment`` without replacing
it.  The lower layer continues to own raw-packet/guard/Mahony/band/stillness/WPE
ordering and the dual WPE/tau compiler histories.  This product adds persistent
sigma_applied and RS_applied histories and requires their one pending bit to be
identical to the exact frontend's pending control state.

At each post-Cold candidate, the exact lower candidate is joined to explicit
binary32 sigma operands, then the entire tau/sigma/R_S machine candidate step is
composed.  At the following IMU boundary the exact frontend boundary and the
machine boundary execute together and clear the same pending bit.  Cold samples
are literal identity on all three machine TuneState scalars while WPE, band
and statistics histories still advance. Their carried band covariance produces
boundary floors; each new band/statistics output produces the sigma operands.

Machine sigma/frontend correspondence and all libm witnesses remain explicit
open obligations; no synthetic machine TuneState is inserted at TunerReady or
Live entry.
"""
from __future__ import annotations
from dataclasses import dataclass, field

from tools.stability.ou3_alt_contraction import finite_guarded_tuner_prefix as FRONT
from tools.stability.ou3_alt_contraction import finite_guarded_tuner_wpe_tau_deployment as LOWER
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as C
from tools.stability.ou3_alt_contraction import finite_tuner_projection_bridge as PROJECTION
from tools.stability.ou3_alt_contraction import finite_tuner_deployment_config as D
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_binary32 as SIGM
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_machine_real_join as SIGJOIN
from tools.stability.ou3_alt_contraction import finite_tuner_machine_tunestate_product as PRODUCT
from tools.stability.ou3_alt_contraction import finite_tuner_machine_candidate_step as CANDSTEP
from tools.stability.ou3_alt_contraction import finite_tuner_boundary_commit as EXACTBOUND
from tools.stability.ou3_alt_contraction import finite_tuner_machine_boundary as MACHBOUND
from tools.stability.ou3_alt_contraction import finite_tuner_commit as COMMIT

from tools.stability.ou3_alt_contraction import finite_machine_frontend_sigma_source as MF

QUALIFICATION='OU3_ALT_GUARDED_MACHINE_TUNESTATE_DEPLOYMENT_V1'


@dataclass(frozen=True)
class State:
    lower:LOWER.State
    machine:PRODUCT.State
    frontends:MF.Pair=field(default_factory=MF.initial_pair)
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.lower,LOWER.State) or not isinstance(self.machine,PRODUCT.State):
            raise TypeError('guarded WPE/tau state and whole machine TuneState required')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong guarded machine TuneState qualification')
        if self.machine.tau!=self.lower.tau:
            raise ValueError('whole machine tau history detached from guarded WPE/tau history')
        if self.machine.pending!=self.lower.frontend.tuner.pending:
            raise ValueError('machine and exact frontend pending bits differ')
        if not isinstance(self.frontends,MF.Pair) or self.frontends.samples!=self.lower.frontend.tuner.sample_index:
            raise ValueError('startup machine frontend sample count detached')


@dataclass(frozen=True)
class StepResult:
    state:State
    lower:LOWER.Result
    machine:CANDSTEP.Result|None
    separate_sigma_join:SIGJOIN.Join|None
    fma_sigma_join:SIGJOIN.Join|None
    separate_frontend:MF.Result
    fma_frontend:MF.Result
    def __post_init__(self):
        if not isinstance(self.state,State) or not isinstance(self.lower,LOWER.Result):
            raise TypeError('guarded whole-TuneState step result malformed')
        if not isinstance(self.separate_frontend,MF.Result) or not isinstance(self.fma_frontend,MF.Result):
            raise TypeError('both executed machine frontend steps required')
        if self.state.lower!=self.lower.state or self.state.frontends!=MF.Pair(self.separate_frontend.state,self.fma_frontend.state):
            raise ValueError('startup successor detached from executed frontend steps')
        has=self.lower.tau_step is not None
        if has != (self.machine is not None):
            raise ValueError('machine TuneState candidate presence detached from lower candidate')
        if has:
            if not isinstance(self.separate_sigma_join,SIGJOIN.Join) or not isinstance(self.fma_sigma_join,SIGJOIN.Join):
                raise TypeError('post-Cold machine candidate requires both sigma joins')
            MF.require_sigma(self.separate_frontend,self.separate_sigma_join.machine)
            MF.require_sigma(self.fma_frontend,self.fma_sigma_join.machine)
        elif self.separate_sigma_join is not None or self.fma_sigma_join is not None:
            raise ValueError('Cold/noncandidate result cannot retain sigma joins')


@dataclass(frozen=True)
class BoundaryResult:
    state:State
    exact:EXACTBOUND.Result
    machine:MACHBOUND.Boundary
    noise_floors:tuple
    def __post_init__(self):
        if not isinstance(self.state,State) or not isinstance(self.exact,EXACTBOUND.Result) or not isinstance(self.machine,MACHBOUND.Boundary):
            raise TypeError('combined exact/machine boundary result required')
        if self.exact.state.pending or self.machine.state.pending:
            raise ValueError('combined pending boundary failed to clear both representations')
        MF.require_boundary(self.state.frontends,self.noise_floors,self.machine)


def initial(frontend:FRONT.State):
    lower=LOWER.initial(frontend)
    machine=PRODUCT.initial()
    return State(lower,machine)


def boundary(state:State,cfg:COMMIT.CommitConfig,*,bench_noise_sigma,
             exact_noise_sqrt=None,exact_rs_sqrt_scale=None,
             separate_noise_sqrt_gain=None,fma_noise_sqrt_gain=None,
             sync_covariance=False,rs_scale=1):
    """Execute one next-IMU pending transaction on exact and machine products."""
    if not isinstance(state,State): raise TypeError('guarded whole machine TuneState State required')
    exact=EXACTBOUND.apply(state.lower.frontend.tuner,cfg,
        bench_noise_sigma=bench_noise_sigma,noise_sqrt=exact_noise_sqrt,
        rs_sqrt_scale=exact_rs_sqrt_scale,sync_covariance=sync_covariance,rs_scale=rs_scale)
    live=(state.lower.frontend.tuner.stage=='Live')
    if rs_scale!=1: raise ValueError('shipping pending boundary uses literal R_S scale one')
    floors=MF.boundary_floors(state.frontends,bench_noise_sigma=bench_noise_sigma,required=state.machine.pending,
        separate_sqrt_gain=separate_noise_sqrt_gain,fma_sqrt_gain=fma_noise_sqrt_gain)
    machine=MACHBOUND.imu_boundary(state.machine,cfg,live=live,
        separate_band_noise_floor_sigma=None if floors[0] is None else floors[0].noise_sigma,
        fma_band_noise_floor_sigma=None if floors[1] is None else floors[1].noise_sigma,
        sync_covariance=sync_covariance)
    next_front=FRONT.State(state.lower.frontend.guard,exact.state)
    next_lower=LOWER.State(next_front,state.lower.tau,state.lower.wpe)
    return BoundaryResult(State(next_lower,machine.state,state.frontends),exact,machine,floors)


def step(state:State,raw,*,dt,deployment_cfg:D.DeploymentConfig,
         separate_frontend:MF.Result,fma_frontend:MF.Result,
         separate_sigma_machine:SIGM.Target|None=None,
         fma_sigma_machine:SIGM.Target|None=None,
         separate_spectral_pow=None,separate_spectral_sqrt=None,
         fma_spectral_pow=None,fma_spectral_sqrt=None,
         separate_rs_exp_decay=None,fma_rs_exp_decay=None,
         **kwargs):
    if not isinstance(state,State): raise TypeError('guarded whole machine TuneState State required')
    if state.machine.pending:
        raise ValueError('pending whole TuneState must execute combined boundary before next sample')
    out=LOWER.step(state.lower,raw,dt=dt,**kwargs)
    sep=MF.bind_step(state.frontends.separate,separate_frontend,frequency=out.separate_frequency,
        band_cfg=kwargs['band_cfg'],stats_cfg=kwargs['stats_cfg'],dt=dt,bench_noise_sigma=kwargs['bench_noise_sigma'])
    fma=MF.bind_step(state.frontends.fma,fma_frontend,frequency=out.fma_frequency,
        band_cfg=kwargs['band_cfg'],stats_cfg=kwargs['stats_cfg'],dt=dt,bench_noise_sigma=kwargs['bench_noise_sigma'])
    pair=MF.Pair(sep.state,fma.state)
    cand=out.frontend.tuner.candidate
    if cand is None:
        supplied=(separate_sigma_machine,fma_sigma_machine,separate_spectral_pow,separate_spectral_sqrt,
                  fma_spectral_pow,fma_spectral_sqrt,separate_rs_exp_decay,fma_rs_exp_decay)
        if any(x is not None for x in supplied):
            raise ValueError('Cold/noncandidate startup sample consumes no sigma/R_S machine witnesses')
        machine=PRODUCT.hold(state.machine)
        return StepResult(State(out.state,machine,pair),out,None,None,None,sep,fma)

    c=kwargs.get('candidate_cfg'); sigma_wave_sqrt=kwargs.get('sigma_wave_sqrt')
    if not isinstance(c,C.CandidateConfig) or sigma_wave_sqrt is None:
        raise TypeError('post-Cold startup machine product requires exact candidate config/sigma witness')
    if not isinstance(deployment_cfg,D.DeploymentConfig): raise TypeError('DeploymentConfig required')
    if not isinstance(separate_sigma_machine,SIGM.Target) or not isinstance(fma_sigma_machine,SIGM.Target):
        raise TypeError('post-Cold startup requires both compiler-mode binary32 sigma targets')
    MF.require_sigma(sep,separate_sigma_machine); MF.require_sigma(fma,fma_sigma_machine)
    exact_sample=PROJECTION.sample_from_projection(out.frontend.tuner.band,out.frontend.tuner.stillness,
                                       sigma_wave_sqrt=sigma_wave_sqrt)
    sj=SIGJOIN.join(exact_sample,c,deployment_cfg,separate_sigma_machine)
    fj=SIGJOIN.join(exact_sample,c,deployment_cfg,fma_sigma_machine)
    m=CANDSTEP.step(state.machine,cand,c,deployment_cfg,out.tau_step,dt=dt,
        separate_sigma_join=sj,fma_sigma_join=fj,
        separate_spectral_pow=separate_spectral_pow,separate_spectral_sqrt=separate_spectral_sqrt,
        fma_spectral_pow=fma_spectral_pow,fma_spectral_sqrt=fma_spectral_sqrt,
        separate_rs_exp_decay=separate_rs_exp_decay,fma_rs_exp_decay=fma_rs_exp_decay,
        uniform_supplies=bool(kwargs.get('machine_wpe_entry') is not None
                              and kwargs['machine_wpe_entry'].bounded_profile))
    return StepResult(State(out.state,m.product.state,pair),out,m,sj,fj,sep,fma)


def readiness():
    lower=LOWER.readiness(); m=CANDSTEP.readiness(); b=MACHBOUND.readiness()
    return {
      'guarded_startup_WPE_tau_order_reused_without_reimplementation':lower['postCold_sample_entry_WPE_frequency_precedes_tau_candidate_and_current_WPE_update'],
      'construction_roots_tau_sigma_RS_machine_product_before_sample_one':True,
      'Cold_samples_hold_whole_machine_TuneState_while_WPE_may_advance':True,
      'postCold_exact_candidate_and_machine_sigma_inputs_joined_on_same_frontend_sample':True,
      'postCold_tau_sigma_RS_machine_candidate_step_composed':m['tau_sigma_RS_scalar_ledgers_advance_once_before_pending_bit_is_attached'],
      'machine_and_exact_frontend_pending_bits_identical_in_product_state':True,
      'combined_next_IMU_boundary_clears_exact_and_machine_pending_together':b['next_boundary_common_machine_commit_attached'],
      'startup_pending_boundary_consumes_binary32_commit_graph':True,
      'machine_band_stats_persist_including_Cold_samples':True,
      'boundary_floor_derived_from_carried_pre_sample_machine_band':True,
      'machine_sigma_variance_and_noise_from_same_frontend_successor':True,
      'startup_frontend_machine_TuneState_product_attached':True,
      'upstream_frontend_binary32_correspondence_closed':False,
      'all_target_libm_correspondence_closed':False,
      'source_uniform_machine_supply_bounds_closed':False,
      'every_admitted_startup_history_reaches_TunerReady_with_this_product':False,
      'goLive_carries_whole_machine_TuneState_product':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_STARTUP_PASS':False,
      'ALT_LIVE_PASS':False,
    }
