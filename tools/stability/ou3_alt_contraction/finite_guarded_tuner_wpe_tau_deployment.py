"""Guarded startup frontend with coherent dual WPE-log and tau machine tracks.

This wraps the literal guarded frontend rather than the older equality-based tau
wrapper.  The exact frontend is the shadow/source graph.  Two global machine
histories are carried from construction:

  separate WPE log -> separate tuner frequency -> separate tau EMA
  FMA WPE log      -> FMA tuner frequency      -> FMA tau EMA.

Cold samples execute the guarded frontend and may advance WPE, but consume no
tau candidate operands.  Post-Cold samples use the WPE state carried into that
sample, advance each matching tau track when the exact frontend executes a
candidate, and only then advance the current sample's WPE log state.  This is
the shipping order.

The layer retains exact-vs-machine frequency/target/decay discrepancies per
mode.  It does not bound WPE/libm supplies and does not prove startup
reachability; those gates remain fail-closed.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_guarded_tuner_prefix as FRONT
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as CAND
from tools.stability.ou3_alt_contraction import finite_tuner_tau_deployment_ledger as TAU
from tools.stability.ou3_alt_contraction import finite_wpe_log_binary32 as WPELOG
from tools.stability.ou3_alt_contraction import finite_wpe_frequency_binary32 as WPEF


@dataclass(frozen=True)
class State:
    frontend:FRONT.State
    tau:TAU.State
    wpe:WPELOG.State
    def __post_init__(self):
        if not isinstance(self.frontend,FRONT.State) or not isinstance(self.tau,TAU.State) or not isinstance(self.wpe,WPELOG.State):
            raise TypeError('guarded frontend, tau ledger and WPE log ledger required')
        exact=self.frontend.tuner.wpe; initialized=exact.log_period is not None
        if initialized != (self.wpe.separate.log_period is not None) or initialized != (self.wpe.fma.log_period is not None):
            raise ValueError('startup machine WPE log initialization detached from exact frontend WPE state')
        if self.tau.updates>self.frontend.tuner.sample_index:
            raise ValueError('startup tau updates cannot exceed physical frontend samples')
        if self.wpe.samples>self.frontend.tuner.sample_index:
            raise ValueError('startup WPE machine samples cannot exceed physical frontend samples')


@dataclass(frozen=True)
class ModeSupply:
    frequency:F
    target_input:F
    decay_input:F
    def __post_init__(self):
        for n in ('frequency','target_input','decay_input'): object.__setattr__(self,n,F(getattr(self,n)))


@dataclass(frozen=True)
class Result:
    state:State
    frontend:FRONT.Result
    tau_step:TAU.StepResult|None
    wpe_step:WPELOG.StepResult
    separate_supply:ModeSupply|None
    fma_supply:ModeSupply|None


def initial(frontend:FRONT.State):
    if not isinstance(frontend,FRONT.State): raise TypeError('guarded frontend State required')
    if F(frontend.tuner.tune.tau_applied)!=TAU.INITIAL:
        raise ValueError('startup exact tuner seed detached from compiled 1.1f tau seed')
    if frontend.tuner.sample_index!=0:
        raise ValueError('startup machine ledgers must root before frontend samples')
    exact=frontend.tuner.wpe
    if exact.log_period is not None:
        raise ValueError('construction-rooted WPE must be uninitialized')
    return State(frontend,TAU.initial(),WPELOG.initial())


def _frequency_sources(state:State,cfg:CAND.CandidateConfig,exact_frequency,*,separate_getter,fma_getter):
    exact=state.frontend.tuner.wpe; lo,hi=F(cfg.min_freq),F(cfg.max_freq)
    if exact.usable_period:
        if exact_frequency is None: raise TypeError('usable startup WPE requires exact shadow frequency')
        if state.wpe.separate.log_period is None or state.wpe.fma.log_period is None:
            raise ValueError('usable startup WPE requires initialized machine log tracks')
        sl=WPEF.StoredLogPeriod(F(exact.log_period),state.wpe.separate.log_period,state.wpe.separate.log_period-F(exact.log_period))
        fl=WPEF.StoredLogPeriod(F(exact.log_period),state.wpe.fma.log_period,state.wpe.fma.log_period-F(exact.log_period))
        if not isinstance(separate_getter,WPEF.GetterResult) or separate_getter.log!=sl:
            raise ValueError('startup separate getter detached from separate WPE log track')
        if not isinstance(fma_getter,WPEF.GetterResult) or fma_getter.log!=fl:
            raise ValueError('startup FMA getter detached from FMA WPE log track')
        return (WPEF.tuner_frequency(exact,min_hz=lo,max_hz=hi,getter=separate_getter,shadow_frequency=exact_frequency),
                WPEF.tuner_frequency(exact,min_hz=lo,max_hz=hi,getter=fma_getter,shadow_frequency=exact_frequency))
    if separate_getter is not None or fma_getter is not None or exact_frequency is not None:
        raise ValueError('preusable startup WPE consumes no getter/shadow-frequency witness')
    return (WPEF.tuner_frequency(exact,min_hz=lo,max_hz=hi),WPEF.tuner_frequency(exact,min_hz=lo,max_hz=hi))


def _advance_wpe(state:State,out:FRONT.Result,*,separate_log_witness,fma_log_witness):
    w=out.tuner.wpe
    if not w.produced_period:
        if separate_log_witness is not None or fma_log_witness is not None:
            raise ValueError('nonproducing startup WPE branch consumes no machine log witness')
        return WPELOG.hold_sample(state.wpe)
    if state.wpe.separate.log_period is None:
        if not isinstance(separate_log_witness,WPELOG.InitWitness) or not isinstance(fma_log_witness,WPELOG.InitWitness):
            raise TypeError('first valid startup WPE period requires mode-specific std::log witnesses')
        return WPELOG.first_valid(state.wpe,separate=separate_log_witness,fma=fma_log_witness)
    if not isinstance(separate_log_witness,WPELOG.SmoothWitness) or not isinstance(fma_log_witness,WPELOG.SmoothWitness):
        raise TypeError('smoothed startup WPE period requires mode-specific log/exp witnesses')
    return WPELOG.smooth_valid(state.wpe,separate=separate_log_witness,fma=fma_log_witness)


def step(state:State,raw,*,dt,separate_getter=None,fma_getter=None,
         separate_tau_exp_decay=None,fma_tau_exp_decay=None,
         separate_log_witness=None,fma_log_witness=None,**kwargs):
    if not isinstance(state,State): raise TypeError('startup WPE/tau State required')
    cfg=kwargs.get('candidate_cfg')
    # Execute the exact guarded frontend first; its candidate is nevertheless
    # sourced from the WPE entry state because the lower recurrence enforces the
    # shipping preupdate ordering.
    out=FRONT.step(state.frontend,raw,dt=dt,**kwargs); cand=out.tuner.candidate
    if cand is None:
        if any(x is not None for x in (separate_getter,fma_getter,separate_tau_exp_decay,fma_tau_exp_decay)):
            raise ValueError('Cold/noncandidate startup branch consumes no tau deployment witnesses')
        tau=state.tau; ss=fs=None; tstep=None
    else:
        if not isinstance(cfg,CAND.CandidateConfig): raise TypeError('post-Cold startup candidate requires CandidateConfig')
        TAU.qualify_config(cfg)
        ema=kwargs.get('ema')
        if not isinstance(ema,CAND.EmaWitness): raise TypeError('post-Cold startup candidate requires exact EmaWitness')
        exact_freq=F(cand.frequency)
        sf,ff=_frequency_sources(state,cfg,exact_freq if state.frontend.tuner.wpe.usable_period else None,
                                 separate_getter=separate_getter,fma_getter=fma_getter)
        if exact_freq!=sf.exact_clamped_frequency or exact_freq!=ff.exact_clamped_frequency:
            raise ValueError('startup exact candidate detached from exact WPE/prior frequency')
        if separate_tau_exp_decay is None or fma_tau_exp_decay is None:
            raise TypeError('post-Cold startup candidate requires mode-specific tau exp witnesses')
        tstep=TAU.step_tracks(state.tau,separate_frequency=sf.stored,fma_frequency=ff.stored,cfg=cfg,dt=dt,
                              separate_exp_decay=separate_tau_exp_decay,fma_exp_decay=fma_tau_exp_decay)
        tau=tstep.state
        ss=ModeSupply(sf.machine_minus_shadow,tstep.separate_target.exact_target-F(cand.tau_target),F(separate_tau_exp_decay)-F(ema.decay_tau_sigma))
        fs=ModeSupply(ff.machine_minus_shadow,tstep.fma_target.exact_target-F(cand.tau_target),F(fma_tau_exp_decay)-F(ema.decay_tau_sigma))
    wstep=_advance_wpe(state,out,separate_log_witness=separate_log_witness,fma_log_witness=fma_log_witness)
    return Result(State(out.state,tau,wstep.state),out,tstep,wstep,ss,fs)


def readiness():
    w=WPELOG.readiness(); t=TAU.readiness(); f=FRONT.readiness()
    return {
      'guarded_startup_frontend_reused_without_reimplementation':f['raw_accel_guard_state_persists_across_samples'],
      'construction_seed_roots_both_tau_and_WPE_machine_ledgers':True,
      'Cold_samples_preserve_tau_but_advance_WPE_machine_sample_history':True,
      'postCold_sample_entry_WPE_frequency_precedes_tau_candidate_and_current_WPE_update':True,
      'startup_separate_WPE_frequency_advances_only_separate_tau_track':True,
      'startup_FMA_WPE_frequency_advances_only_FMA_tau_track':True,
      'startup_global_compiler_track_coherence_includes_WPE_log_and_tau_states': bool(
          w['global_compiler_track_coherence_includes_WPE_log_state'] and
          t['global_compiler_tracks_accept_mode_coherent_distinct_WPE_frequencies']),
      'startup_exact_vs_machine_frequency_target_decay_supplies_retained_per_mode':True,
      'WPE_log_std_log_target_libm_correspondence_closed':False,
      'WPE_log_exp_target_libm_correspondence_closed':False,
      'tuner_exp_libm_binary32_correspondence_closed':False,
      'every_admitted_startup_history_reaches_TunerReady_with_this_product':False,
      'goLive_carries_this_exact_WPE_tau_product':False,
      'storage_search_allowed':False,
      'ALT_STARTUP_PASS':False,'ALT_LIVE_PASS':False,
    }
