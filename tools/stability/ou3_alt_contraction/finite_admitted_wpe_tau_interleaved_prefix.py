"""Admitted Live product with coherent WPE-log -> frequency -> tau compiler tracks.

This is the first layer that keeps the global compiler-mode history coherent
through BOTH adaptive states.  The exact admitted BRMM/BIAS/ISS word remains the
shadow/source graph.  In parallel, two machine histories are carried:

  separate: WPE log EMA separate-mul/add -> WPE frequency -> tau separate EMA
  fma:      WPE log EMA contracted FMA   -> WPE frequency -> tau FMA

One IMU event is ordered exactly as shipping:
  1. read the WPE state carried into the sample;
  2. obtain the tuner frequency (0.2f prior or exp(-logT)) per compiler track;
  3. execute the exact admitted Live event/candidate using the preupdate WPE;
  4. advance the matching tau machine tracks;
  5. advance the current sample's WPE log machine tracks for the next sample.

MAG and HOLD preserve both machine ledgers.  Target-libm correctness and a
source-uniform bound on WPE log/frequency deployment supplies remain open.
Startup provenance for the WPE machine ledger is also still open, so this layer
must not promote ALT gates or storage search.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_admitted_tau_interleaved_prefix as TAUJOIN
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as CAND
from tools.stability.ou3_alt_contraction import finite_tuner_tau_deployment_ledger as TAULEDGER
from tools.stability.ou3_alt_contraction import finite_shipping_tau_target_binary32 as TARGET
from tools.stability.ou3_alt_contraction import finite_wpe_frequency_binary32 as WPEF
from tools.stability.ou3_alt_contraction import finite_wpe_log_binary32 as WPELOG
from tools.stability.ou3_alt_contraction import finite_source_continuation as SOURCE


@dataclass(frozen=True)
class State:
    base:TAUJOIN.State
    wpe:WPELOG.State
    live_entry_wpe_samples:int
    def __post_init__(self):
        if not isinstance(self.base,TAUJOIN.State) or not isinstance(self.wpe,WPELOG.State):
            raise TypeError('admitted tau product and binary32 WPE log ledger required')
        if not isinstance(self.live_entry_wpe_samples,int) or self.live_entry_wpe_samples<0:
            raise ValueError('nonnegative Live-entry WPE sample count required')
        if self.wpe.samples < self.live_entry_wpe_samples:
            raise ValueError('WPE machine ledger sample count moved backward')
        exact=TAUJOIN._entry_wpe(self.base)
        initialized=exact.log_period is not None
        if initialized != (self.wpe.separate.log_period is not None) or initialized != (self.wpe.fma.log_period is not None):
            raise ValueError('machine WPE log initialization detached from exact carried WPE state')


@dataclass(frozen=True)
class ModeSupply:
    frequency:F
    target_input:F
    decay_input:F
    def __post_init__(self):
        object.__setattr__(self,'frequency',F(self.frequency)); object.__setattr__(self,'target_input',F(self.target_input)); object.__setattr__(self,'decay_input',F(self.decay_input))


@dataclass(frozen=True)
class ImuResult:
    state:State
    event:object
    tau_step:TAULEDGER.StepResult
    wpe_step:WPELOG.StepResult
    separate_supply:ModeSupply
    fma_supply:ModeSupply


@dataclass(frozen=True)
class CompleteWord:
    state:State
    lower:TAUJOIN.CompleteWord
    def __post_init__(self):
        if not isinstance(self.state,State) or not isinstance(self.lower,TAUJOIN.CompleteWord):
            raise TypeError('WPE/tau complete word requires strong state and lower complete word')
        if self.lower.state != self.state.base:
            raise ValueError('WPE/tau complete word detached from lower tau word')
        if self.state.wpe.samples-self.state.live_entry_wpe_samples != SOURCE.TRANSITIONS:
            raise ValueError('WPE machine ledger did not consume exactly 600 Live IMU samples')


def begin(base:TAUJOIN.State,wpe:WPELOG.State):
    """Component constructor; startup provenance of ``wpe`` is not yet proved."""
    if not isinstance(wpe,WPELOG.State): raise TypeError('binary32 WPE log State required')
    return State(base,wpe,wpe.samples)


def _suffix(lower_out):
    try: return lower_out.event.event.event.live.tuner_suffix
    except AttributeError as exc: raise TypeError('admitted IMU result lost tuner/WPE suffix') from exc


def _frequency_sources(state:State,*,separate_getter,fma_getter,shadow_frequency):
    exact=TAUJOIN._entry_wpe(state.base); runtime=state.base.prefix.prefix.live.live_word.runtime
    lo,hi=F(runtime.candidate_cfg.min_freq),F(runtime.candidate_cfg.max_freq)
    if exact.usable_period:
        if shadow_frequency is None: raise TypeError('usable WPE branch requires exact shadow frequency')
        if state.wpe.separate.log_period is None or state.wpe.fma.log_period is None:
            raise ValueError('usable WPE branch requires initialized machine log tracks')
        sl=WPEF.StoredLogPeriod(F(exact.log_period),state.wpe.separate.log_period,state.wpe.separate.log_period-F(exact.log_period))
        fl=WPEF.StoredLogPeriod(F(exact.log_period),state.wpe.fma.log_period,state.wpe.fma.log_period-F(exact.log_period))
        if not isinstance(separate_getter,WPEF.GetterResult) or separate_getter.log!=sl:
            raise ValueError('separate WPE frequency getter detached from separate compiler log track')
        if not isinstance(fma_getter,WPEF.GetterResult) or fma_getter.log!=fl:
            raise ValueError('FMA WPE frequency getter detached from FMA compiler log track')
        sf=WPEF.tuner_frequency(exact,min_hz=lo,max_hz=hi,getter=separate_getter,shadow_frequency=shadow_frequency)
        ff=WPEF.tuner_frequency(exact,min_hz=lo,max_hz=hi,getter=fma_getter,shadow_frequency=shadow_frequency)
        return sf,ff
    if separate_getter is not None or fma_getter is not None or shadow_frequency is not None:
        raise ValueError('preusable WPE branch consumes no getter/shadow frequency')
    return (WPEF.tuner_frequency(exact,min_hz=lo,max_hz=hi),
            WPEF.tuner_frequency(exact,min_hz=lo,max_hz=hi))


def _advance_wpe_machine(state:State,suffix,*,separate_log_witness,fma_log_witness):
    w=suffix.wpe
    if not w.produced_period:
        if separate_log_witness is not None or fma_log_witness is not None:
            raise ValueError('nonproducing WPE branch consumes no machine log witness')
        return WPELOG.hold_sample(state.wpe)
    if state.wpe.separate.log_period is None:
        if not isinstance(separate_log_witness,WPELOG.InitWitness) or not isinstance(fma_log_witness,WPELOG.InitWitness):
            raise TypeError('first valid WPE period requires mode-specific std::log witnesses')
        return WPELOG.first_valid(state.wpe,separate=separate_log_witness,fma=fma_log_witness)
    if not isinstance(separate_log_witness,WPELOG.SmoothWitness) or not isinstance(fma_log_witness,WPELOG.SmoothWitness):
        raise TypeError('smoothed WPE period requires mode-specific log/exp witnesses')
    return WPELOG.smooth_valid(state.wpe,separate=separate_log_witness,fma=fma_log_witness)


def imu_step(state:State,*,separate_getter=None,fma_getter=None,shadow_frequency=None,
             separate_tau_exp_decay,fma_tau_exp_decay,
             separate_log_witness=None,fma_log_witness=None,**kwargs):
    if not isinstance(state,State): raise TypeError('admitted WPE/tau State required')
    ema=kwargs.get('ema')
    if not isinstance(ema,CAND.EmaWitness): raise TypeError('same-event exact tuner EmaWitness required')
    sf,ff=_frequency_sources(state,separate_getter=separate_getter,fma_getter=fma_getter,shadow_frequency=shadow_frequency)
    nxt_prefix,out=TAUJOIN.BASE.imu_step(state.base.prefix,**kwargs)
    cand=TAUJOIN._candidate(out)
    if F(cand.frequency)!=sf.exact_clamped_frequency or F(cand.frequency)!=ff.exact_clamped_frequency:
        raise ValueError('exact tuner candidate detached from exact WPE shadow frequency')
    cfg=state.base.prefix.prefix.live.live_word.runtime.candidate_cfg
    tau=TAULEDGER.step_tracks(state.base.tau,separate_frequency=sf.stored,fma_frequency=ff.stored,cfg=cfg,dt=TAULEDGER.DT,
                              separate_exp_decay=separate_tau_exp_decay,fma_exp_decay=fma_tau_exp_decay)
    ss=ModeSupply(sf.machine_minus_shadow,tau.separate_target.exact_target-F(cand.tau_target),F(separate_tau_exp_decay)-F(ema.decay_tau_sigma))
    fs=ModeSupply(ff.machine_minus_shadow,tau.fma_target.exact_target-F(cand.tau_target),F(fma_tau_exp_decay)-F(ema.decay_tau_sigma))
    suffix=_suffix(out)
    wstep=_advance_wpe_machine(state,suffix,separate_log_witness=separate_log_witness,fma_log_witness=fma_log_witness)
    base_next=TAUJOIN.State(nxt_prefix,tau.state,state.base.live_entry_tau_updates)
    return ImuResult(State(base_next,wstep.state,state.live_entry_wpe_samples),out,tau,wstep,ss,fs)


def mag_step(state:State,**kwargs):
    if not isinstance(state,State): raise TypeError('admitted WPE/tau State required')
    base,event=TAUJOIN.mag_step(state.base,**kwargs)
    return State(base,state.wpe,state.live_entry_wpe_samples),event


def set_hold(state:State,*,hold):
    if not isinstance(state,State): raise TypeError('admitted WPE/tau State required')
    base,event=TAUJOIN.set_hold(state.base,hold=hold)
    return State(base,state.wpe,state.live_entry_wpe_samples),event


def complete(state:State): return CompleteWord(state,TAUJOIN.complete(state.base))


def readiness():
    wl=WPELOG.readiness(); tl=TAULEDGER.readiness(); lower=TAUJOIN.readiness()
    return {
      'admitted_BRMM_BIAS_ISS_tau_product_consumed':True,
      'global_separate_and_FMA_WPE_log_tracks_carried_in_Live_product':True,
      'sample_entry_WPE_frequency_precedes_current_sample_WPE_update':True,
      'separate_WPE_frequency_advances_only_separate_tau_track':True,
      'FMA_WPE_frequency_advances_only_FMA_tau_track':True,
      'global_compiler_track_coherence_includes_WPE_log_and_tau_states': bool(
          wl['global_compiler_track_coherence_includes_WPE_log_state'] and
          tl['global_compiler_tracks_accept_mode_coherent_distinct_WPE_frequencies']),
      'exact_vs_machine_frequency_target_decay_supplies_retained_per_mode':True,
      'MAG_and_HOLD_preserve_both_WPE_and_tau_ledgers':True,
      'complete_word_requires_600_WPE_samples_and_600_tau_updates':True,
      'startup_WPE_machine_ledger_provenance_closed':False,
      'WPE_log_std_log_target_libm_correspondence_closed':wl['WPE_log_std_log_target_libm_correspondence_closed'],
      'WPE_log_exp_target_libm_correspondence_closed':wl['WPE_log_exp_target_libm_correspondence_closed'],
      'source_uniform_WPE_frequency_supply_bound_closed':False,
      'tuner_exp_libm_binary32_correspondence_closed':lower['tuner_exp_libm_binary32_correspondence_closed'],
      'all_event_arithmetic_witnesses_source_uniformly_qualified':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,'ALT_STARTUP_PASS':False,'ALT_END_TO_END_PASS':False,
    }
