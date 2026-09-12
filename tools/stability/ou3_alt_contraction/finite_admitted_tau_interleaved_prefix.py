"""Admitted BRMM/BIAS/ISS Live interleaver with persistent binary32 tau ledger.

This is the theorem-facing Live attachment for the tau deployment recurrence.
The lower prefix already owns admitted BRMM, BIAS, bounded IMU ISS, one t_L
origin and IMU/MAG/HOLD ordering. This layer adds a persistent dual-compiler
tau ledger.

Every Live IMU edge must provide the actual stored tuner-frequency object and
the binary32 std::exp result used by the tau EMA. After executing the existing
same-history event, we require the exact tuner candidate frequency, target and
EMA decay to be the same values consumed by the deployment ledger. The separate
and FMA machine tracks then advance from their own persistent predecessors. MAG
and HOLD preserve the ledger exactly.

``imu_step_from_wpe`` is the strongest frequency-source entry. It requires the
binary32 tuner frequency to be derived from the exact WPE state carried into
THIS sample: literal 0.2f prior before the usable latch, or the separate
``std::exp(-log_period)`` getter result afterwards. Thus a free StoredFrequency
can remain available for component algebra without counting as theorem-source
provenance.

``begin_from_goLive`` is the strong construction path. It accepts the tau ledger
carried across the startup goLive edge only when the lower admitted product
contains the exact same ``finite_live_imu_prefix.State``. The generic ``begin``
remains available for component algebra, but readiness does not treat that path
as startup provenance.

Binary32 WPE log-state production and target-libm correctness remain open, so
this is not deployment closure or storage authorization.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_admitted_iss_interleaved_prefix as BASE
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as CAND
from tools.stability.ou3_alt_contraction import finite_tuner_frequency_binary32 as FREQ
from tools.stability.ou3_alt_contraction import finite_tuner_tau_deployment_ledger as LEDGER
from tools.stability.ou3_alt_contraction import finite_shipping_tau_target_binary32 as TARGET
from tools.stability.ou3_alt_contraction import finite_complete_word_tau_qualification as QUAL
from tools.stability.ou3_alt_contraction import finite_startup_live_tau_bridge as GO
from tools.stability.ou3_alt_contraction import finite_source_continuation as SOURCE
from tools.stability.ou3_alt_contraction import finite_wpe_frequency_binary32 as WPEF


@dataclass(frozen=True)
class State:
    prefix:BASE.State
    tau:LEDGER.State
    live_entry_tau_updates:int
    def __post_init__(self):
        if not isinstance(self.prefix,BASE.State) or not isinstance(self.tau,LEDGER.State):
            raise TypeError('admitted ISS prefix and tau deployment ledger required')
        if not isinstance(self.live_entry_tau_updates,int) or self.live_entry_tau_updates<0:
            raise ValueError('nonnegative Live-entry tau update count required')
        if self.tau.updates < self.live_entry_tau_updates:
            raise ValueError('tau ledger update count moved backward after Live entry')
        if self.live_entry_tau_updates > LEDGER.MAX_UPDATES-SOURCE.TRANSITIONS:
            raise ValueError('Live-entry tau ledger leaves no certified room for 600-step word')
        QUAL.qualify_runtime(self.prefix.prefix.live.live_word.runtime)

    @property
    def imu_steps(self): return self.prefix.imu_steps


@dataclass(frozen=True)
class ImuResult:
    state:State
    event:object
    tau_step:LEDGER.StepResult


@dataclass(frozen=True)
class CompleteWord:
    state:State
    lower:BASE.CompleteWord
    def __post_init__(self):
        if not isinstance(self.state,State) or not isinstance(self.lower,BASE.CompleteWord):
            raise TypeError('tau complete word requires strong state and lower complete word')
        if self.lower.state != self.state.prefix:
            raise ValueError('tau complete word detached from lower admitted complete word')
        if self.state.tau.updates-self.state.live_entry_tau_updates != SOURCE.TRANSITIONS:
            raise ValueError('tau complete word did not advance exactly once per 600 Live IMU events')


def begin(prefix:BASE.State,tau:LEDGER.State):
    """Generic component constructor; not startup provenance evidence."""
    if not isinstance(tau,LEDGER.State): raise TypeError('startup-carried tau ledger required')
    return State(prefix,tau,tau.updates)


def begin_from_goLive(prefix:BASE.State,go:GO.Result):
    """Strong Live constructor consuming the exact startup/goLive tau product."""
    if not isinstance(prefix,BASE.State) or not isinstance(go,GO.Result):
        raise TypeError('admitted ISS prefix and goLive tau bridge result required')
    try:
        embedded=prefix.prefix.live.live_word.live.live.live
    except AttributeError as exc:
        raise TypeError('admitted Live product lost embedded goLive IMU state') from exc
    if embedded != go.live.state:
        raise ValueError('admitted Live product detached from exact goLive frontend/filter state')
    return State(prefix,go.tau,go.tau.updates)


def _candidate(lower_out):
    try:
        cand=lower_out.event.event.event.live.tuner_suffix.candidate
    except AttributeError as exc:
        raise TypeError('strong admitted IMU result lost executed Live tuner candidate') from exc
    if not isinstance(cand,CAND.CandidateResult):
        raise ValueError('Live IMU edge did not execute tuner candidate recurrence')
    return cand


def _entry_wpe(state:State):
    try:
        return state.prefix.prefix.live.live_word.live.live.live.tuner.wpe
    except AttributeError as exc:
        raise TypeError('strong admitted tau product lost sample-entry WPE state') from exc


def imu_step(state:State,*,stored_frequency:FREQ.StoredFrequency,tau_exp_decay,**kwargs):
    """Component entry with an already-derived StoredFrequency."""
    if not isinstance(state,State): raise TypeError('admitted tau interleaved State required')
    if not isinstance(stored_frequency,FREQ.StoredFrequency): raise TypeError('actual StoredFrequency required')
    ema=kwargs.get('ema')
    if not isinstance(ema,CAND.EmaWitness): raise TypeError('same-event exact tuner EmaWitness required')
    if F(ema.decay_tau_sigma)!=F(tau_exp_decay):
        raise ValueError('binary32 tau exp witness detached from exact candidate EMA decay')
    nxt_prefix,out=BASE.imu_step(state.prefix,**kwargs)
    cand=_candidate(out); target=TARGET.evaluate(stored_frequency.stored_hz)
    if F(cand.frequency)!=F(stored_frequency.stored_hz):
        raise ValueError('Live candidate frequency detached from stored binary32 tuner frequency')
    if F(cand.tau_target)!=F(target.exact_target):
        raise ValueError('Live candidate tau target detached from source-locked machine-operand target')
    cfg=state.prefix.prefix.live.live_word.runtime.candidate_cfg
    tau_step=LEDGER.step(state.tau,stored_frequency,cfg,dt=LEDGER.DT,exp_decay=tau_exp_decay)
    return ImuResult(State(nxt_prefix,tau_step.state,state.live_entry_tau_updates),out,tau_step)


def imu_step_from_wpe(state:State,*,tuner_frequency:WPEF.TunerFrequencyResult,
                      tau_exp_decay,**kwargs):
    """Strong theorem entry deriving the tau frequency from THIS sample's WPE entry state."""
    if not isinstance(state,State) or not isinstance(tuner_frequency,WPEF.TunerFrequencyResult):
        raise TypeError('admitted tau state and WPE tuner-frequency source required')
    if tuner_frequency.shadow != _entry_wpe(state):
        raise ValueError('binary32 tuner frequency detached from carried sample-entry WPE state')
    return imu_step(state,stored_frequency=tuner_frequency.stored,
                    tau_exp_decay=tau_exp_decay,**kwargs)


def mag_step(state:State,**kwargs):
    if not isinstance(state,State): raise TypeError('admitted tau interleaved State required')
    nxt,event=BASE.mag_step(state.prefix,**kwargs)
    return State(nxt,state.tau,state.live_entry_tau_updates),event


def set_hold(state:State,*,hold):
    if not isinstance(state,State): raise TypeError('admitted tau interleaved State required')
    nxt,event=BASE.set_hold(state.prefix,hold=hold)
    return State(nxt,state.tau,state.live_entry_tau_updates),event


def complete(state:State):
    return CompleteWord(state,BASE.complete(state.prefix))


def readiness():
    ledger=LEDGER.readiness(); qual=QUAL.readiness(); go=GO.readiness(); wf=WPEF.readiness()
    return {
      'admitted_BRMM_BIAS_ISS_interleaver_consumed':True,
      'tau_relevant_persistent_runtime_config_source_qualified':qual['tau_relevant_persistent_runtime_config_source_qualified'],
      'Live_product_carries_persistent_dual_compiler_tau_ledger':True,
      'strong_Live_constructor_requires_exact_goLive_filter_frontend_state_and_tau_ledger':True,
      'goLive_tau_ledger_identity_bridge_available':go['goLive_preserves_binary32_tau_ledger_by_identity'],
      'each_Live_IMU_requires_same_candidate_frequency_target_and_decay_as_tau_ledger':True,
      'strong_Live_IMU_frequency_source_bound_to_sample_entry_WPE_or_prior':True,
      'WPE_getter_to_tuner_store_topology_available':wf['WPE_frequency_getter_to_tuner_binary32_store_topology_closed'],
      'each_Live_IMU_retains_both_binary32_tau_roundoff_certificates':True,
      'MAG_and_HOLD_preserve_tau_ledger_exactly':True,
      'complete_word_requires_exactly_600_tau_updates_after_Live_entry':True,
      'admitted_startup_reachability_with_tau_ledger_closed':False,
      'WPE_binary32_log_period_production_closed':wf['WPE_binary32_log_period_production_closed'],
      'upstream_WPE_to_StoredFrequency_binary32_correspondence_closed':False,
      'tuner_exp_libm_binary32_correspondence_closed':ledger['tuner_exp_libm_binary32_correspondence_closed'],
      'all_event_arithmetic_witnesses_source_uniformly_qualified':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
      'ALT_STARTUP_PASS':False,
      'ALT_END_TO_END_PASS':False,
    }
