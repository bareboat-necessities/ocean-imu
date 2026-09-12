"""Persistent binary32 tau deployment ledger for startup and Live ALT words.

The exact-real tuner state is useful for algebra but is not the machine state.
This ledger carries the actual binary32 ``tau_applied`` recurrence in parallel.
Because compiler FP contraction is not yet qualified, it carries two coherent
GLOBAL tracks rather than selecting one per sample:

  * ``separate`` always uses rounded multiply then rounded add;
  * ``fma`` always uses the contracted FMA shape.

The compiler choice is global, so two tracks are sufficient; there is no 2^N
branch explosion.  Importantly, the two global histories are allowed to consume
different stored frequencies and exp results.  Once WPE log-period smoothing is
modeled in binary32, its state may itself depend on the global contraction mode;
forcing both tau tracks to share one WPE frequency would splice two incompatible
machine histories.

``step_tracks`` is therefore the theorem-facing primitive.  ``step`` is only a
common-input convenience for component tests and delegates to ``step_tracks``.
Each mode derives its own source-locked tau target and retains its own uniform
roundoff certificate.  Cold/other branches that execute no ``adapt_mekf`` call
use ``hold`` and consume no arithmetic witness.

WPE binary32 log-state production and target-libm correctness remain open.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as CAND
from tools.stability.ou3_alt_contraction import finite_tuner_frequency_binary32 as FREQ
from tools.stability.ou3_alt_contraction import finite_tuner_tau_binary32 as TAU
from tools.stability.ou3_alt_contraction import finite_shipping_tau_target_binary32 as TARGET
from tools.stability.ou3_alt_contraction import finite_tuner_tau_roundoff_bound as ROUND

DT=B.rn32(F(1,200)); INITIAL=B.rn32(F(11,10)); MAX_UPDATES=30600
ADAPT_SEC=B.rn32(F(9,5)); ADAPT_PERIODS=B.rn32(F(2,5))


@dataclass(frozen=True)
class State:
    separate:F=INITIAL
    fma:F=INITIAL
    updates:int=0
    def __post_init__(self):
        a,b=F(self.separate),F(self.fma)
        object.__setattr__(self,'separate',a); object.__setattr__(self,'fma',b)
        if not B.is_binary32(a) or not B.is_binary32(b) or a<=0 or b<=0:
            raise ValueError('tau deployment ledger stores positive binary32 values only')
        if not isinstance(self.updates,int) or not 0<=self.updates<=MAX_UPDATES:
            raise ValueError('tau deployment update count outside bounded startup+word horizon')


@dataclass(frozen=True)
class StepResult:
    state:State
    separate_target:TARGET.TargetPair
    fma_target:TARGET.TargetPair
    separate_step:TAU.TauStep
    fma_step:TAU.TauStep
    separate_certificate:ROUND.BoundCertificate
    fma_certificate:ROUND.BoundCertificate
    @property
    def target(self):
        """Compatibility view only when both mode histories used one target."""
        if self.separate_target != self.fma_target:
            raise ValueError('compiler-mode tau targets differ; no common target exists')
        return self.separate_target


def initial(): return State()


def qualify_config(cfg:CAND.CandidateConfig):
    if not isinstance(cfg,CAND.CandidateConfig): raise TypeError('CandidateConfig required')
    required={'min_freq':TARGET.FLOOR,'max_freq':TARGET.CEIL,'tau_coeff':B.rn32(1),
              'min_tau':TARGET.TAU_MIN,'max_tau':TARGET.TAU_MAX,
              'adapt_tau_sec':ADAPT_SEC,'adapt_tau_sea_periods':ADAPT_PERIODS}
    for name,value in required.items():
        if F(getattr(cfg,name))!=F(value):
            raise ValueError(f'tau ledger {name} detached from shipping binary32 default')
    if cfg.clamp_enabled is not True: raise ValueError('tau ledger requires shipping clamp-enabled branch')
    return cfg


def hold(state:State):
    if not isinstance(state,State): raise TypeError('tau deployment State required')
    return state


def step_tracks(state:State,*,separate_frequency:FREQ.StoredFrequency,
                fma_frequency:FREQ.StoredFrequency,cfg:CAND.CandidateConfig,dt,
                separate_exp_decay,fma_exp_decay):
    """Advance the two coherent global compiler histories independently."""
    if not isinstance(state,State): raise TypeError('tau deployment State required')
    if state.updates>=MAX_UPDATES: raise ValueError('tau deployment ledger exceeded bounded startup+word horizon')
    qualify_config(cfg)
    if B.rn32(dt)!=DT: raise ValueError('tau deployment ledger requires canonical shipping 5 ms update')
    if not isinstance(separate_frequency,FREQ.StoredFrequency) or not isinstance(fma_frequency,FREQ.StoredFrequency):
        raise TypeError('actual stored tuner frequency required for each compiler track')
    st=TARGET.evaluate(separate_frequency.stored_hz)
    ft=TARGET.evaluate(fma_frequency.stored_hz)
    sep=TAU.step_from_stored_frequency(state.separate,separate_frequency,cfg,dt=dt,
                                       exp_decay=separate_exp_decay)
    fma=TAU.step_from_stored_frequency(state.fma,fma_frequency,cfg,dt=dt,
                                       exp_decay=fma_exp_decay)
    cs=ROUND.certify_source_target(sep,st); cf=ROUND.certify_source_target(fma,ft)
    nxt=State(sep.next_separate,fma.next_fma,state.updates+1)
    return StepResult(nxt,st,ft,sep,fma,cs,cf)


def step(state:State,stored:FREQ.StoredFrequency,cfg:CAND.CandidateConfig,*,dt,exp_decay):
    """Common-input component helper; not evidence that global WPE tracks coincide."""
    return step_tracks(state,separate_frequency=stored,fma_frequency=stored,cfg=cfg,dt=dt,
                       separate_exp_decay=exp_decay,fma_exp_decay=exp_decay)


def readiness():
    return {
      'shipping_tau_initial_binary32_seed_materialized':True,
      'global_separate_and_FMA_compiler_tracks_persist_without_branch_explosion':True,
      'global_compiler_tracks_accept_mode_coherent_distinct_WPE_frequencies':True,
      'global_compiler_tracks_accept_mode_coherent_distinct_exp_results':True,
      'common_frequency_across_compiler_tracks_assumed':False,
      'cold_or_nonadapting_identity_branch_materialized':True,
      'each_update_requires_actual_StoredFrequency_object_per_compiler_track':True,
      'each_update_retains_source_locked_target_and_roundoff_certificate_per_track':True,
      'bounded_startup_plus_600_update_cap_enforced':True,
      'global_compiler_track_coherence_includes_WPE_log_state':False,
      'upstream_WPE_to_StoredFrequency_binary32_correspondence_closed':False,
      'tuner_exp_libm_binary32_correspondence_closed':False,
      'startup_master_product_carries_tau_ledger':False,
      'Live_master_product_carries_tau_ledger':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
