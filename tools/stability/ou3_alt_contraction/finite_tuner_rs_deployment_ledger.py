"""Persistent dual-compiler deployment ledger for shipping ``RS_applied``.

The global ALT deployment word already carries two coherent compiler histories
for WPE/tau.  The SpectralMSE R_S channel must not collapse them back together:
different upstream histories may produce different tau/sigma targets, different
pow/sqrt inputs/results, and therefore different alpha_RS values.

This ledger carries TWO persistent stored binary32 R_S values:

* ``separate`` advances only with the separate multiply/add EMA successor;
* ``fma`` advances only with the contracted-FMA EMA successor.

Each update consumes its OWN same-source SpectralMSE machine/real target join
and its OWN same-source alpha join.  Both must share the same deployment config,
but equality of target or alpha values across compiler histories is deliberately
not required.

Shipping seed is source-locked at ``TuneState::RS_applied = 0.5f``.  A hold
sample (Cold branch/no tuner candidate) is literal identity.  Pending commit
snapshotting is intentionally separate: this object proves persistent candidate
state, not yet the next-sample MEKF commit transaction.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_deployment_config as D
from tools.stability.ou3_alt_contraction import finite_tuner_rs_alpha_machine_real_join as A
from tools.stability.ou3_alt_contraction import finite_tuner_rs_ema_compiler_modes as EMA
from tools.stability.ou3_alt_contraction import finite_tuner_spectral_machine_real_join as T

SOURCE=Path(__file__).resolve().parents[3]/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h'
SEED_RS=B.rn32(F(1,2))
MAX_UPDATES=4096
QUALIFICATION='OU3_ALT_RS_DEPLOYMENT_LEDGER_V1'


@dataclass(frozen=True)
class State:
    separate:F=SEED_RS
    fma:F=SEED_RS
    updates:int=0
    qualification:str=QUALIFICATION
    def __post_init__(self):
        object.__setattr__(self,'separate',F(self.separate)); object.__setattr__(self,'fma',F(self.fma))
        if self.qualification!=QUALIFICATION: raise ValueError('wrong RS ledger qualification')
        if not B.is_binary32(self.separate) or not B.is_binary32(self.fma):
            raise ValueError('persistent RS ledger values must be actual binary32')
        if not isinstance(self.updates,int) or not 0<=self.updates<=MAX_UPDATES:
            raise ValueError('invalid RS ledger update count')


@dataclass(frozen=True)
class TrackStep:
    target:T.Join
    alpha:A.Join
    ema:EMA.Step
    selected:F
    mode:str
    def __post_init__(self):
        if self.mode not in ('separate','fma'): raise ValueError('invalid RS compiler mode')
        if not isinstance(self.target,T.Join) or not isinstance(self.alpha,A.Join) or not isinstance(self.ema,EMA.Step):
            raise TypeError('target, alpha and compiler-mode EMA step required')
        object.__setattr__(self,'selected',F(self.selected))
        if self.alpha.cfg!=self.target.cfg or self.alpha.tau_target!=self.target.target.tau_target:
            raise ValueError('RS alpha detached from same target/config')
        if self.ema.target!=self.target.machine_target_RS or self.ema.alpha_source!=self.alpha.machine:
            raise ValueError('RS EMA detached from same target/alpha')
        expected=self.ema.next_separate if self.mode=='separate' else self.ema.next_fma
        if self.selected!=expected: raise ValueError('selected RS successor detached from compiler mode')


@dataclass(frozen=True)
class Result:
    state:State
    separate:TrackStep
    fma:TrackStep
    def __post_init__(self):
        if not isinstance(self.state,State) or not isinstance(self.separate,TrackStep) or not isinstance(self.fma,TrackStep):
            raise TypeError('RS ledger result components required')
        if self.separate.mode!='separate' or self.fma.mode!='fma':
            raise ValueError('global RS compiler histories crossed')
        if self.state.separate!=self.separate.selected or self.state.fma!=self.fma.selected:
            raise ValueError('RS ledger successor detached from selected compiler histories')
        if self.separate.target.cfg!=self.fma.target.cfg:
            raise ValueError('RS compiler histories must share one deployment config')


def initial(): return State()


def _track(previous,target,alpha,mode):
    if not isinstance(target,T.Join) or not isinstance(alpha,A.Join):
        raise TypeError('same-source RS target and alpha joins required')
    if alpha.cfg!=target.cfg or alpha.tau_target!=target.target.tau_target:
        raise ValueError('RS alpha detached from same target/config')
    machine=EMA.step(previous,target.machine_target_RS,alpha.machine)
    selected=machine.next_separate if mode=='separate' else machine.next_fma
    return TrackStep(target,alpha,machine,selected,mode)


def step(state:State,*,separate_target:T.Join,separate_alpha:A.Join,
         fma_target:T.Join,fma_alpha:A.Join):
    if not isinstance(state,State): raise TypeError('RS deployment ledger State required')
    if state.updates>=MAX_UPDATES: raise ValueError('RS ledger exhausted certified finite update budget')
    if separate_target.cfg!=fma_target.cfg:
        raise ValueError('RS compiler histories detached from common deployment config')
    s=_track(state.separate,separate_target,separate_alpha,'separate')
    f=_track(state.fma,fma_target,fma_alpha,'fma')
    return Result(State(s.selected,f.selected,state.updates+1),s,f)


def hold(state:State):
    """Cold/no-candidate sample: persistent R_S machine state is identity."""
    if not isinstance(state,State): raise TypeError('RS deployment ledger State required')
    return state


def _source_shape_matches():
    s=SOURCE.read_text()
    return all(x in s for x in (
      'float RS_applied    = 0.5f;',
      'tune_.RS_applied    += alpha_RS * (RS_t    - tune_.RS_applied);'))


def readiness():
    return {
      'qualification':QUALIFICATION,
      'shipping_RS_seed_and_update_source_shape_matches':_source_shape_matches(),
      'source_locked_binary32_RS_seed_0p5':True,
      'persistent_separate_and_FMA_RS_histories_carried':True,
      'compiler_histories_may_consume_distinct_SpectralMSE_targets':True,
      'compiler_histories_may_consume_distinct_alpha_RS_values':True,
      'each_track_alpha_bound_to_its_same_target_tau_and_config':True,
      'Cold_or_no_candidate_is_literal_RS_identity':True,
      'pending_next_sample_commit_snapshot_attached':False,
      'startup_frontend_RS_machine_history_attached':False,
      'Live_600_step_RS_machine_history_attached':False,
      'target_libm_correspondence_closed':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
