"""Persistent dual-compiler binary32 ledger for shipping ``sigma_applied``.

Shipping seeds ``TuneState::sigma_applied = 1e-2f`` and updates it with the SAME
``alpha`` used by tau.  The theorem-facing sigma ledger therefore consumes a
qualified common-alpha object, not a free coefficient.  As with tau and R_S,
compiler contraction is not guessed: one coherent global track always uses
separate multiply/add and the other always uses FMA.

The two tracks may consume different mode-coherent sigma targets and different
common-alpha witnesses inherited from distinct WPE/tau histories.  They must
still share one deployment configuration.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_common_alpha_qualification as A
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_binary32 as S

SOURCE=Path(__file__).resolve().parents[3]/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h'
INITIAL=B.rn32(F(1,100)); MAX_UPDATES=30600
QUALIFICATION='OU3_ALT_SIGMA_DEPLOYMENT_LEDGER_V1'


@dataclass(frozen=True)
class State:
    separate:F=INITIAL
    fma:F=INITIAL
    updates:int=0
    qualification:str=QUALIFICATION
    def __post_init__(self):
        a,b=F(self.separate),F(self.fma)
        object.__setattr__(self,'separate',a); object.__setattr__(self,'fma',b)
        if self.qualification!=QUALIFICATION: raise ValueError('wrong sigma ledger qualification')
        if not B.is_binary32(a) or not B.is_binary32(b) or a<0 or b<0:
            raise ValueError('sigma deployment ledger stores nonnegative binary32 only')
        if not isinstance(self.updates,int) or not 0<=self.updates<=MAX_UPDATES:
            raise ValueError('sigma deployment update count outside bounded startup+word horizon')


@dataclass(frozen=True)
class TrackStep:
    target:S.Target
    alpha:A.Qualified
    previous:F
    selected:F
    mode:str
    def __post_init__(self):
        if self.mode not in ('separate','fma'): raise ValueError('invalid sigma compiler mode')
        if not isinstance(self.target,S.Target) or not isinstance(self.alpha,A.Qualified):
            raise TypeError('sigma target and qualified common alpha required')
        object.__setattr__(self,'previous',F(self.previous)); object.__setattr__(self,'selected',F(self.selected))
        if not B.is_binary32(self.previous) or not B.is_binary32(self.selected):
            raise ValueError('sigma ledger track stores binary32 values only')
        if self.target.cfg!=self.alpha.deployment_cfg:
            raise ValueError('sigma target detached from qualified common-alpha deployment config')
        expected=B.ema(self.previous,self.target.sigma_target,self.alpha.alpha,
                       contracted=(self.mode=='fma'))
        if self.selected!=expected: raise ValueError('sigma successor detached from compiler mode/common alpha')


@dataclass(frozen=True)
class Result:
    state:State
    separate:TrackStep
    fma:TrackStep
    def __post_init__(self):
        if not isinstance(self.state,State) or not isinstance(self.separate,TrackStep) or not isinstance(self.fma,TrackStep):
            raise TypeError('sigma ledger result components required')
        if self.separate.mode!='separate' or self.fma.mode!='fma': raise ValueError('sigma compiler histories crossed')
        if self.state.separate!=self.separate.selected or self.state.fma!=self.fma.selected:
            raise ValueError('sigma ledger successor detached from selected histories')
        if self.separate.target.cfg!=self.fma.target.cfg:
            raise ValueError('sigma histories must share one deployment config')


def initial(): return State()

def hold(state:State):
    if not isinstance(state,State): raise TypeError('sigma deployment State required')
    return state


def _track(previous,target,alpha,mode):
    if not isinstance(target,S.Target) or not isinstance(alpha,A.Qualified):
        raise TypeError('sigma target and qualified common alpha required')
    if target.cfg!=alpha.deployment_cfg:
        raise ValueError('sigma target detached from qualified common-alpha deployment config')
    selected=B.ema(previous,target.sigma_target,alpha.alpha,contracted=(mode=='fma'))
    return TrackStep(target,alpha,previous,selected,mode)


def step(state:State,*,separate_target:S.Target,separate_alpha:A.Qualified,
         fma_target:S.Target,fma_alpha:A.Qualified):
    if not isinstance(state,State): raise TypeError('sigma deployment State required')
    if state.updates>=MAX_UPDATES: raise ValueError('sigma deployment ledger exceeded bounded startup+word horizon')
    if separate_target.cfg!=fma_target.cfg:
        raise ValueError('sigma compiler histories detached from common deployment config')
    s=_track(state.separate,separate_target,separate_alpha,'separate')
    f=_track(state.fma,fma_target,fma_alpha,'fma')
    return Result(State(s.selected,f.selected,state.updates+1),s,f)


def _source_shape_matches():
    s=SOURCE.read_text()
    return all(x in s for x in (
      'float sigma_applied = 1e-2f;',
      'tune_.sigma_applied += alpha    * (sigma_t - tune_.sigma_applied);'))


def readiness():
    return {
      'qualification':QUALIFICATION,
      'shipping_sigma_seed_and_update_source_shape_matches':_source_shape_matches(),
      'source_locked_binary32_sigma_seed_0p01':True,
      'persistent_separate_and_FMA_sigma_histories_carried':True,
      'sigma_each_track_consumes_qualified_common_tau_sigma_alpha':True,
      'compiler_histories_may_consume_distinct_sigma_targets_and_alphas':True,
      'cold_or_nonadapting_sigma_identity_branch_materialized':True,
      'pending_common_TuneState_boundary_attached':False,
      'startup_frontend_sigma_machine_history_attached':False,
      'Live_600_step_sigma_machine_history_attached':False,
      'sigma_sqrt_and_still_exp_libm_correspondence_closed':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
