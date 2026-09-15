"""Binary32 sigma-relevant projection of shipping StillnessAdapter.

The sigma channel reads only ``last_is_still`` and ``still_time_sec`` from the
shared StillnessAdapter.  Those two states are determined by the energy EMA and
do not depend on the adapter's tracker-frequency relaxation state.  This module
therefore proves the narrow projection actually consumed by the tuner instead
of introducing an irrelevant frequency-side state into the sigma theorem.

Represented source order:

  a_norm = a_vert_up_lp / gravity_std
  inst = a_norm * a_norm
  energy = (1-energy_alpha)*old + energy_alpha*inst
  is_still = energy < energy_thresh
  still_time = min(old+dt, 60) if still else 0
  attenuation = clamp(exp(-still_time/1),0,1) if still else 1

The energy multiply-add carries every local no-reassociation contraction choice.
The actual stored energy successor is a witness constrained to that finite set;
the still predicate is then evaluated from that actual successor.  The later
attenuation exp is tied to the SAME stored still-time by a tight real enclosure
and binary32 RNE cell. Platform expf/compiler correspondence remains open.
"""
from __future__ import annotations

from tools.stability.ou3_alt_contraction import finite_wrapper_clock_binary32 as HORIZON
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_stillness_runtime as R
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_binary32 as E

COMMON=Path(__file__).resolve().parents[3]/'src/tuner/SeaStateFusionTunerCommon.h'
FILTER=Path(__file__).resolve().parents[3]/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h'
ZERO=B.rn32(0); ONE=B.rn32(1); SIXTY=B.rn32(60)
QUALIFICATION='OU3_ALT_STILLNESS_SIGMA_BINARY32_V2'; MAX_SAMPLES=HORIZON.MAX_STEPS


def _q(x,name):
    q=F(x)
    if not B.is_binary32(q): raise ValueError(f'{name} must be actual binary32')
    return q

def _uniq(values): return tuple(sorted(set(F(v) for v in values)))

def _sum_products(a,b,c,d):
    ab=B.mul(a,b); cd=B.mul(c,d)
    return _uniq((B.add(ab,cd),B.fma(a,b,cd),B.fma(c,d,ab)))


@dataclass(frozen=True)
class State:
    energy:F=ZERO
    still_time:F=ZERO
    is_still:bool=False
    samples:int=0
    qualification:str=QUALIFICATION
    def __post_init__(self):
        object.__setattr__(self,'energy',_q(self.energy,'stillness energy'))
        object.__setattr__(self,'still_time',_q(self.still_time,'stillness time'))
        if not isinstance(self.is_still,bool): raise TypeError('literal stillness flag required')
        if not isinstance(self.samples,int) or not 0<=self.samples<=MAX_SAMPLES: raise ValueError('stillness sample count outside bounded startup+word horizon')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong stillness sigma qualification')
        if self.energy<0 or self.still_time<0: raise ValueError('nonnegative stillness state required')


@dataclass(frozen=True)
class Step:
    before:State
    state:State
    vertical_lp:F
    dt:F
    threshold:F
    a_norm:F
    inst_energy:F
    energy_values:tuple
    attenuation:F
    exp_result:F|None
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.before,State) or not isinstance(self.state,State): raise TypeError('stillness predecessor/successor required')
        for n in ('vertical_lp','dt','threshold','a_norm','inst_energy','attenuation'):
            object.__setattr__(self,n,_q(getattr(self,n),n))
        vals=_uniq(self.energy_values); object.__setattr__(self,'energy_values',vals)
        if self.exp_result is not None: object.__setattr__(self,'exp_result',_q(self.exp_result,'stillness attenuation exp'))
        if self.qualification!=QUALIFICATION: raise ValueError('wrong stillness step qualification')
        if self.state.energy not in vals: raise ValueError('actual stillness energy outside local contraction set')
        if self.state.samples!=self.before.samples+1: raise ValueError('stillness sample count did not advance once')
        if self.state.is_still!=(self.state.energy < self.threshold): raise ValueError('still predicate detached from actual machine energy/threshold')
        expected_time=min(B.add(self.before.still_time,self.dt),SIXTY) if self.state.is_still else ZERO
        if self.state.still_time!=expected_time: raise ValueError('still time detached from actual machine predicate and source update')
        if self.state.is_still:
            if self.exp_result is None or self.attenuation!=min(max(self.exp_result,ZERO),ONE): raise ValueError('still branch attenuation detached from exp result')
        else:
            if self.exp_result is not None or self.attenuation!=ONE: raise ValueError('moving branch must use attenuation one and no exp witness')


def step(state:State,cfg:R.Config,*,vertical_lp,dt,energy_successor,attenuation_exp=None):
    if not isinstance(state,State) or not isinstance(cfg,R.Config): raise TypeError('stillness State and Config required')
    if state.samples>=MAX_SAMPLES: raise ValueError('stillness projection exceeded bounded startup+word horizon')
    a=_q(vertical_lp,'stillness vertical LP'); h=_q(dt,'stillness dt')
    if h<=0: raise ValueError('positive stillness dt required')
    g=B.rn32(cfg.gravity); alpha=B.rn32(cfg.energy_alpha); threshold=B.rn32(cfg.energy_thresh)
    if g<=0 or not 0<alpha<=1: raise ValueError('invalid compiled stillness config')
    an=B.div(a,g); inst=B.mul(an,an); decay=B.sub(ONE,alpha)
    values=_sum_products(decay,state.energy,alpha,inst)
    en=_q(energy_successor,'actual stillness energy successor')
    if en not in values: raise ValueError('actual stillness energy successor outside local contraction set')
    still=en<threshold
    if still:
        st=min(B.add(state.still_time,h),SIXTY)
        if attenuation_exp is None: raise ValueError('still branch requires attenuation exp result')
        e=_q(attenuation_exp,'stillness attenuation exp')
        lo,hi=E.exp_minus_enclosure(st)
        if not E._interval_hits_rne_cell(lo,hi,e): raise ValueError('stillness attenuation exp detached from SAME stored still-time RNE cell')
        atten=min(max(e,ZERO),ONE)
    else:
        if attenuation_exp is not None: raise ValueError('moving branch consumes no attenuation exp result')
        st=ZERO; atten=ONE; e=None
    nxt=State(en,st,still,state.samples+1)
    return Step(state,nxt,a,h,threshold,an,inst,values,atten,e)


def _source_shape_matches():
    s=COMMON.read_text(); f=FILTER.read_text()
    return all(x in s for x in (
      'const float a_norm      = a_vert_up_lp / gravity_std_;',
      'const float inst_energy = a_norm * a_norm;',
      'energy_ema = (1.0f - energy_alpha) * energy_ema + energy_alpha * inst_energy;',
      'const bool is_still = (energy_ema < energy_thresh);','last_is_still = is_still;',
      'still_time_sec += dt;','if (still_time_sec > 60.0f) {','still_time_sec = 60.0f;',
      'still_time_sec = 0.0f;')) and all(x in f for x in (
      'if (freq_stillness_.isStill()) {','const float still_t = std::max(0.0f, freq_stillness_.getStillTime());',
      'float atten = std::exp(-still_t / STILL_VAR_DECAY_SEC);','var_wave *= atten;'))


def readiness():
    return {
      'qualification':QUALIFICATION,
      'shipping_sigma_relevant_stillness_source_shape_matches':_source_shape_matches(),
      'sigma_projection_independent_of_tracker_frequency_relaxation_state':True,
      'energy_normalization_square_and_EMA_binary32_graph_materialized':True,
      'actual_energy_successor_bound_to_all_local_contraction_choices':True,
      'still_predicate_and_capped_time_derived_from_actual_machine_energy':True,
      'attenuation_exp_bound_to_same_stored_machine_still_time_by_RNE_cell':True,
      'moving_branch_attenuation_is_literal_one':True,
      'target_exp_libm_correspondence_closed':False,
      'upstream_vertical_LP_machine_production_closed':False,
      'target_compiler_contraction_membership_closed':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
