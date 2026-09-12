"""Shared measurement-only frontend composition for the OU-III ALT word.

Shipping computes one private-Mahony vertical acceleration and reuses it for the
adaptive sigma band, tracker-input LPF, and—later in the same sample—the WPE.
The ordering matters: sigma/tuner reads a read-only view of the WPE state carried
into the sample, while ``wave_period_.update`` consumes the current vertical
sample only after tuning and direction work.  This module makes both the shared
vertical alias and that causal split structural.

The dominant-frequency tracker itself remains an explicit runtime output witness
only for the full StillnessAdapter/direction helper.  The stability-side tuner
projection removes it exactly elsewhere.  Frontend transcendental binary32
evaluations remain later obligations.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_prediction_graph as P
from tools.stability.ou3_alt_contraction import finite_vertical_complementary_runtime as V
from tools.stability.ou3_alt_contraction import finite_wpe_runtime as W
from tools.stability.ou3_alt_contraction import finite_band_variance_runtime as B
from tools.stability.ou3_alt_contraction import finite_stillness_runtime as S


def R(x): return P.rational(x)


@dataclass(frozen=True)
class LPFState:
    state:F=F(0)
    initialized:bool=False
    def __post_init__(self):
        object.__setattr__(self,'state',R(self.state))
        if not isinstance(self.initialized,bool): raise TypeError('literal LPF initialization flag required')


@dataclass(frozen=True)
class LPFDecayWitness:
    decay:F
    def __post_init__(self):
        d=R(self.decay)
        if not 0<d<=1: raise ValueError('LPF exp decay must lie in (0,1]')
        object.__setattr__(self,'decay',d)


@dataclass(frozen=True)
class LPFResult:
    state:LPFState
    output:F


@dataclass(frozen=True)
class TrackerOutputWitness:
    frequency:F
    def __post_init__(self):
        f=R(self.frequency)
        if f<=0: raise ValueError('finite positive tracker frequency required for represented branch')
        object.__setattr__(self,'frequency',f)


def tracker_lpf_step(s:LPFState,vertical:V.Result,*,decay:LPFDecayWitness):
    if not isinstance(vertical,V.Result): raise TypeError('private vertical successor required')
    x=vertical.vertical_accel
    if not s.initialized:
        return LPFResult(LPFState(x,True),x)
    d=decay.decay
    y=(1-d)*x+d*s.state
    return LPFResult(LPFState(y,True),y)


def wpe_view(wpe_state:W.WPEState,*,current_period=None,current_frequency=None):
    """Tuner-visible WPE view BEFORE this sample advances WPE."""
    if not isinstance(wpe_state,W.WPEState): raise TypeError('carried WPE state required')
    if wpe_state.usable_period:
        return B.wpe_frequency_view(wpe_state,period=current_period,frequency=current_frequency)
    return B.wpe_frequency_view(wpe_state)


def wpe_step_from_vertical(wpe_state:W.WPEState,wpe_cfg:W.WPEConfig,vertical:V.Result,**kwargs):
    if not isinstance(vertical,V.Result): raise TypeError('private vertical successor required')
    if 'vertical_accel' in kwargs: raise TypeError('vertical acceleration is owned by private observer')
    return W.update(wpe_state,wpe_cfg,vertical_accel=vertical.vertical_accel,**kwargs)


def band_step_from_vertical_view(band:B.BandState,stats:B.StatsState,view:B.WPEFrequencyView,vertical:V.Result,**kwargs):
    if not isinstance(vertical,V.Result): raise TypeError('private vertical successor required')
    if not isinstance(view,B.WPEFrequencyView): raise TypeError('pre-update WPE tuner view required')
    if 'vertical_accel' in kwargs: raise TypeError('vertical acceleration is owned by private observer')
    return B.frontend_step_from_view(band,stats,view,vertical_accel=vertical.vertical_accel,**kwargs)


def band_step_from_vertical(band:B.BandState,stats:B.StatsState,wpe:W.UpdateResult,vertical:V.Result,**kwargs):
    """Legacy isolated helper; shipping temporal word uses *_from_vertical_view."""
    if not isinstance(vertical,V.Result): raise TypeError('private vertical successor required')
    if 'vertical_accel' in kwargs: raise TypeError('vertical acceleration is owned by private observer')
    return B.frontend_step(band,stats,wpe,vertical_accel=vertical.vertical_accel,**kwargs)


def stillness_step_from_vertical(still:S.State,still_cfg:S.Config,lpf:LPFResult,tracker:TrackerOutputWitness,**kwargs):
    if not isinstance(lpf,LPFResult) or not isinstance(tracker,TrackerOutputWitness):
        raise TypeError('finite tracker LPF successor and tracker output witness required')
    if 'a_vert_up_lp' in kwargs or 'tracker_frequency' in kwargs:
        raise TypeError('stillness inputs are owned by tracker frontend')
    return S.step(still,still_cfg,a_vert_up_lp=lpf.output,tracker_frequency=tracker.frequency,**kwargs)


def readiness():
    return {
      'one_private_vertical_output_shared_by_tuner_LPF_and_later_WPE':True,
      'tuner_uses_preupdate_WPE_state':True,
      'current_vertical_advances_WPE_only_after_tuner_entry':True,
      'same_private_vertical_output_drives_tracker_input_LPF':True,
      'tracker_input_LPF_materialized':True,
      'tracker_output_algorithm_attached':False,
      'LPF_exp_binary32_attached':False,
      'private_vertical_seed_and_fast_inv_sqrt_attached':False,
      'raw_sensor_BRMM_disturbance_relation_attached':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
