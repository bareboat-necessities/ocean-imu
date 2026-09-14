"""Persistent binary32 private-Mahony -> tracker-LPF -> stillness source.

This closes the missing upstream source relation for the machine sigma path.
Shipping computes one private-Mahony vertical acceleration, sends that SAME
vertical value directly to the adaptive wave band, and also sends it through
FreqInputLPF before StillnessAdapter.  The tuner sigma projection consumes only
StillnessAdapter's energy/still-time state, so the acceleration-band tracker is
irrelevant to this narrow relation.

Represented machine order:

  guarded body IMU -> first-sample or initialized private Mahony -> vertical_accel
  vertical_accel -> FreqInputLPF -> a_vert_up_lp
  a_vert_up_lp -> stillness energy/predicate/time/attenuation

The LPF source expression is literal binary32 shipping order
``exp(-2*pi_f*fc*dt)`` followed by
``(1-alpha)*x + alpha*state``.  Both legal no-reassociation contraction outcomes
for the final sum are retained.  The exp result is tied to the same rounded
argument by a rigorous real enclosure and its binary32 RNE cell.  Target libm
and compiler contraction-profile correspondence remain open.  The mutable
shipping LPF cutoff is carried by this local state but is not yet rooted in the
theorem RuntimeConfig in this local component. The joined startup/Live wrappers
fix the default cutoff; mutable setter ancestry remains open. The ordinary
first-sample seed is attached here; nearly antiparallel SVD and nonfinite
branches are still unqualified, not removed from the admitted source.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_binary32_mahony as MAHONY
from tools.stability.ou3_alt_contraction import finite_binary32_mahony_startup as STARTUP
from tools.stability.ou3_alt_contraction import finite_vertical_complementary_runtime as VERT
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as SENSOR
from tools.stability.ou3_alt_contraction import finite_stillness_runtime as STILL_CFG
from tools.stability.ou3_alt_contraction import finite_stillness_sigma_binary32 as STILL
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_binary32 as EXP

COMMON=Path(__file__).resolve().parents[3]/'src/tuner/SeaStateFusionTunerCommon.h'
FILTER=Path(__file__).resolve().parents[3]/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h'
QUALIFICATION='OU3_ALT_MACHINE_VERTICAL_STILLNESS_SOURCE_V1'
ZERO=B.rn32(0); ONE=B.rn32(1); TWO=B.rn32(2)
# static_cast<float>(M_PI) under the declared IEEE binary32 arithmetic profile.
PI_F=MAHONY.value(0x40490FDB)
MAX_SAMPLES=STILL.MAX_SAMPLES


def _q(x,name):
    q=F(x)
    if not B.is_binary32(q): raise ValueError(f'{name} must be actual binary32')
    return q


def _uniq(values): return tuple(sorted(set(F(v) for v in values)))


@dataclass(frozen=True)
class LPFState:
    value:F=ZERO
    cutoff_hz:F=ONE
    initialized:bool=False
    samples:int=0
    qualification:str=QUALIFICATION
    def __post_init__(self):
        object.__setattr__(self,'value',_q(self.value,'LPF state'))
        object.__setattr__(self,'cutoff_hz',_q(self.cutoff_hz,'LPF cutoff'))
        if self.cutoff_hz<=0: raise ValueError('positive FreqInputLPF cutoff required')
        if not isinstance(self.initialized,bool): raise TypeError('literal LPF initialized flag required')
        if not isinstance(self.samples,int) or not 0<=self.samples<=MAX_SAMPLES:
            raise ValueError('LPF sample count outside bounded startup+word horizon')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong machine LPF qualification')


@dataclass(frozen=True)
class LPFStep:
    before:LPFState
    state:LPFState
    input:F
    dt:F
    exp_argument:F|None
    alpha:F|None
    successor_values:tuple
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.before,LPFState) or not isinstance(self.state,LPFState):
            raise TypeError('machine LPF predecessor/successor required')
        object.__setattr__(self,'input',_q(self.input,'LPF input'))
        object.__setattr__(self,'dt',_q(self.dt,'LPF dt'))
        if self.dt<=0: raise ValueError('positive LPF dt required')
        vals=_uniq(self.successor_values); object.__setattr__(self,'successor_values',vals)
        if self.exp_argument is not None: object.__setattr__(self,'exp_argument',_q(self.exp_argument,'LPF exp argument'))
        if self.alpha is not None: object.__setattr__(self,'alpha',_q(self.alpha,'LPF alpha'))
        if self.state.value not in vals: raise ValueError('LPF successor outside same-expression contraction set')
        if self.state.cutoff_hz!=self.before.cutoff_hz or self.state.samples!=self.before.samples+1:
            raise ValueError('LPF persistent configuration/count detached')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong machine LPF step qualification')
        if not self.before.initialized:
            if self.exp_argument is not None or self.alpha is not None or vals!=(self.input,):
                raise ValueError('first FreqInputLPF sample is exact seed and consumes no exp')
            if not self.state.initialized or self.state.value!=self.input:
                raise ValueError('first FreqInputLPF sample did not seed exact input')
        elif self.exp_argument is None or self.alpha is None:
            raise ValueError('initialized FreqInputLPF step requires its exp relation')


def lpf_step(state:LPFState,*,x,dt,alpha_exp=None,successor=None):
    if not isinstance(state,LPFState): raise TypeError('machine LPF State required')
    if state.samples>=MAX_SAMPLES: raise ValueError('machine LPF exceeded bounded startup+word horizon')
    xv=_q(x,'LPF input'); h=_q(dt,'LPF dt')
    if h<=0: raise ValueError('positive LPF dt required')
    if not state.initialized:
        if alpha_exp is not None or successor is not None:
            raise ValueError('first FreqInputLPF sample consumes no exp/successor witness')
        nxt=LPFState(xv,state.cutoff_hz,True,state.samples+1)
        return LPFStep(state,nxt,xv,h,None,None,(xv,))
    # C++: -2.0f * static_cast<float>(M_PI) * fc_hz * dt, left-associated.
    mag=B.mul(B.mul(B.mul(TWO,PI_F),state.cutoff_hz),h)
    if mag<0 or mag>60: raise ValueError('LPF exp argument outside retained proof enclosure')
    if alpha_exp is None: raise TypeError('initialized FreqInputLPF requires exp result')
    alpha=_q(alpha_exp,'LPF exp result')
    lo,hi=EXP.exp_minus_enclosure(mag)
    if not EXP._interval_hits_rne_cell(lo,hi,alpha):
        raise ValueError('LPF exp result detached from SAME rounded argument RNE cell')
    one_minus=B.sub(ONE,alpha)
    ax=B.mul(one_minus,xv); ap=B.mul(alpha,state.value)
    vals=_uniq((B.add(ax,ap),B.fma(one_minus,xv,ap),B.fma(alpha,state.value,ax)))
    if successor is None: raise TypeError('initialized LPF requires actual stored successor')
    y=_q(successor,'LPF stored successor')
    if y not in vals: raise ValueError('LPF stored successor outside local contraction set')
    nxt=LPFState(y,state.cutoff_hz,True,state.samples+1)
    return LPFStep(state,nxt,xv,h,mag,alpha,vals)


@dataclass(frozen=True)
class State:
    vertical:VERT.State
    lpf:LPFState
    stillness:STILL.State
    samples:int=0
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.vertical,VERT.State) or not isinstance(self.lpf,LPFState) or not isinstance(self.stillness,STILL.State):
            raise TypeError('machine Mahony, LPF and stillness states required')
        if not self.vertical.initialized and self.vertical != VERT.State():
            raise ValueError('uninitialized machine observer lost its literal reset state')
        if self.lpf.samples!=self.stillness.samples or self.lpf.samples!=self.samples:
            raise ValueError('machine vertical frontend sample counts detached')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong machine vertical/stillness qualification')


@dataclass(frozen=True)
class Result:
    before:State
    state:State
    mahony:MAHONY.StepResult
    lpf:LPFStep
    stillness:STILL.Step
    band_input:F
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.before,State) or not isinstance(self.state,State): raise TypeError('machine vertical source states required')
        if not isinstance(self.mahony,MAHONY.StepResult) or not isinstance(self.lpf,LPFStep) or not isinstance(self.stillness,STILL.Step):
            raise TypeError('machine Mahony/LPF/stillness results required')
        object.__setattr__(self,'band_input',_q(self.band_input,'machine band input'))
        if self.mahony.vertical.state!=self.state.vertical or self.band_input!=self.mahony.vertical.vertical_accel:
            raise ValueError('machine band input detached from SAME Mahony vertical successor')
        if self.lpf.before!=self.before.lpf or self.lpf.state!=self.state.lpf or self.lpf.input!=self.band_input:
            raise ValueError('machine LPF detached from SAME Mahony vertical output')
        if self.stillness.before!=self.before.stillness or self.stillness.state!=self.state.stillness:
            raise ValueError('machine stillness detached from persistent predecessor/successor')
        if self.stillness.vertical_lp!=self.lpf.state.value:
            raise ValueError('machine stillness input detached from SAME stored LPF output')
        if self.state.samples!=self.before.samples+1:
            raise ValueError('machine vertical/stillness source did not advance one sample')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong machine vertical source result qualification')


def begin(vertical:VERT.State,tracker_lpf_value,tracker_lpf_initialized,stillness:STILL.State,*,cutoff_hz,samples=0):
    if not isinstance(vertical,VERT.State) or not isinstance(stillness,STILL.State): raise TypeError('carried vertical/stillness states required')
    if stillness.samples!=samples: raise ValueError('carried stillness/sample ordinal detached')
    lpf=LPFState(tracker_lpf_value,cutoff_hz,tracker_lpf_initialized,samples)
    return State(vertical,lpf,stillness,samples)


def step(state:State,guarded:SENSOR.GuardedImuSample,vertical_cfg:VERT.Config,still_cfg:STILL_CFG.Config,*,dt,
         lpf_alpha_exp=None,lpf_successor=None,still_energy_successor,still_attenuation_exp=None,
         profile=MAHONY.PROFILE):
    if not isinstance(state,State) or not isinstance(guarded,SENSOR.GuardedImuSample):
        raise TypeError('persistent machine source and SAME guarded sample required')
    h=_q(dt,'machine frontend dt')
    mah=STARTUP.step(state.vertical,vertical_cfg,dt=h,
        gyro=guarded.raw_gyro_body,acc=guarded.conditioned_accel_body,profile=profile)
    band_input=_q(mah.vertical.vertical_accel,'machine Mahony vertical')
    lp=lpf_step(state.lpf,x=band_input,dt=h,alpha_exp=lpf_alpha_exp,successor=lpf_successor)
    st=STILL.step(state.stillness,still_cfg,vertical_lp=lp.state.value,dt=h,
        energy_successor=still_energy_successor,attenuation_exp=still_attenuation_exp)
    nxt=State(mah.vertical.state,lp.state,st.state,state.samples+1)
    return Result(state,nxt,mah,lp,st,band_input)


def require_frontend_input(source:Result,frontend_result):
    """Bind AdaptiveWaveBandPass input to this exact machine vertical sample."""
    from tools.stability.ou3_alt_contraction import finite_machine_frontend_sigma_source as FRONTEND
    if not isinstance(source,Result) or not isinstance(frontend_result,FRONTEND.Result):
        raise TypeError('machine vertical source and frontend result required')
    if frontend_result.band.envelope.x!=source.band_input:
        raise ValueError('machine adaptive-band input detached from SAME Mahony vertical output')
    return frontend_result


def require_sigma_stillness(source:Result,target):
    """Bind sigma stillness operands to this exact machine stillness successor."""
    from tools.stability.ou3_alt_contraction import finite_tuner_sigma_binary32 as SIGMA
    if not isinstance(source,Result) or not isinstance(target,SIGMA.Target):
        raise TypeError('machine vertical source and sigma target required')
    st=source.stillness
    if (target.still,target.still_time,target.attenuation)!=(st.state.is_still,st.state.still_time,st.attenuation):
        raise ValueError('sigma stillness operands detached from SAME machine stillness successor')
    return target


def _source_shape_matches():
    c=COMMON.read_text(); f=FILTER.read_text()
    return all(x in c for x in (
      'const float alpha = std::exp(-2.0f * static_cast<float>(M_PI) * fc_hz * dt);',
      'state = (1.0f - alpha) * x + alpha * state;',
      'const float a_norm      = a_vert_up_lp / gravity_std_;',
      'energy_ema = (1.0f - energy_alpha) * energy_ema + energy_alpha * inst_energy;')) and all(x in f for x in (
      'const float a_vert_measurement = vertical_accel_comp_.verticalAccelUpMs2();',
      'const float a_vert_lp = freq_input_lpf_.step(a_vert_measurement, dt);',
      'const float f_after_still = freq_stillness_.step(a_vert_lp, dt, f_tracker);'))


def readiness():
    return {
      'qualification':QUALIFICATION,
      'shipping_vertical_LPF_stillness_source_shape_matches':_source_shape_matches(),
      'initialized_private_Mahony_binary32_graph_consumed':True,
      'same_machine_Mahony_vertical_drives_band_and_tracker_LPF':True,
      'FreqInputLPF_binary32_recurrence_materialized':True,
      'FreqInputLPF_first_sample_exact_seed_materialized':True,
      'FreqInputLPF_exp_bound_to_same_rounded_argument_RNE_cell':True,
      'FreqInputLPF_contraction_choices_retained':True,
      'same_stored_LPF_output_drives_binary32_stillness_projection':True,
      'sigma_stillness_operands_bindable_to_same_machine_successor':True,
      'machine_band_input_bindable_to_same_Mahony_successor':True,
      'tracker_LPF_cutoff_runtime_ancestry_closed':False,
      'target_libm_and_compiler_profile_correspondence_closed':False,
      'guard_binary32_machine_history_attached':False,
      'startup_machine_history_attached':False,
      'Live_600_step_machine_history_attached':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,'ALT_STARTUP_PASS':False,'ALT_END_TO_END_PASS':False,
    }
