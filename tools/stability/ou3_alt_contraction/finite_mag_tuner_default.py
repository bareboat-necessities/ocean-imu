"""Exact default MagAutoTuner finite state machine for ALT.

This module composes the deployed default (unweighted, no startup hard-iron
solve) MagAutoTuner path around ``finite_mag_accumulator`` and
``finite_mag_gauge_fix``.  It preserves shipping ordering:

  ready latch -> input/norm gates -> running-norm outlier gate -> unit-weight
  accumulation -> count/window/effective-weight finalize gates -> mean norm /
  horizontal-field gates -> gauge-fixed reference.

Binary32 allFinite/sqrt/quaternion-normalization correspondence and the physical
source ancestry of the startup magnetic packet remain explicit open obligations.
The optional quality-weighted and hard-iron branches are not silently replaced
by this default branch.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as SENSOR
from tools.stability.ou3_alt_contraction import finite_mag_accumulator as ACC
from tools.stability.ou3_alt_contraction import finite_mag_gauge_fix as GAUGE


def R(x): return M.rational(x)


@dataclass(frozen=True)
class SqrtWitness:
    radicand:F
    value:F
    def __post_init__(self):
        a,v=R(self.radicand),R(self.value)
        if a<0 or v<0 or v*v!=a:
            raise ValueError('exact nonnegative sqrt witness required')
        object.__setattr__(self,'radicand',a); object.__setattr__(self,'value',v)


@dataclass(frozen=True)
class Config:
    mag_norm_min:F=F(1,1000)
    max_norm_ratio:F=F(35,100)
    min_sample_weight:F=F(3,100)
    min_horizontal_fraction:F=F(5,100)
    min_effective_weight:F=F(0)
    min_samples:int=250
    min_window:F=F(10)
    max_window:F=F(0)
    sample_dt:F=F(1,200)
    quality_weighting:bool=False
    estimate_hard_iron:bool=False
    def __post_init__(self):
        for n in ('mag_norm_min','max_norm_ratio','min_sample_weight',
                  'min_horizontal_fraction','min_effective_weight','min_window',
                  'max_window','sample_dt'):
            v=R(getattr(self,n))
            if v<0: raise ValueError('nonnegative MagAutoTuner configuration required')
            object.__setattr__(self,n,v)
        if not isinstance(self.min_samples,int) or self.min_samples<0:
            raise ValueError('nonnegative integer min_samples required')
        if not isinstance(self.quality_weighting,bool) or not isinstance(self.estimate_hard_iron,bool):
            raise TypeError('literal MagAutoTuner branch flags required')
    def accumulator_config(self):
        return ACC.Config(self.min_samples,self.min_window,self.max_window,
                          self.sample_dt,self.quality_weighting)


@dataclass(frozen=True)
class State:
    accumulator:ACC.State=ACC.State()
    rejected_count:int=0
    ready:bool=False
    mean:tuple|None=None
    world_reference:tuple|None=None
    last_sample_weight:F=F(0)
    last_world_sample:tuple=(F(0),F(0),F(0))
    def __post_init__(self):
        if not isinstance(self.accumulator,ACC.State): raise TypeError('Mag accumulator state required')
        if not isinstance(self.rejected_count,int) or self.rejected_count<0:
            raise ValueError('nonnegative rejected count required')
        if not isinstance(self.ready,bool): raise TypeError('literal ready latch required')
        w=R(self.last_sample_weight)
        if w<0: raise ValueError('nonnegative last sample weight required')
        object.__setattr__(self,'last_sample_weight',w)
        object.__setattr__(self,'last_world_sample',tuple(M.vec(self.last_world_sample,3)))
        if self.ready:
            if self.mean is None or self.world_reference is None:
                raise ValueError('ready tuner must retain mean and world reference')
            object.__setattr__(self,'mean',tuple(M.vec(self.mean,3)))
            object.__setattr__(self,'world_reference',tuple(M.vec(self.world_reference,3)))
        elif self.mean is not None or self.world_reference is not None:
            raise ValueError('unready tuner cannot carry finalized magnetic reference')


@dataclass(frozen=True)
class StepResult:
    state:State
    returned_ready:bool
    accepted:bool
    rejection:str|None


def _reject(state:State,why:str):
    return StepResult(State(state.accumulator,state.rejected_count+1,False,None,None,
                            F(0),(F(0),F(0),F(0))),False,False,why)


def step(state:State,cfg:Config,*,q_tilt_bw,mag_body,dt,
         q_norm:SqrtWitness|None=None,mag_norm:SqrtWitness|None=None,
         mean_norm:SqrtWitness|None=None,horizontal_sqrt:GAUGE.HorizontalSqrt|None=None):
    if not isinstance(state,State) or not isinstance(cfg,Config):
        raise TypeError('default MagAutoTuner state/config required')
    q=tuple(M.vec(q_tilt_bw,4)); mag=tuple(M.vec(mag_body,3)); d=R(dt)

    # Shipping returns immediately before inspecting any new sample once ready.
    if state.ready:
        if any(x is not None for x in (q_norm,mag_norm,mean_norm,horizontal_sqrt)):
            raise ValueError('ready MagAutoTuner latch consumes no new arithmetic witnesses')
        return StepResult(state,True,False,None)

    if cfg.quality_weighting:
        raise NotImplementedError('quality-weighted MagAutoTuner is not the deployed default proof branch')
    if cfg.estimate_hard_iron:
        raise NotImplementedError('startup hard-iron solve is not the deployed default proof branch')

    q2=M.dot(q,q)
    if not isinstance(q_norm,SqrtWitness) or q_norm.radicand!=q2:
        raise ValueError('tilt quaternion norm witness detached from same input')
    if q_norm.value<=F(1,10**6):
        return _reject(state,'quaternion_norm')
    qn=tuple(x/q_norm.value for x in q)
    if M.dot(qn,qn)!=1: raise AssertionError('normalized tilt quaternion is not unit')

    m2=M.dot(mag,mag)
    if not isinstance(mag_norm,SqrtWitness) or mag_norm.radicand!=m2:
        raise ValueError('magnetometer norm witness detached from same input')
    if mag_norm.value<=cfg.mag_norm_min:
        return _reject(state,'mag_norm')

    world=tuple(SENSOR.q_rotate(qn,mag))
    # Unit quaternion rotation preserves the norm exactly.  Shipping computes a
    # second norm after rotation; keeping this identity prevents a detached
    # world-norm operand from entering the proof.
    if M.dot(world,world)!=m2:
        raise AssertionError('unit quaternion rotation lost magnetic norm')
    world_n=mag_norm.value

    a=state.accumulator
    if a.accepted_count>0 and a.weight_sum>F(1,10**6):
        running_mean_n=a.norm_sum/a.weight_sum
        if running_mean_n>cfg.mag_norm_min and cfg.max_norm_ratio>0:
            rel=abs(world_n-running_mean_n)/running_mean_n
            if rel>cfg.max_norm_ratio:
                return _reject(state,'running_norm_outlier')

    # Default weighting is exactly w=1, then clamp01(1)=1.
    if F(1)<cfg.min_sample_weight:
        return _reject(state,'min_sample_weight')

    accepted=ACC.accepted_sample(a,cfg.accumulator_config(),world_sample=world,
                                 world_norm=world_n,dt=d)
    mid=State(accepted.state,state.rejected_count,False,None,None,F(1),world)
    if not accepted.finalize_due:
        if mean_norm is not None or horizontal_sqrt is not None:
            raise ValueError('pre-finalize sample consumes no finalization sqrt witnesses')
        return StepResult(mid,False,True,None)

    timed_out=cfg.max_window>0 and accepted.state.accepted_window>=cfg.max_window
    if (cfg.min_effective_weight>0 and
        accepted.state.weight_sum<cfg.min_effective_weight and not timed_out):
        if mean_norm is not None or horizontal_sqrt is not None:
            raise ValueError('effective-weight-blocked finalize consumes no mean witnesses')
        return StepResult(mid,False,True,None)

    mean=ACC.mean(accepted.state)
    mean2=M.dot(mean,mean)
    if not isinstance(mean_norm,SqrtWitness) or mean_norm.radicand!=mean2:
        raise ValueError('mean norm witness detached from same accepted accumulator')
    if mean_norm.value<=cfg.mag_norm_min:
        return StepResult(mid,False,True,'mean_norm')

    h2=mean[0]*mean[0]+mean[1]*mean[1]
    if not isinstance(horizontal_sqrt,GAUGE.HorizontalSqrt) or horizontal_sqrt.radicand!=h2:
        raise ValueError('horizontal sqrt witness detached from same accepted mean')
    h=horizontal_sqrt.value
    if h<=cfg.mag_norm_min:
        return StepResult(mid,False,True,'horizontal_norm')
    if cfg.min_horizontal_fraction>0 and h/mean_norm.value<cfg.min_horizontal_fraction:
        return StepResult(mid,False,True,'horizontal_fraction')

    gauge=GAUGE.gauge_fix(mean,mag_norm_min=cfg.mag_norm_min,sqrt=horizontal_sqrt)
    if not gauge.ready:
        return StepResult(mid,False,True,'reference_norm')
    out=State(accepted.state,state.rejected_count,True,mean,gauge.world_reference,F(1),world)
    return StepResult(out,True,True,None)


def readiness():
    return {
      'default_ready_latch_precedes_new_sample_arithmetic':True,
      'tilt_quaternion_normalization_relation_materialized':True,
      'raw_and_rotated_mag_norm_gate_materialized':True,
      'running_mean_norm_outlier_gate_materialized':True,
      'default_unit_weight_and_min_weight_gate_materialized':True,
      'accepted_accumulator_and_finalize_gate_composed':True,
      'mean_norm_horizontal_fraction_and_gauge_fix_composed':True,
      'raw_sample_acceptance_boolean_removed':True,
      'input_allFinite_binary32_attached':False,
      'sqrt_and_quaternion_normalization_binary32_attached':False,
      'startup_tilt_quaternion_runtime_ancestry_attached':False,
      'startup_mag_sensor_source_ancestry_attached':False,
      'optional_quality_weighted_branch_attached':False,
      'optional_hard_iron_branch_attached':False,
      'complete_word_finite_identity':False,
      'ALT_STARTUP_PASS':False,
      'ALT_LIVE_PASS':False,
    }
