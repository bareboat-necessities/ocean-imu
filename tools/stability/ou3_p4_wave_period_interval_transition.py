#!/usr/bin/env python3
"""Outward WavePeriodEstimator state transition for the P4 BRMM source cover.

This is the executable state recurrence missing from the earlier steady-transfer
frontend lemma.  A caller supplies one SAME-SOURCE interval for the private
complementary vertical-acceleration input.  The deployed two high-pass stages,
leaky velocity/elevation states and EW moment states are then propagated with
outward interval arithmetic.

The moment-start predicate is branch preserving: a cell wholly before 3/lambda
has no moment update, a cell wholly after it applies the EW update, and a cell
straddling the threshold returns both branches.  The moment horizon is retained
inside the deployed [20,180] s guards rather than frozen to a replay value.

This module intentionally stops before converting moment boxes to a valid
period.  An independent interval box for velocity/elevation moments generally
contains combinations with zero variance, so taking their ratio would invent a
false singularity.  The missing theorem step must retain enough SAME-HISTORY
correlation to prove positive elevation/velocity variance and enclose the raw
moment-ratio/log-period update.  Likewise this module does not prove that the
physical BRMM motion plus private Mahony observer produces a particular input
interval; that source attachment remains separate.
"""
from __future__ import annotations

from dataclasses import dataclass
import argparse
import json
import math
from pathlib import Path
import re

from ou3_interval import Interval
import ou3_validated_transcendentals as VT
import ou3_brmm_wave_period_frontend as FRONTEND

REPO=Path(__file__).resolve().parents[2]
ESTIMATOR=REPO/"src"/"tuner"/"WavePeriodEstimator.h"
SCHEMA=2
QUALIFICATION="OU3_P4_WAVE_PERIOD_ESTIMATOR_INTERVAL_MOMENT_TRANSITION_V2"


def I(x:float)->Interval:return Interval.point(float(x))
def finite(x:Interval)->bool:return isinstance(x,Interval) and math.isfinite(x.lo) and math.isfinite(x.hi) and x.lo<=x.hi

def _ctor_default(text:str,name:str)->float:
    m=re.search(rf"\b{name}\s*=\s*([0-9.eE+-]+)f",text)
    if not m:raise RuntimeError("cannot extract WavePeriodEstimator default "+name)
    return float(m.group(1))


def constants()->dict:
    text=ESTIMATOR.read_text(encoding="utf-8")
    hp=_ctor_default(text,"high_pass_hz")
    mh=_ctor_default(text,"moment_horizon_periods")
    lh=_ctor_default(text,"log_smoothing_periods")
    hmin=_ctor_default(text,"min_horizon_sec")
    hmax=_ctor_default(text,"max_horizon_sec")
    pi=Interval.outward_bounds(3.141592653589793,3.141592653589794)
    lam=I(2.0)*pi*I(max(1e-4,hp))
    return {"high_pass_hz":hp,"moment_horizon_periods":mh,"log_smoothing_periods":lh,
            "min_horizon_sec":hmin,"max_horizon_sec":hmax,"lambda":lam}


@dataclass(frozen=True)
class WavePeriodState:
    accel_prev:Interval
    high_pass_1:Interval
    high_pass_1_prev:Interval
    high_pass_2:Interval
    velocity:Interval
    elevation:Interval
    velocity_mean:Interval
    velocity_sq:Interval
    elevation_mean:Interval
    elevation_sq:Interval
    weight:Interval
    elapsed_sec:Interval


def zero_state()->WavePeriodState:
    z=I(0.0)
    return WavePeriodState(z,z,z,z,z,z,z,z,z,z,z,z)


def validate_state(s:WavePeriodState)->list[str]:
    f=[]
    for name,x in s.__dict__.items():
        if not finite(x):f.append(name+" nonfinite")
    if finite(s.weight) and (s.weight.lo < -1e-300 or s.weight.hi>1.0+1e-12):f.append("weight outside [0,1] enclosure")
    if finite(s.elapsed_sec) and s.elapsed_sec.lo < -1e-300:f.append("negative elapsed")
    return f


def _decay_gain(dt:Interval,lam:Interval)->tuple[Interval,Interval]:
    if not finite(dt) or dt.lo<=0 or not finite(lam) or lam.lo<=0:raise ValueError("invalid dt/lambda")
    decay=VT.exp_interval(-(lam*dt))
    gain=(I(1.0)-decay)/lam
    return decay,gain


def linear_frontend_step(s:WavePeriodState,dt:Interval,input_accel:Interval)->WavePeriodState:
    failures=validate_state(s)
    if failures:raise ValueError("invalid predecessor state: "+repr(failures))
    if not finite(input_accel):raise ValueError("same-source input interval required")
    c=constants();decay,gain=_decay_gain(dt,c["lambda"])
    stage1=decay*(s.high_pass_1+input_accel-s.accel_prev)
    stage2=decay*(s.high_pass_2+stage1-s.high_pass_1_prev)
    velocity=decay*s.velocity+gain*stage2
    elevation=decay*s.elevation+gain*velocity
    return WavePeriodState(
        input_accel,stage1,stage1,stage2,velocity,elevation,
        s.velocity_mean,s.velocity_sq,s.elevation_mean,s.elevation_sq,s.weight,
        s.elapsed_sec+dt)


def _alpha_for_horizon(dt:Interval,horizon:Interval)->Interval:
    if not finite(horizon) or horizon.lo<=0:raise ValueError("invalid horizon")
    return I(1.0)-VT.exp_interval(-(dt/horizon))


def moment_update(s:WavePeriodState,dt:Interval,horizon:Interval)->WavePeriodState:
    failures=validate_state(s)
    if failures:raise ValueError("invalid moment predecessor: "+repr(failures))
    alpha=_alpha_for_horizon(dt,horizon)
    one_minus=I(1.0)-alpha
    weight=one_minus*s.weight+alpha
    vm=one_minus*s.velocity_mean+alpha*s.velocity
    vs=one_minus*s.velocity_sq+alpha*s.velocity.square()
    em=one_minus*s.elevation_mean+alpha*s.elevation
    es=one_minus*s.elevation_sq+alpha*s.elevation.square()
    return WavePeriodState(s.accel_prev,s.high_pass_1,s.high_pass_1_prev,s.high_pass_2,
                           s.velocity,s.elevation,vm,vs,em,es,weight,s.elapsed_sec)


def step_branches(s:WavePeriodState,dt:Interval,input_accel:Interval)->list[tuple[str,WavePeriodState]]:
    lin=linear_frontend_step(s,dt,input_accel)
    c=constants(); start=I(3.0)/c["lambda"]
    definitely_before=lin.elapsed_sec.hi < start.lo
    definitely_after=lin.elapsed_sec.lo >= start.hi
    if definitely_before:return [("pre_moment",lin)]
    horizon=Interval(c["min_horizon_sec"],c["max_horizon_sec"])
    updated=moment_update(lin,dt,horizon)
    if definitely_after:return [("moment",updated)]
    return [("pre_moment",lin),("moment",updated)]


def build()->dict:
    frontend=FRONTEND.build(REPO);ff=FRONTEND.validate(frontend)
    if ff:raise RuntimeError("frontend prerequisite failed: "+repr(ff))
    c=constants(); dt=Interval(*map(float,frontend["declared_inputs"]["imu_dt_s"]))
    z=zero_state()
    one=step_branches(z,dt,Interval.outward_bounds(-1.0,1.0))
    smoke=len(one)==1 and one[0][0]=="pre_moment" and not validate_state(one[0][1])
    start=I(3.0)/c["lambda"]
    e=WavePeriodState(z.accel_prev,z.high_pass_1,z.high_pass_1_prev,z.high_pass_2,z.velocity,z.elevation,
                      z.velocity_mean,z.velocity_sq,z.elevation_mean,z.elevation_sq,z.weight,
                      Interval(max(0.0,start.lo-dt.hi*1.5),start.hi))
    split={name for name,_ in step_branches(e,dt,I(0.0))}=={"pre_moment","moment"}
    return {
        "schema":SCHEMA,"qualification":QUALIFICATION,
        "canonical_source":"COMPLETE_BRMM_NORMAL_LIVE_WORD",
        "shipping_two_high_pass_and_leaky_integrator_transition_materialized":True,
        "EW_first_second_moment_transition_materialized":True,
        "moment_horizon_full_shipping_guard_interval_retained":True,
        "moment_start_branch_preserved":True,
        "moment_start_straddling_requires_source_split":split,
        "same_source_input_interval_required":True,
        "exact_structural_zero_not_outward_widened":True,
        "physical_BRMM_to_private_complementary_input_attached_here":False,
        "positive_correlated_moment_variance_proved_here":False,
        "raw_period_ratio_transition_materialized_here":False,
        "log_period_transition_materialized_here":False,
        "usable_period_latch_transition_materialized_here":False,
        "independent_moment_boxes_may_be_divided_to_make_period":False,
        "point_smoke_pass":smoke,
        "complete_wave_period_source_transition_closed_here":False,
        "P4_promoted_here":False,
        "next_obligation":(
            "attach the private complementary vertical-acceleration input to the admitted BRMM physical state, then retain same-history moment correlation strongly enough to prove positive velocity/elevation variances and enclose the raw-period/log-period/usable-latch branches"
        ),
    }


def validate(d:dict)->list[str]:
    f=[]
    if d.get("schema")!=SCHEMA or d.get("qualification")!=QUALIFICATION:f.append("schema/qualification mismatch")
    for k in ("shipping_two_high_pass_and_leaky_integrator_transition_materialized","EW_first_second_moment_transition_materialized","moment_horizon_full_shipping_guard_interval_retained","moment_start_branch_preserved","moment_start_straddling_requires_source_split","same_source_input_interval_required","exact_structural_zero_not_outward_widened","point_smoke_pass"):
        if d.get(k) is not True:f.append(k+" not true")
    for k in ("physical_BRMM_to_private_complementary_input_attached_here","positive_correlated_moment_variance_proved_here","raw_period_ratio_transition_materialized_here","log_period_transition_materialized_here","usable_period_latch_transition_materialized_here","independent_moment_boxes_may_be_divided_to_make_period","complete_wave_period_source_transition_closed_here","P4_promoted_here"):
        if d.get(k) is not False:f.append(k+" not false")
    return f


def main()->int:
    ap=argparse.ArgumentParser();ap.add_argument("--output",type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d["validation_pass"]=not f;d["validation_failures"]=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"linear":d["shipping_two_high_pass_and_leaky_integrator_transition_materialized"],"moments":d["EW_first_second_moment_transition_materialized"],"ratio":d["raw_period_ratio_transition_materialized_here"],"split":d["moment_start_straddling_requires_source_split"],"failures":f},sort_keys=True))
    return int(bool(f))
if __name__=="__main__":raise SystemExit(main())
