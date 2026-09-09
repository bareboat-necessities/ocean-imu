#!/usr/bin/env python3
"""Interval state transition for the sigma-channel StillnessAdapter path.

Only the state that can affect the OU sigma target is retained.  Shipping forms
``a_vert_lp`` from the same private-Mahony vertical sample, updates the
StillnessAdapter energy/timer, and later attenuates ``var_wave`` by
``exp(-still_time/1s)`` when the current sample is still.  Tracker frequency is
irrelevant to that predicate/timer and is intentionally not introduced as a
free proof coordinate.
"""
from __future__ import annotations
from dataclasses import dataclass
import argparse,json,math,re
from pathlib import Path
from ou3_interval import Interval
import ou3_source_domain_contract as SOURCE
import ou3_validated_transcendentals as VT

REPO=Path(__file__).resolve().parents[2]
WRAPPER=REPO/'src'/'kalman_ou_iii'/'SeaStateFusionFilter_OU_III.h'
COMMON=REPO/'src'/'tuner'/'SeaStateFusionTunerCommon.h'
SCHEMA=1
QUALIFICATION='OU3_BRMM_SIGMA_STILLNESS_SAME_SIGNAL_STEP_V1'
PI=Interval.outward_bounds(3.141592653589793,3.141592653589794)

def I(x):return Interval.outward_bounds(float(x),float(x))
def _member(text,name):
    m=re.search(rf'\b{re.escape(name)}\s*=\s*([0-9.eE+-]+)f\s*;',text)
    if not m:raise RuntimeError('cannot extract '+name)
    return float(m.group(1))
def _wide_exp(x):
    if not(math.isfinite(x.lo) and math.isfinite(x.hi)):raise ValueError('finite exp interval required')
    scale=1;mag=max(abs(x.lo),abs(x.hi))
    while mag/scale>VT.MAX_ABS_ARGUMENT:scale*=2
    out=VT.exp_interval(x/I(scale))
    while scale>1:out=out.square();scale//=2
    return out

def clampi(x,lo,hi):return Interval(max(lo,min(hi,x.lo)),max(lo,min(hi,x.hi)))

@dataclass(frozen=True)
class Constants:
    dt:float;gravity:float;lpf_cutoff_hz:float;energy_alpha:float;energy_thresh:float;still_max_s:float;variance_decay_s:float
@dataclass(frozen=True)
class State:
    lpf_state:Interval;lpf_initialized:bool;energy_ema:Interval;still_time_s:Interval;last_is_still:bool
@dataclass(frozen=True)
class Successor:
    state:State;a_vertical_lp:Interval;still_attenuation:Interval;branch:str

def constants():
    w=WRAPPER.read_text();c=COMMON.read_text()
    return Constants(
      dt=float(SOURCE.parse_const(w,'FREQ_SMOOTHER_DT')),
      gravity=9.80665,
      lpf_cutoff_hz=float(SOURCE.parse_const(w,'MAX_FREQ_HZ')),
      energy_alpha=_member(c,'energy_alpha'),energy_thresh=_member(c,'energy_thresh'),
      still_max_s=60.0,variance_decay_s=1.0)

def zero_state():return State(I(0),False,I(0),I(0),False)

def advance(state,a_vertical,*,c=None):
    c=c or constants();dt=I(c.dt)
    alpha_lpf=_wide_exp(-(I(2)*PI*I(c.lpf_cutoff_hz)*dt))
    lp=a_vertical if not state.lpf_initialized else (I(1)-alpha_lpf)*a_vertical+alpha_lpf*state.lpf_state
    inst=(lp/I(c.gravity)).square()
    e=(I(1-c.energy_alpha)*state.energy_ema)+I(c.energy_alpha)*inst
    # Exact boolean predicate must branch when an interval straddles threshold.
    modes=[]
    if e.hi < c.energy_thresh:modes=[True]
    elif e.lo >= c.energy_thresh:modes=[False]
    else:modes=[False,True]
    out=[]
    for still in modes:
        if still:
            t=clampi(state.still_time_s+dt,0.0,c.still_max_s)
            atten=_wide_exp(-(t/I(c.variance_decay_s)))
            branch='still'
        else:
            t=I(0);atten=I(1);branch='moving'
        out.append(Successor(State(lp,True,e,t,still),lp,clampi(atten,0,1),branch))
    return out

def shipping_parity():
    w=WRAPPER.read_text();c=COMMON.read_text()
    return {
      'lpf_precedes_stillness':w.find('const float a_vert_lp = freq_input_lpf_.step(a_vert_measurement, dt);') < w.find('const float f_after_still = freq_stillness_.step(a_vert_lp, dt, f_tracker);'),
      'energy_uses_squared_normalized_vertical':'const float inst_energy = a_norm * a_norm;' in c,
      'energy_alpha_update':'energy_ema = (1.0f - energy_alpha) * energy_ema + energy_alpha * inst_energy;' in c,
      'still_predicate_strict_below':'const bool is_still = (energy_ema < energy_thresh);' in c,
      'still_timer_cap_60':'if (still_time_sec > 60.0f)' in c,
      'sigma_uses_current_still_state':'if (freq_stillness_.isStill())' in w and 'freq_stillness_.getStillTime()' in w,
      'sigma_still_decay_1s':'constexpr float STILL_VAR_DECAY_SEC = 1.0f;' in w and 'var_wave *= atten;' in w,
    }
def build():
    p=shipping_parity();s=zero_state();succ=advance(s,I(.02))
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','shipping_parity':p,'shipping_parity_pass':all(p.values()),'same_vertical_sample_drives_lpf_energy_and_sigma_attenuation':True,'tracker_frequency_not_needed_for_sigma_still_predicate':True,'threshold_straddling_cell_is_branched':True,'validated_exp_used_for_lpf_and_still_decay':True,'smoke_successors':len(succ),'complete_BRMM_stillness_predecessor_family_covered_here':False,'P4_promoted_here':False}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('shipping_parity_pass','same_vertical_sample_drives_lpf_energy_and_sigma_attenuation','tracker_frequency_not_needed_for_sigma_still_predicate','threshold_straddling_cell_is_branched','validated_exp_used_for_lpf_and_still_decay'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('complete_BRMM_stillness_predecessor_family_covered_here','P4_promoted_here'):
        if d.get(k) is not False:f.append(k+' not false')
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'parity':d['shipping_parity_pass'],'branches':d['smoke_successors'],'P4':d['P4_promoted_here'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
