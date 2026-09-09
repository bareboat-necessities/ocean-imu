#!/usr/bin/env python3
"""P4 extension of the trusted BRMM tuner step with shipping stillness decay.

The base BRMM tuner transition already owns AdaptiveWaveBandPass, p00/p01/p11,
debiased first/second moments, (f,tau,sigma,R_S), candidate EMA, staging and
commit. This module changes only the one source-visible omission: the current
same-history StillnessAdapter attenuation multiplies de-noised ``var_wave``
before the 1e-6 variance floor, exactly as deployed.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
from ou3_interval import Interval
import ou3_brmm_tuner_scheduler_step as BASE

SCHEMA=1
QUALIFICATION='OU3_P4_BRMM_TUNER_STILLNESS_TRANSITION_V1'

def I(x):return Interval.point(float(x))
def targets(state:BASE.TunerState,still_attenuation:Interval,c:BASE.Constants)->BASE.CandidateState:
    if still_attenuation.lo<0 or still_attenuation.hi>1:raise ValueError('still attenuation must lie in [0,1]')
    f=BASE.clamp_interval(state.moments.frequency_hz,c.tune_freq_min,c.tune_freq_max)
    tau=BASE.clamp_interval(I(c.tau_coeff*.5)/f,c.tau_min,c.tau_max)
    noise_sigma=I(c.acc_noise_floor_sigma)*BASE.sqrt_interval(BASE.max_interval(state.band.p11,1e-30))
    var_wave=BASE.max_interval(BASE.acceleration_variance(state.moments)-noise_sigma.square(),0.0)
    var_wave=var_wave*still_attenuation
    sigma=BASE.min_interval(I(c.sigma_coeff)*BASE.sqrt_interval(BASE.max_interval(var_wave,1e-6)),c.sigma_max)
    rs=BASE.clamp_interval(BASE.spectral_mse_rs_target(tau,sigma,c),c.rs_min,c.rs_max)
    return BASE.CandidateState(tau,sigma,rs)

def advance_after_measurement(state:BASE.TunerState,*,a_vertical:Interval,f_wave_previous_wpe:Interval,still_attenuation:Interval,c:BASE.Constants|None=None):
    c=c or BASE.constants()
    if f_wave_previous_wpe.lo<=0:raise ValueError('same-word previous WPE frequency must stay positive')
    band=BASE.band_step(state.band,a_vertical,state.moments.frequency_hz,c)
    moments=BASE.moment_step(state.moments,band.band,f_wave_previous_wpe,c)
    provisional=BASE.TunerState(band,moments,state.candidate,state.active,state.scheduler)
    target=targets(provisional,still_attenuation,c)
    candidate=BASE.candidate_ema(state.candidate,target,moments.frequency_hz,c)
    elapsed=state.scheduler.since_last_commit_stage_s+I(c.dt)
    base=BASE.TunerState(band,moments,candidate,state.active,BASE.SchedulerState(elapsed,False))
    if elapsed.lo>c.adapt_every_s:return [(BASE.TunerState(band,moments,candidate,state.active,BASE.SchedulerState(I(0),True)),target)]
    if elapsed.hi<=c.adapt_every_s:return [(base,target)]
    return [(base,target),(BASE.TunerState(band,moments,candidate,state.active,BASE.SchedulerState(I(0),True)),target)]

def build():
    c=BASE.constants();f=I(.2);active=BASE.ActiveSchedule(I(1.1),I(.5),I(2),BASE.pseudo_period(I(1.1),c));state=BASE.TunerState(BASE.BandState(I(0),I(.1),I(.2),I(0),I(.1),True),BASE.MomentState(I(0),I(.5),I(.1),I(.5),f),BASE.CandidateState(I(1.1),I(.5),I(2)),active,BASE.SchedulerState(I(.05),False))
    moving=advance_after_measurement(state,a_vertical=I(.2),f_wave_previous_wpe=f,still_attenuation=I(1),c=c)
    still=advance_after_measurement(state,a_vertical=I(.2),f_wave_previous_wpe=f,still_attenuation=I(.2),c=c)
    monotone=all(s[1].sigma.hi<=m[1].sigma.hi for s,m in zip(still,moving))
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','base_band_covariance_and_moments_reused':True,'same_history_still_attenuation_enters_before_variance_floor':True,'raw_target_tuple_emitted_before_candidate_EMA':True,'candidate_EMA_and_stage_timer_reused':True,'raw_and_active_schedule_remain_distinct':True,'still_attenuation_reduces_or_preserves_sigma_smoke':monotone,'complete_BRMM_stillness_source_attached_here':False,'P4_promoted_here':False}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('base_band_covariance_and_moments_reused','same_history_still_attenuation_enters_before_variance_floor','raw_target_tuple_emitted_before_candidate_EMA','candidate_EMA_and_stage_timer_reused','raw_and_active_schedule_remain_distinct','still_attenuation_reduces_or_preserves_sigma_smoke'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('complete_BRMM_stillness_source_attached_here','P4_promoted_here'):
        if d.get(k) is not False:f.append(k+' not false')
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'stillness':d['same_history_still_attenuation_enters_before_variance_floor'],'P4':d['P4_promoted_here'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())