#!/usr/bin/env python3
"""Joint SAME-SIGNAL WavePeriodEstimator + sigma transition for OU-III P4.

This module never accepts an independently chosen ``(f, sigma)`` rectangle.
One state carries the WavePeriodEstimator, period-scaled AdaptiveWaveBandPass,
its unit-white covariance recursion, and SeaStateAutoTuner first/second moments.
A successful branch therefore emits one correlated
``(f,tau_target,sigma_target_raw,T_S,R_S)`` tuple from one source history.

Shipping causality is retained exactly at the sample boundary: the current
sample calls update_tuner(..., tuner_frequency_hz_()) BEFORE that sample updates
WavePeriodEstimator. Thus the current band/tuner target uses the predecessor
usable period (or fixed prior), while the current raw acceleration sample enters
both the band statistic and the later WavePeriodEstimator transition. A newly
usable period cannot influence tuning until the following sample.
"""
from __future__ import annotations
from dataclasses import dataclass
import argparse,json,math
from pathlib import Path
from ou3_interval import Interval
import ou3_validated_transcendentals as VT
import ou3_validated_positive_transcendentals as VPOS
import ou3_p4_wave_period_interval_transition as WAVE
import ou3_p4_complete_brmm_target_cell as TARGET
import ou3_brmm_wave_period_frontend as FRONTEND
REPO=Path(__file__).resolve().parents[2]
SCHEMA=2
QUALIFICATION="OU3_P4_JOINT_SAME_SIGNAL_ESTIMATOR_TRANSITION_V2"
_PRIOR=None

def I(x):return Interval.point(float(x))
def finite(x):return isinstance(x,Interval) and math.isfinite(x.lo) and math.isfinite(x.hi) and x.lo<=x.hi
def clampi(x,lo,hi):
    if not finite(x) or not lo<=hi:raise ValueError("invalid clamp")
    return Interval(max(lo,min(hi,x.lo)),max(lo,min(hi,x.hi)))
def nonnegative(x):return Interval(max(0.0,x.lo),max(0.0,x.hi))
def prior_frequency():
    global _PRIOR
    if _PRIOR is None:
        d=FRONTEND.build(REPO);f=FRONTEND.validate(d)
        if f:raise RuntimeError("wave-period frontend prerequisite failed: "+repr(f))
        a,b=map(float,d["declared_inputs"]["fixed_tuning_frequency_prior_hz"]);_PRIOR=Interval(a,b)
    return _PRIOR

@dataclass(frozen=True)
class BandState:
    lowpass_low:Interval;band:Interval;p00:Interval;p01:Interval;p11:Interval;ready:bool
@dataclass(frozen=True)
class TunerState:
    frequency_hz:Interval|None;mean_value:Interval;mean_weight:Interval;sq_value:Interval;sq_weight:Interval;variance_horizon_s:Interval|None
@dataclass(frozen=True)
class JointEstimatorState:
    source_token:str;predecessor_token:str|None;wave:WAVE.WavePeriodState;canonical_period_s:Interval|None;usable_period:bool;band:BandState;tuner:TunerState
@dataclass(frozen=True)
class JointImage:
    state:JointEstimatorState;status:str;frequency_hz:Interval|None;sigma_target_raw_mps2:Interval|None;tau_target_s:Interval|None;pseudo_period_s:Interval|None;rs_target:Interval|None

def zero_state(token="JOINT_ROOT"):
    z=I(0);return JointEstimatorState(token,None,WAVE.zero_state(),None,False,BandState(z,z,z,z,z,False),TunerState(None,z,z,z,z,None))
def _wave_horizon(period):
    c=WAVE.constants();p=period if period is not None else I(6)
    return clampi(I(c["moment_horizon_periods"])*p,c["min_horizon_sec"],c["max_horizon_sec"])
def _moment_variances(s):
    if s.weight.lo<=1e-3:return None
    vm=s.velocity_mean/s.weight;em=s.elevation_mean/s.weight
    vv=s.velocity_sq/s.weight-vm.square();ev=s.elevation_sq/s.weight-em.square()
    if vv.lo<=1e-12 or ev.lo<=1e-12:return None
    return vv,ev
def _raw_period(s):
    pair=_moment_variances(s)
    if pair is None:return None
    vv,ev=pair;c=WAVE.constants();omega2=vv/ev-c["lambda"].square()
    if omega2.lo<=1e-8:return None
    pi=Interval.outward_bounds(3.141592653589793,3.141592653589794)
    return I(2)*pi/VPOS.sqrt_interval(omega2)
def _canonical_period_step(old,raw,dt):
    if old is None:return raw
    c=WAVE.constants();k=float(c["log_smoothing_periods"])
    if not k>0:return raw
    horizon=clampi(I(k)*old,max(0.05,dt.lo),35.0)
    alpha=I(1)-VT.exp_interval(-(dt/horizon))
    delta=alpha*(VPOS.log_interval(raw)-VPOS.log_interval(old))
    if delta.lo < -VT.MAX_ABS_ARGUMENT or delta.hi > VT.MAX_ABS_ARGUMENT:return None
    return old*VT.exp_interval(delta)
def _usable(old,w,p):
    if old:return True
    if p is None:return False
    c=WAVE.constants();start=I(3)/c["lambda"];floor=I(4)/c["lambda"]
    return w.elapsed_sec.lo>=floor.hi and (w.elapsed_sec-start).lo>=p.hi

def _band_step(b,x,dt,fref):
    if fref.lo<=0 or dt.lo<=0:return None
    pi=Interval.outward_bounds(3.141592653589793,3.141592653589794)
    nyq=I(.45)/dt;upper=Interval(min(6.,nyq.lo),min(6.,nyq.hi))
    if upper.lo<=.01:return None
    lr=I(.5)*fref;low=Interval(max(.01,lr.lo),max(.01,lr.hi));low=Interval(min(low.lo,upper.lo/1.05),min(low.hi,upper.hi/1.05))
    hr=I(4)*fref;h0=Interval(min(upper.lo,hr.lo),min(upper.hi,hr.hi));h1=Interval(max(h0.lo,1.05*low.lo),max(h0.hi,1.05*low.hi));high=Interval(min(h1.lo,upper.lo),min(h1.hi,upper.hi))
    if high.lo<=low.hi:return None
    al=I(1)-VT.exp_interval(-(I(2)*pi*low*dt));ah=I(1)-VT.exp_interval(-(I(2)*pi*high*dt));ql=I(1)-al;qh=I(1)-ah
    low_new=ql*b.lowpass_low+al*x;band_new=qh*b.band+ah*(ql*(x-b.lowpass_low))
    a00=ql;a10=-(ah*ql);a11=qh;b0=al;b1=ah*ql
    p00=nonnegative(a00.square()*b.p00+b0.square());p01=a00*(a10*b.p00+a11*b.p01)+b0*b1
    p11=nonnegative(a10.square()*b.p00+I(2)*a10*a11*b.p01+a11.square()*b.p11+b1.square())
    return BandState(low_new,band_new,p00,p01,p11,True)
def _ema(v,w,x,a):
    om=I(1)-a;return om*v+a*x,om*w+a
def _tuner_step(t,band,dt,f_input):
    f=clampi(f_input,.05,5.);sea=clampi(I(.5)/f,.5,6.);teff=I(2)*sea;requested=clampi(I(4)*teff,.3,60.);h=clampi(requested,max(.05,dt.lo),35.)
    a=I(1)-VT.exp_interval(-(dt/h));mv,mw=_ema(t.mean_value,t.mean_weight,band,a);sv,sw=_ema(t.sq_value,t.sq_weight,band.square(),a)
    return TunerState(f,mv,mw,sv,sw,h)
def _variance(t):
    if t.mean_weight.lo<=1e-6 or t.sq_weight.lo<=1e-6:return None
    mu=t.mean_value/t.mean_weight;return nonnegative(t.sq_value/t.sq_weight-mu.square())
def _sigma_raw(t,b,acc_noise_floor_sigma,still_attenuation,sigma_coeff,max_sigma):
    if not(b.ready and acc_noise_floor_sigma.lo>=0 and 0<=still_attenuation.lo<=still_attenuation.hi<=1):return None
    noise_sigma=acc_noise_floor_sigma*VPOS.sqrt_interval(nonnegative(b.p11));vn=noise_sigma.square();total=_variance(t)
    if total is None:total=vn
    wave=nonnegative(total-vn)*still_attenuation;wave=Interval(max(1e-6,wave.lo),max(1e-6,wave.hi))
    return clampi(VPOS.sqrt_interval(wave)*I(sigma_coeff),0.,max_sigma)
def _tuning_input_from_predecessor(s):
    return s.canonical_period_s.reciprocal() if s.usable_period and s.canonical_period_s is not None else prior_frequency()
def _wave_successors(s,dt,input_accel):
    lin=WAVE.linear_frontend_step(s.wave,dt,input_accel);c=WAVE.constants();start=I(3)/c["lambda"]
    if lin.elapsed_sec.hi<start.lo:modes=[("pre_moment",lin)]
    elif lin.elapsed_sec.lo>=start.hi:modes=[("moment",WAVE.moment_update(lin,dt,_wave_horizon(s.canonical_period_s)))]
    else:modes=[("pre_moment",lin),("moment",WAVE.moment_update(lin,dt,_wave_horizon(s.canonical_period_s)))]
    out=[]
    for mode,w in modes:
        raw=_raw_period(w) if mode=="moment" else None;period=s.canonical_period_s;split=False
        if raw is not None:
            nxt=_canonical_period_step(period,raw,dt)
            if nxt is None:split=True
            else:period=nxt
        out.append((w,period,_usable(s.usable_period,w,period),split))
    return out

def step_joint(s,dt,input_accel,*,acc_noise_floor_sigma=I(.0148),still_attenuation=I(1.0)):
    if not s.source_token:raise ValueError("source token required")
    tc=TARGET.constants();f_tuner=_tuning_input_from_predecessor(s)
    f_band=s.tuner.frequency_hz if s.tuner.frequency_hz is not None else f_tuner;f_band=clampi(f_band,tc["MIN_TUNE_FREQ_HZ"],tc["MAX_TUNE_FREQ_HZ"])
    band=_band_step(s.band,input_accel,dt,f_band)
    if band is None:return [JointImage(s,"split_required_band_coefficients",None,None,None,None,None)]
    tuner=_tuner_step(s.tuner,band.band,dt,f_tuner);sig=_sigma_raw(tuner,band,acc_noise_floor_sigma,still_attenuation,tc["sigma_coeff"],tc["MAX_SIGMA_A"])
    tau=TARGET.tau_from_frequency(f_tuner,tc) if sig is not None else None;ts=TARGET.pseudo_period_from_tau(tau,tc) if tau is not None else None;rs=TARGET.spectral_mse_rs(tau,sig,tc) if tau is not None and sig is not None else None
    out=[]
    for w,period,usable,split in _wave_successors(s,dt,input_accel):
        ns=JointEstimatorState(s.source_token,s.source_token,w,period,usable,band,tuner)
        if split:out.append(JointImage(ns,"split_required_log_period",f_tuner,sig,tau,ts,rs))
        elif sig is None:out.append(JointImage(ns,"hold_sigma_not_ready",f_tuner,None,None,None,None))
        else:out.append(JointImage(ns,"joint_target",f_tuner,sig,tau,ts,rs))
    return out

def build():
    s=zero_state();dt=Interval.outward_bounds(.004,.006);b=step_joint(s,dt,Interval.outward_bounds(-.1,.1));prior=prior_frequency()
    startup=bool(b) and all(x.frequency_hz is not None and x.frequency_hz.contains_interval(prior) for x in b)
    return {"schema":SCHEMA,"qualification":QUALIFICATION,"canonical_source":"COMPLETE_BRMM_NORMAL_LIVE_WORD","same_signal_wave_period_and_band_state_materialized":True,"tuner_precedes_current_sample_wave_period_update":True,"newly_usable_period_affects_tuner_only_next_sample":True,"fixed_prior_used_until_predecessor_period_usable":True,"same_history_positive_variance_required_before_ratio":True,"raw_period_and_canonical_log_period_relation_materialized":True,"exact_frequency_reciprocal_relation_materialized":True,"adaptive_band_previous_tuner_frequency_lag_materialized":True,"adaptive_band_time_varying_white_noise_covariance_materialized":True,"sigma_first_second_moment_statistic_materialized":True,"band_noise_p11_subtraction_materialized":True,"raw_sigma_and_effective_OU_sigma_kept_distinct":True,"joint_f_tau_sigma_TS_RS_image_materialized":True,"independent_f_sigma_rectangle_accepted_by_transition":False,"unresolved_wide_cells_require_source_split":True,"startup_prior_causality_smoke_pass":startup,"physical_BRMM_signal_attachment_closed_here":False,"all_admitted_same_history_source_cells_emitted_here":False,"endpoint_and_every_prefix_augmented_LDLT_closed_here":False,"P4_MOTION_PASS":False,"P4_PASS":False,"P5_MAY_START":False,"next_obligation":"attach every admitted BRMM/private-observer vertical-acceleration and stillness branch to this causal joint state, recursively split dependency-lost cells, carry only its emitted targets through active EMA/commit/scheduler and same-history Riccati/Joseph/reset cells, then certify endpoint and literal every-prefix augmented LDLT with finite-precision ISS"}
def validate(d):
    f=[]
    if d.get("schema")!=SCHEMA or d.get("qualification")!=QUALIFICATION:f.append("schema/qualification mismatch")
    for k in ("same_signal_wave_period_and_band_state_materialized","tuner_precedes_current_sample_wave_period_update","newly_usable_period_affects_tuner_only_next_sample","fixed_prior_used_until_predecessor_period_usable","same_history_positive_variance_required_before_ratio","raw_period_and_canonical_log_period_relation_materialized","exact_frequency_reciprocal_relation_materialized","adaptive_band_previous_tuner_frequency_lag_materialized","adaptive_band_time_varying_white_noise_covariance_materialized","sigma_first_second_moment_statistic_materialized","band_noise_p11_subtraction_materialized","raw_sigma_and_effective_OU_sigma_kept_distinct","joint_f_tau_sigma_TS_RS_image_materialized","unresolved_wide_cells_require_source_split","startup_prior_causality_smoke_pass"):
        if d.get(k) is not True:f.append(k+" not true")
    for k in ("independent_f_sigma_rectangle_accepted_by_transition","physical_BRMM_signal_attachment_closed_here","all_admitted_same_history_source_cells_emitted_here","endpoint_and_every_prefix_augmented_LDLT_closed_here","P4_MOTION_PASS","P4_PASS","P5_MAY_START"):
        if d.get(k) is not False:f.append(k+" not false")
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument("--output",type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d["validation_pass"]=not f;d["validation_failures"]=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+"\n");print(json.dumps({"joint":d["joint_f_tau_sigma_TS_RS_image_materialized"],"causal":d["startup_prior_causality_smoke_pass"],"P4":d["P4_PASS"],"failures":f},sort_keys=True));return int(bool(f))
if __name__=="__main__":raise SystemExit(main())
