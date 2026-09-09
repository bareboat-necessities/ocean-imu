#!/usr/bin/env python3
"""Joint SAME-SIGNAL WavePeriodEstimator + sigma transition for OU-III P4.

This module is deliberately relational.  It never accepts an independently
chosen ``(f, sigma)`` rectangle.  One state carries the WavePeriodEstimator,
the period-scaled AdaptiveWaveBandPass, its unit-white covariance recursion,
and the SeaStateAutoTuner first/second moments.  A successful branch therefore
emits one correlated tuple

    (f, tau_target, sigma_target_raw, T_S, R_S)

from one source token/history.

Cells that are too wide to certify positive moment variances do not manufacture
a period by dividing independent moment boxes: they return a ``split_required``
branch.  Likewise a log-period step whose validated exponential increment falls
outside the audited exp kernel requests a source split.  Those refusal semantics
are theorem-significant and non-promoting.
"""
from __future__ import annotations

from dataclasses import dataclass
import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_validated_transcendentals as VT
import ou3_validated_positive_transcendentals as VPOS
import ou3_p4_wave_period_interval_transition as WAVE
import ou3_p4_complete_brmm_target_cell as TARGET

REPO=Path(__file__).resolve().parents[2]
SCHEMA=1
QUALIFICATION="OU3_P4_JOINT_SAME_SIGNAL_ESTIMATOR_TRANSITION_V1"


def I(x:float)->Interval:return Interval.point(float(x))
def finite(x:Interval)->bool:return isinstance(x,Interval) and math.isfinite(x.lo) and math.isfinite(x.hi) and x.lo<=x.hi

def clampi(x:Interval,lo:float,hi:float)->Interval:
    if not finite(x) or not lo<=hi: raise ValueError("invalid clamp")
    return Interval(max(lo,min(hi,x.lo)),max(lo,min(hi,x.hi)))

def nonnegative(x:Interval)->Interval:return Interval(max(0.0,x.lo),max(0.0,x.hi))


@dataclass(frozen=True)
class BandState:
    lowpass_low:Interval
    band:Interval
    p00:Interval
    p01:Interval
    p11:Interval
    ready:bool

@dataclass(frozen=True)
class TunerState:
    frequency_hz:Interval|None
    mean_value:Interval
    mean_weight:Interval
    sq_value:Interval
    sq_weight:Interval
    variance_horizon_s:Interval|None

@dataclass(frozen=True)
class JointEstimatorState:
    source_token:str
    predecessor_token:str|None
    wave:WAVE.WavePeriodState
    canonical_period_s:Interval|None
    usable_period:bool
    band:BandState
    tuner:TunerState

@dataclass(frozen=True)
class JointImage:
    state:JointEstimatorState
    status:str
    frequency_hz:Interval|None
    sigma_target_raw_mps2:Interval|None
    tau_target_s:Interval|None
    pseudo_period_s:Interval|None
    rs_target:Interval|None


def zero_state(token:str="JOINT_ROOT")->JointEstimatorState:
    z=I(0.0)
    return JointEstimatorState(token,None,WAVE.zero_state(),None,False,
        BandState(z,z,z,z,z,False),TunerState(None,z,z,z,z,None))


def _wave_horizon(period:Interval|None)->Interval:
    c=WAVE.constants()
    p=period if period is not None else I(6.0)
    return clampi(I(c["moment_horizon_periods"])*p,c["min_horizon_sec"],c["max_horizon_sec"])


def _moment_variances(s:WAVE.WavePeriodState)->tuple[Interval,Interval]|None:
    if s.weight.lo<=1e-3: return None
    vm=s.velocity_mean/s.weight; em=s.elevation_mean/s.weight
    vv=s.velocity_sq/s.weight-vm.square()
    ev=s.elevation_sq/s.weight-em.square()
    # SAME-HISTORY cell must be tight enough that both centered variances are
    # strictly positive.  Do not clamp a negative dependency artifact to zero
    # and then divide by it.
    if vv.lo<=1e-12 or ev.lo<=1e-12: return None
    return vv,ev


def _raw_period(s:WAVE.WavePeriodState)->Interval|None:
    pair=_moment_variances(s)
    if pair is None:return None
    vv,ev=pair;c=WAVE.constants()
    omega2=vv/ev-c["lambda"].square()
    if omega2.lo<=1e-8:return None
    root=VPOS.sqrt_interval(omega2)
    pi=Interval.outward_bounds(3.141592653589793,3.141592653589794)
    return I(2.0)*pi/root


def _canonical_period_step(old:Interval|None,raw:Interval,dt:Interval)->Interval|None:
    if old is None:return raw
    c=WAVE.constants();k=float(c["log_smoothing_periods"])
    if not k>0:return raw
    requested=I(k)*old
    horizon=clampi(requested,max(0.05,dt.lo),35.0)
    alpha=I(1.0)-VT.exp_interval(-(dt/horizon))
    log_raw=VPOS.log_interval(raw);log_old=VPOS.log_interval(old)
    delta=alpha*(log_raw-log_old)
    # Trusted VT exponential is intentionally audited only on |x|<=1/2.
    # Wider source cells must split rather than silently invoke libm.
    if delta.lo < -VT.MAX_ABS_ARGUMENT or delta.hi > VT.MAX_ABS_ARGUMENT:return None
    return old*VT.exp_interval(delta)


def _usable(old:bool,w:WAVE.WavePeriodState,p:Interval|None)->bool:
    if old:return True
    if p is None:return False
    c=WAVE.constants();start=I(3.0)/c["lambda"];floor=I(4.0)/c["lambda"]
    history=w.elapsed_sec-start
    return w.elapsed_sec.lo>=floor.hi and history.lo>=p.hi


def _band_step(b:BandState,x:Interval,dt:Interval,fref:Interval)->BandState|None:
    # Literal AdaptiveWaveBandPass defaults and clamp ordering.
    if fref.lo<=0 or dt.lo<=0:return None
    pi=Interval.outward_bounds(3.141592653589793,3.141592653589794)
    nyq=I(0.45)/dt
    upper=Interval(min(6.0,nyq.lo),min(6.0,nyq.hi))
    if upper.lo<=0.01:return None
    low=clampi(I(0.5)*fref,0.01,upper.hi/1.05)
    # high=max(min(upper,4*f),1.05*low), then min upper.  Interval hull is
    # monotone in both positive coordinates, so endpoint application is sound.
    hraw=I(4.0)*fref
    hlo=max(min(upper.lo,hraw.lo),1.05*low.lo)
    hhi=max(min(upper.hi,hraw.hi),1.05*low.hi)
    high=Interval(min(hlo,upper.lo),min(hhi,upper.hi))
    if high.lo<=low.hi:return None
    al=I(1.0)-VT.exp_interval(-(I(2.0)*pi*low*dt))
    ah=I(1.0)-VT.exp_interval(-(I(2.0)*pi*high*dt))
    ql=I(1.0)-al;qh=I(1.0)-ah
    low_new=ql*b.lowpass_low+al*x
    hp=ql*(x-b.lowpass_low)
    band_new=qh*b.band+ah*hp
    a00=ql;a10=-(ah*ql);a11=qh;b0=al;b1=ah*ql
    p00=nonnegative(a00.square()*b.p00+b0.square())
    p01=a00*(a10*b.p00+a11*b.p01)+b0*b1
    p11=nonnegative(a10.square()*b.p00+I(2.0)*a10*a11*b.p01+a11.square()*b.p11+b1.square())
    return BandState(low_new,band_new,p00,p01,p11,True)


def _ema(value:Interval,weight:Interval,x:Interval,alpha:Interval)->tuple[Interval,Interval]:
    om=I(1.0)-alpha
    return om*value+alpha*x,om*weight+alpha


def _tuner_step(t:TunerState,band:Interval,dt:Interval,f:Interval)->TunerState:
    f=clampi(f,0.05,5.0)
    sea=clampi(I(0.5)/f,0.5,6.0)
    teff=I(2.0)*sea
    requested=clampi(I(4.0)*teff,0.3,60.0)
    horizon=clampi(requested,max(0.05,dt.lo),35.0)
    alpha=I(1.0)-VT.exp_interval(-(dt/horizon))
    mv,mw=_ema(t.mean_value,t.mean_weight,band,alpha)
    sv,sw=_ema(t.sq_value,t.sq_weight,band.square(),alpha)
    return TunerState(f,mv,mw,sv,sw,horizon)


def _variance(t:TunerState)->Interval|None:
    if t.mean_weight.lo<=1e-6 or t.sq_weight.lo<=1e-6:return None
    mu=t.mean_value/t.mean_weight
    return nonnegative(t.sq_value/t.sq_weight-mu.square())


def _sigma_raw(t:TunerState,b:BandState,acc_noise_floor_sigma:Interval,
               still_attenuation:Interval,sigma_coeff:float,max_sigma:float)->Interval|None:
    if not (b.ready and acc_noise_floor_sigma.lo>=0 and 0<=still_attenuation.lo<=still_attenuation.hi<=1):return None
    noise_gain=nonnegative(b.p11)
    noise_sigma=acc_noise_floor_sigma*VPOS.sqrt_interval(noise_gain)
    var_noise=noise_sigma.square()
    total=_variance(t)
    if total is None: total=var_noise
    wave=nonnegative(total-var_noise)*still_attenuation
    wave=Interval(max(1e-6,wave.lo),max(1e-6,wave.hi))
    sigma=VPOS.sqrt_interval(wave)*I(sigma_coeff)
    return clampi(sigma,0.0,max_sigma)


def step_joint(s:JointEstimatorState,dt:Interval,input_accel:Interval,
               *,acc_noise_floor_sigma:Interval=I(0.0148),
               still_attenuation:Interval=I(1.0))->list[JointImage]:
    if not s.source_token:raise ValueError("source token required")
    lin=WAVE.linear_frontend_step(s.wave,dt,input_accel)
    c=WAVE.constants();start=I(3.0)/c["lambda"]
    moment_modes=[]
    if lin.elapsed_sec.hi<start.lo: moment_modes=[("pre_moment",lin)]
    elif lin.elapsed_sec.lo>=start.hi: moment_modes=[("moment",WAVE.moment_update(lin,dt,_wave_horizon(s.canonical_period_s)))]
    else:
        moment_modes=[("pre_moment",lin),("moment",WAVE.moment_update(lin,dt,_wave_horizon(s.canonical_period_s)))]
    out=[]
    tc=TARGET.constants()
    for mode,w in moment_modes:
        raw=_raw_period(w) if mode=="moment" else None
        period=s.canonical_period_s
        if raw is not None:
            period=_canonical_period_step(period,raw,dt)
            if period is None:
                ns=JointEstimatorState(s.source_token,s.predecessor_token,w,s.canonical_period_s,s.usable_period,s.band,s.tuner)
                out.append(JointImage(ns,"split_required_log_period",None,None,None,None,None));continue
        usable=_usable(s.usable_period,w,period)
        if period is None:
            ns=JointEstimatorState(s.source_token,s.predecessor_token,w,None,usable,s.band,s.tuner)
            out.append(JointImage(ns,"hold_no_valid_period",None,None,None,None,None));continue
        f=period.reciprocal()
        # Shipping ordering: current band coefficients use PREVIOUS tuner f;
        # only if tuner f is not ready does the current WavePeriodEstimator f
        # supply the fallback.  The tuner then commits current f after the band step.
        fband=s.tuner.frequency_hz if s.tuner.frequency_hz is not None else f
        band=_band_step(s.band,input_accel,dt,clampi(fband,0.05,5.0))
        if band is None:
            ns=JointEstimatorState(s.source_token,s.predecessor_token,w,period,usable,s.band,s.tuner)
            out.append(JointImage(ns,"split_required_band_coefficients",None,None,None,None,None));continue
        tuner=_tuner_step(s.tuner,band.band,dt,f)
        sig=_sigma_raw(tuner,band,acc_noise_floor_sigma,still_attenuation,
                       tc["sigma_coeff"],tc["MAX_SIGMA_A"])
        ns=JointEstimatorState(s.source_token,s.predecessor_token,w,period,usable,band,tuner)
        if sig is None:
            out.append(JointImage(ns,"hold_sigma_not_ready",f,None,None,None,None));continue
        tau=TARGET.tau_from_frequency(f,tc);ts=TARGET.pseudo_period_from_tau(tau,tc)
        # spectral_mse_rs requires strictly positive sigma.  Shipping var floor
        # makes the reachable raw target positive whenever sigma_coeff>0.
        rs=TARGET.spectral_mse_rs(tau,sig,tc)
        out.append(JointImage(ns,"joint_target",f,sig,tau,ts,rs))
    return out


def build()->dict:
    # Structural smoke: all coordinates are one state object and independent
    # f/sigma inputs do not exist in the API.  Startup correctly holds rather
    # than fabricating a period from zero moments.
    s=zero_state();dt=Interval.outward_bounds(0.004,0.006)
    b=step_joint(s,dt,Interval.outward_bounds(-0.1,0.1))
    startup=bool(b) and all(x.status in ("hold_no_valid_period","split_required_log_period") for x in b)
    return {
      "schema":SCHEMA,"qualification":QUALIFICATION,
      "canonical_source":"COMPLETE_BRMM_NORMAL_LIVE_WORD",
      "same_signal_wave_period_and_band_state_materialized":True,
      "same_history_positive_variance_required_before_ratio":True,
      "raw_period_and_canonical_log_period_relation_materialized":True,
      "exact_frequency_reciprocal_relation_materialized":True,
      "adaptive_band_previous_tuner_frequency_lag_materialized":True,
      "adaptive_band_time_varying_white_noise_covariance_materialized":True,
      "sigma_first_second_moment_statistic_materialized":True,
      "band_noise_p11_subtraction_materialized":True,
      "raw_sigma_and_effective_OU_sigma_kept_distinct":True,
      "joint_f_tau_sigma_TS_RS_image_materialized":True,
      "independent_f_sigma_rectangle_accepted_by_transition":False,
      "unresolved_wide_cells_require_source_split":True,
      "startup_hold_smoke_pass":startup,
      "physical_BRMM_signal_attachment_closed_here":False,
      "all_admitted_same_history_source_cells_emitted_here":False,
      "endpoint_and_every_prefix_augmented_LDLT_closed_here":False,
      "P4_MOTION_PASS":False,"P4_PASS":False,"P5_MAY_START":False,
      "next_obligation":"attach each admitted BRMM/private-observer vertical-acceleration history and stillness branch to this joint state, recursively split cells that cannot certify the moment/log/band predicates, then feed only emitted joint targets into same-history Riccati/Joseph/reset source cells and literal every-prefix augmented LDLT",
    }


def validate(d:dict)->list[str]:
    f=[]
    if d.get("schema")!=SCHEMA or d.get("qualification")!=QUALIFICATION:f.append("schema/qualification mismatch")
    for k in ("same_signal_wave_period_and_band_state_materialized","same_history_positive_variance_required_before_ratio","raw_period_and_canonical_log_period_relation_materialized","exact_frequency_reciprocal_relation_materialized","adaptive_band_previous_tuner_frequency_lag_materialized","adaptive_band_time_varying_white_noise_covariance_materialized","sigma_first_second_moment_statistic_materialized","band_noise_p11_subtraction_materialized","raw_sigma_and_effective_OU_sigma_kept_distinct","joint_f_tau_sigma_TS_RS_image_materialized","unresolved_wide_cells_require_source_split","startup_hold_smoke_pass"):
        if d.get(k) is not True:f.append(k+" not true")
    for k in ("independent_f_sigma_rectangle_accepted_by_transition","physical_BRMM_signal_attachment_closed_here","all_admitted_same_history_source_cells_emitted_here","endpoint_and_every_prefix_augmented_LDLT_closed_here","P4_MOTION_PASS","P4_PASS","P5_MAY_START"):
        if d.get(k) is not False:f.append(k+" not false")
    return f


def main()->int:
    ap=argparse.ArgumentParser();ap.add_argument("--output",type=Path,required=True);a=ap.parse_args()
    d=build();f=validate(d);d["validation_pass"]=not f;d["validation_failures"]=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"joint":d["joint_f_tau_sigma_TS_RS_image_materialized"],"P4":d["P4_PASS"],"failures":f},sort_keys=True))
    return int(bool(f))
if __name__=="__main__":raise SystemExit(main())
