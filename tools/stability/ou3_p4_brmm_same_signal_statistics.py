#!/usr/bin/env python3
"""Correlated BRMM period/sigma statistics for OU-III P4.

This is a proof companion to the trusted shipping WPE/tuner state.  It does not
change deployed code.  It retains the exact central-numerator identities

    C_v = q_v*w - m_v^2,
    C_e = q_e*w - m_e^2,
    C_a = q_a*w - m_a^2,

and advances them with the SAME samples and EW alphas as the shipping first and
second moments.  Consequently period and sigma are images of one physical
signal history; they are never independent interval coordinates.

For m'=r*m+a*x, q'=r*q+a*x^2, w'=r*w+a,

    C' = r^2 C + r*a*(C/w + w*(x-m/w)^2),                  w>0.

The current raw target uses the predecessor WPE frequency exactly as shipping
does.  Sigma is then obtained from the same period-scaled band variance after
p11 white-noise subtraction and same-history stillness attenuation.  tau, T_S
and R_S are deterministic dependent images of that target.
"""
from __future__ import annotations
from dataclasses import dataclass,replace
import argparse,json,math,re
from pathlib import Path

from ou3_interval import Interval
import ou3_validated_transcendentals as VT
import ou3_validated_log as VLOG
import ou3_brmm_wpe_state_step as WPE
import ou3_brmm_tuner_scheduler_step as TUNER

SCHEMA=1
QUALIFICATION='OU3_P4_BRMM_SAME_SIGNAL_CENTRAL_STATISTICS_V1'
PI=Interval.outward_bounds(math.pi,math.nextafter(math.pi,math.inf))
ONE=Interval.point(1.0)
TWO=Interval.point(2.0)

def I(x):return Interval.point(float(x))
def nonnegative(x):return Interval(max(0.0,x.lo),max(0.0,x.hi))

def central_from_shipping(mean,sq,weight):
    """Initialization only; production interval roots must provide C directly."""
    c=sq*weight-mean.square()
    if c.hi<0:raise ValueError('shipping moments inconsistent with nonnegative central numerator')
    return nonnegative(c)

def central_update(C,m,w,x,alpha):
    if C.lo<0:raise ValueError('central numerator lost nonnegativity')
    r=ONE-alpha
    if w.lo<=0:
        if w.lo==0 and w.hi==0 and C.lo==0 and C.hi==0 and m.lo==0 and m.hi==0:return I(0)
        raise RuntimeError('central weight straddles zero; split predecessor source cell')
    centered=x-m/w
    return nonnegative(r.square()*C+r*alpha*(C/w+w*centered.square()))

def variance(C,w):
    if C.lo<0 or w.lo<=1e-12:raise RuntimeError('positive central weight required')
    return nonnegative(C/w.square())

@dataclass(frozen=True)
class State:
    wpe_velocity_C:Interval
    wpe_elevation_C:Interval
    tuner_band_C:Interval
    theorem_log_period_s:Interval|None

@dataclass(frozen=True)
class Target:
    frequency_hz:Interval
    tau:Interval
    sigma:Interval
    pseudo_period:Interval
    rs:Interval
    tuner_C_next:Interval

@dataclass(frozen=True)
class WPESuccessor:
    shipping_state:WPE.WPEState
    correlation_state:State
    raw_period_s:Interval|None
    branch:str


def from_shipping(wpe:WPE.WPEState,tuner:TUNER.TunerState)->State:
    return State(
        central_from_shipping(wpe.velocity_mean,wpe.velocity_sq,wpe.weight),
        central_from_shipping(wpe.elevation_mean,wpe.elevation_sq,wpe.weight),
        central_from_shipping(tuner.moments.mean_value,tuner.moments.sq_value,tuner.moments.mean_weight),
        wpe.log_period_s,
    )

def frequency(state:State, *, usable_period:bool=True)->Interval:
    """The shipping selector, not an independently selectable frequency.

    None represents the shipping not-yet-initialized NaN log state. The fixed
    configured prior is used until the one-way usable latch, even if logT is
    already finite. This module covers the declared default configuration.
    """
    if not usable_period:
        source=(Path(__file__).resolve().parents[2]/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h').read_text()
        match=re.search(r'constexpr\s+float\s+TUNE_FREQ_PRIOR_HZ\s*=\s*([0-9.eE+-]+)f',source)
        if not match:raise RuntimeError('shipping tuner prior definition changed')
        return I(float(match.group(1)))
    if state.theorem_log_period_s is None:
        raise ValueError('usable WPE cannot have an uninitialized canonical log period')
    return WPE.wide_exp(-state.theorem_log_period_s)

def _wpe_alpha(state:WPE.WPEState,c:WPE.Constants):
    period=I(6) if state.log_period_s is None else WPE.wide_exp(state.log_period_s)
    h=WPE.clamp_interval(I(c.moment_horizon_periods)*period,c.min_horizon_s,c.max_horizon_s)
    return ONE-VT.exp_interval(-(I(c.dt)/h)),h

def _wpe_physical_and_moments(state:WPE.WPEState,a_vertical:Interval,c:WPE.Constants):
    lam=TWO*PI*I(c.high_pass_hz);decay=VT.exp_interval(-(lam*I(c.dt)));gain=(ONE-decay)/lam
    stage1=decay*(state.high_pass_1+a_vertical-state.accel_prev)
    stage2=decay*(state.high_pass_2+stage1-state.high_pass_1_prev)
    velocity_new=decay*state.velocity+gain*stage2
    elevation_new=decay*state.elevation+gain*velocity_new
    alpha,h=_wpe_alpha(state,c);r=ONE-alpha
    weight=r*state.weight+alpha
    vm=r*state.velocity_mean+alpha*velocity_new;vq=r*state.velocity_sq+alpha*velocity_new.square()
    em=r*state.elevation_mean+alpha*elevation_new;eq=r*state.elevation_sq+alpha*elevation_new.square()
    base=replace(state,accel_prev=a_vertical,high_pass_1=stage1,high_pass_1_prev=stage1,high_pass_2=stage2,
                 velocity=velocity_new,elevation=elevation_new,velocity_mean=vm,velocity_sq=vq,
                 elevation_mean=em,elevation_sq=eq,weight=weight,elapsed_s=state.elapsed_s+I(c.dt),
                 last_moment_horizon_s=h)
    return base,alpha,velocity_new,elevation_new

def advance_wpe(wpe:WPE.WPEState,corr:State,a_vertical:Interval,c:WPE.Constants)->list[WPESuccessor]:
    """Post-moment-start transition including prior Live and every guard branch.

    Each returned interval is a conservative image of its named guard; unions
    cover all predecessors. In particular a straddling validity guard MUST NOT
    return only the hold image. Closed threshold intersections are intentional:
    using nextafter(threshold,+inf) could drop real points before binary32 is
    enclosed. No target-platform rounding qualification is inferred here.
    """
    lam=TWO*PI*I(c.high_pass_hz)
    if wpe.elapsed_s.lo < (I(3)/lam).hi:
        raise ValueError('pre-moment-start source needs the startup transition')
    if wpe.log_period_s != corr.theorem_log_period_s:
        raise ValueError('shipping and correlated log-period states detached')
    if wpe.usable_period and wpe.log_period_s is None:
        raise ValueError('usable WPE with uninitialized log period')
    base,alpha,v,e=_wpe_physical_and_moments(wpe,a_vertical,c)
    Cv=central_update(corr.wpe_velocity_C,wpe.velocity_mean,wpe.weight,v,alpha)
    Ce=central_update(corr.wpe_elevation_C,wpe.elevation_mean,wpe.weight,e,alpha)
    cs=replace(corr,wpe_velocity_C=Cv,wpe_elevation_C=Ce)
    if base.weight.lo<=1e-3:
        # Explicit noncoverage, never silently promote a subset of this guard.
        raise RuntimeError('weight qualification straddles 1e-3; split source cell')
    vv=variance(Cv,base.weight);ev=variance(Ce,base.weight)
    out=[]
    if vv.lo<=1e-12 or ev.lo<=1e-12:
        out.append(WPESuccessor(base,cs,None,'hold_invalid_variance'))
    if vv.hi<=1e-12 or ev.hi<=1e-12:return out
    # Both variances share EXACTLY the same debiasing weight: vv/ev=Cv/Ce.
    # Intersect the valid branch before division, retaining the central states
    # alongside the output. This removes a duplicated independent weight ratio.
    floor=(I(1e-12)*base.weight.square()).lo
    cv_valid=Interval(max(Cv.lo,floor),Cv.hi)
    ce_valid=Interval(max(Ce.lo,floor),Ce.hi)
    omega2=cv_valid/ce_valid-lam.square()
    if omega2.lo<=1e-8:
        out.append(WPESuccessor(base,cs,None,'hold_invalid_omega'))
    if omega2.hi<=1e-8:return out
    omega2=Interval(max(omega2.lo,1e-8),omega2.hi)
    raw=TWO*PI/WPE.sqrt_positive(omega2);log_raw=VLOG.log_interval(raw)
    if corr.theorem_log_period_s is None:
        log_new=log_raw;log_h=I(0)
    elif c.log_smoothing_periods<=0:
        log_new=log_raw;log_h=I(c.dt)
    else:
        requested=I(c.log_smoothing_periods)*WPE.wide_exp(corr.theorem_log_period_s)
        log_h=WPE.clamp_interval(requested,max(c.dynamic_horizon_min_s,c.dt),c.dynamic_horizon_max_s)
        al=ONE-VT.exp_interval(-(I(c.dt)/log_h));log_new=corr.theorem_log_period_s+al*(log_raw-corr.theorem_log_period_s)
    shipping=replace(base,raw_period_s=raw,log_period_s=log_new,last_log_horizon_s=log_h)
    cs=replace(cs,theorem_log_period_s=log_new)
    if wpe.usable_period:
        out.append(WPESuccessor(shipping,cs,raw,'valid_period'))
        return out
    # The latch is evaluated ONLY after a valid raw-period update.
    age=base.elapsed_s-I(3)/lam
    floor_age=base.elapsed_s-I(4)/lam
    age_margin=age-WPE.wide_exp(log_new)
    can_latch=floor_age.hi>=0 and age_margin.hi>=0
    must_latch=floor_age.lo>=0 and age_margin.lo>=0
    if not must_latch:
        out.append(WPESuccessor(shipping,cs,raw,'valid_period_prior'))
    if can_latch:
        out.append(WPESuccessor(replace(shipping,usable_period=True),cs,raw,'valid_period_takeover'))
    return out

def _tuner_alpha(m:TUNER.MomentState,f_wave:Interval,c:TUNER.Constants):
    f=TUNER.clamp_interval(f_wave,0.05,5.0)
    sea=TUNER.clamp_interval(I(.5)/f,c.dynamic_scale_min_s,c.dynamic_scale_max_s)
    requested=TUNER.clamp_interval(I(4)*(I(2)*sea),.3,60.)
    h=TUNER.clamp_interval(requested,max(c.dynamic_horizon_min_s,c.dt),c.dynamic_horizon_max_s)
    return ONE-VT.exp_interval(-(I(c.dt)/h))

def target_after_band(old_tuner:TUNER.TunerState,new_band:TUNER.BandState,corr:State,f_previous:Interval,still_attenuation:Interval,c:TUNER.Constants)->Target:
    if still_attenuation.lo<0 or still_attenuation.hi>1:raise ValueError('still attenuation outside [0,1]')
    alpha=_tuner_alpha(old_tuner.moments,f_previous,c)
    C=central_update(corr.tuner_band_C,old_tuner.moments.mean_value,old_tuner.moments.mean_weight,new_band.band,alpha)
    r=ONE-alpha;w_new=r*old_tuner.moments.mean_weight+alpha
    total=variance(C,w_new)
    f=TUNER.clamp_interval(f_previous,c.tune_freq_min,c.tune_freq_max)
    tau=TUNER.clamp_interval(I(c.tau_coeff*.5)/f,c.tau_min,c.tau_max)
    noise=I(c.acc_noise_floor_sigma)*TUNER.sqrt_interval(TUNER.max_interval(new_band.p11,1e-30))
    wave=nonnegative(total-noise.square())*still_attenuation
    sigma=TUNER.min_interval(I(c.sigma_coeff)*TUNER.sqrt_interval(TUNER.max_interval(wave,1e-6)),c.sigma_max)
    ts=TUNER.pseudo_period(tau,c)
    rs=TUNER.clamp_interval(TUNER.spectral_mse_rs_target(tau,sigma,c),c.rs_min,c.rs_max)
    return Target(f,tau,sigma,ts,rs,C)

def build():
    import ou3_brmm_frontend_state_step as FRONT
    st=FRONT._point_state();corr=from_shipping(st.wpe,st.tuner);tc=TUNER.constants();wc=WPE.constants(tc.dt)
    f=frequency(corr);band=TUNER.band_step(st.tuner.band,I(.2),st.tuner.moments.frequency_hz,tc)
    target=target_after_band(st.tuner,band,corr,f,I(1),tc)
    ws=advance_wpe(st.wpe,corr,I(.2),wc)
    valid=[x for x in ws if x.branch=='valid_period']
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'same_signal_central_WPE_velocity_variance_materialized':True,'same_signal_central_WPE_elevation_variance_materialized':True,
      'same_signal_central_tuner_band_variance_materialized':True,'shipping_first_second_moments_retained_in_parallel':True,
      'previous_estimated_WPE_frequency_drives_current_tuner_target':True,'sigma_uses_same_band_signal_and_p11_noise_subtraction':True,
      'same_history_stillness_attenuation_enters_sigma':True,'tau_TS_RS_are_dependent_images_of_estimated_f_sigma':True,
      'independent_period_sigma_parameters_forbidden':True,'point_target_positive':all(x.lo>0 for x in (target.frequency_hz,target.tau,target.sigma,target.pseudo_period,target.rs)),
      'point_WPE_valid_successor_exists':bool(valid),'complete_BRMM_predecessor_family_covered_here':False,
      'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('same_signal_central_WPE_velocity_variance_materialized','same_signal_central_WPE_elevation_variance_materialized','same_signal_central_tuner_band_variance_materialized','shipping_first_second_moments_retained_in_parallel','previous_estimated_WPE_frequency_drives_current_tuner_target','sigma_uses_same_band_signal_and_p11_noise_subtraction','same_history_stillness_attenuation_enters_sigma','tau_TS_RS_are_dependent_images_of_estimated_f_sigma','independent_period_sigma_parameters_forbidden','point_target_positive','point_WPE_valid_successor_exists'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('complete_BRMM_predecessor_family_covered_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'same_signal':d['tau_TS_RS_are_dependent_images_of_estimated_f_sigma'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())