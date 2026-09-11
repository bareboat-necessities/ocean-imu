#!/usr/bin/env python3
"""Universal physical joint24 local-event outer using inverse-free gain bounds.

This supersedes the earlier experimental covariance-box local outer. Actual
same-cell Kalman gains are enclosed by ``inverse_free_gain_outer`` from a proved
event-local covariance norm and S>=R; no endpoint covariance is reused locally.
The state domain comes from ALT's conditional regional Normal-Live binding and
does not execute or assume startup reachability.

For each H/A x BIAS family this module emits prediction, physical S,
accelerometer, and magnetometer homogeneous outers. The two IMU-core patterns
are pred->acc and pred->S->acc. Magnetometer remains a separate event outer for
arbitrary-finite star composition.
"""
from __future__ import annotations
import json,math
from pathlib import Path
from ou3_interval import Interval,matrix_identity,matrix_mul
import ou3_p4_complete_brmm_differential_prediction as PRED
import ou3_p4_complete_brmm_differential_events as EVENTS
import ou3_brmm_complete_window_execution_kernel as KERNEL
import ou3_brmm_tuner_scheduler_step as TUNER
from tools.stability.ou3_alt_contraction import joint24_events as J24
from tools.stability.ou3_alt_contraction import bias_families as BIAS
from tools.stability.ou3_alt_contraction import physical_numeric_bounds as PHYS
from tools.stability.ou3_alt_contraction import inverse_free_gain_outer as GAIN
from tools.stability.ou3_alt_contraction import regional_error_domain as REGION

REPO=Path(__file__).resolve().parents[3]
DOMAIN=REPO/'tools/stability/ou3_proof_operating_domain.json'
QUALIFICATION='OU3_ALT_INVERSE_FREE_UNIVERSAL_LOCAL_PHYSICAL_OUTER_V2'
def I(x):return Interval.point(float(x))
def B(x):x=float(x);return Interval.outward_bounds(-x,x)
def _zero(r,c):return [[I(0) for _ in range(c)] for _ in range(r)]
def _finite(A):return len(A)==24 and all(len(r)==24 for r in A) and all(math.isfinite(x.lo) and math.isfinite(x.hi) for r in A for x in r)
def _compose(*maps):
    A=matrix_identity(24)
    for J in maps:A=matrix_mul(J,A)
    return A
def _state_box(mode,radii):
    out=[]
    for name in ('attitude_cayley_norm','gyro_bias_norm_rad_s','velocity_norm_mps','position_norm_m','integral_displacement_norm_m_s','latent_acceleration_norm_mps2'):out.extend([B(radii[name])]*3)
    if mode=='A':out.extend([B(radii['accelerometer_bias_error_norm_mps2'])]*3)
    return out
def _bias_box(contract,projection_radius):
    beta=[B(contract.true_bias_component_bound) for _ in range(3)]
    held=[B(contract.true_bias_component_bound+projection_radius) for _ in range(3)]
    return held,beta
def _a_prediction_lift(J21,phi_hat,phi_true):
    A=_zero(24,24)
    for i in range(21):
        for j in range(21):A[i][j]=J21[i][j]
    for i in range(3):A[18+i][21+i]=phi_true-phi_hat;A[21+i][21+i]=phi_true
    return A
def _prediction(mode,state,omega,kc,tau,contract):
    p=PRED.prediction_event(mode,state,omega,kc.h,tau,tau_ba=kc.accel_bias_tau_s if mode=='A' else None)
    return J24.h18_prediction_lift(p['J_state'],contract.phi_true)[0] if mode=='H' else _a_prediction_lift(p['J_state'],p['phi_ba'],contract.phi_true)
def _h_event(kind,state,K,held,beta,*,f=None,Rhat=None,m=None,Sphys=None):
    extra=list(Sphys or ());total=24+len(extra);u=EVENTS.AD.independent_vector(list(state)+list(held)+list(beta)+extra,n=total,offset=0);z=u[:18];eba=u[18:21]
    if kind=='accelerometer':y=J24._h18_accel_residual(z,eba,f,Rhat)
    elif kind=='magnetometer':y=EVENTS.residual_magnetometer(z,m)
    elif kind=='S_zero':y=[z[EVENTS.OFF_S+i]-u[24+i] for i in range(3)]
    else:raise ValueError(kind)
    out=list(EVENTS._apply_physical_correction(z,K,y))+list(u[18:24]);return [r[:24] for r in EVENTS.AD.jacobian(out)]
def _a_event(kind,state,K,beta,radius,*,f=None,Rhat=None,m=None,Sphys=None):
    extra=list(Sphys or ());total=24+len(extra);u=EVENTS.AD.independent_vector(list(state)+list(beta)+extra,n=total,offset=0);z=u[:21];b=u[21:24]
    if kind=='accelerometer':y=EVENTS.residual_accelerometer(z,f,Rhat)
    elif kind=='magnetometer':y=EVENTS.residual_magnetometer(z,m)
    elif kind=='S_zero':y=[z[EVENTS.OFF_S+i]-u[24+i] for i in range(3)]
    else:raise ValueError(kind)
    out21=EVENTS._apply_physical_correction(z,K,y);J21,_,_=J24._project_joint(out21,b,total,radius);J=[list(r) for r in J21]
    for i in range(3):row=[I(0) for _ in range(total)];row[21+i]=I(1);J.append(row)
    return [r[:24] for r in J]
def build():
    gain=GAIN.build();gf=GAIN.validate(gain);phys=PHYS.build();pf=PHYS.validate(phys);region=REGION.build();rf=REGION.validate(region)
    if gf or pf or rf:raise RuntimeError(f'local outer prerequisites failed gain={gf} phys={pf} region={rf}')
    d=json.loads(DOMAIN.read_text())['normal_live'];kc=KERNEL._process_constants(DOMAIN);tc=TUNER.constants();tau=Interval.outward_bounds(tc.tau_min,tc.tau_max);omega=[B(math.radians(float(phys['body_rate_norm_upper_deg_s'])))]*3
    f=[B(float(phys['specific_force_norm_upper_mps2']))]*3;Rhat=[[Interval(-1.0,1.0) for _ in range(3)] for _ in range(3)];m=[B(float(d['magnetic_vector_norm_upper_uT']))]*3;Sphys=[B(float(phys['centered_S_norm_upper_m_s']))]*3;radii=region['coordinate_radii'];radius=float(d['active_accelerometer_bias_projection_limit_mps2']);reports={}
    for mode in ('H','A'):
        state=_state_box(mode,radii);reports[mode]={}
        for contract in BIAS.contracts():
            held,beta=_bias_box(contract,radius);Krow=gain['reports'][mode];row={'failures':[],'maps':{}}
            def attempt(name,fn):
                try:A=fn();row['maps'][name]=A;row[name+'_finite']=_finite(A);row['failures']+=[] if row[name+'_finite'] else [name+': nonfinite']
                except Exception as exc:row['maps'][name]=None;row[name+'_finite']=False;row['failures'].append(name+': '+type(exc).__name__+': '+str(exc))
            attempt('prediction',lambda:_prediction(mode,state,omega,kc,tau,contract))
            if mode=='H':
                attempt('S_zero',lambda:_h_event('S_zero',state,Krow['S_zero']['K'],held,beta,Sphys=Sphys));attempt('accelerometer',lambda:_h_event('accelerometer',state,Krow['accelerometer']['K'],held,beta,f=f,Rhat=Rhat));attempt('magnetometer',lambda:_h_event('magnetometer',state,Krow['magnetometer']['K'],held,beta,m=m))
            else:
                attempt('S_zero',lambda:_a_event('S_zero',state,Krow['S_zero']['K'],beta,radius,Sphys=Sphys));attempt('accelerometer',lambda:_a_event('accelerometer',state,Krow['accelerometer']['K'],beta,radius,f=f,Rhat=Rhat));attempt('magnetometer',lambda:_a_event('magnetometer',state,Krow['magnetometer']['K'],beta,radius,m=m))
            core=all(row.get(k+'_finite') for k in ('prediction','S_zero','accelerometer'));mag=bool(row.get('magnetometer_finite'));row['imu_core_family']=[]
            if core:
                P=row['maps']['prediction'];S=row['maps']['S_zero'];A=row['maps']['accelerometer'];Id=matrix_identity(24);row['imu_core_family']=[_compose(P,Id,A),_compose(P,S,A)]
            row['imu_core_family_finite']=core and all(_finite(x) for x in row['imu_core_family']);row['magnetometer_event_outer_finite']=mag;reports[mode][contract.name]=row
    core_closed=all(r['imu_core_family_finite'] for rows in reports.values() for r in rows.values());mag_closed=all(r['magnetometer_event_outer_finite'] for rows in reports.values() for r in rows.values())
    return {'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','inverse_free_gain_outer_consumed':True,'regional_error_domain_consumed_without_startup':True,'fresh_startup_reachability_consumed':False,'endpoint_Pbar_used_for_local_gain':False,'independent_covariance_box_used':False,'physical_D_S_bound_m_s':phys['centered_S_norm_upper_m_s'],'reports':reports,'all_IMU_core_local_outers_finite':core_closed,'all_single_magnetometer_local_outers_finite':mag_closed,'local_physical_outer_materialized':core_closed and mag_closed,'magnetometer_count_assumption_used':False,'ALT_LIVE_PASS':False,'next_obligation':'induct IMU-core family over 600 steps; common-M certificate must also make every separate magnetometer event outer nonexpansive so arbitrary finite async tuples close'}
def summary(d):return {m:{f:{'failures':r['failures'],'imu_core':r['imu_core_family_finite'],'mag':r['magnetometer_event_outer_finite']} for f,r in rows.items()} for m,rows in d['reports'].items()}
def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    for k in ('inverse_free_gain_outer_consumed','regional_error_domain_consumed_without_startup'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('fresh_startup_reachability_consumed','endpoint_Pbar_used_for_local_gain','independent_covariance_box_used','magnetometer_count_assumption_used','ALT_LIVE_PASS'):
        if d.get(k) is not False:f.append(k+' not false')
    if float(d.get('physical_D_S_bound_m_s',0))!=1100.0:f.append('physical D_S drifted')
    return f
