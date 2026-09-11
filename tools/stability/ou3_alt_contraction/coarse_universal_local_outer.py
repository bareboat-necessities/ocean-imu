#!/usr/bin/env python3
"""First numerical universal local-map outer attempt for ALT.

This deliberately conservative construction is NOT a source-history generator.
The exact same-cell relation has already been defined by exact_endpoint_relation.
Here we project that relation onto separately certified coordinate enclosures and
evaluate the exact shipping differential formulas over their Cartesian product.
Because R subset product(proj_i R), the resulting interval Jacobians are valid
outer enclosures of every actual same-cell Jacobian.  The product may contain
unrealizable combinations; that only makes the certificate conservative.

The purpose of this module is to obtain the first concrete numerical local map
family that can be fed to endpoint_family_induction.  If the broad projection
product causes interval inversion/chart failure or destroys contraction, the
failure classifies this representation as too coarse and does NOT falsify the
physical theorem.
"""
from __future__ import annotations

import argparse,json,math
from pathlib import Path

from ou3_interval import Interval,matrix_identity,matrix_mul
import ou3_brmm_riccati_tube_factored as TUBE
import ou3_brmm_complete_window_execution_kernel as KERNEL
import ou3_brmm_tuner_scheduler_step as TUNER
import ou3_brmm_full_normal_live_word as WORD
import ou3_p4_complete_brmm_differential_prediction as PRED
import ou3_p4_complete_brmm_differential_events as EVENTS
import ou3_p4_hard_entry_set as ENTRY
from tools.stability.ou3_alt_contraction import bias_families as BIAS
from tools.stability.ou3_alt_contraction import joint24_events as J24
from tools.stability.ou3_alt_contraction import physical_numeric_bounds as PHYS
from tools.stability.ou3_alt_contraction import exact_endpoint_relation as EXACT

REPO=Path(__file__).resolve().parents[3]
DOMAIN=REPO/'tools/stability/ou3_proof_operating_domain.json'
QUALIFICATION='OU3_ALT_COARSE_UNIVERSAL_LOCAL_MAP_OUTER_V1'

def I(x):return Interval.point(float(x))
def B(r):r=float(r);return Interval.outward_bounds(-r,r)
def _zero(r,c):return [[I(0) for _ in range(c)] for _ in range(r)]
def _shape(A):return len(A),len(A[0]) if A else 0

def psd_covariance_box_from_diagonal_upper(diag):
    """Outer box for PSD P with P_ii<=d_i; never a reachable-P generator."""
    d=[float(x) for x in diag]
    if not d or any(not math.isfinite(x) or x<=0 for x in d):raise ValueError('positive finite covariance diagonal ceiling required')
    n=len(d);out=[]
    for i in range(n):
        row=[]
        for j in range(n):
            if i==j:row.append(Interval.outward_bounds(0.0,d[i]))
            else:
                q=math.nextafter(math.sqrt(d[i]*d[j]),math.inf);row.append(Interval(-q,q))
        out.append(row)
    return out

def error_state_box(mode,radii):
    groups=(('attitude_cayley_norm',3),('gyro_bias_norm_rad_s',3),('velocity_norm_mps',3),('position_norm_m',3),('integral_displacement_norm_m_s',3),('latent_acceleration_norm_mps2',3))
    out=[]
    for name,n in groups:out.extend([B(radii[name]) for _ in range(n)])
    if mode=='A':out.extend([B(radii['accelerometer_bias_error_norm_mps2']) for _ in range(3)])
    if len(out)!=(18 if mode=='H' else 21):raise RuntimeError('state box dimension drift')
    return out

def _bias_domains(contract,projection_radius):
    beta=[B(contract.true_bias_component_bound) for _ in range(3)]
    ecomp=float(contract.true_bias_component_bound)+float(projection_radius)
    return [B(ecomp) for _ in range(3)],beta

def _a_prediction_lift(J21,phi_hat,phi_true):
    A=_zero(24,24)
    for i in range(21):
        for j in range(21):A[i][j]=J21[i][j]
    for i in range(3):A[18+i][21+i]=phi_true-phi_hat;A[21+i][21+i]=phi_true
    return A

def prediction_outer(mode,state,omega,h,tau,tau_ba,contract):
    p=PRED.prediction_event(mode,state,omega,h,tau,tau_ba=tau_ba if mode=='A' else None)
    return J24.h18_prediction_lift(p['J_state'],contract.phi_true)[0] if mode=='H' else _a_prediction_lift(p['J_state'],p['phi_ba'],contract.phi_true)

def _h_joseph_outer(kind,state,P,R,held,beta,*,f=None,Rhat=None,m=None,Sphys=None):
    extra=list(Sphys or ());total=24+len(extra);u=EVENTS.AD.independent_vector(list(state)+list(held)+list(beta)+extra,n=total,offset=0);z=u[:18];eba=u[18:21]
    H=EVENTS._event_H('H',kind,f_hat=f,R_hat=Rhat,m_body=m);K,_=EVENTS.source_joseph_gain(P,H,R)
    if kind=='accelerometer':y=J24._h18_accel_residual(z,eba,f,Rhat)
    elif kind=='magnetometer':y=EVENTS.residual_magnetometer(z,m)
    elif kind=='S_zero':y=[z[EVENTS.OFF_S+i]-u[24+i] for i in range(3)]
    else:raise ValueError(kind)
    out=list(EVENTS._apply_physical_correction(z,K,y))+list(u[18:24]);return [row[:24] for row in EVENTS.AD.jacobian(out)]
def _a_joseph_outer(kind,state,P,R,beta,radius,*,f=None,Rhat=None,m=None,Sphys=None):
    extra=list(Sphys or ());total=24+len(extra);u=EVENTS.AD.independent_vector(list(state)+list(beta)+extra,n=total,offset=0);z=u[:21];b=u[21:24]
    H=EVENTS._event_H('A',kind,f_hat=f,R_hat=Rhat,m_body=m);K,_=EVENTS.source_joseph_gain(P,H,R)
    if kind=='accelerometer':y=EVENTS.residual_accelerometer(z,f,Rhat)
    elif kind=='magnetometer':y=EVENTS.residual_magnetometer(z,m)
    elif kind=='S_zero':y=[z[EVENTS.OFF_S+i]-u[24+i] for i in range(3)]
    else:raise ValueError(kind)
    out21=EVENTS._apply_physical_correction(z,K,y);J21,_,_=J24._project_joint(out21,b,total,float(radius));J=[list(row) for row in J21]
    for i in range(3):row=[I(0) for _ in range(total)];row[21+i]=I(1);J.append(row)
    return [row[:24] for row in J]
def _compose(*maps):
    A=matrix_identity(24)
    for J in maps:A=matrix_mul(J,A)
    return A
def _finite_matrix(A):return _shape(A)==(24,24) and all(math.isfinite(x.lo) and math.isfinite(x.hi) for row in A for x in row)

def build():
    exact=EXACT.build();xf=EXACT.validate(exact);phys=PHYS.build();pf=PHYS.validate(phys);entry=ENTRY.build();ef=ENTRY.validate(entry);tube=TUBE.build();tf=TUBE.validate_covariance_ceiling(tube)
    if xf or pf or ef or tf:raise RuntimeError(f'coarse outer prerequisites failed exact={xf} phys={pf} entry={ef} tube={tf}')
    domain=json.loads(DOMAIN.read_text())['normal_live'];kc=KERNEL._process_constants(DOMAIN);tc=TUNER.constants();tau=Interval.outward_bounds(tc.tau_min,tc.tau_max)
    omega_cap=math.radians(float(phys['body_rate_norm_upper_deg_s']));omega=[B(omega_cap) for _ in range(3)];fcap=float(phys['specific_force_norm_upper_mps2']);f=[B(fcap) for _ in range(3)];Rhat=[[Interval(-1.0,1.0) for _ in range(3)] for _ in range(3)];mcap=float(domain['magnetic_vector_norm_upper_uT']);m=[B(mcap) for _ in range(3)];Sphys=[B(float(phys['centered_S_norm_upper_m_s'])) for _ in range(3)]
    active=TUNER.ActiveSchedule(tau,Interval.outward_bounds(max(tc.acc_noise_floor_sigma,1e-6),tc.sigma_max),Interval.outward_bounds(tc.rs_min,tc.rs_max),Interval.outward_bounds(tc.pseudo_min_s,tc.pseudo_max_s));RS=WORD.R_S_zero(TUNER.active_rs_std_xyz(active,tc));radii=entry['coordinate_radii'];radius=float(domain['active_accelerometer_bias_projection_limit_mps2']);reports={}
    for mode in ('H','A'):
        state=error_state_box(mode,radii);P=psd_covariance_box_from_diagonal_upper(tube['modes'][mode]['Pbar_diagonal_variance_upper']);mode_rows={}
        for contract in BIAS.contracts():
            held,beta=_bias_domains(contract,radius);row={'bias_family':contract.name,'maps':{},'failures':[]}
            def attempt(name,fn):
                try:
                    A=fn();row['maps'][name]=A;row[name+'_finite']=_finite_matrix(A)
                    if not row[name+'_finite']:row['failures'].append(name+': nonfinite interval map')
                except Exception as exc:row['maps'][name]=None;row[name+'_finite']=False;row['failures'].append(name+': '+type(exc).__name__+': '+str(exc))
            attempt('prediction',lambda:prediction_outer(mode,state,omega,kc.h,tau,kc.accel_bias_tau_s,contract))
            if mode=='H':
                attempt('S_zero',lambda:_h_joseph_outer('S_zero',state,P,RS,held,beta,Sphys=Sphys));attempt('accelerometer',lambda:_h_joseph_outer('accelerometer',state,P,kc.Racc,held,beta,f=f,Rhat=Rhat));attempt('magnetometer',lambda:_h_joseph_outer('magnetometer',state,P,kc.Rmag,held,beta,m=m))
            else:
                attempt('S_zero',lambda:_a_joseph_outer('S_zero',state,P,RS,beta,radius,Sphys=Sphys));attempt('accelerometer',lambda:_a_joseph_outer('accelerometer',state,P,kc.Racc,beta,radius,f=f,Rhat=Rhat));attempt('magnetometer',lambda:_a_joseph_outer('magnetometer',state,P,kc.Rmag,beta,radius,m=m))
            if all(row.get(k+'_finite') for k in ('prediction','S_zero','accelerometer','magnetometer')):
                Pm=row['maps']['prediction'];Sm=row['maps']['S_zero'];Am=row['maps']['accelerometer'];Mm=row['maps']['magnetometer'];Id=matrix_identity(24);row['sample_family']=[_compose(Pm,s,Am,mx) for s in (Id,Sm) for mx in (Id,Mm)];row['sample_family_finite']=all(_finite_matrix(x) for x in row['sample_family'])
            else:row['sample_family']=[];row['sample_family_finite']=False
            mode_rows[contract.name]=row
        reports[mode]=mode_rows
    all_local=all(r['sample_family_finite'] for modes in reports.values() for r in modes.values())
    return {'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','exact_endpoint_relation_consumed':exact['exact_universal_endpoint_map_relation_defined'],'outer_projection_product_is_history_generator':False,'same_cell_formula_evaluated_after_certified_coordinate_projection':True,'PSD_covariance_cross_bounds_derived_from_rigorous_diagonal_ceiling':True,'physical_centered_S_bound_m_s':phys['centered_S_norm_upper_m_s'],'legacy_300m_s_used_for_physical_S':False,'replay_or_finite_source_sample_used':False,'reports':reports,'all_local_event_outer_maps_finite':all_local,'coarse_universal_sample_family_materialized':all_local,'600_step_endpoint_outer_materialized_here':False,'common_storage_certified_here':False,'ALT_LIVE_PASS':False,'next_obligation':('feed the finite per-sample outer families into partitioned 600-step endpoint induction, then try one common M/rho with outward projected LDLT' if all_local else 'split the certified coordinate projection responsible for each failed local interval map; do not shrink the theorem source or replace it by replay')}
def summary(d):return {mode:{fam:{'failures':r['failures'],'sample_family_finite':r['sample_family_finite']} for fam,r in rows.items()} for mode,rows in d['reports'].items()}
def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION or d.get('canonical_source')!='COMPLETE_BRMM_NORMAL_LIVE_WORD':f.append('qualification/source mismatch')
    for k in ('exact_endpoint_relation_consumed','same_cell_formula_evaluated_after_certified_coordinate_projection','PSD_covariance_cross_bounds_derived_from_rigorous_diagonal_ceiling'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('outer_projection_product_is_history_generator','legacy_300m_s_used_for_physical_S','replay_or_finite_source_sample_used','600_step_endpoint_outer_materialized_here','common_storage_certified_here','ALT_LIVE_PASS'):
        if d.get(k) is not False:f.append(k+' not false')
    if float(d.get('physical_centered_S_bound_m_s',0))!=1100.0:f.append('qualified physical S bound changed')
    # Finiteness is an experiment result, not a validation prerequisite: this
    # module is allowed to fail closed and report which projection needs split.
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);out={k:v for k,v in d.items() if k!='reports'};out['local_summary']=summary(d);out['validation_pass']=not f;out['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps({'all_local':d['all_local_event_outer_maps_finite'],'summary':summary(d),'validation_failures':f},indent=2,sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
