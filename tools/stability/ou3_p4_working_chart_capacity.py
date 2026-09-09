#!/usr/bin/env python3
"""Universal P4 working-chart capacity beyond the 45-degree entry set.

The deterministic hard entry set is an entry hypothesis, not an invariant box.
For a bootstrap/prefix proof we need a larger finite-angle chart in which the
nonlinear vector/reset sectors remain valid.  This producer repeats the same
source-uniform H18 prior-free completion at a fixed 0.90 rad outer attitude
angle (~51.57 deg), then inherits the existing A21 first-active bias completion.

Nothing about the declared 45-degree entry set is changed.  A pass proves only
that the *differential backbone* remains available throughout the wider chart;
it does not prove that every prefix remains inside it.  That is a separate
bootstrap/working-tube obligation.
"""
from __future__ import annotations
import argparse,copy,json,math
from pathlib import Path

import ou3_p4_cayley_sector_certificate as CAYLEY
import ou3_brmm_h18_information_composition as HINFO
import ou3_brmm_h18_prior_free_completion as HPF
import ou3_brmm_dynamic_source_certificate as DYNAMIC
import ou3_full_process_ucc as PROCESS
import ou3_brmm_full_word_event_algebra as EVENT
import ou3_brmm_live_covariance_seed as LIVE
import ou3_brmm_a21_prior_free_completion as ABASE
import ou3_brmm_windowed_vector_pe as PE

QUALIFICATION='OU3_P4_UNIVERSAL_WORKING_CHART_CAPACITY_V1'
WORKING_ANGLE_RAD=0.90
DELTA=1e-18

def down(x):return math.nextafter(float(x),-math.inf)
def up(x):return math.nextafter(float(x),math.inf)

def _h18(angle):
    cay=CAYLEY.build(outer_angle_rad=angle);base=HINFO.build();dyn=DYNAMIC.build();proc=PROCESS.build();event=EVENT.build()
    bad={'cayley':CAYLEY.validate(cay),'hinfo':HINFO.validate(base),'dynamic':DYNAMIC.validate(dyn),'process':PROCESS.validate(proc),'event':EVENT.validate(event)};bad={k:v for k,v in bad.items() if v}
    if bad:raise RuntimeError('working-chart H18 prerequisites failed: '+repr(bad))
    k=float(cay['exact_vector_information_retention_factor_lower']);alpha0=float(base['eta6_information_lower']);alpha=down(k*alpha0)
    comp=base['triangular_information_composition'];daw=float(comp['aw_direction_information_lower']);cross=float(comp['C_aw_spectral_norm_squared_upper']);trace=up(alpha+daw+cross);coupled=down(alpha*daw/trace);nonaw=float(comp['non_aw_translation_lambda_min_lower']);D=down(min(coupled,nonaw))
    if min(alpha,coupled,D)<=0:raise RuntimeError('working-chart information lost positivity')
    h=copy.deepcopy(base);h['eta6_information_lower']=alpha;h['H18_information_useful_gate_pass']=D>=DELTA
    h['triangular_information_composition']['A_transpose_A_lower']=alpha;h['triangular_information_composition']['coupled_eta6_aw_scalar_2x2_trace_upper']=trace;h['triangular_information_composition']['coupled_eta6_aw_scalar_2x2_determinant_lower']=down(alpha*daw);h['triangular_information_composition']['coupled_eta6_aw_lambda_min_lower']=coupled;h['triangular_information_composition']['D_H18_lambda_min_lower']=D
    pbar=HPF._same_word_covariance_upper(Path(HPF.DEFAULT_DOMAIN).resolve(),dyn,proc,h);fnorm=HPF._prediction_norm_sq_upper(proc);penalty=up((DELTA**2/4.0)*fnorm*float(pbar['Pbar_trace_upper']))
    rows=[];fail=[];worst=math.inf
    for x in HPF._x_cover(dyn):
        ok,row=HPF._full_H18_cell(x,process=proc,dynamic=dyn,penalty_physical=penalty);rows.append(row)
        if ok:worst=min(worst,float(row['pivot_lower']))
        else:fail.append(row)
    preserve=event['full_matrix_margin_preservation'];suffix=all(bool(preserve[n]) for n in ('covers_prediction','covers_every_due_S_update','covers_every_Normal_Live_accelerometer_update','covers_asynchronous_magnetometer_update','covers_immediate_left_error_reset','covers_aw_covariance_floor','covers_not_due_or_rejected_identity_branches'))
    closed=bool(rows) and not fail and math.isfinite(worst) and worst>0 and suffix
    return {'outer_angle_rad':angle,'outer_angle_deg':angle*180/math.pi,'cayley_radius_upper':cay['cayley_radius_upper'],'entry_covered':cay['declared_filter_entrance_covered'],'vector_information_retention_lower':k,'eta6_information_lower':alpha,'D_H18_information_lower':D,'worst_LDLT_pivot_lower':worst if closed else None,'x_cells':len(rows),'failures':fail,'suffix_event_algebra':suffix,'closed':closed}

def build():
    h=_h18(WORKING_ANGLE_RAD);proc=PROCESS.build();live=LIVE.build();pe=PE.build();bad={'process':PROCESS.validate(proc),'live':LIVE.validate(live),'pe':PE.validate(pe)};bad={k:v for k,v in bad.items() if v}
    if bad:raise RuntimeError('working-chart A21 prerequisites failed: '+repr(bad))
    held=live['held_ba'];p0=float(held['seed_variance']);q=float(proc['active_accelerometer_bias']['Q_accel_bias_lambda_min_lower']);margin=down((1-DELTA)*q-DELTA*p0);parity=ABASE._prediction_source_parity();pf=[k for k,v in parity.items() if not v];route=pe['A_mode_bias_route'];phi=float(route['homogeneous_bias_contraction_upper_over_word']);biasgap=down(1-up(phi*phi));a21=bool(h['closed'] and margin>0 and not pf)
    entry_deg=45.0;headroom=h['outer_angle_deg']-entry_deg
    return {'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','declared_entry_attitude_deg':entry_deg,'working_chart_outer_angle_rad':WORKING_ANGLE_RAD,'working_chart_outer_angle_deg':h['outer_angle_deg'],'attitude_headroom_deg':headroom,'declared_entry_set_changed':False,'declared_entry_set_shrunk':False,'working_tube_is_larger_than_entry_set':headroom>0,'H18':h,'A21':{'first_active_ba_margin_lower':margin,'finite_bias_homogeneous_energy_gap_lower':biasgap,'prediction_source_parity_failures':pf,'closed':a21},'universal_wider_working_chart_differential_backbone_closed':bool(h['closed'] and a21),'every_prefix_membership_in_working_chart_closed_here':False,'P4_MOTION_PASS':False,'P4_PASS':False}
def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    for k in ('working_tube_is_larger_than_entry_set','universal_wider_working_chart_differential_backbone_closed'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('declared_entry_set_changed','declared_entry_set_shrunk','every_prefix_membership_in_working_chart_closed_here','P4_MOTION_PASS','P4_PASS'):
        if d.get(k) is not False:f.append(k+' not false')
    h=d['H18'];a=d['A21']
    if h.get('entry_covered') is not True or h.get('closed') is not True:f.append('H18 working chart not closed')
    if a.get('closed') is not True:f.append('A21 working chart not closed')
    if not float(d.get('attitude_headroom_deg',0))>0:f.append('no chart headroom')
    if not 0.64<float(h.get('vector_information_retention_lower',0))<1:f.append('working-chart information retention invalid')
    if h.get('failures'):f.append('working-chart x-cell failure')
    return f

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'angle_deg':d['working_chart_outer_angle_deg'],'headroom_deg':d['attitude_headroom_deg'],'k':d['H18']['vector_information_retention_lower'],'pivot':d['H18']['worst_LDLT_pivot_lower'],'A21_ba':d['A21']['first_active_ba_margin_lower'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
