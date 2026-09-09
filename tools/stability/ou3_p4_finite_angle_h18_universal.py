#!/usr/bin/env python3
"""Universal H18 finite-angle backbone for P4 on complete BRMM.

P3 already proves its H18 prior-free full-matrix condition for every admitted
COMPLETE_BRMM_NORMAL_LIVE_WORD without enumerating a finite source family.  This
producer repeats that quantitative chain after replacing the tangent vector-PE
information by the source-independent finite-Cayley differential information
retained on the FULL declared 45-degree entry sector.

For a finite Cayley state c, vector residual differential is

    J_v(c) = -[R v]x A(c),
    sigma_min(A(c)) >= a_min.

Thus every tangent vector-information lower is reduced by at most a_min^2.
The existing H18 information composition already keeps the accelerometer a_w
cross as a matrix cross term and the four-S translation information separately;
we recompute that same 2x2 completion with the reduced eta6 information rather
than charging a nonlinear eta radius.

This is a differential/mean-value backbone.  Exact chord IQCs and reset graph
sectors are still required to turn it into the finite-map endpoint certificate.
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

QUALIFICATION='OU3_P4_UNIVERSAL_H18_FINITE_ANGLE_BACKBONE_V1'


def build():
    # 0.80 rad is the existing outward sector that strictly contains the 45deg entry.
    cay=CAYLEY.build(outer_angle_rad=0.80); base=HINFO.build(); dyn=DYNAMIC.build(); proc=PROCESS.build(); event=EVENT.build()
    bad={'cayley':CAYLEY.validate(cay),'hinfo':HINFO.validate(base),'dynamic':DYNAMIC.validate(dyn),'process':PROCESS.validate(proc),'event':EVENT.validate(event)}
    bad={k:v for k,v in bad.items() if v}
    if bad: raise RuntimeError('finite-angle H18 prerequisites failed: '+repr(bad))

    k=float(cay['exact_vector_information_retention_factor_lower'])
    alpha0=float(base['eta6_information_lower'])
    alpha=math.nextafter(k*alpha0,-math.inf)
    comp=base['triangular_information_composition']
    d_aw=float(comp['aw_direction_information_lower'])
    cross=float(comp['C_aw_spectral_norm_squared_upper'])
    trace=math.nextafter(alpha+d_aw+cross,math.inf)
    coupled=math.nextafter(alpha*d_aw/trace,-math.inf)
    non_aw=float(comp['non_aw_translation_lambda_min_lower'])
    D=math.nextafter(min(coupled,non_aw),-math.inf)
    if not (alpha>0 and coupled>0 and D>0): raise RuntimeError('finite-angle information lost positivity')

    h=copy.deepcopy(base)
    h['eta6_information_lower']=alpha
    h['H18_information_useful_gate_pass']=D>=1e-18
    h['triangular_information_composition']['A_transpose_A_lower']=alpha
    h['triangular_information_composition']['coupled_eta6_aw_scalar_2x2_trace_upper']=trace
    h['triangular_information_composition']['coupled_eta6_aw_scalar_2x2_determinant_lower']=math.nextafter(alpha*d_aw,-math.inf)
    h['triangular_information_composition']['coupled_eta6_aw_lambda_min_lower']=coupled
    h['triangular_information_composition']['D_H18_lambda_min_lower']=D

    pbar=HPF._same_word_covariance_upper(Path(HPF.DEFAULT_DOMAIN).resolve(),dyn,proc,h)
    fnorm=HPF._prediction_norm_sq_upper(proc)
    penalty=math.nextafter((HPF.USEFUL_GATE**2/4.0)*fnorm*float(pbar['Pbar_trace_upper']),math.inf)
    rows=[];fail=[];worst=math.inf
    for x in HPF._x_cover(dyn):
        ok,row=HPF._full_H18_cell(x,process=proc,dynamic=dyn,penalty_physical=penalty)
        rows.append(row)
        if ok: worst=min(worst,float(row['pivot_lower']))
        else: fail.append(row)
    closed=bool(rows) and not fail and math.isfinite(worst) and worst>0

    preserve=event['full_matrix_margin_preservation']
    suffix=all(bool(preserve[k]) for k in ('covers_prediction','covers_every_due_S_update','covers_every_Normal_Live_accelerometer_update','covers_asynchronous_magnetometer_update','covers_immediate_left_error_reset','covers_aw_covariance_floor','covers_not_due_or_rejected_identity_branches'))
    return {
      'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'full_declared_45deg_entry_covered':cay['declared_filter_entrance_covered'],
      'finite_angle_vector_information_retention_lower':k,
      'tangent_eta6_information_lower':alpha0,'finite_angle_eta6_information_lower':alpha,
      'finite_angle_H18_information_lower':D,'finite_angle_information_above_P3_delta':D>=1e-18,
      'same_word_Pbar_trace_upper':float(pbar['Pbar_trace_upper']),
      'delta_squared_completion_penalty':penalty,'x_cells_certified':len(rows),'x_cell_failures':fail,
      'worst_full_H18_LDLT_pivot_lower':worst if closed else None,
      'universal_H18_finite_angle_prior_free_LDLT_closed':closed,
      'linear_suffix_event_algebra_available':suffix,
      'exact_chord_finite_reset_graph_still_required':True,
      'finite_map_endpoint_closed_here':False,'every_prefix_closed_here':False,
      'source_enumeration_used':False,'trajectory_replay_used':False,'domain_shrunk':False,'P3_delta':1e-18,
      'P4_MOTION_PASS':False,'P4_PASS':False,
    }

def validate(d):
    f=[]
    for k in ('full_declared_45deg_entry_covered','finite_angle_information_above_P3_delta','universal_H18_finite_angle_prior_free_LDLT_closed','linear_suffix_event_algebra_available','exact_chord_finite_reset_graph_still_required'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('finite_map_endpoint_closed_here','every_prefix_closed_here','source_enumeration_used','trajectory_replay_used','domain_shrunk','P4_MOTION_PASS','P4_PASS'):
        if d.get(k) is not False:f.append(k+' not false')
    if float(d.get('P3_delta',0))!=1e-18:f.append('P3 delta changed')
    if not 0.64<float(d.get('finite_angle_vector_information_retention_lower',0))<1:f.append('retention invalid')
    if d.get('x_cell_failures'):f.append('finite-angle H18 x-cell failure')
    return f

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'k':d['finite_angle_vector_information_retention_lower'],'alpha6':d['finite_angle_eta6_information_lower'],'D18':d['finite_angle_H18_information_lower'],'pivot':d['worst_full_H18_LDLT_pivot_lower'],'closed':d['universal_H18_finite_angle_prior_free_LDLT_closed'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())