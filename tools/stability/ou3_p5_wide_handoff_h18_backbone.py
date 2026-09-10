#!/usr/bin/env python3
"""P5 wide-handoff H18 finite-angle differential backbone.

The shipping handoff theorem places the tilt quotient below pi/2+0.02 rad.
For P5 transient capture we recompute the same prior-free H18 information/LDLT
chain on a 1.60-rad sector, strictly containing the handoff section. The P4-only
0.64 information-retention convenience floor is not a geometric singularity and
is not used here. This is a differential/detectability backbone, not yet the
finite-map capture theorem.
"""
from __future__ import annotations
import copy,json,math
from pathlib import Path
import ou3_p4_cayley_sector_certificate as CAYLEY
import ou3_brmm_h18_information_composition as HINFO
import ou3_brmm_h18_prior_free_completion as HPF
import ou3_brmm_dynamic_source_certificate as DYNAMIC
import ou3_full_process_ucc as PROCESS
import ou3_brmm_full_word_event_algebra as EVENT
import ou3_startup_handoff_tilt_hemisphere as HANDOFF
QUALIFICATION='OU3_P5_WIDE_HANDOFF_H18_BACKBONE_V1';OUTER_ANGLE_RAD=1.60

def build():
 hand=HANDOFF.build()
 if HANDOFF.validate(hand):raise RuntimeError('handoff tilt prerequisite invalid')
 hu=float(hand['true_gravity_quotient_tilt_strict_upper_rad'])
 if not hu<OUTER_ANGLE_RAD<math.pi:raise RuntimeError('wide sector does not contain handoff')
 cay=CAYLEY.build(outer_angle_rad=OUTER_ANGLE_RAD)
 if not(cay['source_generated_not_trajectory_fit'] and cay['chart_antipode_excluded'] and float(cay['chart_sigma_min_lower'])>0 and float(cay['exact_vector_information_retention_factor_lower'])>0):raise RuntimeError('wide Cayley geometry invalid')
 base=HINFO.build();dyn=DYNAMIC.build();proc=PROCESS.build();event=EVENT.build();bad={'hinfo':HINFO.validate(base),'dynamic':DYNAMIC.validate(dyn),'process':PROCESS.validate(proc),'event':EVENT.validate(event)};bad={k:v for k,v in bad.items() if v}
 if bad:raise RuntimeError('wide H18 prerequisites failed: '+repr(bad))
 k=float(cay['exact_vector_information_retention_factor_lower']);alpha=math.nextafter(k*float(base['eta6_information_lower']),-math.inf);comp=base['triangular_information_composition'];daw=float(comp['aw_direction_information_lower']);cross=float(comp['C_aw_spectral_norm_squared_upper']);trace=math.nextafter(alpha+daw+cross,math.inf);coupled=math.nextafter(alpha*daw/trace,-math.inf);D=math.nextafter(min(coupled,float(comp['non_aw_translation_lambda_min_lower'])),-math.inf)
 if min(alpha,coupled,D)<=0:raise RuntimeError('wide information lost positivity')
 h=copy.deepcopy(base);h['eta6_information_lower']=alpha;h['H18_information_useful_gate_pass']=D>=1e-18;h['triangular_information_composition']['A_transpose_A_lower']=alpha;h['triangular_information_composition']['coupled_eta6_aw_scalar_2x2_trace_upper']=trace;h['triangular_information_composition']['coupled_eta6_aw_scalar_2x2_determinant_lower']=math.nextafter(alpha*daw,-math.inf);h['triangular_information_composition']['coupled_eta6_aw_lambda_min_lower']=coupled;h['triangular_information_composition']['D_H18_lambda_min_lower']=D
 pbar=HPF._same_word_covariance_upper(Path(HPF.DEFAULT_DOMAIN).resolve(),dyn,proc,h);fnorm=HPF._prediction_norm_sq_upper(proc);penalty=math.nextafter((HPF.USEFUL_GATE**2/4.0)*fnorm*float(pbar['Pbar_trace_upper']),math.inf)
 rows=[];fails=[];worst=math.inf
 for x in HPF._x_cover(dyn):
  ok,row=HPF._full_H18_cell(x,process=proc,dynamic=dyn,penalty_physical=penalty);rows.append(row)
  if ok:worst=min(worst,float(row['pivot_lower']))
  else:fails.append(row)
 closed=bool(rows) and not fails and math.isfinite(worst) and worst>0
 preserve=event['full_matrix_margin_preservation'];suffix=all(bool(preserve[q]) for q in ('covers_prediction','covers_every_due_S_update','covers_every_Normal_Live_accelerometer_update','covers_asynchronous_magnetometer_update','covers_immediate_left_error_reset','covers_aw_covariance_floor','covers_not_due_or_rejected_identity_branches'))
 return {'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','filter_changed':False,'quality_gates_changed':False,'domain_shrunk':False,'source_enumeration_used':False,'trajectory_replay_used':False,'shipping_handoff_tilt_strict_upper_rad':hu,'shipping_handoff_tilt_strict_upper_deg':hand['true_gravity_quotient_tilt_strict_upper_deg'],'transient_outer_angle_rad':OUTER_ANGLE_RAD,'transient_outer_angle_deg':math.degrees(OUTER_ANGLE_RAD),'strictly_contains_complete_handoff_tilt_section':hu<OUTER_ANGLE_RAD,'P4_usable_floor_is_not_consumed_by_P5_backbone':True,'wide_chart_sigma_min_lower':float(cay['chart_sigma_min_lower']),'wide_vector_information_retention_lower':k,'wide_finite_angle_H18_information_lower':D,'wide_information_above_P3_delta':D>=1e-18,'x_cells_certified':len(rows),'x_cell_failures':fails,'worst_full_H18_LDLT_pivot_lower':worst if closed else None,'WIDE_HANDOFF_H18_PRIOR_FREE_LDLT_CLOSED':closed,'linear_suffix_event_algebra_available':suffix,'exact_chord_Joseph_reset_capture_word_still_required':True,'finite_map_contraction_closed_here':False,'finite_capture_time_closed_here':False,'P4_PASS':False,'P5_PASS':False}
def validate(d):
 f=[]
 if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
 for q in ('strictly_contains_complete_handoff_tilt_section','P4_usable_floor_is_not_consumed_by_P5_backbone','wide_information_above_P3_delta','WIDE_HANDOFF_H18_PRIOR_FREE_LDLT_CLOSED','linear_suffix_event_algebra_available','exact_chord_Joseph_reset_capture_word_still_required'):
  if d.get(q) is not True:f.append(q+' not true')
 for q in ('filter_changed','quality_gates_changed','domain_shrunk','source_enumeration_used','trajectory_replay_used','finite_map_contraction_closed_here','finite_capture_time_closed_here','P4_PASS','P5_PASS'):
  if d.get(q) is not False:f.append(q+' not false')
 if d.get('x_cell_failures'):f.append('wide H18 x-cell failures')
 if not float(d.get('worst_full_H18_LDLT_pivot_lower') or 0)>0:f.append('nonpositive H18 pivot')
 return f
if __name__=='__main__':
 d=build();f=validate(d);print(json.dumps({**d,'validation_pass':not f,'validation_failures':f},indent=2,sort_keys=True));raise SystemExit(bool(f))
