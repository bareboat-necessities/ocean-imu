#!/usr/bin/env python3
"""Tight covariance/spacing read of the complete-BRMM four-S information lemma.

This is not a second source, word, estimator or proof architecture.  It invokes
``ou3_brmm_four_s_translation_information`` as the source/parity prerequisite,
then sharpens two conservative choices while staying inside the SAME canonical
3 s word.

First, configured S measurement noise is independent between selected updates,
so its 4x4 covariance contributes max_i R_{S,i}^2 rather than four times that
value.  Process nuisance may remain correlated and keeps its PSD trace bound.

Second, the base lemma selected the first four convenient scheduler witnesses
in windows [0,g],[2g,3g],[4g,5g],[6g,7g].  The progress certificate guarantees
an actual S firing in every g-window.  We may therefore select, from the same
literal word, the better-conditioned windows

    [0,g], [4g,5g], [8g,9g], [12g,13g],

while retaining every other shipping S update.  Since 13g <= 1.95 s < 3 s,
these witnesses remain a subset of the canonical word.  All OU decay and
process-noise penalties are recomputed over the longer selected subword.

No gate, filter constant, source family or runtime operation is changed.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path

import ou3_brmm_four_s_translation_information as BASE

DEFAULT_DOMAIN = BASE.DEFAULT_DOMAIN
SCHEMA = 2
QUALIFICATION = "OU3_COMPLETE_BRMM_FOUR_S_PHYSICAL_INFORMATION_TIGHT_COVARIANCE_SPACING"
SPACING_M = 4
WINDOWS = ((0.0,1.0),(4.0,5.0),(8.0,9.0),(12.0,13.0))
CANONICAL_WORD_S = 3.0


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def _newton_rows(m: int) -> tuple[list[float], dict]:
    d01=Fraction(m-1); d02=Fraction(2*m-1); d03=Fraction(3*m-1)
    d12=Fraction(m-1); d13=Fraction(2*m-1); d23=Fraction(m-1)
    r0=Fraction(1)
    r1=Fraction(2,1)/d01
    r2=sum((Fraction(1,1)/(d01*d02),Fraction(1,1)/(d01*d12),Fraction(1,1)/(d02*d12)),Fraction(0))
    r3=sum((Fraction(1,1)/(d01*d02*d03),Fraction(1,1)/(d01*d12*d13),Fraction(1,1)/(d02*d12*d23),Fraction(1,1)/(d03*d13*d23)),Fraction(0))
    return [up(float(x)) for x in (r0,r1,r2,r3)], {
        'scaled_pair_gap_lower': {'u1-u0':float(d01),'u2-u0':float(d02),'u3-u0':float(d03),'u2-u1':float(d12),'u3-u1':float(d13),'u3-u2':float(d23)},
        'divided_difference_raw_row_l1_exact':[str(x) for x in (r0,r1,r2,r3)],
        'divided_difference_raw_row_l1_upper':[up(float(x)) for x in (r0,r1,r2,r3)],
        'exact_rational_arithmetic_used_before_outward_float_conversion':True,
    }


def _physical_rows(m: int, third: float) -> tuple[dict[str,float], dict]:
    q, newton = _newton_rows(m); q0,q1,q2,q3=q
    c0=up(1.0/6.0)
    c01=up(((m+1.0)**2)/2.0)
    c012=up((2.0*m+1.0)/2.0)
    halfsum=up((m+2.0)/2.0)
    rowA=up(q3/third)
    rowV=up(2.0*up(q2+up(c012*rowA)))
    rowP=up(q1+up(halfsum*rowV)+up(c01*rowA))
    rowS=up(q0+rowP+up(0.5*rowV)+up(c0*rowA))
    rows={'S':rowS,'g*p':rowP,'g^2*v':rowV,'g^3*a_w':rowA}
    frob=0.0
    for r in rows.values(): frob=up(frob+up(r*r))
    mtm=down(1.0/frob)
    physical={
        'state_order':['S','g*p','g^2*v','g^3*a_w'],
        'third_divided_difference_lower':third,
        'third_divided_difference_enters_quantitative_inverse':True,
        'derivative_upper_bounds':{
            'c_u0':c0,'c_first_divided_difference':c01,'c_second_divided_difference':c012,
            'basis_half_u0_plus_u1':halfsum,'basis_u0':1.0,'basis_half_u0_squared':0.5},
        'physical_state_inverse_raw_record_row_l1_upper':rows,
        'physical_state_inverse_frobenius_squared_upper':frob,
        'physical_observation_MtM_lambda_min_lower':mtm,
        'inverse_matrix_frobenius_norm_bound_used':True,
        'determinant_used_for_quantitative_bound':False,
        'raw_to_Newton_condition_alone_used_as_physical_bound':False,
    }
    return rows, {**newton, 'physical_state_recovery':physical}


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    d = BASE.build(Path(domain_path).resolve())
    failures = BASE.validate(d)
    if failures:
        raise RuntimeError(f"base four-S physical information invalid: {failures}")

    g=float(d['uniform_S_gap_s_upper'])
    T=up(13.0*g)
    if not (T < CANONICAL_WORD_S):
        raise RuntimeError(f"spread four-S witnesses leave canonical word: {T}")
    tau_lo=float(d['tau_applied_s'][0])
    decay=up(T/tau_lo)
    aresp, exp_cert=BASE._validated_exp_negative_lower(decay)
    third=down(aresp/6.0)
    rows, newton = _physical_rows(SPACING_M, third)
    physical=newton['physical_state_recovery']

    base_noise=d['selected_S_record_noise']
    measurement_lambda=up(max(map(float,base_noise['measurement_variance_axis_upper'])))
    qc=float(base_noise['OU_driving_intensity_upper'])
    process_per_record=up(qc*up(T**7)/252.0)
    process_stack_lambda=up(4.0*process_per_record)
    total_lambda=up(measurement_lambda+process_stack_lambda)
    inv_lower=down(1.0/total_lambda)

    mtm=float(physical['physical_observation_MtM_lambda_min_lower'])
    info=down(inv_lower*mtm)
    if not (math.isfinite(info) and info>0.0):
        raise RuntimeError('spread four-S physical information lower is not strict')

    old_noise=dict(base_noise)
    old_lambda=float(old_noise['four_record_covariance_lambda_max_upper'])
    old_info=float(d['newton_coordinate_information']['D_S_physical_lambda_min_lower'])
    noise=dict(old_noise)
    noise.update({
        'legacy_loose_four_record_covariance_lambda_max_upper_diagnostic':old_lambda,
        'measurement_covariance_is_diagonal_across_selected_updates':True,
        'measurement_covariance_lambda_max_upper':measurement_lambda,
        'process_covariance_may_be_correlated_across_selected_updates':True,
        'process_S_variance_per_record_upper':process_per_record,
        'process_four_record_covariance_lambda_max_trace_upper':process_stack_lambda,
        'four_record_covariance_lambda_max_upper':total_lambda,
        'Sigma_Y_inverse_scalar_lower':inv_lower,
        'measurement_noise_not_multiplied_by_record_count':True,
        'process_cross_record_correlation_still_covered_by_trace_bound':True,
    })
    ni={**d['newton_coordinate_information'], **newton}
    ni.update({
        'raw_record_inverse_covariance_scalar_lower':inv_lower,
        'third_divided_difference_lower_enters_quantitative_bound':True,
        'D_S_physical_lambda_min_lower':info,
        'D_S_newton_lambda_min_lower':info,
        'D_S_newton_matrix_lower':f"D_S,z >= {info:.17g} * I_4",
        'tight_covariance_bound_consumed':True,
        'spread_scheduler_witnesses_consumed':True,
        'base_corrected_physical_information_lower_diagnostic':old_info,
    })

    out=dict(d)
    out.update({
        'schema':SCHEMA,
        'qualification':QUALIFICATION,
        'scaled_u_windows':[[a,b] for a,b in WINDOWS],
        'four_S_windows_s':[[down(a*g),up(b*g)] for a,b in WINDOWS],
        'selected_window_horizon_s_upper':T,
        'validated_exponential_lower_certificate':exp_cert,
        'aw_homogeneous_response_lower':aresp,
        'aw_scaled_third_divided_difference_lower':third,
        'selected_S_record_noise':noise,
        'newton_coordinate_information':ni,
        'same_complete_BRMM_component_sharpened':True,
        'spread_scheduler_witnesses_inside_same_word':True,
        'selected_spacing_multiple':SPACING_M,
        'canonical_word_horizon_s':CANONICAL_WORD_S,
        'source_family_replaced':False,
        'P3_architecture_replaced':False,
        'P3_promoted':False,
    })
    return out


def validate(d: dict) -> list[str]:
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    if d.get('canonical_source')!='COMPLETE_BRMM_NORMAL_LIVE_WORD':f.append('canonical source changed')
    for key in ('component_of_complete_BRMM_full_word','actual_applied_SpectralMSE_R_S_consumed','all_due_S_updates_remain_in_literal_word','same_complete_BRMM_component_sharpened','spread_scheduler_witnesses_inside_same_word'):
        if d.get(key) is not True:f.append(f'{key} is not true')
    for key in ('source_family_replaced','P3_architecture_replaced','P3_promoted'):
        if d.get(key) is not False:f.append(f'{key} is not false')
    if d.get('selected_spacing_multiple')!=4:f.append('four-S spacing changed')
    if not float(d.get('selected_window_horizon_s_upper',math.inf))<float(d.get('canonical_word_horizon_s',0)):f.append('selected S witnesses do not fit same word')
    noise=d.get('selected_S_record_noise',{})
    for key in ('measurement_covariance_is_diagonal_across_selected_updates','process_covariance_may_be_correlated_across_selected_updates','measurement_noise_not_multiplied_by_record_count','process_cross_record_correlation_still_covered_by_trace_bound'):
        if noise.get(key) is not True:f.append('noise structure missing '+key)
    ni=d.get('newton_coordinate_information',{})
    if ni.get('third_divided_difference_lower_enters_quantitative_bound') is not True:f.append('physical divided-difference correction was lost')
    if ni.get('tight_covariance_bound_consumed') is not True:f.append('tight covariance bound not consumed')
    if ni.get('spread_scheduler_witnesses_consumed') is not True:f.append('spread scheduler witnesses not consumed')
    rows=ni.get('physical_state_recovery',{}).get('physical_state_inverse_raw_record_row_l1_upper',{})
    for k in ('S','g*p','g^2*v','g^3*a_w'):
        if not (isinstance(rows.get(k),(int,float)) and float(rows[k])>0):f.append('invalid physical inverse row '+k)
    info=ni.get('D_S_physical_lambda_min_lower')
    if not isinstance(info,(int,float)) or not (math.isfinite(float(info)) and float(info)>0):f.append('tight physical information lower is not strict')
    return list(dict.fromkeys(f))


def main() -> int:
    ap=argparse.ArgumentParser();ap.add_argument('--domain',type=Path,default=DEFAULT_DOMAIN);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    d=build(args.domain);failures=validate(d);d['validation_pass']=not failures;d['validation_failures']=failures
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(d,indent=2,sort_keys=True),encoding='utf-8')
    n=d['selected_S_record_noise'];ni=d['newton_coordinate_information']
    print(json.dumps({'spacing':d['selected_spacing_multiple'],'horizon_s':d['selected_window_horizon_s_upper'],'measurement_cov_lambda_max':n['measurement_covariance_lambda_max_upper'],'process_stack_cov_lambda_max':n['process_four_record_covariance_lambda_max_trace_upper'],'total_cov_lambda_max':n['four_record_covariance_lambda_max_upper'],'physical_rows':ni['physical_state_recovery']['physical_state_inverse_raw_record_row_l1_upper'],'tight_physical_information':ni['D_S_physical_lambda_min_lower'],'failures':failures},indent=2,sort_keys=True))
    return 0 if not failures else 2

if __name__=='__main__':raise SystemExit(main())
