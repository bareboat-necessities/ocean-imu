#!/usr/bin/env python3
"""Source-uniform binary32 ISS enclosure for OU-III Kalman/reset arithmetic.

Explicit scalar-loop roundoff uses conservative source-uniform magnitude
envelopes from the canonical endpoint-referenced covariance ceiling.  The
independent rowwise K*residual correction bound may be enormous; that is
acceptable for avoiding overflow/roundoff omissions, but it is NOT a
mathematical reset-domain or storage bound.  The exact reset sector obtains its
correction domain from the same augmented graph instead.

The abandoned scalar moving-Riccati relative-contraction margin is not consumed
by this arithmetic enclosure.  Eigen LDLT and libm/time update are consumed
through explicit conditional platform postconditions. Roundoff remains an
additive practical-ISS input.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
import ou3_p4_rowwise_coefficient_enclosure_fast as COEFF
import ou3_brmm_riccati_tube_factored as TUBE
import ou3_p4_hard_entry_set as ENTRY
import ou3_p4_projection_binary32_enclosure as PROJFP
import ou3_p4_innovation_binary32_bounds as INNOV
import ou3_p4_binary32_platform_arithmetic_contract as PLATFORM
U=2.0**-24;NX={'H18':18,'A21':21}
GROUPS_H=('attitude_cayley_norm','gyro_bias_norm_rad_s','velocity_norm_mps','position_norm_m','integral_displacement_norm_m_s','latent_acceleration_norm_mps2');GROUPS_A=GROUPS_H+('accelerometer_bias_error_norm_mps2',)
def gamma(k):
    if k<0 or k*U>=1:raise ValueError('invalid gamma index')
    return math.nextafter((k*U)/(1-k*U),math.inf)
def up(x):return math.nextafter(float(x),math.inf)
def build():
    coeff=COEFF.build();tube=TUBE.build();entry=ENTRY.build();proj=PROJFP.build();innov=INNOV.build();platform=PLATFORM.build();bad={'coeff':COEFF.validate(coeff),'endpoint_covariance_envelope':TUBE.validate_covariance_ceiling(tube),'entry':ENTRY.validate(entry),'projection':PROJFP.validate(proj),'innovation':INNOV.validate(innov),'platform':PLATFORM.validate(platform)};bad={k:v for k,v in bad.items() if v}
    if bad:raise RuntimeError('binary32 ISS prerequisites failed: '+repr(bad))
    er=entry['coordinate_radii'];gain_eps=float(platform['eigen_ldlt3_gain_row_forward_relative_upper']);time_eps=float(platform['time_update_state_forward_relative_upper']);modes={}
    for mode,tkey,groups in (('H18','H',GROUPS_H),('A21','A',GROUPS_A)):
        n=NX[mode];mm=coeff['modes'][mode];pdiag=[float(x) for x in tube['modes'][tkey]['Pbar_diagonal_variance_upper']];pabs=up(max(pdiag))
        residual=max(float(v) for v in mm['hard_entry_residual_norm_upper_diagnostic'].values());kmax=0.0
        for row in mm['rows'].values():
            for key in ('K_row_norm_upper_accelerometer','K_row_norm_upper_magnetometer','K_row_norm_upper_S_zero'):kmax=max(kmax,max(map(float,row[key])))
        correction_outer=float(mm['independent_rowbox_attitude_correction_norm_upper_diagnostic'])
        correction_component_round=up(gamma(5)*3.0*kmax*residual);correction_norm_round=up(math.sqrt(n)*correction_component_round)
        kcp_abs=up(3.0*kmax*pabs);kcp_round=up(gamma(5)*3.0*kmax*pabs);ksk_round=max(float(x['KSK_roundoff_abs_upper_per_covariance_entry']) for x in innov['modes'][mode]['events'].values());symmetry_relative=gamma(2)
        # Arithmetic-only outer reset magnitude. This does NOT certify the
        # reset Cayley chart and is not consumed by exact-real storage.
        reset_matvec_component_round=up(gamma(5)*(1.0+correction_outer)*max(1.0,residual));reset_matvec_norm_round=up(math.sqrt(3)*reset_matvec_component_round)
        ldlt_component=up(gain_eps*max(1.0,kmax)*residual);ldlt_state=up(math.sqrt(n)*ldlt_component);state_radius=up(math.sqrt(sum(float(er[g])**2 for g in groups)));time_state=up(time_eps*max(1.0,state_radius));one_event_state=up(correction_norm_round+ldlt_state+reset_matvec_norm_round+float(proj['single_evaluation_projection_value_error_norm_upper_mps2']))
        modes[mode]={'state_dimension':n,'hard_entry_full_state_norm_upper':state_radius,'P_entry_abs_upper_from_PSD_diagonal_ceiling':pabs,'K_scalar_abs_upper_from_row_norms':kmax,'physical_residual_norm_upper_arithmetic_only':residual,'independent_rowbox_correction_norm_upper_arithmetic_only':correction_outer,'independent_rowbox_correction_not_reset_domain':True,'correction_matvec_roundoff_norm_upper':correction_norm_round,'KCP_term_abs_upper':kcp_abs,'KCP_roundoff_abs_upper_per_covariance_entry':kcp_round,'KSK_roundoff_abs_upper_per_covariance_entry':ksk_round,'symmetry_average_relative_roundoff_upper':symmetry_relative,'reset_state_matvec_roundoff_norm_upper':reset_matvec_norm_round,'projection_roundoff_norm_upper':float(proj['single_evaluation_projection_value_error_norm_upper_mps2']),'Eigen_LDLT_gain_solve_additive_state_error_norm_upper_per_accepted_update':ldlt_state,'time_update_libm_and_coefficient_additive_state_error_norm_upper_per_prediction':time_state,'explicit_measurement_reset_state_roundoff_norm_upper_per_event':one_event_state,'explicit_scalar_loop_roundoff_channels_finite':all(math.isfinite(x) for x in (correction_norm_round,kcp_round,ksk_round,symmetry_relative,reset_matvec_norm_round)),'KSK_roundoff_closed':bool(innov['modes'][mode]['all_event_innovation_and_KSK_bounds_closed']),'Eigen_LDLT_gain_solve_roundoff_closed_conditionally':math.isfinite(ldlt_state),'time_update_libm_coefficient_roundoff_closed_conditionally':math.isfinite(time_state),'all_additive_state_roundoff_channels_finite':all(math.isfinite(x) for x in (one_event_state,time_state))}
    explicit=all(m['explicit_scalar_loop_roundoff_channels_finite'] and m['KSK_roundoff_closed'] for m in modes.values());conditional=explicit and all(m['Eigen_LDLT_gain_solve_roundoff_closed_conditionally'] and m['time_update_libm_coefficient_roundoff_closed_conditionally'] and m['all_additive_state_roundoff_channels_finite'] for m in modes.values())
    return {'qualification':'OU3_P4_BINARY32_KALMAN_RESET_ISS_V5','runtime_scalar_format':'IEEE754_binary32','unit_roundoff':U,
      'canonical_endpoint_referenced_covariance_envelope_consumed':True,'failed_moving_Riccati_relative_margin_consumed':False,
      'source_uniform_magnitude_envelopes_only_not_word_composition':True,'rowwise_K_reset_domain_used':False,'rowwise_K_used_for_roundoff_magnitude_only':True,'PSD_cross_covariance_bound_used':'|Pij|<=sqrt(Pii*Pjj)','FMA_contraction_safe_by_separate_operation_overcount':True,'innovation_and_KSK_binary32_enclosure_consumed':True,'explicit_Joseph_reset_scalar_loop_roundoff_enclosed':explicit,'projection_binary32_enclosure_consumed':True,'platform_execution_contract_consumed':True,'platform_execution_contract_qualification':platform['qualification'],'platform_contract_target_qualified_here':False,'modes':modes,'remaining_arithmetic_prerequisites':[],'full_shipping_Kalman_reset_finite_precision_enclosure_closed_conditionally':conditional,'additive_ISS_channel_complete_for_conditional_P4':conditional,'deployment_finite_precision_qualification_closed':False,'roundoff_absorbed_into_contraction_margin':False,'roundoff_is_explicit_additive_ISS_input':True,'P4_MOTION_PASS':False,'P4_PASS':False,'P4_DEPLOYMENT_PASS':False}
def validate(d):
    f=[]
    if d.get('qualification')!='OU3_P4_BINARY32_KALMAN_RESET_ISS_V5':f.append('qualification mismatch')
    for k in ('canonical_endpoint_referenced_covariance_envelope_consumed','source_uniform_magnitude_envelopes_only_not_word_composition','rowwise_K_used_for_roundoff_magnitude_only','FMA_contraction_safe_by_separate_operation_overcount','innovation_and_KSK_binary32_enclosure_consumed','explicit_Joseph_reset_scalar_loop_roundoff_enclosed','projection_binary32_enclosure_consumed','platform_execution_contract_consumed','full_shipping_Kalman_reset_finite_precision_enclosure_closed_conditionally','additive_ISS_channel_complete_for_conditional_P4','roundoff_is_explicit_additive_ISS_input'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('failed_moving_Riccati_relative_margin_consumed','rowwise_K_reset_domain_used','platform_contract_target_qualified_here','deployment_finite_precision_qualification_closed','roundoff_absorbed_into_contraction_margin','P4_MOTION_PASS','P4_PASS','P4_DEPLOYMENT_PASS'):
        if d.get(k) is not False:f.append(k+' not false')
    if d.get('remaining_arithmetic_prerequisites')!=[]:f.append('conditional arithmetic blocker remains')
    for mode,m in d.get('modes',{}).items():
        if m.get('independent_rowbox_correction_not_reset_domain') is not True:f.append(mode+' arithmetic correction mislabeled')
        for k in ('explicit_scalar_loop_roundoff_channels_finite','KSK_roundoff_closed','Eigen_LDLT_gain_solve_roundoff_closed_conditionally','time_update_libm_coefficient_roundoff_closed_conditionally','all_additive_state_roundoff_channels_finite'):
            if m.get(k) is not True:f.append(mode+' '+k+' not true')
        for k in ('Eigen_LDLT_gain_solve_additive_state_error_norm_upper_per_accepted_update','time_update_libm_and_coefficient_additive_state_error_norm_upper_per_prediction','explicit_measurement_reset_state_roundoff_norm_upper_per_event'):
            if not(math.isfinite(float(m.get(k,math.nan))) and float(m[k])>=0):f.append(mode+' '+k+' invalid')
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'conditional_full':d['full_shipping_Kalman_reset_finite_precision_enclosure_closed_conditionally'],'deployment':d['deployment_finite_precision_qualification_closed'],'H18':d['modes']['H18'],'A21':d['modes']['A21'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())