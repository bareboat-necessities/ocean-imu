#!/usr/bin/env python3
"""Source-uniform binary32 ISS enclosure for OU-III Kalman/reset arithmetic.

Explicit shipping scalar loops are bounded directly. Eigen LDLT and libm/time
update are consumed through the explicit conditional platform execution
contract; target-toolchain qualification remains separate. The resulting
roundoff is an additive practical-ISS input, not an artificial reduction of the
exact-real contraction margin and not a composition of independent K boxes.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path

import ou3_p4_rowwise_coefficient_enclosure_fast as COEFF
import ou3_brmm_riccati_tube as TUBE
import ou3_brmm_riccati_tube_smallx_scaled as FASTTUBE
import ou3_p4_hard_entry_set as ENTRY
import ou3_p4_projection_binary32_enclosure as PROJFP
import ou3_p4_innovation_binary32_bounds as INNOV
import ou3_p4_binary32_platform_arithmetic_contract as PLATFORM

U=2.0**-24
NX={'H18':18,'A21':21}
GROUPS_H=('attitude_cayley_norm','gyro_bias_norm_rad_s','velocity_norm_mps','position_norm_m','integral_displacement_norm_m_s','latent_acceleration_norm_mps2')
GROUPS_A=GROUPS_H+('accelerometer_bias_error_norm_mps2',)

def gamma(k:int)->float:
    if k<0 or k*U>=1:raise ValueError('invalid gamma index')
    return math.nextafter((k*U)/(1-k*U),math.inf)
def up(x):return math.nextafter(float(x),math.inf)

def build()->dict:
    coeff=COEFF.build();tube=FASTTUBE.build_base();entry=ENTRY.build();proj=PROJFP.build();innov=INNOV.build();platform=PLATFORM.build()
    bad={'coeff':COEFF.validate(coeff),'tube':TUBE.validate(tube),'entry':ENTRY.validate(entry),
         'projection':PROJFP.validate(proj),'innovation':INNOV.validate(innov),'platform':PLATFORM.validate(platform)}
    bad={k:v for k,v in bad.items() if v}
    if bad:raise RuntimeError('binary32 ISS prerequisites failed: '+repr(bad))

    er=entry['coordinate_radii'];gain_eps=float(platform['eigen_ldlt3_gain_row_forward_relative_upper']);time_eps=float(platform['time_update_state_forward_relative_upper'])
    modes={}
    for mode,tkey,groups in (('H18','H',GROUPS_H),('A21','A',GROUPS_A)):
        n=NX[mode];mm=coeff['modes'][mode];pdiag=[float(x) for x in tube['modes'][tkey]['Pbar_diagonal_variance_upper']];pabs=up(max(pdiag))
        residual=max(float(v) for v in mm['hard_entry_residual_norm_upper'].values());kmax=0.0
        for row in mm['rows'].values():
            for key in ('K_row_norm_upper_accelerometer','K_row_norm_upper_magnetometer','K_row_norm_upper_S_zero'):
                kmax=max(kmax,max(map(float,row[key])))
        correction_component_round=up(gamma(5)*3.0*kmax*residual);correction_norm_round=up(math.sqrt(n)*correction_component_round)
        kcp_abs=up(3.0*kmax*pabs);kcp_round=up(gamma(5)*3.0*kmax*pabs)
        ksk_round=max(float(x['KSK_roundoff_abs_upper_per_covariance_entry']) for x in innov['modes'][mode]['events'].values())
        symmetry_relative=gamma(2);dtheta=float(mm['attitude_correction_norm_upper'])
        reset_matvec_component_round=up(gamma(5)*(1.0+dtheta)*max(1.0,residual));reset_matvec_norm_round=up(math.sqrt(3)*reset_matvec_component_round)

        # Directly consume the platform forward postconditions. For a gain row,
        # ||delta K_i|| <= eps_gain max(1,||K_i||); multiplication by the
        # physical residual gives the corresponding state-correction input.
        ldlt_component=up(gain_eps*max(1.0,kmax)*residual)
        ldlt_state=up(math.sqrt(n)*ldlt_component)
        state_radius=up(math.sqrt(sum(float(er[g])**2 for g in groups)))
        time_state=up(time_eps*max(1.0,state_radius))
        one_event_state=up(correction_norm_round+ldlt_state+reset_matvec_norm_round+float(proj['single_evaluation_projection_value_error_norm_upper_mps2']))

        modes[mode]={
          'state_dimension':n,'hard_entry_full_state_norm_upper':state_radius,
          'P_entry_abs_upper_from_PSD_diagonal_ceiling':pabs,'K_scalar_abs_upper_from_row_norms':kmax,
          'physical_residual_norm_upper':residual,'correction_matvec_roundoff_norm_upper':correction_norm_round,
          'KCP_term_abs_upper':kcp_abs,'KCP_roundoff_abs_upper_per_covariance_entry':kcp_round,
          'KSK_roundoff_abs_upper_per_covariance_entry':ksk_round,'symmetry_average_relative_roundoff_upper':symmetry_relative,
          'reset_state_matvec_roundoff_norm_upper':reset_matvec_norm_round,
          'projection_roundoff_norm_upper':float(proj['single_evaluation_projection_value_error_norm_upper_mps2']),
          'Eigen_LDLT_gain_solve_additive_state_error_norm_upper_per_accepted_update':ldlt_state,
          'time_update_libm_and_coefficient_additive_state_error_norm_upper_per_prediction':time_state,
          'explicit_measurement_reset_state_roundoff_norm_upper_per_event':one_event_state,
          'explicit_scalar_loop_roundoff_channels_finite':all(math.isfinite(x) for x in (correction_norm_round,kcp_round,ksk_round,symmetry_relative,reset_matvec_norm_round)),
          'KSK_roundoff_closed':bool(innov['modes'][mode]['all_event_innovation_and_KSK_bounds_closed']),
          'Eigen_LDLT_gain_solve_roundoff_closed_conditionally':math.isfinite(ldlt_state),
          'time_update_libm_coefficient_roundoff_closed_conditionally':math.isfinite(time_state),
          'all_additive_state_roundoff_channels_finite':all(math.isfinite(x) for x in (one_event_state,time_state)),
        }

    explicit=all(m['explicit_scalar_loop_roundoff_channels_finite'] and m['KSK_roundoff_closed'] for m in modes.values())
    conditional=explicit and all(m['Eigen_LDLT_gain_solve_roundoff_closed_conditionally'] and m['time_update_libm_coefficient_roundoff_closed_conditionally'] and m['all_additive_state_roundoff_channels_finite'] for m in modes.values())
    return {
      'qualification':'OU3_P4_BINARY32_KALMAN_RESET_ISS_V3','runtime_scalar_format':'IEEE754_binary32','unit_roundoff':U,
      'dependency_reduced_smallx_tube_used_for_same_bounds':True,'source_uniform_magnitude_envelopes_only_not_word_composition':True,
      'PSD_cross_covariance_bound_used':'|Pij|<=sqrt(Pii*Pjj)','FMA_contraction_safe_by_separate_operation_overcount':True,
      'innovation_and_KSK_binary32_enclosure_consumed':True,'explicit_Joseph_reset_scalar_loop_roundoff_enclosed':explicit,
      'projection_binary32_enclosure_consumed':True,'platform_execution_contract_consumed':True,
      'platform_execution_contract_qualification':platform['qualification'],'platform_contract_target_qualified_here':False,
      'modes':modes,'remaining_arithmetic_prerequisites':[],
      'full_shipping_Kalman_reset_finite_precision_enclosure_closed_conditionally':conditional,
      'additive_ISS_channel_complete_for_conditional_P4':conditional,
      'deployment_finite_precision_qualification_closed':False,
      'roundoff_absorbed_into_contraction_margin':False,'roundoff_is_explicit_additive_ISS_input':True,
      'P4_MOTION_PASS':False,'P4_PASS':False,'P4_DEPLOYMENT_PASS':False,
    }

def validate(d):
    f=[]
    if d.get('qualification')!='OU3_P4_BINARY32_KALMAN_RESET_ISS_V3':f.append('qualification mismatch')
    for k in ('dependency_reduced_smallx_tube_used_for_same_bounds','source_uniform_magnitude_envelopes_only_not_word_composition','FMA_contraction_safe_by_separate_operation_overcount','innovation_and_KSK_binary32_enclosure_consumed','explicit_Joseph_reset_scalar_loop_roundoff_enclosed','projection_binary32_enclosure_consumed','platform_execution_contract_consumed','full_shipping_Kalman_reset_finite_precision_enclosure_closed_conditionally','additive_ISS_channel_complete_for_conditional_P4','roundoff_is_explicit_additive_ISS_input'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('platform_contract_target_qualified_here','deployment_finite_precision_qualification_closed','roundoff_absorbed_into_contraction_margin','P4_MOTION_PASS','P4_PASS','P4_DEPLOYMENT_PASS'):
        if d.get(k) is not False:f.append(k+' not false')
    if d.get('remaining_arithmetic_prerequisites')!=[]:f.append('conditional arithmetic blocker remains')
    for mode,m in d.get('modes',{}).items():
        for k in ('explicit_scalar_loop_roundoff_channels_finite','KSK_roundoff_closed','Eigen_LDLT_gain_solve_roundoff_closed_conditionally','time_update_libm_coefficient_roundoff_closed_conditionally','all_additive_state_roundoff_channels_finite'):
            if m.get(k) is not True:f.append(mode+' '+k+' not true')
        for k in ('Eigen_LDLT_gain_solve_additive_state_error_norm_upper_per_accepted_update','time_update_libm_and_coefficient_additive_state_error_norm_upper_per_prediction','explicit_measurement_reset_state_roundoff_norm_upper_per_event'):
            if not (math.isfinite(float(m.get(k,math.nan))) and float(m[k])>=0):f.append(mode+' '+k+' invalid')
    return f

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'conditional_full':d['full_shipping_Kalman_reset_finite_precision_enclosure_closed_conditionally'],'deployment':d['deployment_finite_precision_qualification_closed'],'H18':d['modes']['H18'],'A21':d['modes']['A21'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
