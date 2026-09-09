#!/usr/bin/env python3
"""Source-uniform binary32 ISS enclosure for explicit OU-III Joseph/reset loops.

This producer bounds rounding in the scalar arithmetic that is explicit in the
shipping implementation: residual/correction mat-vecs, the hand-written Joseph
covariance loops, symmetry averaging, and finite reset matrix-vector products.
It uses the source-uniform Riccati diagonal ceiling and rowwise K bounds only as
magnitude ceilings; those boxes are never composed as a stability word.

The source-correlated innovation magnitude certificate closes the K S K^T loop.
The remaining non-explicit-library obligations are the 3x3 Eigen LDLT solve and
libm-generated time-update coefficients.  Until their implementation contracts
are supplied, full shipping finite precision remains fail-closed.
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

U=2.0**-24
NX={'H18':18,'A21':21}


def gamma(k:int)->float:
    if k<0 or k*U>=1: raise ValueError('invalid gamma index')
    return math.nextafter((k*U)/(1-k*U),math.inf)
def up(x): return math.nextafter(float(x),math.inf)


def build()->dict:
    coeff=COEFF.build(); tube=FASTTUBE.build_base(); entry=ENTRY.build(); proj=PROJFP.build(); innov=INNOV.build()
    bad={'coeff':COEFF.validate(coeff),'tube':TUBE.validate(tube),'entry':ENTRY.validate(entry),
         'projection':PROJFP.validate(proj),'innovation':INNOV.validate(innov)}
    bad={k:v for k,v in bad.items() if v}
    if bad: raise RuntimeError('binary32 ISS prerequisites failed: '+repr(bad))

    modes={}
    for mode,tkey in (('H18','H'),('A21','A')):
        n=NX[mode]; mm=coeff['modes'][mode]
        pdiag=[float(x) for x in tube['modes'][tkey]['Pbar_diagonal_variance_upper']]
        pabs=up(max(pdiag))
        residual=max(float(v) for v in mm['hard_entry_residual_norm_upper'].values())
        kmax=0.0
        for row in mm['rows'].values():
            for key in ('K_row_norm_upper_accelerometer','K_row_norm_upper_magnetometer','K_row_norm_upper_S_zero'):
                kmax=max(kmax,max(map(float,row[key])))

        correction_component_round=up(gamma(5)*3.0*kmax*residual)
        correction_norm_round=up(math.sqrt(n)*correction_component_round)
        kcp_abs=up(3.0*kmax*pabs)
        kcp_round=up(gamma(5)*3.0*kmax*pabs)
        ksk_round=max(float(x['KSK_roundoff_abs_upper_per_covariance_entry']) for x in innov['modes'][mode]['events'].values())
        symmetry_relative=gamma(2)
        dtheta=float(mm['attitude_correction_norm_upper'])
        reset_matvec_component_round=up(gamma(5)*(1.0+dtheta)*max(1.0,residual))
        reset_matvec_norm_round=up(math.sqrt(3)*reset_matvec_component_round)

        modes[mode]={
          'state_dimension':n,'P_entry_abs_upper_from_PSD_diagonal_ceiling':pabs,
          'K_scalar_abs_upper_from_row_norms':kmax,'physical_residual_norm_upper':residual,
          'correction_matvec_roundoff_norm_upper':correction_norm_round,
          'KCP_term_abs_upper':kcp_abs,'KCP_roundoff_abs_upper_per_covariance_entry':kcp_round,
          'KSK_roundoff_abs_upper_per_covariance_entry':ksk_round,
          'symmetry_average_relative_roundoff_upper':symmetry_relative,
          'reset_state_matvec_roundoff_norm_upper':reset_matvec_norm_round,
          'projection_roundoff_norm_upper':float(proj['single_evaluation_projection_value_error_norm_upper_mps2']),
          'explicit_scalar_loop_roundoff_channels_finite':all(math.isfinite(x) for x in (
             correction_norm_round,kcp_round,ksk_round,symmetry_relative,reset_matvec_norm_round)),
          'KSK_roundoff_closed':bool(innov['modes'][mode]['all_event_innovation_and_KSK_bounds_closed']),
          'Eigen_LDLT_gain_solve_roundoff_closed':False,
          'time_update_libm_coefficient_roundoff_closed':False,
        }

    explicit=all(m['explicit_scalar_loop_roundoff_channels_finite'] and m['KSK_roundoff_closed'] for m in modes.values())
    full=False
    return {
      'qualification':'OU3_P4_BINARY32_KALMAN_RESET_ISS_V2',
      'runtime_scalar_format':'IEEE754_binary32','unit_roundoff':U,
      'dependency_reduced_smallx_tube_used_for_same_bounds':True,
      'source_uniform_magnitude_envelopes_only_not_word_composition':True,
      'PSD_cross_covariance_bound_used':'|Pij|<=sqrt(Pii*Pjj)',
      'FMA_contraction_safe_by_separate_operation_overcount':True,
      'innovation_and_KSK_binary32_enclosure_consumed':True,
      'explicit_Joseph_reset_scalar_loop_roundoff_enclosed':explicit,
      'projection_binary32_enclosure_consumed':True,
      'modes':modes,
      'remaining_arithmetic_prerequisites':[
        'Eigen 3x3 LDLT solve backward-error contract including safe_ldlt3 bump branch',
        'libm sin/cos/exp/expm1/sqrt ulp contracts for time-update/tuner coefficient generation'],
      'full_shipping_Kalman_reset_finite_precision_enclosure_closed':full,
      'additive_ISS_channel_complete_for_P4':False,
      'P4_MOTION_PASS':False,'P4_PASS':False,
    }


def validate(d):
    f=[]
    if d.get('qualification')!='OU3_P4_BINARY32_KALMAN_RESET_ISS_V2':f.append('qualification mismatch')
    for k in ('dependency_reduced_smallx_tube_used_for_same_bounds','source_uniform_magnitude_envelopes_only_not_word_composition',
              'FMA_contraction_safe_by_separate_operation_overcount','innovation_and_KSK_binary32_enclosure_consumed',
              'explicit_Joseph_reset_scalar_loop_roundoff_enclosed','projection_binary32_enclosure_consumed'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('full_shipping_Kalman_reset_finite_precision_enclosure_closed','additive_ISS_channel_complete_for_P4','P4_MOTION_PASS','P4_PASS'):
        if d.get(k) is not False:f.append(k+' not false')
    if len(d.get('remaining_arithmetic_prerequisites',[]))!=2:f.append('arithmetic blockers changed')
    for mode,m in d.get('modes',{}).items():
        if m.get('explicit_scalar_loop_roundoff_channels_finite') is not True:f.append(mode+' explicit roundoff not finite')
        if m.get('KSK_roundoff_closed') is not True:f.append(mode+' KSK roundoff not closed')
        for k in ('Eigen_LDLT_gain_solve_roundoff_closed','time_update_libm_coefficient_roundoff_closed'):
            if m.get(k) is not False:f.append(mode+' '+k+' falsely closed')
    return f


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'explicit':d['explicit_Joseph_reset_scalar_loop_roundoff_enclosed'],
      'full':d['full_shipping_Kalman_reset_finite_precision_enclosure_closed'],
      'blockers':d['remaining_arithmetic_prerequisites'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
