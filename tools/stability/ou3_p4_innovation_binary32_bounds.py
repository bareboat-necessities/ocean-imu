#!/usr/bin/env python3
"""Source-uniform innovation magnitude and hand-written KSK roundoff bounds.

For each shipping 3-vector measurement, S=H P H^T+R is SPD.  The canonical
endpoint-referenced factored Riccati covariance envelope gives a diagonal Pbar
ceiling, hence lambda_max(P)<=trace(Pbar).  This producer consumes only that
rigorous covariance-magnitude contract; the abandoned scalar moving-Riccati
relative-contraction margin is not an arithmetic premise.

Using the literal H block structure gives

  S-zero: ||H||_2^2 = 1,
  mag:    ||-[m]x||_2^2 <= ||m||^2,
  accel H18: ||[-f]x, R||_2^2 <= ||f||^2+1,
  accel A21: ||[-f]x, R, I||_2^2 <= ||f||^2+2.

Thus lambda_max(S)<=trace(Pbar)||H||^2+lambda_max(R), while
lambda_min(S)>=lambda_min(R).  The deployed S-zero anisotropy is read out of the
shipping header rather than restated here, and the horizontal pair is not
isotropic.  R_S=diag((rho_x r_S)^2,(rho_y r_S)^2,r_S^2), so both R eigenvalue
bounds are monotone in the axis factors: lambda_min falls with the smallest
deployed factor and lambda_max rises with the largest.  The pair therefore
enters as min/max over the three axes, and which axis binds each side is
reported rather than assumed -- on the deployed pair the smaller horizontal
factor binds Rmin and the vertical factor binds Rmax.

The shipping Joseph loop forms K S K^T with nine scalar triple products.  An
absolute binary32 error bound follows from the source-uniform scalar K ceiling,
|S_ab|<=lambda_max(S), and an overcounted separate multiply/add model.  FMA
contraction can only reduce that rounding count.  These bounds are arithmetic
magnitude enclosures only; they do not compose independent K/S boxes as a word
or supply a reset correction domain.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path

import ou3_brmm_riccati_tube_factored as TUBE
import ou3_brmm_dynamic_source_certificate as DYNAMIC
import ou3_brmm_complete_source as SOURCE
import ou3_p4_rowwise_coefficient_enclosure_fast as COEFF

U=2.0**-24

def deployed_axis_std_factors():
    """Deployed [rho_x,rho_y,1] R_S std factors; fail closed on an unreadable pair."""
    axis,reason=SOURCE.deployed_axis_std_factors()
    if axis is None:raise RuntimeError('deployed R_S axis factors unreadable: '+reason)
    return axis

def up(x):return math.nextafter(float(x),math.inf)
def gamma(k):
    if k<0 or k*U>=1:raise ValueError('invalid gamma')
    return up((k*U)/(1-k*U))

def build():
    tube=TUBE.build(); dyn=DYNAMIC.build(); coeff=COEFF.build()
    bad={'endpoint_covariance_envelope':TUBE.validate_covariance_ceiling(tube),'dynamic':DYNAMIC.validate(dyn),'coeff':COEFF.validate(coeff)}
    bad={k:v for k,v in bad.items() if v}
    if bad:raise RuntimeError('innovation prerequisites failed: '+repr(bad))
    domain=json.loads(TUBE.DEFAULT_DOMAIN.read_text());live=domain['normal_live']; fmax=float(live['specific_force_norm_upper_mps2']);mmax=float(live['magnetic_vector_norm_upper_uT'])
    rslo,rshi=map(float,dyn['dynamic_invariant']['R_S_applied'])
    axis=deployed_axis_std_factors();axis_min=min(axis);axis_max=max(axis)
    # Runtime configured vector measurement variances on this proof branch.  The
    # anisotropic S-zero axes enter through the monotone min/max of the deployed
    # factors, never through one horizontal scalar.
    r={'accelerometer':(0.2**2,0.2**2),'magnetometer':(0.3**2,0.3**2),
       'S_zero':((axis_min*rslo)**2,(axis_max*rshi)**2)}
    modes={}
    for mode,tkey in (('H18','H'),('A21','A')):
        pdiag=list(map(float,tube['modes'][tkey]['Pbar_diagonal_variance_upper']));ptrace=up(sum(pdiag))
        h2={'accelerometer':up(fmax*fmax+(1.0 if mode=='H18' else 2.0)),
            'magnetometer':up(mmax*mmax),'S_zero':1.0}
        kmax={}
        for ev,key in (('accelerometer','K_row_norm_upper_accelerometer'),('magnetometer','K_row_norm_upper_magnetometer'),('S_zero','K_row_norm_upper_S_zero')):
            kmax[ev]=max(max(map(float,row[key])) for row in coeff['modes'][mode]['rows'].values())
        events={}
        for ev in ('accelerometer','magnetometer','S_zero'):
            rmin,rmax=r[ev]
            smax=up(ptrace*h2[ev]+rmax);smin=math.nextafter(rmin,-math.inf)
            cond=up(smax/smin)
            # 9 terms. Each term K*S*K: two rounded multiplies; accumulation of
            # nine terms costs eight additions. gamma_26 safely overcounts all.
            abs_sum=up(9.0*kmax[ev]*kmax[ev]*smax)
            rnd=up(gamma(26)*abs_sum)
            events[ev]={'H_operator_norm_squared_upper':h2[ev],
              'R_lambda_min_lower':smin,'R_lambda_max_upper':rmax,
              'S_lambda_min_lower':smin,'S_lambda_max_upper':smax,
              'S_condition_number_upper':cond,'K_scalar_abs_upper':kmax[ev],
              'KSK_absolute_term_sum_upper_per_entry':abs_sum,
              'KSK_roundoff_abs_upper_per_covariance_entry':rnd,
              'KSK_roundoff_closed':math.isfinite(rnd)}
        modes[mode]={'Pbar_trace_upper':ptrace,'events':events,
          'all_event_innovation_and_KSK_bounds_closed':all(x['KSK_roundoff_closed'] for x in events.values())}
    return {'qualification':'OU3_P4_INNOVATION_AND_KSK_BINARY32_BOUNDS_V2','runtime_scalar_format':'IEEE754_binary32','unit_roundoff':U,
      'canonical_endpoint_referenced_covariance_envelope_consumed':True,
      'failed_moving_Riccati_relative_margin_consumed':False,
      'same_source_R_S_range_consumed':True,'RS_axis_factors_read_from_deployed_source':True,
      'actual_RS_axis_std_factors':axis,'actual_RS_horizontal_factors':axis[:2],
      'S_zero_R_lambda_min_binding_axis_factor':axis_min,'S_zero_R_lambda_max_binding_axis_factor':axis_max,
      'PSD_innovation_entry_bound_from_lambda_max':True,'FMA_safe_by_rounding_overcount':True,
      'independent_K_S_word_boxes_used':False,'rowwise_K_reset_domain_used':False,'modes':modes,
      'all_innovation_magnitude_bounds_closed':all(m['all_event_innovation_and_KSK_bounds_closed'] for m in modes.values()),
      'Eigen_LDLT_solve_closed_here':False,'full_P4_finite_precision_closed_here':False}

def validate(d):
    f=[]
    if d.get('qualification')!='OU3_P4_INNOVATION_AND_KSK_BINARY32_BOUNDS_V2':f.append('qualification mismatch')
    for k in ('canonical_endpoint_referenced_covariance_envelope_consumed','same_source_R_S_range_consumed','RS_axis_factors_read_from_deployed_source','PSD_innovation_entry_bound_from_lambda_max','FMA_safe_by_rounding_overcount','all_innovation_magnitude_bounds_closed'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('failed_moving_Riccati_relative_margin_consumed','independent_K_S_word_boxes_used','rowwise_K_reset_domain_used','Eigen_LDLT_solve_closed_here','full_P4_finite_precision_closed_here'):
        if d.get(k) is not False:f.append(k+' not false')
    axis,reason=SOURCE.deployed_axis_std_factors()
    if axis is None:f.append('deployed R_S axis factors unreadable: '+reason)
    else:
        if [float(x) for x in d.get('actual_RS_axis_std_factors',[])]!=axis:f.append('RS axis factors no longer match the deployed source')
        if [float(x) for x in d.get('actual_RS_horizontal_factors',[])]!=axis[:2]:f.append('RS horizontal pair no longer matches the deployed source')
        if float(d.get('S_zero_R_lambda_min_binding_axis_factor',0))!=min(axis) or float(d.get('S_zero_R_lambda_max_binding_axis_factor',0))!=max(axis):f.append('RS binding axis factors are not the deployed min/max')
    for mode,m in d.get('modes',{}).items():
        for ev,x in m['events'].items():
            if not (x['S_lambda_min_lower']>0 and x['S_lambda_max_upper']>=x['S_lambda_min_lower'] and x['KSK_roundoff_closed']):f.append(mode+' '+ev+' invalid')
    return f

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'closed':d['all_innovation_magnitude_bounds_closed'],'modes':d['modes'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
