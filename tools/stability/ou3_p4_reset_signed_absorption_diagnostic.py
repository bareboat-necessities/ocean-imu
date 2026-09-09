#!/usr/bin/env python3
"""Homogeneous exact-reset sector for the signed P4 Joseph ledger.

For one accepted shipping correction let d be the injected attitude correction,
a the Cayley vector of the *deployed* normalized correction quaternion, and c
the pre-injection true-minus-estimated Cayley error.  Exact composition gives

  c+ = (c-a+.5 a x c)/(1+.25 a'c)

while shipping transports covariance with G(d)=I+.5[d]x.  Therefore

  e_reset = G(d)t + rho,             ||G(d)^-1||_2 = 1.

The useful property for the augmented LDLT is not an absolute radius for rho,
but a homogeneous graph sector rho=O(d).  On ||c||<=q and ||d||<=delta_bar,
write

  ||a|| <= alpha ||d||,    ||a-d|| <= eps_a ||d||.

The polynomial small-correction branch and the axis-angle branch both satisfy
these inequalities with the constants constructed below.  Substitution in the
exact Cayley composition bound yields

  ||rho_theta|| <= mu_R ||d||,

  mu_R = [eps_a(1+q/2)
          + alpha q/4 (q+delta_bar+delta_bar q/2)]
         / [1-alpha delta_bar q/4].

Thus the exact reset defect can be retained in the SAME Joseph augmented
coordinate with d=K y.  It is not an independent packet disturbance and no
event-count multiplier is introduced.  Rowwise K bounds are consumed only to
obtain the source-uniform correction ceiling delta_bar; production storage must
still use the same-cell K/S/Joseph directions.

This producer is deliberately non-promoting: it supplies the reset IQC needed
by the endpoint/every-prefix augmented LDLT, but does not assert that those
matrices have closed.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path

import ou3_p4_rowwise_coefficient_enclosure_fast as COEFF
import ou3_p4_hard_entry_set as ENTRY
import ou3_p4_exact_reset_transport as RESET

QUALIFICATION='OU3_P4_RESET_HOMOGENEOUS_SAME_CELL_SECTOR_V2'


def up(x): return math.nextafter(float(x), math.inf)
def down(x): return math.nextafter(float(x), -math.inf)


def correction_ratio_bounds(delta_bar: float) -> tuple[float,float]:
    """Return global alpha and eps_a on 0<=||d||<=delta_bar.

    Above the 1e-2 branch threshold the normalized shipping correction is the
    axis-angle quaternion, so ||a||/||d||=2 tan(r/2)/r and
    ||a-d||/||d||=that-1, both increasing for r in [0,pi).  Below the threshold
    the source polynomial bounds in RESET provide a direct coefficient hull.
    Taking the maximum of the two branch endpoint ratios covers the branch
    switch without assuming which side executes.
    """
    d=float(delta_bar)
    if not (math.isfinite(d) and 0.0 <= d <= RESET.CAYLEY_MONOTONE_NORM_MAX):
        raise ValueError('correction ceiling outside validated reset utility domain')
    ds=min(d,RESET.SERIES_BRANCH_NORM)
    if ds>0:
        _,_,series_coeff_hi=RESET._series_cayley_scalar_bounds(ds)
        # Recompute the lower coefficient exactly as the source helper does so
        # the relative |a-d| sector covers both sides of one.
        d2=up(ds*ds);d4=up(d2*d2)
        w_lo=down(1.0-up(d2/8.0));w_hi=up(1.0+up(d4/384.0))
        k_lo=down(0.5-up(d2/48.0));k_hi=up(0.5+up(d4/3840.0))
        coeff_lo=down((2.0*k_lo)/w_hi)
        coeff_hi=up((2.0*k_hi)/w_lo)
        alpha_series=max(1.0,float(series_coeff_hi),coeff_hi)
        eps_series=max(abs(coeff_lo-1.0),abs(coeff_hi-1.0))
    else:
        alpha_series=1.0;eps_series=0.0
    if d >= RESET.SERIES_BRANCH_NORM:
        a_hi,_=RESET._axis_cayley_scalar_bounds(d)
        alpha_axis=up(a_hi/d) if d>0 else 1.0
        eps_axis=up(max(0.0,alpha_axis-1.0))
    else:
        alpha_axis=1.0;eps_axis=0.0
    return up(max(alpha_series,alpha_axis)),up(max(eps_series,eps_axis))


def homogeneous_sector(q: float, delta_bar: float) -> dict:
    q=float(q);d=float(delta_bar)
    alpha,eps=correction_ratio_bounds(d)
    denom=down(1.0-up(0.25*alpha*d*q))
    if not denom>0.0:
        raise RuntimeError('uniform reset sector reaches Cayley composition antipode')
    target=up(q+d+0.5*d*q)
    num=up(eps*up(1.0+0.5*q) + up(0.25*alpha*q*target))
    mu=up(num/denom)
    # Endpoint absolute primitive is a useful consistency check only.  Since
    # the homogeneous derivation is a direct substitution into the same exact
    # formula, its mu*d envelope must dominate the primitive endpoint bound.
    absolute=RESET.reset_defect_bound(q,d)
    rho=float(absolute['reset_attitude_defect_norm_upper'])
    consistency=up(mu*d) >= rho if d>0 else rho==0.0
    return {
      'correction_norm_upper':d,'alpha_cayley_over_rotation_vector_upper':alpha,
      'relative_cayley_minus_rotation_vector_upper':eps,
      'cayley_composition_denominator_lower':denom,
      'reset_defect_over_correction_norm_upper':mu,
      'homogeneous_reset_sector':'||rho_theta|| <= mu_R ||d_theta||; d_theta=K_theta*y in the same Joseph cell',
      'endpoint_absolute_reset_defect_norm_upper':rho,
      'homogeneous_sector_dominates_endpoint_absolute_bound':consistency,
      'reset_inverse_operator_norm_upper':1.0,
      'chart_safe':bool(absolute['chart_safe']),
    }


def build():
    c=COEFF.build();e=ENTRY.build();bad={'coeff':COEFF.validate(c),'entry':ENTRY.validate(e)};bad={k:v for k,v in bad.items() if v}
    if bad:raise RuntimeError('reset sector prerequisites failed: '+repr(bad))
    q=float(e['coordinate_radii']['attitude_cayley_norm']);modes={}
    for mode in ('H18','A21'):
        delta=float(c['modes'][mode]['attitude_correction_norm_upper'])
        row=homogeneous_sector(q,delta)
        row['limiting_correction_event']=c['modes'][mode]['limiting_correction_event']
        modes[mode]=row
    return {'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','full_declared_entry_cayley_norm':q,
      'same_source_uniform_correction_magnitude_envelope_consumed':True,'exact_reset_transport_consumed':True,
      'reset_defect_is_homogeneous_in_same_cell_correction':True,'same_cell_d_equals_Ky_must_be_retained_downstream':True,
      'rowwise_K_used_only_for_correction_ceiling_not_storage':True,'reset_inverse_operator_norm_upper':1.0,'modes':modes,
      'covariance_inverse_scalarized_here':False,'packet_count_multiplier_used':False,'independent_reset_disturbance_used':False,
      'endpoint_augmented_LDLT_closed_here':False,'every_prefix_augmented_LDLT_closed_here':False,'P4_promoted_here':False}


def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    for k in ('same_source_uniform_correction_magnitude_envelope_consumed','exact_reset_transport_consumed',
              'reset_defect_is_homogeneous_in_same_cell_correction','same_cell_d_equals_Ky_must_be_retained_downstream',
              'rowwise_K_used_only_for_correction_ceiling_not_storage'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('covariance_inverse_scalarized_here','packet_count_multiplier_used','independent_reset_disturbance_used',
              'endpoint_augmented_LDLT_closed_here','every_prefix_augmented_LDLT_closed_here','P4_promoted_here'):
        if d.get(k) is not False:f.append(k+' not false')
    if d.get('reset_inverse_operator_norm_upper')!=1.0:f.append('reset inverse norm changed')
    for mode,m in d.get('modes',{}).items():
        if m.get('chart_safe') is not True:f.append(mode+' chart unsafe')
        if m.get('homogeneous_sector_dominates_endpoint_absolute_bound') is not True:f.append(mode+' homogeneous/absolute inconsistency')
        for k in ('correction_norm_upper','alpha_cayley_over_rotation_vector_upper','relative_cayley_minus_rotation_vector_upper',
                  'reset_defect_over_correction_norm_upper','cayley_composition_denominator_lower'):
            if not (math.isfinite(float(m.get(k,math.nan))) and float(m[k])>=0):f.append(mode+' '+k+' invalid')
        if not float(m['cayley_composition_denominator_lower'])>0:f.append(mode+' denominator nonpositive')
    return f


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'modes':d['modes'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
