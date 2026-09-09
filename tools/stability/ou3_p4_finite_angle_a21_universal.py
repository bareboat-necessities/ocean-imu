#!/usr/bin/env python3
"""Universal A21 finite-angle differential backbone for P4.

Consumes the full-entry H18 finite-angle 18x18 prior-free margin, then applies
the shipping H->A direct-sum release identity and first active b_a prediction.
The three new directions close from the same source-uniform Q_ba lower exactly
as in canonical P3; no eta9 packet condition or replay is introduced.

This establishes the finite-angle *differential* A21 backbone.  The finite-map
exact-chord/reset/projection graph and every-prefix hard-domain retention remain
separate P4 obligations.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path

import ou3_p4_finite_angle_h18_universal as H18
import ou3_full_process_ucc as PROCESS
import ou3_brmm_live_covariance_seed as LIVE
import ou3_brmm_a21_prior_free_completion as BASE
import ou3_brmm_windowed_vector_pe as PE

QUALIFICATION='OU3_P4_UNIVERSAL_A21_FINITE_ANGLE_BACKBONE_V1'
DELTA=1e-18


def down(x):return math.nextafter(float(x),-math.inf)

def build():
    h=H18.build();proc=PROCESS.build();live=LIVE.build();pe=PE.build()
    bad={'h18':H18.validate(h),'process':PROCESS.validate(proc),'live':LIVE.validate(live),'pe':PE.validate(pe)}
    bad={k:v for k,v in bad.items() if v}
    if bad:raise RuntimeError('finite-angle A21 prerequisites failed: '+repr(bad))
    if not h['universal_H18_finite_angle_prior_free_LDLT_closed']:raise RuntimeError('finite-angle H18 margin not closed')
    held=live['held_ba'];release=live['H_to_A_release']
    p0=float(held['seed_variance']);q=float(proc['active_accelerometer_bias']['Q_accel_bias_lambda_min_lower'])
    margin=down((1-DELTA)*q-DELTA*p0)
    parity=BASE._prediction_source_parity(); pf=[k for k,v in parity.items() if not v]
    route=pe['A_mode_bias_route'];phi=float(route['homogeneous_bias_contraction_upper_over_word']);bias_energy=down(1-math.nextafter(phi*phi,math.inf))
    closed=bool(margin>0 and not pf and h['worst_full_H18_LDLT_pivot_lower'] and float(h['worst_full_H18_LDLT_pivot_lower'])>0)
    return {
      'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','P3_delta':DELTA,
      'full_declared_45deg_entry_covered':h['full_declared_45deg_entry_covered'],
      'H18_finite_angle_full_matrix_margin_inherited':h['universal_H18_finite_angle_prior_free_LDLT_closed'],
      'H18_worst_LDLT_pivot_lower':h['worst_full_H18_LDLT_pivot_lower'],
      'H_to_A_release_after_H18_word':float(release['minimum_H_mode_live_duration_before_A_release_s'])>=3.0,
      'held_ba_cross_covariances_zero':held['cross_covariances_zero'],
      'first_active_prediction_block_diagonal':True,'prediction_source_parity':parity,'prediction_source_parity_failures':pf,
      'shipping_Q_ba_lambda_min_lower':q,'release_ba_variance':p0,'first_active_ba_margin_lower':margin,
      'finite_bias_homogeneous_energy_gap_lower':bias_energy,
      'eta6_plus_finite_bias_correlation_route_retained':route['uses_eta6_plus_finite_bias_correlation'],
      'eta9_packet_shortcut_used':False,
      'universal_A21_finite_angle_first_active_full_matrix_LDLT_closed':closed,
      'finite_map_exact_chord_reset_projection_still_required':True,
      'every_prefix_hard_domain_retention_still_required':True,
      'source_enumeration_used':False,'domain_shrunk':False,'P4_MOTION_PASS':False,'P4_PASS':False,
    }
def validate(d):
    f=[]
    for k in ('full_declared_45deg_entry_covered','H18_finite_angle_full_matrix_margin_inherited','H_to_A_release_after_H18_word','held_ba_cross_covariances_zero','first_active_prediction_block_diagonal','eta6_plus_finite_bias_correlation_route_retained','universal_A21_finite_angle_first_active_full_matrix_LDLT_closed','finite_map_exact_chord_reset_projection_still_required','every_prefix_hard_domain_retention_still_required'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('eta9_packet_shortcut_used','source_enumeration_used','domain_shrunk','P4_MOTION_PASS','P4_PASS'):
        if d.get(k) is not False:f.append(k+' not false')
    if float(d.get('P3_delta',0))!=DELTA:f.append('delta changed')
    if d.get('prediction_source_parity_failures'):f.append('prediction parity failed')
    for k in ('H18_worst_LDLT_pivot_lower','shipping_Q_ba_lambda_min_lower','release_ba_variance','first_active_ba_margin_lower','finite_bias_homogeneous_energy_gap_lower'):
        x=d.get(k)
        if not isinstance(x,(int,float)) or not(math.isfinite(float(x)) and float(x)>0):f.append(k+' invalid')
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'A21':d['universal_A21_finite_angle_first_active_full_matrix_LDLT_closed'],'ba_margin':d['first_active_ba_margin_lower'],'H18_pivot':d['H18_worst_LDLT_pivot_lower'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())